"""Minimal capability ablations for bacterial surface-colonization history.

Exploratory agents. Parameters are NOT fitted to Lee et al. 2018 raw microscopy.
The literature reports an exponential-lag definition different from the
first-passage threshold reported by this simulator. The comparison must remain
qualitative unless the observation and likelihood models are matched.

Agent state: (one continuous intracellular surface-exposure memory m, lineage
founder id, generation). A stable, optional, inherited adhesiveness z is a
comparator. A dynamic, optional, local extracellular signal S is a comparator.
The simulation represents only persistent visible surface occupants; unstable
contacts have explicitly modeled short dwell times used for the transient-visit
check. Each independent realization is a simulated experimental field.
"""
from dataclasses import dataclass,replace,asdict
import numpy as np

@dataclass(frozen=True)
class Parameters:
    horizon_h: float=40.0
    dt_h: float=0.10
    inflow_per_h: float=25.0
    carrying_capacity: int=600
    threshold_cells: int=100
    division_rate_h: float=0.80
    detach_naive_h: float=0.14
    detach_memory_h: float=0.035
    stable_attach_naive: float=0.0035
    stable_attach_memory: float=0.085
    daughter_ret_naive: float=0.65
    daughter_ret_memory: float=0.98
    contact_acquisition: float=0.55
    on_surface_acquisition_h: float=0.75
    free_memory_tau_h: float=10.0
    memory_midpoint: float=0.115
    memory_hill: float=7.0
    donor_contacts: int=4
    transient_dwell_mean_seconds: float=10.0
    stop_at_lag_and_8h: bool=True
    # Uncalibrated comparator rates
    stable_selection_prevalence: float=0.80
    selection_background_prevalence: float=0.015
    selection_strength: float=0.085
    surface_deposition_h: float=0.12
    surface_decay_h: float=0.10
    surface_attachment_strength: float=0.15

@dataclass(frozen=True)
class Capabilities:
    contact_sensing: bool=True
    internal_persistence: bool=True
    adhesion_feedback: bool=True
    daughter_inheritance: bool=True
    immutable_selection: bool=False
    environment: bool=False
    ongoing_sensing: bool=True
    binary_memory: bool=False

CAPABILITIES={
    'passive':Capabilities(False,False,False,False,False,False,False),
    'selection':Capabilities(False,False,False,False,True,False,False),
    'cell_no_inheritance':Capabilities(True,True,True,False,False,False,True),
    'cell_no_ongoing_sensing':Capabilities(True,True,True,True,False,False,False),
    'cell_memory':Capabilities(True,True,True,True,False,False,True),
    'cell_one_bit':Capabilities(True,True,True,True,False,False,True,True),
    'cell_plus_environment':Capabilities(True,True,True,True,False,True,True),
    'environment':Capabilities(False,False,False,False,False,True,False),
    'cell_no_decay':Capabilities(True,False,True,True,False,False,True),
    'cell_linear_feedback':Capabilities(True,True,True,True,False,False,True),
}

def sigmoid_hill(x,p,linear=False):
    if linear:
        return np.clip(np.asarray(x,dtype=float),0,1)
    x=np.clip(np.asarray(x,dtype=float),0,1)
    # Symmetric Hill response; h(0)=0 even with low midpoint.
    z=np.power(x,p.memory_hill)
    return z/(z+p.memory_midpoint**p.memory_hill)

def initial_memory(hold_h,p,c):
    if not(c.contact_sensing and c.adhesion_feedback):
        return 0.
    fresh_memory=1.-(1.-p.contact_acquisition)**p.donor_contacts
    if c.internal_persistence:
        return float(fresh_memory*np.exp(-hold_h/p.free_memory_tau_h))
    return float(fresh_memory)

def run_field(p:Parameters,c:Capabilities,hold_h:float,seed:int,transfer:bool=True,
              fresh_surface:bool=True,store_curve:bool=False,linear_feedback:bool=False):
    rng=np.random.default_rng(seed)
    # The same sequence of time-varying arrival conditions for naive and exposed.
    # Transfer is a cell property; all main assays use a virgin challenge surface.
    m_init=initial_memory(hold_h,p,c) if transfer else 0.
    selected=(c.immutable_selection and transfer)
    m=np.zeros(0,dtype=float)
    gen=np.zeros(0,dtype=np.int16)
    founder=np.zeros(0,dtype=np.int32)
    z=np.zeros(0,dtype=float)
    next_id=0
    time_lag=float('nan')
    n_events=0
    both_retained=0
    transient_total=0
    transient_lt30=0
    first_hour_present=[]
    n8=0; deep8=0.; maxgen8=0
    surface=0.70 if c.environment and not fresh_surface else 0.
    step_n=int(round(p.horizon_h/p.dt_h))
    curve=np.empty(step_n+1,dtype=float)
    curve[0]=0.
    for step in range(step_n):
        if c.environment:
            surface*=np.exp(-p.surface_decay_h*p.dt_h)
            surface=1-(1-surface)*np.exp(-p.surface_deposition_h*len(m)/p.carrying_capacity*p.dt_h)
        if len(m):
            h=sigmoid_hill(m,p,linear=linear_feedback) if c.adhesion_feedback else np.zeros_like(m)
            h=np.maximum(h,z) if c.immutable_selection else h
            pdeath=1-np.exp(-(p.detach_naive_h*(1-h)+p.detach_memory_h*h)*p.dt_h)
            keep=rng.random(len(m))>=pdeath
            m,gen,founder,z=m[keep],gen[keep],founder[keep],z[keep]
            if c.contact_sensing and c.ongoing_sensing and len(m):
                if c.binary_memory:
                    flip=(rng.random(len(m))<1-np.exp(-p.on_surface_acquisition_h*p.dt_h))
                    m=np.maximum(m,flip.astype(float))
                else:
                    m=1-(1-m)*np.exp(-p.on_surface_acquisition_h*p.dt_h)
            if len(m):
                n=len(m)
                divide=rng.random(n)<1-np.exp(-p.division_rate_h*p.dt_h)
                idx=np.flatnonzero(divide)
                n_events+=len(idx)
                if len(idx):
                    # Daughter stays according to mother's present state, while
                    # the mother remains. This records two-legged retention.
                    hh=sigmoid_hill(m[idx],p,linear=linear_feedback) if c.adhesion_feedback else np.zeros(len(idx))
                    hh=np.maximum(hh,z[idx]) if c.immutable_selection else hh
                    retention=p.daughter_ret_naive+(p.daughter_ret_memory-p.daughter_ret_naive)*hh
                    good=rng.random(len(idx))<retention
                    idx=idx[good]
                    both_retained+=len(idx)
                    if len(idx):
                        child_m=m[idx].copy() if c.daughter_inheritance else np.zeros(len(idx))
                        child_z=z[idx].copy() if c.immutable_selection else np.zeros(len(idx))
                        m=np.concatenate((m,child_m))
                        gen=np.concatenate((gen,gen[idx]+1))
                        founder=np.concatenate((founder,founder[idx]))
                        z=np.concatenate((z,child_z))
        # Age-independent inflow supplies founding agents, selected or not.
        narr=rng.poisson(p.inflow_per_h*p.dt_h)
        if narr:
            # Individual donor agents experience different numbers of prior contacts.
            incoming_m=np.full(narr,m_init)
            if transfer and c.contact_sensing and c.adhesion_feedback:
                visits=rng.poisson(p.donor_contacts,narr)
                incoming_m=1-np.power(1-p.contact_acquisition,visits)
                if c.internal_persistence:incoming_m*=np.exp(-hold_h/p.free_memory_tau_h)
                if c.binary_memory:incoming_m=(rng.random(narr)<incoming_m).astype(float)
            incoming_z=np.zeros(narr)
            if c.immutable_selection:
                pi=p.stable_selection_prevalence if selected else p.selection_background_prevalence
                incoming_z=(rng.random(narr)<pi).astype(float)
            h=sigmoid_hill(incoming_m,p,linear=linear_feedback) if c.adhesion_feedback else np.zeros(narr)
            h=np.maximum(h,incoming_z) if c.immutable_selection else h
            attach=p.stable_attach_naive+(p.stable_attach_memory-p.stable_attach_naive)*h
            if c.environment and (not fresh_surface or surface>0):
                attach=attach+(1-attach)*p.surface_attachment_strength*surface
            accepted=rng.random(narr)<attach
            if accepted.any():
                ids=np.arange(next_id,next_id+int(accepted.sum()),dtype=np.int32)
                next_id+=len(ids)
                m=np.concatenate((m,incoming_m[accepted]))
                gen=np.concatenate((gen,np.zeros(len(ids),dtype=np.int16)))
                founder=np.concatenate((founder,ids))
                z=np.concatenate((z,incoming_z[accepted]))
            n_short=int((~accepted).sum())
            transient_total+=n_short
            transient_lt30+=int((rng.exponential(p.transient_dwell_mean_seconds,n_short)<30).sum())
        if len(m)>p.carrying_capacity:
            keep=rng.choice(len(m),p.carrying_capacity,replace=False)
            m,gen,founder,z=m[keep],gen[keep],founder[keep],z[keep]
        curve[step+1]=len(m)
        t=(step+1)*p.dt_h
        if t<=1+1e-9:first_hour_present.append(len(m))
        if abs(t-8.)<p.dt_h/3:
            n8=len(m);deep8=float(np.mean(gen>=3)) if len(gen) else 0.;maxgen8=int(gen.max()) if len(gen) else 0
        if np.isnan(time_lag) and len(m)>=p.threshold_cells:
            time_lag=t
        if p.stop_at_lag_and_8h and np.isfinite(time_lag) and t>=8.-1e-6:
            break
    deep=float(np.mean(gen>=3)) if len(gen) else 0.
    maxgen=int(np.max(gen)) if len(gen) else 0
    out=dict(lag_threshold_h=float(time_lag),
             restricted_threshold_h=float(time_lag) if np.isfinite(time_lag) else float(p.horizon_h),
             threshold_reached=bool(np.isfinite(time_lag)),
             short_fraction_all_contacts=(transient_lt30/(transient_total+next_id) if transient_total+next_id else float('nan')),
             final_count=int(len(m)),
             early_count=float(np.mean(first_hour_present)) if first_hour_present else 0.,
             deep_fraction=deep, max_generation=maxgen,
             two_leg_retention=(both_retained/n_events if n_events else float('nan')),
             transient_lt30=(transient_lt30/transient_total if transient_total else float('nan')),
             n_divisions=n_events, founder_count=next_id, hold_h=float(hold_h),
             n_at_8h=n8,deep_fraction_8h=deep8,max_generation_8h=maxgen8)
    if store_curve:
        out['curve']=curve[:step+2]
    return out


def daughter_transfer_probability(p,c,parent_memory=0.8,linear_feedback=False):
    """Predicted immediate attachment on a fresh surface BEFORE reacquisition.

    For an acquired-state model, this is an ideal selective daughter assay;
    it is not observed in the published microscopy archive.
    """
    child_state=parent_memory if c.daughter_inheritance else 0.0
    h=float(sigmoid_hill(child_state,p,linear=linear_feedback)) if c.adhesion_feedback else 0.0
    return p.stable_attach_naive+(p.stable_attach_memory-p.stable_attach_naive)*h

def replicated(p,c,holds=(0.,14.,37.),reps=30,seed=1000,model='cell_memory',store_curve=False):
    scenarios=[('fresh',0.,False)]+[(f'transfer_{int(t)}h',t,True) for t in holds]
    result={}
    for j,(label,t,tr) in enumerate(scenarios):
        runs=[run_field(p,c,t,seed+100003*j+17*k,transfer=tr,store_curve=store_curve,
                        linear_feedback=model=='cell_linear_feedback') for k in range(reps)]
        result[label]=runs
    return result

def summarize(res):
    out={}
    for key,runs in res.items():
        d={}
        for var in ('lag_threshold_h','restricted_threshold_h','threshold_reached','deep_fraction','early_count','short_fraction_all_contacts','transient_lt30','two_leg_retention','final_count','n_at_8h','deep_fraction_8h','max_generation_8h'):
            a=np.array([v[var] for v in runs],float)
            observed=a[np.isfinite(a)]
            d[var]={'mean':float(observed.mean()) if len(observed) else None,
                    'median':float(np.median(observed)) if len(observed) else None,
                    'q025':float(np.quantile(observed,.025)) if len(observed) else None,
                    'q975':float(np.quantile(observed,.975)) if len(observed) else None,
                    'n_observed':len(observed),'n_total':len(runs)}
        out[key]=d
    return out
