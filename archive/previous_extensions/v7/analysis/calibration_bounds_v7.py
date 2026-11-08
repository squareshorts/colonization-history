"""Calibration-aware partial identification of four preparation-carried effects.

This is an additive, transformed-endpoint *synthetic benchmark*, not a bacterial
model fit or a verified intervention protocol. The four-arm transfer identifies
A+S and E+R. Two additional selective readouts have bounded sensitivities.
The method returns simultaneous-feasible sets via linear programming, and is
compared to a miscalibrated point-inversion procedure.

Run: python analysis/calibration_bounds_v7.py --replicates 300
"""
from __future__ import annotations
import argparse
import csv
import json
from dataclasses import dataclass, asdict
from pathlib import Path
import numpy as np
from functools import lru_cache
from scipy.stats import t

ALPHA = 0.05
N_MEANS = 4  # C, F, cellular assay and surface assay; joint Bonferroni CI

@dataclass(frozen=True)
class Scenario:
    name: str
    A: float
    S: float
    E: float
    R: float
    I: float
    qA: float
    qS: float
    qE: float
    qR: float
    outcome_AxS: float = 0.0
    outcome_ExR: float = 0.0
    assay_AxS: float = 0.0
    assay_ExR: float = 0.0

# Values are illustrative, dimensionless and deliberately include mechanisms
# that are absent; not inferred from any of the 2018 microscopy recordings.
SCENARIOS = (
    Scenario("absent_acquired_and_environment", 0., .8, 0., .7, -.12, .78, .18, .14, .67),
    Scenario("mixed_four_carriers", .24, .56, .18, .42, -.12, .78, .18, .14, .67),
    Scenario("all_carriers_absent", 0., 0., 0., 0., 0., .78, .18, .14, .67),
    Scenario("calibration_invalid_negative_control", 0., .8, 0., .7, -.12, .78, .28, .14, .48),
    Scenario("nonlinear_misspecification", .24, .56, .18, .42, -.12, .78, .18, .14, .67,
             outcome_AxS=1.2, outcome_ExR=1.5, assay_AxS=1.2, assay_ExR=1.8),
)

# Conservative *stipulated* calibration intervals:  no biological calibration
# study is available. One scenario intentionally violates these intervals.
BOUNDS = {"A": (.70, .85), "S": (0., .20),
          "E": (0., .18), "R": (.60, .90)}
# What an inadequately calibrated point-analysis could assume.
NOMINAL = {"A": .78, "S": .075, "E": .065, "R": .78}


@lru_cache(maxsize=20)
def _critical(n, multiple):
    return t.ppf(1 - ALPHA/(2*N_MEANS) if multiple else 1 - ALPHA/2, n-1)


def mean_ci(sample, multiple=True):
    """Student-t marginal intervals at simultaneous Bonferroni coverage."""
    x = np.asarray(sample, float)
    if len(x) < 3:
        raise ValueError("At least three independent fields are required")
    sd = x.std(ddof=1)
    half = _critical(len(x), multiple)*sd/np.sqrt(len(x))
    return float(x.mean()-half), float(x.mean()+half)


def feasible_pair(total_ci, assay_ci, coefficients_1, coefficients_2, upper=3.):
    """Exact 2-D convex polygon clipping for linear partial identification.

    Given T=x1+x2 and y=q1*x1+q2*x2, qj in [lo,hi],
    union over allowed sensitivities has the linear constraints
        q1_lo*x1+q2_lo*x2 <= y_hi
        q1_hi*x1+q2_hi*x2 >= y_lo.
    The positivity and max effect restrictions are explicit. The polygon
    vertices provide sharp finite-dimensional bounds without optimization.
    """
    tlo, thi = total_ci
    ylo, yhi = assay_ci
    a1, b1 = coefficients_1
    a2, b2 = coefficients_2
    poly = [(0.,0.),(upper,0.),(upper,upper),(0.,upper)]
    constraints = ((1.,1.,thi),(-1.,-1.,-tlo),
                   (a1,a2,yhi),(-b1,-b2,-ylo))
    for a, b, bound in constraints:
        if not poly:
            return None
        result=[]
        old=poly[-1]
        vold=a*old[0]+b*old[1]-bound
        for new in poly:
            vnew=a*new[0]+b*new[1]-bound
            was_in=vold <= 1e-10
            is_in=vnew <= 1e-10
            if is_in != was_in:
                ratio=vold/(vold-vnew)
                inter=(old[0]+ratio*(new[0]-old[0]),
                       old[1]+ratio*(new[1]-old[1]))
                result.append(inter)
            if is_in:
                result.append(new)
            old,vold=new,vnew
        poly=result
    if not poly:
        return None
    return [(float(min(p[j] for p in poly)),float(max(p[j] for p in poly)))
            for j in range(2)]


def point_estimate(total, assay, coeff1, coeff2):
    q1,q2=coeff1,coeff2
    if abs(q1-q2)<1e-9:
        raise ValueError("Assay selectivity is insufficient")
    return (assay-q2*total)/(q1-q2), (q1*total-assay)/(q1-q2)


def point_interval(total_samples, assay_samples, coeff1, coeff2):
    x=np.asarray(total_samples,float)
    y=np.asarray(assay_samples,float)
    q1,q2=coeff1,coeff2
    est,_=point_estimate(x.mean(),y.mean(),q1,q2)
    # In this synthetic benchmark the assays are generated independently of the
    # four arm outcomes; for paired assay data this SE needs the covariance.
    v=(y.var(ddof=1)/len(y)+q2*q2*x.var(ddof=1)/len(x))/(q1-q2)**2
    half=_critical(min(len(x),len(y)), False)*np.sqrt(v)
    return est,float(est-half),float(est+half)


def simulate_once(s:Scenario,n:int,rng,field_sd=.26,assay_sd=.12):
    """Blockwise observations; shared effects cancel in within-block contrasts."""
    shared=rng.normal(0.,.28,n)
    err=rng.normal(0.,field_sd,(n,4))
    base=.2+shared+err[:,0]
    cell_total=s.A+s.S+s.outcome_AxS*s.A*s.S
    surface_total=s.E+s.R+s.outcome_ExR*s.E*s.R
    c=.2+shared+cell_total+err[:,1]
    f=.2+shared+surface_total+err[:,2]
    both=.2+shared+cell_total+surface_total+s.I+err[:,3]
    total_c=c-base
    total_f=f-base
    interaction=both-c-f+base
    assay_c=rng.normal(s.qA*s.A+s.qS*s.S+s.assay_AxS*s.A*s.S,assay_sd,n)
    assay_f=rng.normal(s.qE*s.E+s.qR*s.R+s.assay_ExR*s.E*s.R,assay_sd,n)
    return total_c,total_f,interaction,assay_c,assay_f


def experiment(s,n,rng,calibration=BOUNDS):
    C,F,I,LC,LF=simulate_once(s,n,rng)
    cell=feasible_pair(mean_ci(C),mean_ci(LC),calibration["A"],calibration["S"])
    surf=feasible_pair(mean_ci(F),mean_ci(LF),calibration["E"],calibration["R"])
    pc=point_interval(C,LC,NOMINAL["A"],NOMINAL["S"])
    pe=point_interval(F,LF,NOMINAL["E"],NOMINAL["R"])
    return {"set_A": cell[0] if cell else None, "set_S": cell[1] if cell else None,
            "set_E": surf[0] if surf else None,"set_R": surf[1] if surf else None,
            "point_A":pc,"point_E":pe,"I_mean":float(I.mean()),
            "nonempty":cell is not None and surf is not None}


def aggregate(s,n,replicates,rng):
    stats={k:[] for k in ("set_A","set_S","set_E","set_R")}
    point={"A":[],"E":[]};feasible=0
    for _ in range(replicates):
        e=experiment(s,n,rng)
        feasible+=int(e["nonempty"])
        for k in stats:
            v=e[k]
            stats[k].append(v)
        for k in point:
            point[k].append(e["point_"+k])
    out={"scenario":s.name,"fields_per_arm":n,"replicates":replicates,
         "feasible_rate":feasible/replicates}
    for k,x in (("A",s.A),("S",s.S),("E",s.E),("R",s.R)):
        arr=[v for v in stats["set_"+k] if v is not None]
        out[f"set_{k}_coverage_conditional_on_feasibility"]=float(np.mean([lo-1e-9<=x<=hi+1e-9 for lo,hi in arr])) if arr else None
        out[f"set_{k}_positive_fraction"]=float(np.mean([v is not None and v[0]>0 for v in stats["set_"+k]]))
        out[f"set_{k}_mean_width"]=float(np.mean([hi-lo for lo,hi in arr])) if arr else None
    for k,x in (("A",s.A),("E",s.E)):
        out[f"point_{k}_positive_fraction"]=float(np.mean([v[1]>0 for v in point[k]]))
        out[f"point_{k}_coverage"]=float(np.mean([v[1]<=x<=v[2] for v in point[k]]))
        out[f"point_{k}_bias"]=float(np.mean([v[0]-x for v in point[k]]))
    return out


def run(outdir=Path("results/calibration_bounds_v7"),replicates=300,ns=(16,64,160),seed=20261106):
    outdir=Path(outdir);outdir.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(seed)
    result=[aggregate(s,n,replicates,rng) for s in SCENARIOS for n in ns]
    with (outdir/"calibration_bounds_v7.json").open("w") as f:
        json.dump({"seed":seed,"bounded_coefficients":BOUNDS,"nominal_coefficients":NOMINAL,
                   "scenarios":[asdict(s) for s in SCENARIOS],
                   "independent_fields_per_arm":list(ns),"replicates":replicates,
                   "intervals":"Bonferroni simultaneous Student-t for four measured means; assay coefficient bounds treated as certain",
                   "result":result},f,indent=2)
    with (outdir/"calibration_bounds_v7.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax=plt.subplots(figsize=(9,5))
    q=[r for r in result if r["scenario"]=="absent_acquired_and_environment"]
    xs=[r["fields_per_arm"] for r in q]
    for name,key,ls in (("Point calibration: false A", "point_A_positive_fraction", "--"),
                        ("Bounded calibration: false A", "set_A_positive_fraction", "-"),
                        ("Point calibration: false E", "point_E_positive_fraction", ":"),
                        ("Bounded calibration: false E", "set_E_positive_fraction", "-.")):
        ax.plot(xs,[r[key] for r in q],label=name,linestyle=ls,marker="o")
    ax.set(xlabel="Independent fields per arm / assay",ylabel="Fraction of false positive 95% lower bounds",
           ylim=(-.025,1.025),title="Illustrative calibration error can dominate replication")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(outdir/"false_attribution_calibration.png",dpi=200)
    plt.close(fig)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--n",type=int,nargs="+",default=[16,64,160])
    p.add_argument("--seed",type=int,default=20261106)
    p.add_argument("--outdir",type=Path,default=Path("results/calibration_bounds_v7"))
    a=p.parse_args()
    result=run(a.outdir,a.replicates,a.n,a.seed)
    for r in result:
        print(f"{r['scenario']:37s} n={r['fields_per_arm']:3d} ",
              f"point_positive_A={r['point_A_positive_fraction']:.3f} bounds_positive_A={r['set_A_positive_fraction']:.3f}",
              f"point_positive_E={r['point_E_positive_fraction']:.3f} bounds_positive_E={r['set_E_positive_fraction']:.3f}")
