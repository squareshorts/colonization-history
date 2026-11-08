import sys,time
sys.path.insert(0,'src')
from minimal_agents import *
if __name__=='__main__':
 t=time.perf_counter()
 p=Parameters()
 for name in CAPABILITIES:
  r=replicated(p,CAPABILITIES[name],reps=16,seed=817,model=name)
  print(name,[(k,round(v['lag_threshold_h']['mean'] or -1,2),v['lag_threshold_h']['n_observed'],round(v['deep_fraction']['mean'],2)) for k,v in summarize(r).items()],flush=True)
 print('seconds',time.perf_counter()-t)
