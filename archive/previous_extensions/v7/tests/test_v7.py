from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"empirical"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/"analysis"))
from field_level_v7 import km_at, summarize, exact_sign_flip
from calibration_bounds_v7 import feasible_pair, point_estimate, experiment, SCENARIOS, BOUNDS, mean_ci
from import_archival_field_v7 import P
import numpy as np


def test_kaplan_meier_event_and_censoring():
    tracks=[(0,2,False),(0,3,True),(0,200,False)]
    r=km_at(tracks,threshold=2,residence_s=300)
    assert r["n"]==3
    assert abs(r["survival"]-2/3)<1e-12


def test_operational_attachment_threshold():
    tracks=[(0,2,False),(0,3,False),(0,100,True)]
    assert abs(km_at(tracks,threshold=2,residence_s=6)["survival"] - 2/3) < 1e-12
    assert km_at(tracks,threshold=3,residence_s=6)["survival"] == 1.


def test_exact_signflip_p_and_interval():
    values=[.13618714812657722,.03846889952153122,-.12626262626262624]
    assert exact_sign_flip(values)==.75
    s=summarize(values)
    assert abs(s["mean"]-.01613)<.0001
    assert s["ci95_t"][0]<0<s["ci95_t"][1]


def test_matrix_inversion():
    true=(.23,.77)
    total=sum(true)
    y=.78*true[0]+.08*true[1]
    x=point_estimate(total,y,.78,.08)
    assert np.allclose(true,x)


def test_partial_identification_includes_truth():
    q=SCENARIOS[1]
    cell=feasible_pair((q.A+q.S,q.A+q.S),(q.qA*q.A+q.qS*q.S,)*2,BOUNDS["A"],BOUNDS["S"])
    surf=feasible_pair((q.E+q.R,q.E+q.R),(q.qE*q.E+q.qR*q.R,)*2,BOUNDS["E"],BOUNDS["R"])
    assert cell is not None and surf is not None
    assert cell[0][0]<=q.A<=cell[0][1]
    assert surf[0][0]<=q.E<=surf[0][1]


def test_true_absent_mechanisms_not_excluded_by_correct_calibration():
    q=SCENARIOS[0]
    cell=feasible_pair((q.A+q.S,)*2,(q.qA*q.A+q.qS*q.S,)*2,BOUNDS["A"],BOUNDS["S"])
    surf=feasible_pair((q.E+q.R,)*2,(q.qE*q.E+q.qR*q.R,)*2,BOUNDS["E"],BOUNDS["R"])
    assert cell is not None and surf is not None
    assert abs(cell[0][0])<1e-9
    assert abs(surf[0][0])<1e-9


def test_invalid_calibration_can_force_wrong_carrier():
    q=SCENARIOS[3]
    # A=0, but assay signal is too large for stipulated maximal cross-reactivity.
    cell=feasible_pair((q.A+q.S,)*2,(q.qA*q.A+q.qS*q.S,)*2,BOUNDS["A"],BOUNDS["S"])
    assert cell is not None
    assert cell[0][0] > 0


def test_field_aggregate_against_connector_values():
    d = [P[x][2][3]-P[x][2][1] for x in P]
    assert abs(sum(d)/len(d)-.01613114)<1e-6


def test_polygon_bounds_against_linear_programming():
    from scipy.optimize import linprog
    rng=np.random.default_rng(9123)
    for _ in range(50):
        tlo,thi=sorted(rng.uniform(-.1,2.0,2))
        ylo,yhi=sorted(rng.uniform(-.1,1.5,2))
        a1,a2=(.7,.1);b1,b2=(.85,.2)
        r=feasible_pair((tlo,thi),(ylo,yhi),(a1,b1),(a2,b2))
        aa=[[1,1],[-1,-1],[a1,a2],[-b1,-b2]]
        bb=[thi,-tlo,yhi,-ylo]
        results=[]
        for j in range(2):
            c=np.eye(2)[j]
            lo=linprog(c,A_ub=aa,b_ub=bb,bounds=(0,3),method='highs')
            hi=linprog(-c,A_ub=aa,b_ub=bb,bounds=(0,3),method='highs')
            results.append((lo,hi))
        if results[0][0].success:
            assert r is not None
            for j,(lo,hi) in enumerate(results):
                assert abs(r[j][0]-lo.fun)<1e-7
                assert abs(r[j][1]+hi.fun)<1e-7
        else:
            assert r is None
