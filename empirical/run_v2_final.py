"""Final v2 empirical run, resumable: each part saves its own JSON; 'combine' merges them.

    python run_v2_final.py <part> <B>     part in P1_first3h, P1_extended, P2_P3, S2,
                                          time_matched, occupancy_adjusted
    python run_v2_final.py combine
"""
import json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_v2 as V

OUT = Path(__file__).resolve().parent.parent / "results" / "empirical_v2"
PARTS = {"P1_first3h": lambda B: V.p1("first3h", B), "P1_extended": lambda B: V.p1("extended", B),
         "P2_P3": lambda B: V.p2p3(min(B, 100)), "S2": lambda B: V.s2(),
         "time_matched": lambda B: V.time_matched(B, "time"),
         "occupancy_adjusted": lambda B: V.time_matched(B, "occupancy")}

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    if sys.argv[1] == "combine":
        res = {"plan": "ANALYSIS_PLAN_v2.md"}
        for p in PARTS:
            res[p] = json.loads((OUT / f"part_{p}.json").read_text())
        (OUT / "fits_v2.json").write_text(json.dumps(res, indent=1))
        print("combined")
    else:
        part, B = sys.argv[1], int(sys.argv[2])
        f = OUT / f"part_{part}.json"
        if not f.exists():
            f.write_text(json.dumps(PARTS[part](B), indent=1))
        print("done", part)
