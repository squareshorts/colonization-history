"""Full-record tracking with per-cell cAMP readings (analysis plan v2).

Extends extract.py: tracks every frame from inoculation until the field holds more
than STOP_COMPONENTS components (tracks still active then are right-censored), and,
at every fluorescence frame, records background-subtracted YFP and CFP inside each
tracked cell's mask. Fluorescence frame k is registered to binary frame
300(k+1)-1 by FFT cross-correlation at one frame chosen per dataset.
"""
import argparse, csv, io, json, re, sys, time, zipfile
from pathlib import Path
import numpy as np
import scipy.io as sio
from scipy import ndimage as ndi
from numpy.fft import fft2, ifft2
sys.path.insert(0, str(Path(__file__).resolve().parent))
from extract import (FRAME_SECONDS, COUNT_STRIDE, MIN_AREA, STRUCT, Tracker, frame_files)

STOP_COMPONENTS = 250
FRAMES_PER_FILE = 300


def label(img):
    lab, n = ndi.label(img > 0, structure=STRUCT)
    if n == 0:
        return lab, np.zeros((0, 2)), np.zeros(0, int), np.zeros(0, int)
    area = np.bincount(lab.ravel(), minlength=n + 1)[1:]
    yy, xx = np.nonzero(lab)
    l = lab[yy, xx] - 1
    cy = np.bincount(l, weights=yy, minlength=n) / np.maximum(area, 1)
    cx = np.bincount(l, weights=xx, minlength=n) / np.maximum(area, 1)
    keep = np.flatnonzero(area >= MIN_AREA)
    return lab, np.column_stack([cy, cx])[keep], area[keep], keep + 1


def register(fluo, mask):
    """Integer shift (dy, dx) that best aligns mask onto fluorescence, and its correlation."""
    F = fluo.astype(float) - np.median(fluo)
    F[F < 0] = 0
    F /= F.std() + 1e-9
    M = mask.astype(float)
    M = (M - M.mean()) / (M.std() + 1e-9)
    cc = np.real(ifft2(fft2(F) * np.conj(fft2(M)))) / F.size
    k = np.unravel_index(cc.argmax(), cc.shape)
    dy = k[0] if k[0] < 512 else k[0] - 1024
    dx = k[1] if k[1] < 512 else k[1] - 1024
    return int(dy), int(dx), float(cc.max())


def intensities(lab, labels, yfp, cfp, shift):
    dy, dx = shift
    lab_s = np.roll(np.roll(lab, dy, 0), dx, 1)
    fg = ndi.binary_dilation(lab_s > 0, iterations=5)
    by, bc = np.median(yfp[~fg]), np.median(cfp[~fg])
    n = lab.max() + 1
    cnt = np.bincount(lab_s.ravel(), minlength=n)
    sy = np.bincount(lab_s.ravel(), weights=yfp.ravel().astype(float), minlength=n)
    sc = np.bincount(lab_s.ravel(), weights=cfp.ravel().astype(float), minlength=n)
    my = sy[labels] / np.maximum(cnt[labels], 1) - by
    mc = sc[labels] / np.maximum(cnt[labels], 1) - bc
    return my, mc


def process(zip_path, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    z = zipfile.ZipFile(zip_path)
    files = frame_files(z)
    fl_name = [n for n in z.namelist() if n.endswith("fluorescence.mat")]
    yfp = cfp = None
    reg = {"dy": 0, "dx": 0, "corr": None, "k": None}
    total = int(re.search(r"-(\d+)\.mat", files[-1]).group(1))
    step = None
    if fl_name:
        d = sio.loadmat(io.BytesIO(z.read(fl_name[0])))
        yfp = d["yfp"]
        cfp = d["cfp"] if "cfp" in d else d["dapi"]
        K = yfp.shape[2]
        step = int(round(total / K))
        ranges = [tuple(int(v) - 1 for v in re.search(r"Frames (\d+)-(\d+)", n).groups()) for n in files]
        def bin_frame(fr):
            j = next(i for i, (a, b) in enumerate(ranges) if a <= fr <= b)
            st = sio.loadmat(io.BytesIO(z.read(files[j])))["imageStack"]
            return st[:, :, fr - ranges[j][0]] > 0
        # Register on the fluorescence frames whose matching binary frame holds the most cells.
        cand = sorted({K - 1, (3 * K) // 4, K // 2, max(1, K // 4), 0})
        load = [(int(bin_frame(min((kk + 1) * step - 1, total - 1)).sum()), kk) for kk in cand]
        ks = [kk for _, kk in sorted(load, reverse=True)[:2]]
        best = None
        for k in ks:
            phases = range(step) if step <= 30 else (step - 1, 0)
            for fr in (k * step + p for p in phases):
                fr = min(fr, total - 1)
                dy, dx, cc = register(yfp[:, :, k], bin_frame(fr))
                if best is None or cc > best[2]:
                    best = (dy, dx, cc, k, fr - k * step)
        reg = {"dy": best[0], "dx": best[1], "corr": best[2], "k": int(best[3]),
               "phase": int(best[4]), "step_frames": step, "K": int(K)}
    use_fluo = yfp is not None and reg["corr"] is not None and reg["corr"] >= 0.3
    tracker = Tracker()
    reading = {}
    counts, frame, stopped = [], 0, None
    t0 = time.perf_counter()
    for fi, name in enumerate(files):
        stack = sio.loadmat(io.BytesIO(z.read(name)))["imageStack"]
        for i in range(stack.shape[2]):
            lab, cents, area, labels = label(stack[:, :, i])
            if frame % COUNT_STRIDE == 0:
                counts.append({"frame": frame, "time_h": frame * FRAME_SECONDS / 3600,
                               "components": int(len(area)), "area_px": int(area.sum())})
            if stopped is None:
                tracker.step(frame, cents, first=(frame == 0))
                kf = (frame - reg["phase"]) // step if use_fluo else -1
                if (use_fluo and (frame - reg["phase"]) % step == 0 and 0 <= kf < yfp.shape[2]
                        and len(labels)):
                    my, mc = intensities(lab, labels, yfp[:, :, kf], cfp[:, :, kf], (reg["dy"], reg["dx"]))
                    # tracker.active maps track -> [centroid, last frame]; match by centroid
                    for tid, (c, last) in tracker.active.items():
                        if last != frame or tid in reading:
                            continue
                        j = np.argmin(np.sum((cents - c) ** 2, axis=1))
                        reading[tid] = (frame, float(my[j]), float(mc[j]), int(area[j]))
                if len(area) > STOP_COMPONENTS:
                    stopped = frame
                    tracker.close(frame)
            frame += 1
    if stopped is None:
        tracker.close(frame - 1)
    rows = []
    for tid, r in tracker.records.items():
        rd = reading.get(tid)
        rows.append({**r, "read_frame": rd[0] if rd else "", "yfp": rd[1] if rd else "",
                     "cfp": rd[2] if rd else "", "read_area": rd[3] if rd else ""})
    with (out_dir / "tracks_full.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    meta = {"zip": Path(zip_path).name, "frames": frame, "tracking_stopped_frame": stopped,
            "stop_components": STOP_COMPONENTS, "registration": reg, "fluorescence_used": bool(use_fluo),
            "n_tracks": len(rows), "n_readings": len(reading), "seconds": time.perf_counter() - t0}
    (out_dir / "meta_full.json").write_text(json.dumps(meta, indent=1))
    return meta


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("zip"); p.add_argument("out")
    a = p.parse_args()
    print(json.dumps(process(a.zip, a.out), indent=1))
