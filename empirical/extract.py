"""Observable extraction from the Lee et al. (2018) binary image deposit.

Input: one dataset zip from figshare 10.6084/m9.figshare.5969335, containing
binary/Frames XXXXX-YYYYY.mat stacks (variable imageStack, [y, x, frame],
bright-field frames every 3 s, 300 frames per file).

Output (per dataset):
  counts.csv     every COUNT_STRIDE frames: time (h), attached components,
                 foreground area (px)
  tracks.csv     full-rate tracks inside tracking windows: start/end frame,
                 origin (arrival or split), censoring flags
  meta.json      frame count, single-cell area estimate, windows, settings

Tracking links connected components between consecutive frames by nearest
centroid within LINK_RADIUS pixels, closing gaps of up to GAP_FRAMES frames
to suppress segmentation flicker. A component without a predecessor is an
arrival; a second component claiming the same predecessor is a split
(division or separation of touching cells). A track ends when no successor
appears within the gap allowance.
"""
import argparse, csv, json, re, sys, time, zipfile, io
from pathlib import Path
import numpy as np
import scipy.io as sio
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

FRAME_SECONDS = 3.0
COUNT_STRIDE = 10                 # 30 s
LINK_RADIUS = 20.0                # px (about 1.3 um at 0.065 um/px)
GAP_FRAMES = 2
MIN_AREA = 40                     # px; smaller components are debris
WINDOW_HOURS = 1.0
WINDOW_STARTS_H = (0.0, 1.0, 2.0) + tuple(float(h) for h in range(4, 49, 4))
STRUCT = np.ones((3, 3), bool)


def frame_files(z):
    names = [n for n in z.namelist() if "/binary/" in n and n.endswith(".mat")]
    key = lambda n: int(re.search(r"Frames (\d+)-", n).group(1))
    return sorted(names, key=key)


def components(img):
    lab, n = ndi.label(img > 0, structure=STRUCT)
    if n == 0:
        return np.zeros((0, 2)), np.zeros(0, int)
    area = np.bincount(lab.ravel(), minlength=n + 1)[1:]
    yy, xx = np.nonzero(lab)
    l = lab[yy, xx] - 1
    cy = np.bincount(l, weights=yy, minlength=n) / np.maximum(area, 1)
    cx = np.bincount(l, weights=xx, minlength=n) / np.maximum(area, 1)
    keep = area >= MIN_AREA
    return np.column_stack([cy, cx])[keep], area[keep]


class Tracker:
    def __init__(self):
        self.active = {}      # track id -> [centroid, last_frame]
        self.records = {}     # track id -> dict
        self.next_id = 0

    def _new(self, frame, origin, parent=-1):
        tid = self.next_id
        self.next_id += 1
        self.records[tid] = {"track": tid, "start": frame, "end": frame, "origin": origin,
                             "parent": parent, "left_censored": False, "right_censored": False}
        return tid

    def step(self, frame, cents, first=False):
        ids = list(self.active)
        assigned = np.full(len(cents), -1)
        if ids and len(cents):
            prev = np.array([self.active[i][0] for i in ids])
            tree = cKDTree(prev)
            d, j = tree.query(cents, distance_upper_bound=LINK_RADIUS)
            order = np.argsort(d)
            claimed = set()
            for k in order:
                if not np.isfinite(d[k]):
                    break
                tid = ids[j[k]]
                if tid not in claimed:
                    assigned[k] = tid
                    claimed.add(tid)
                else:
                    assigned[k] = self._new(frame, "split", parent=tid)
        for k in np.flatnonzero(assigned < 0):
            assigned[k] = self._new(frame, "present" if first else "arrival")
            if first:
                self.records[assigned[k]]["left_censored"] = True
        for k, tid in enumerate(assigned):
            self.active[tid] = [cents[k], frame]
            self.records[tid]["end"] = frame
        for tid in list(self.active):
            if frame - self.active[tid][1] > GAP_FRAMES:
                del self.active[tid]

    def close(self, last_frame):
        for tid in self.active:
            if last_frame - self.records[tid]["end"] <= GAP_FRAMES:
                self.records[tid]["right_censored"] = True
        self.active = {}


def process(zip_path, out_dir, max_frames=None):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    z = zipfile.ZipFile(zip_path)
    files = frame_files(z)
    per_hour = int(3600 / FRAME_SECONDS)
    windows = [(int(h * per_hour), int((h + WINDOW_HOURS) * per_hour)) for h in WINDOW_STARTS_H]
    counts, tracks, single_areas = [], [], []
    tracker, in_window, win_idx = None, False, -1
    frame = 0
    t0 = time.perf_counter()
    for name in files:
        stack = sio.loadmat(io.BytesIO(z.read(name)))["imageStack"]
        for i in range(stack.shape[2]):
            w = next((k for k, (a, b) in enumerate(windows) if a <= frame < b), -1)
            need_track = w >= 0
            if need_track or frame % COUNT_STRIDE == 0:
                cents, area = components(stack[:, :, i])
                if frame % COUNT_STRIDE == 0:
                    counts.append({"frame": frame, "time_h": frame * FRAME_SECONDS / 3600,
                                   "components": int(len(area)), "area_px": int(area.sum())})
                    if frame < 2 * per_hour and len(area):
                        single_areas.extend(area.tolist())
                if need_track:
                    if w != win_idx:
                        if tracker is not None:
                            tracker.close(frame - 1)
                            tracks.extend({**r, "window": win_idx} for r in tracker.records.values())
                        tracker, win_idx = Tracker(), w
                        tracker.step(frame, cents, first=True)
                    else:
                        tracker.step(frame, cents)
            elif tracker is not None:
                tracker.close(frame - 1)
                tracks.extend({**r, "window": win_idx} for r in tracker.records.values())
                tracker, win_idx = None, -1
            frame += 1
            if max_frames and frame >= max_frames:
                break
        if max_frames and frame >= max_frames:
            break
    if tracker is not None:
        tracker.close(frame - 1)
        tracks.extend({**r, "window": win_idx} for r in tracker.records.values())
    for fn, rows in (("counts.csv", counts), ("tracks.csv", tracks)):
        with (out_dir / fn).open("w", newline="") as f:
            if rows:
                w = csv.DictWriter(f, fieldnames=list(rows[0]))
                w.writeheader(); w.writerows(rows)
    sa = np.array(single_areas) if single_areas else np.array([np.nan])
    meta = {"zip": Path(zip_path).name, "frames": frame, "hours": frame * FRAME_SECONDS / 3600,
            "frame_seconds": FRAME_SECONDS, "count_stride": COUNT_STRIDE,
            "link_radius_px": LINK_RADIUS, "gap_frames": GAP_FRAMES, "min_area_px": MIN_AREA,
            "windows_frames": windows, "median_component_area_first2h": float(np.nanmedian(sa)),
            "n_tracks": len(tracks), "seconds": time.perf_counter() - t0}
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("zip")
    p.add_argument("out")
    p.add_argument("--max-frames", type=int, default=None)
    a = p.parse_args()
    print(json.dumps(process(a.zip, a.out, a.max_frames), indent=1))
