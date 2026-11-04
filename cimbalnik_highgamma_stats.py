# -*- coding: utf-8 -*-
"""
Reproducible statistics for the High Gamma event-rate validation (Cimbalnik).

Tests whether the High Gamma spectral-event RATE differs across the three event
classes (IN_FRAG, IN_NB_OUT, OUT_NB), paired by cell-line x DIV block.

Input : summary_3class_Tampere.csv  (one row per cellLine x DIV x class)
        columns: cellLine, DIV, class, n_events_HG, time_s, rate_evHz
Output: Friedman omnibus + pairwise Wilcoxon signed-rank with Holm correction,
        and per-class median event rates.

Run:  python3 cimbalnik_highgamma_stats.py
"""
import pandas as pd
import numpy as np
from scipy.stats import friedmanchisquare, wilcoxon

CSV = "summary_3class_Tampere.csv"   # adjust path if needed
ORDER = ["IN_FRAG", "IN_NB_OUT", "OUT_NB"]

df = pd.read_csv(CSV)

# wide: one row per (cellLine, DIV) block, one column per class = High Gamma event rate (Hz)
wide = df.pivot_table(index=["cellLine", "DIV"], columns="class", values="rate_evHz")[ORDER]
comp = wide.dropna()                  # keep only blocks with all three classes
n = len(comp)
print(f"Complete cell-line x DIV blocks: n = {n}")

a, b, c = comp["IN_FRAG"], comp["IN_NB_OUT"], comp["OUT_NB"]

# --- Friedman omnibus ---
chi, p = friedmanchisquare(a, b, c)
print(f"\nFriedman: chi2(2) = {chi:.3g}, p = {p:.3g}")

# --- pairwise Wilcoxon signed-rank (two-sided) + Holm correction ---
pairs = [("IN_FRAG", "IN_NB_OUT"), ("IN_FRAG", "OUT_NB"), ("IN_NB_OUT", "OUT_NB")]
praw = [wilcoxon(comp[x], comp[y]).pvalue for x, y in pairs]

# Holm step-down
m = len(praw)
o = np.argsort(praw)
holm = [0.0] * m
for rank, idx in enumerate(o):
    holm[idx] = min(1.0, praw[idx] * (m - rank))
for i in range(1, m):
    holm[o[i]] = max(holm[o[i]], holm[o[i - 1]])

print("\nPairwise Wilcoxon signed-rank (Holm-corrected):")
for (x, y), pr, hp in zip(pairs, praw, holm):
    print(f"  {x} vs {y}:  p_raw = {pr:.3g}   p_Holm = {hp:.3g}")

# --- per-class medians ---
print("\nMedian High Gamma event rate (events/s):")
for cl in ORDER:
    print(f"  {cl}: {comp[cl].median():.3g}")
