# Network-Burst Detection & Validation in hiPSC-Derived Neuronal Networks on MEAs

MATLAB pipeline accompanying the MSc thesis
*"Network burst detection and fragmentation analysis in hiPSC-derived neuronal
networks recorded on micro-electrode arrays"* (J. Canestrelli, Tampere
University).

The code (i) detects single-electrode bursts and **network bursts (NBs)** from
MEA spike data using three detectors — **Tampere V1**, **Tampere V2** and
**Doorn (2024)** — (ii) quantifies spontaneous-activity and NB-fragmentation
metrics across four hiPSC cell lines and several DIVs, (iii) runs the
statistical comparisons, and (iv) performs a **spike-independent spectral
validation** of the detected NBs based on the Cimbalník high-frequency-
oscillation framework.

---

## 1. Requirements

- **MATLAB** R2021b or newer.
  - *Signal Processing Toolbox* (`designfilt`, `filtfilt`, `findpeaks`,
    `hilbert`).
  - *Statistics and Machine Learning Toolbox* (`fitlme`, `ranksum`,
    `kruskalwallis`, `friedman`, etc.).
  - *Curve Fitting Toolbox* (`smooth`) for the log-ISI histograms.

No external MATLAB packages are required; all helper functions are included.
The pipeline is now MATLAB-only (the former Python decoder was rewritten in
MATLAB as `decode_validation_mat.m`).

---

## 2. Input data format

The pipeline expects, for each cell line, one `.mat` file (`DD1C.mat`,
`DD5A.mat`, `TUBA.mat`, `C10902.mat`) containing a single nested struct:

```
CellLine.(meaPlate).(DIV).(well) = table with columns:
    Time     spike times in seconds
    Channel  electrode label
```

The four cell lines are two Dravet-syndrome patient lines (DD1C, DD5A) and two
controls (TUBA, C10902). Recordings are 600 s, 16 electrodes per well.

---

## 3. How to run

Two master scripts reproduce the whole thesis, in order. Run them **as scripts**
(not as functions) in the **same MATLAB session**, so the data stays in the base
workspace that the figure/statistics functions read from.

```matlab
% 1) Pre-processing + processing: load data, FR filter, run the 3 detectors
run_preprocessing_processing      % edit the CONFIG block first

% 2) Metrics + figures + statistics
run_analysis_metrics_stats
```

`run_preprocessing_processing.m` — edit the `CONFIG` block: `DATA_ROOT`
(folder with the four cell-line `.mat` files), `RESULTS_ROOT` (output folder,
default `<pwd>/results`) and `PARAM_SET` (1–4; **Set 2 is the reference**, Sets
1/3/4 reproduce the parameter-sensitivity analysis of Chapter 4.4).

All outputs (figures, CSVs) are written under `RESULTS_ROOT`. Set `RESULTS_ROOT`
in the base workspace to override the default.

---

## 4. Parameter sets

`get_parameter_set(setNum)` returns the three detector parameter structs
(`Params` = Tampere NB, `Pdoorn` = Doorn, `ParamsBurst` = single-electrode
bursts) for set `1..4`, reproducing Tables 3.3–3.5 of the thesis. Set 2 is the
reference configuration; Set 1 deliberately keeps the historical Gaussian-kernel
size (6 instead of 6·σ+1) for the sensitivity comparison.

---

## 5. File inventory

### Core detection helpers
| File | Purpose |
|------|---------|
| `gaussian_kernel.m` | Build a 1-D Gaussian smoothing kernel. |
| `smooth_signal.m` | Convolve a signal with a kernel (same-length output). |
| `construct_ISIhist.m` | Smoothed, normalised logarithmic ISI-N histogram. |
| `calcISImax.m` | Find the ISI separating intra- from inter-burst ISIs (Pasquale et al., 2010). |
| `burstDetection.m` | Single-train burst detection, "string" method (Pasquale et al., 2010). |
| `BurstDetectISIn_MOD_AV.m` | ISI-N burst detection (Bakkum et al., 2014), with contributing-channel tracking. |
| `DetectBurstsNeuro.m` | Single-channel burst detection via the log-ISI histogram method. |
| `CutBurstsUniversal.m` | Remove bursts shorter than a minimum duration. |
| `MergeBurstsUniversal.m` | Merge consecutive bursts separated by a small gap. |

### Network-burst detectors
| File | Purpose |
|------|---------|
| `DetectNetworkBurstsNeuro.m` | Tampere **V1** NB detector (threshold = first local minimum after the main peak). |
| `DetectNetworkBurstsNeuroV2.m` | Tampere **V2** NB detector (threshold = global minimum of the valley between the first two prominent peaks). |
| `detectNetworkBursts_Doorn2024.m` | **Doorn** NB detector: locally-adaptive start/stop thresholds on the smoothed global spike-rate signal. |

### Pre-processing & processing
| File | Purpose |
|------|---------|
| `filter_electrodes_by_FR.m` | Three-level firing-rate filter (electrode ≥ 0.2 Hz; well mean FR; min active electrodes). |
| `run_NBdetection_tampere.m` | Run a Tampere detector (V1 or V2 via handle) over all plates/DIVs/wells of one line. |
| `run_NBdetection_doorn.m` | Run the Doorn detector over all plates/DIVs/wells of one line. |
| `get_parameter_set.m` | Return the detector parameters for parameter sets 1–4. |
| `run_preprocessing_processing.m` | **Master script 1** — load → FR filter → V1/V2/Doorn detection → save. |

### Metrics & statistics
| File | Purpose |
|------|---------|
| `compute_all_metrics_csv.m` | Per-well activity-metrics table (MFR, SBR, SBD, NBR, FragPerNB, FragPrev, NBD, CV-IBI). |
| `script3_statistics_LME_new_MFR_aligned_integrated.m` | Full LME statistics (Plate random effect), pairwise post-hoc with FDR, heatmaps and CI line plots. |
| `outputsection4_LME2.m` | Parameter-set (1–4) comparison boxplots + LME/Kruskal-Wallis workflow. |
| `stats_RNDSPIKES_new.m` | Pairwise statistics and heatmaps for the % of random (out-of-burst) spikes. |
| `run_analysis_metrics_stats.m` | **Master script 2** — metrics CSV → every figure/statistics function. |

### Figures
| File | Purpose |
|------|---------|
| `RASTER_NBandSEB.m` | Raster figures for single-electrode and network bursts (Doorn/Tampere). |
| `FR_raster_TvD_modified.m` | Per-well firing-rate + raster figures, Tampere vs Doorn side by side. |
| `ISIhist_threshold_V1vsV2_fix4.m` | Per-well log ISI-N histograms with V1 vs V2 thresholds. |
| `mfr_burstR_NBR.m` | MFR / single-burst-rate / NB-rate line plots across DIVs. |
| `mfr_burstR_NBR_CI_fixed.m` | Same with 95 % CI + IQR outlier removal. |
| `NBrate_boxplot_compare_TvsD.m` | NB-rate boxplots, Tampere vs Doorn per DIV. |
| `NBrate_boxplot_V1vsV2_2fig.m` | NB-rate and NB-count boxplots, Tampere V1 vs V2. |
| `NB_duration_CV_IBI_CI_fixed.m` | Mean NB-duration line plot and CV/LvR inter-burst-interval boxplots. |
| `lineplots_tutti_def_fixed_.m` | Unified NB-rate and fragments/NB line plots (Tampere + Doorn). |
| `boxplot_TvD_fragcomaprison.m` | Fragments-per-NB boxplots (absolute and ratio), Tampere vs Doorn. |
| `venn_fragments_fixed.m` | Venn diagrams of fragmented-NB overlap, Tampere vs Doorn. |
| `mean_SEB_duration_CI.m` | Mean single-electrode burst duration across DIVs (IQR-trimmed). |
| `mean_perc_rnd_spikes_fixed.m` | Mean % of random (out-of-burst) spikes across DIVs. |
| `plot_synchrony_lineplot_.m` | Mean active channels per NB across DIVs (95 % CI). |
| `radar_S_inner_pentagon.m` | Radar/spider chart of per-cell-line activity metrics. |

### Spectral (Cimbalník) validation
| File | Purpose |
|------|---------|
| `NB_Validation_Cimbalnik_v2.m` | HFO-based, spike-independent NB validation with 3-level event stratification (IN-FRAG / IN-NB-OUT / OUT-NB). |
| `NB_Validation_Plots_v4_representative.m` | Validation plots with representative-raster selection. |
| `decode_validation_mat.m` | Load the per-line validation `table` files, assign each event to one of the 3 classes (IN-FRAG / IN-NB-OUT / OUT-NB) and write the tidy `events_3class_Tampere.csv`; called automatically by `run_analysis_metrics_stats.m`. |
| `Validation_Tampere_<line>_v2.mat` | Per-event/per-well validation tables saved by the validation script. |

---

## 6. The Cimbalník spectral validation (summary)

For each electrode the raw extracellular signal is band-pass filtered in three
HFO bands — **High Gamma 65–100 Hz**, **Ripple 100–250 Hz**, **Fast Ripple
250–600 Hz** — with an 8th-order Butterworth IIR filter applied with zero-phase
`filtfilt`. The Hilbert envelope is z-scored against the whole recording per
electrode; the cross-electrode maximum z-score trace defines candidate spectral
events (`z > 3.0`, lasting ≥ 3 cycles at the band's lower edge). Each event is
assigned to exactly one of three mutually exclusive classes by where its peak
falls: inside a fragment (±25 ms window), inside an NB but outside any fragment,
or outside any NB. Events concentrating inside NBs/fragments — measured by peak
z-score, event rate, hit-rate and envelope amplitude — corroborate the
spike-based detection with a physically independent signal. Block-level
(cell-line × DIV) medians are compared with Friedman + Wilcoxon signed-rank
(Holm correction). High Gamma is by far the most informative band.

---

## 7. Notes & known inconsistencies

See `CLEANING_NOTES.md` for the full publication-cleanup log. Highlights:

- The three original `NBdet_*` per-line scripts were collapsed into the two
  reusable runners (`run_NBdetection_tampere.m`, `run_NBdetection_doorn.m`).
- All hard-coded Windows output paths were removed; outputs derive from
  `RESULTS_ROOT`.
- **V1 vs V2 field-name mismatch:** the two Tampere detectors store the same
  fragment/activity fields under slightly different names
  (`is_fragmented`/`isfragmented`, etc.). Both were left as-is so existing
  downstream scripts keep working; unify them before further reuse.
- `script3` intentionally displays DIV 62 as "60" (`DIVlabels`).
- Verbose `fprintf` debug output was removed from the detector functions;
  numerical results are unchanged.

---

## 8. References

- Pasquale, V., Martinoia, S., & Chiappalone, M. (2010). A self-adapting
  approach for the detection of bursts and network bursts in neuronal cultures.
  *Journal of Computational Neuroscience, 29*(1–2), 213–229.
  https://doi.org/10.1007/s10827-009-0175-1
- Bakkum, D. J., et al. (2014). Parameters for burst detection.
  *Frontiers in Computational Neuroscience, 7*, 193.
  https://doi.org/10.3389/fncom.2013.00193
- van der Molen, T., et al. / Doorn et al. (2024). Network-burst detection
  with locally adaptive thresholds (as implemented in
  `detectNetworkBursts_Doorn2024.m`).
- Cimbalník, J., et al. (2018). High-frequency-oscillation framework for the
  spectral characterisation of intracranial signals (adapted here to in-vitro
  MEA recordings).

---

## 9. License

Add a license of your choice (e.g. MIT) before publishing.
