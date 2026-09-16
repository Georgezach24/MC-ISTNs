# MC-ISTNs — Multi-Connectivity Integrated Satellite-Terrestrial Networks

MATLAB simulation of a joint terrestrial (5G NR) and non-terrestrial (LEO satellite) network, developed as Part 1 of a thesis on ML-aided connectivity management for 6G. It models a set of users, terrestrial base stations, and a LEO satellite; computes the full link budget from each user to every candidate node; selects the serving node; and classifies the resulting service state against a standards-derived requirement.

The simulation is the reference generator for a dataset (Part 2, proof-of-concept — see `Model/`) intended to train a machine-learning model that predicts the best connectivity option per user, in place of the SNR-comparison rule used here.

> **Scope note.** This branch implements **single connectivity**: each user is served by exactly one node, or by none (outage). Simultaneous BS+satellite service, handover hysteresis, and fairness-aware load balancing were explored on the abandoned `ml_v1` branch and are **not** part of this code.

## What it does

For a set of geographic positions (base stations, users, one satellite):

1. **Terrestrial link budget** per user–BS pair, using the 3GPP TR 38.901 `UMa`/`UMi` channel model via MATLAB's 5G Toolbox:
   - geodesic 2D distance and *physical* antenna heights (never ENU vertical coordinates, which drift with Earth curvature);
   - an applicability gate (TR 38.901 Table 7.4.1-1: 10 m ≤ d2D ≤ 5 km, 1.5 m ≤ hUT ≤ 22.5 m) — links outside the model's valid range are not computed at all, and are reported as such;
   - a spatially-consistent LOS/NLOS state (§7.4.2 probability, §7.6.3.3 consistency);
   - log-normal shadow fading with the model's own σ (§7.4.1), exponentially correlated across time steps (§7.4.4);
   - small-scale fading: Rician (K = 9 dB, LOS) or Rayleigh (NLOS), flat.
2. **Satellite link budget**: free-space path loss + gaseous attenuation on the slant path (TR 38.811 §6.6.4 / ITU-R P.676), Shadowed-Rician fading (Abdi et al., 2003) parameterized by elevation, gated by a minimum-elevation visibility mask.
3. **SNR** per candidate link from EIRP and thermal noise. Terrestrial EIRP = conducted power + element gain + array gain (TR 38.901 Table 7.8-1 / §7.3); satellite EIRP is derived from the **EIRP density** of TR 38.821 Table 6.1.1.1-1 (LEO-600 S-band) and the channel bandwidth, with RF power kept distinct from amplifier electrical draw.
4. **Node selection**: the candidate with the highest SNR wins. If even the best candidate is below the minimum usable SNR (≈ −7.53 dB, the Shannon-equivalent of MCS 0, TS 38.214 Table 5.1.3.1-2), the user is declared **out of coverage** instead of being attached to an unusable node.
5. **Resource allocation**: the serving node's bandwidth is split equally among its users; per-user capacity is Shannon, capped at the NR maximum spectral efficiency (5.5547 bit/s/Hz, MCS 27).
6. **Service classification**: `Served` / `BelowTarget` / `Outage`, against the 5th-percentile user spectral-efficiency requirement for Dense Urban-eMBB (0.3 bit/s/Hz, TR 37.910 Table 5.4.1.1.1-1). Capacity and *delivered throughput* are reported separately, and the reason a candidate was unusable (below SNR floor / not visible / outside model range) is recorded per link.
7. **Energy**: per-user energy-per-bit (EARTH linear power model for the BS, Auer et al. 2011; linear PA-efficiency model for the satellite) plus a network-level power / bit-per-joule tally, including idle-BS power and a symmetric amplifier-only figure for fair BS-vs-satellite comparison.
8. **Reporting**: a per-user results table, a 3D plot of the topology and serving links, and a versioned run folder with the exact parameters used.

## Requirements

- MATLAB **R2024b** (recorded automatically in each run's `params.txt`) with:
  - **5G Toolbox** (`nrCarrierConfig`, `nrPathLossConfig`, `nrPathLoss`)
  - **Mapping** or **Aerospace Toolbox** (`wgs84Ellipsoid`, `geodetic2enu`, `geodetic2aer`, `distance`, `ecef2lla`)
  - **Phased Array System** / **Communications Toolbox** (`fspl`, `physconst`, `gaspl`)
- Python 3 with `Model/requirements.txt` for Part 2.

## Running

```matlab
addpath('PROD'); runSimulation()
```

or from a shell with MATLAB on `PATH`:

```
matlab -batch "addpath('PROD'); runSimulation()"
```

There is **one** simulation and it runs **once**. It produces every number in the thesis and the ML training set from the same execution. Options are passed as a struct:

```matlab
runSimulation(struct('maxPasses', 10))          % quick check
runSimulation(struct('ciTolerance', 0.02))      % tighter convergence
runSimulation(struct('label', 'v2'))            % tag the versioned run folder
```

### How the simulated time is structured

The run is a loop over **satellite passes**; inside each pass a 1 s time loop advances both the satellite and the users. Every timing choice is taken from a standard rather than picked:

| Quantity | Value | Source |
|---|---|---|
| Elementary window | one pass, 900 s | TR 38.821 Table 4.2-3 NOTE 1 — "a period of time corresponding to the visibility time of the satellite" |
| Time step | 1 s | ITU-R M.2412-0 Annex 1 §5.3.2 (UE displacement below 1 m per step) and TR 38.821 §7.3.2.1.4 Table 7.3.2.1.4-1 (fastest LEO mobility timescale 6.61 s) |
| User speed / direction | 3 km/h, fixed per pass, uniformly random azimuth | ITU-R M.2412-0 §8.4, TABLE 5 b)/c) — "fixed and identical speed of all UEs of the same mobility class, randomly and uniformly distributed direction" |
| Total duration | until the KPIs converge | ITU-R M.2412-0 §7.1 — "a sufficient number of drops … to ensure convergence" plus "the width of confidence intervals"; extended to satellite evaluation by ITU-R M.2514-0 §8.2.4 |

No standard fixes a total duration — TR 38.821's own NTN system-level calibration table (6.1.1.1-5) has no duration field at all. The stopping rule is therefore a convergence criterion: the loop ends when the 95% confidence-interval half-width of the tracked KPIs, computed **across passes**, falls below 5% (relative, for rate and energy metrics) or 2 percentage points (absolute, for fractions). The convergence trace is saved as a result in its own right (`Results/convergence.csv` and `convergence.png`), which is what makes the sample size defensible instead of arbitrary.

Each pass is an independent *drop* in the ITU-R sense: users are re-dropped at new positions and the channel state is reset between passes, while being threaded step-to-step **within** a pass. `PassID` is therefore the unit for the ML train/test split, for confidence intervals, and for any per-repeat aggregation — users inside one pass are not independent observations.

Pass geometry varies: the RAAN is offset per pass so transits range from grazing (peak elevation at the 20° mask) to near-zenith. The largest useful offset is found numerically at startup rather than assumed.

### Moving users, and what is checked every step

Users walk continuously during a pass. That is also what activates the correlated shadow fading and the spatially-consistent LOS state — with static users the per-step displacement was zero and both models were inert.

Because the users move, the run asserts on **every step of every pass** that their state is what it should be: antenna height still exactly 1.5 m, coordinates finite and in range, and every user still holding at least one terrestrial link inside the TR 38.901 Table 7.4.1-1 validity box (10 m ≤ d2D ≤ 5 km). It also checks the outputs for consistency — outage rows carry no capacity, zero throughput and zero load; served rows carry finite positive capacity, finite SNR and load ≥ 1. Any violation aborts the run naming the pass, step and user. Initial radii are sampled so that no walk can leave the validity box, so these assertions are a check rather than a correction.

### Outputs

| File | Content |
|---|---|
| `Dataset/dataset.csv` | one row per pass/step/user — geometry, per-candidate diagnostics, serving decision, capacity, delivered throughput, service state, link state, energy. Feeds both the results and the ML side. Gitignored. |
| `Results/convergence.csv` | per-pass running means and CI half-widths of the tracked KPIs |
| `Results/convergence.png` | the convergence curve against the stopping threshold |
| `Results/temporal_*.png` | time series over the reference pass |
| `Results/network_3d.png` | topology and serving links at closest approach |
| `Results/runs/<timestamp>_runSimulation[_label]/` | versioned copy of all of the above plus the exact parameters |

### Verification scripts

```matlab
geometryValidation()                  % antenna-height & distance regression check
energyModelValidation()               % user-count dependence of the energy metrics
```

Both print explicit pass/fail lines and write CSV/PNG. `geometryValidation` confirms the user height reaching the path-loss model never drifts with distance, that 3D range is consistent with 2D distance and heights, and that the 5 km validity gate fires — and quantifies the error the earlier ENU-based geometry introduced. `energyModelValidation` shows that per-user capacity scales as 1/L while energy-per-bit stays flat (the user count cancels algebraically), and that the network bit/J metric — unlike energy-per-bit — does respond to user composition.

### Result versioning

Every run ends by calling `saveRunVersion`, producing:

```
Results/runs/<YYYYMMDD_HHMMSS>_<script>[_<label>]/
  params.mat      exact parameter restore
  params.txt      readable snapshot, incl. git commit and MATLAB version
  changed.txt     parameter diff vs. the previous run of the same script
  <CSV/PNG output of that run>
```

This is what makes a figure or a table traceable back to the code and parameters that produced it.

### Part 2 (ML)

```
pip install -r Model/requirements.txt
python Model/train_model.py
```

Three passes on the generated dataset — exact SNR (pipeline sanity check), geometry-only, and noisy SNR — predicting the serving type from candidate-level features, split **by pass** (`PassID`) so no pass contributes to both train and test. Because the dataset is a 1 s time series, the training scripts decimate it (one sample every 10 s) so consecutive rows are not near-duplicates. See `Model/README.md`.

## Project structure

```
PROD/
  simulateScenario.m        Core: link budget, fading, node selection, allocation, energy, service state
  runSimulation.m           The entry point — one run: passes, moving users, convergence, dataset
  geometryValidation.m      Regression check on link geometry / antenna heights
  energyModelValidation.m   Check of how the energy metrics depend on user count
  saveRunVersion.m          Versioned run folders (parameters, commit, diff, outputs)
  array.m                   Prints the per-user results table
  visual.m                  3D plot of base stations, users, satellite, serving links
Results/runs/               Versioned outputs of every run (tracked in git)
Dataset/                    CSV output of runSimulation.m (generated, gitignored)
Model/                      Part 2: Python training scripts, metrics and plots
```

## Current limitations / scope

- **Single connectivity.** One serving node per user, chosen by max SNR with an outage floor — not a joint network-wide optimization, and not simultaneous multi-connectivity.
- **No interference.** SNR only; inter-BS co-channel interference is not modeled, and no frequency-reuse scheme is stated. Satellite and terrestrial segments are separated in frequency (2.0 vs 3.5 GHz), so the gap is terrestrial-side.
- **Constant satellite antenna gain.** 30 dBi regardless of off-axis angle; no beam-pointing policy or radiation pattern. UE antenna gain is assumed 0 dBi.
- **Satellite atmosphere partially modeled.** Gaseous attenuation is included; rain/cloud attenuation and ionospheric scintillation are not (low impact at S-band, but not quantified here).
- **Energy is a modeled proxy**, not a measurement, and the two segments have different subsystem scopes: the satellite model's fixed power term is set to zero, which is an *optimistic* assumption for the satellite side. Use the amplifier-only network figure for like-for-like comparison.
- **Energy per bit is invariant to node load** by construction (the user count cancels); network bit/J is the metric that responds to load and user composition.
- **No user mobility.** Users are stationary; only the satellite moves.
- **Handovers are state transitions with an interruption cost**, not a full 3GPP RRC handover procedure with signaling, failure, and re-establishment.

## Roadmap

Following the external assessment of the thesis, in rough priority order:

1. Regenerate the dataset from the current model, realign the Part 2 scripts, and re-run all ML results.
2. Paired, same-realization comparison of explicit policies (terrestrial-only / satellite-only / actual) with proper confidence intervals, across multiple topologies.
3. SINR: inter-BS interference under a stated frequency-reuse assumption — or an explicit justification of orthogonality and its spectral-efficiency cost.
4. Satellite antenna pattern with a per-user off-axis angle and a stated beam-pointing policy.
5. Extended verification: worked numeric link-budget examples, known-behavior checks, and pre-specified edge-case handling.
6. Part 2: a non-trivial ML target — prediction before the decision point (pre-handover), regression/ranking on capacity or energy, evaluated by achieved service rather than classification accuracy.
