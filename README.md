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
run('PROD/test_simulation.m')      % single reference scenario
```

or from a shell with MATLAB on `PATH`:

```
matlab -batch "run('PROD/test_simulation.m')"
```

Prints the per-user table (serving node, distance, path loss, SNR, capacity, satellite elevation, node power, energy per bit), opens the 3D figure, and saves a versioned run folder under `Results/runs/`.

### Generating a dataset (Monte-Carlo driver)

```matlab
monteCarloDriver()                    % 200 randomized scenarios -> Dataset/dataset.csv
monteCarloDriver(1000)                % 1000 scenarios
monteCarloDriver(500, [], 'tag')      % tag the versioned run folder
```

Randomizes BS/user counts and positions, the sub-satellite point (±10° lat/lon), and the UMa/UMi scenario; keeps the radio configuration fixed and identical to `test_simulation.m`. Users are placed uniformly in an annulus around a reference BS, inside the terrestrial models' validity range. One CSV row per user per scenario, with scenario metadata, geometry, per-user metrics, node load, service state, unavailability reasons, and the `CandBS_*`/`CandSat_*` per-candidate diagnostics (so a model can learn the comparison instead of reading off the winner). `Dataset/` is gitignored — it's generated output.

### KPI variance over repeated runs

```matlab
kpiRepeatedRuns()                     % 500 repeats of the same topology -> Results/
```

Re-runs the *same* static topology with a different seed each time, isolating variance caused by channel stochasticity from variance caused by topology. Writes per-run rows, a summary grouped by serving type, and overlaid histograms.

> Caveat: grouping by serving type compares different users and geometries, not the same users under different policies — it does not by itself isolate a terrestrial-vs-satellite effect. A paired, per-policy comparison is still to be implemented.

### Satellite pass over time

```matlab
temporalPassSimulation()              % dt = 5 s -> Results/
temporalPassSimulation(1)             % dt = 1 s
```

Steps the same static BS/user topology through time while the satellite moves along a **circular Keplerian orbit**: the state is propagated in an inertial frame and converted to geodetic coordinates through an explicit Earth-rotation step, so orbital speed and ground-track geometry come from one consistent model. Channel state (LOS + shadow fading) is threaded step to step. Handovers are recorded as explicit link-state transitions with an interruption cost of 2·RTT (TR 38.821 §7.3.2.1.1) deducted from delivered throughput. The run stops when the satellite leaves visibility after a pass. Writes a per-step CSV, per-KPI time-series plots, and per-user handover / outage-event / lost-time counts.

### Verification scripts

```matlab
geometryValidation()                  % antenna-height & distance regression check
energyModelValidation()               % user-count dependence of the energy metrics
```

Both print explicit pass/fail lines and write CSV/PNG. `geometryValidation` confirms the user height reaching the path-loss model never drifts with distance, that 3D range is consistent with 2D distance and heights, and that the 5 km validity gate fires — and quantifies the error the earlier ENU-based geometry introduced. `energyModelValidation` shows that per-user capacity scales as 1/L while energy-per-bit stays flat (the user count cancels algebraically), and that the network bit/J metric — unlike energy-per-bit — does respond to user composition.

### Result versioning

Every script ends by calling `saveRunVersion`, producing:

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

Three passes on the generated dataset — exact SNR (pipeline sanity check), geometry-only, and noisy SNR — predicting the serving type from candidate-level features, split by scenario. See `Model/README.md`.

> **Currently out of sync:** the committed dataset predates the physical-model corrections listed below, and the scripts hardcode an older elevation mask. Regenerate the dataset and realign the scripts before quoting any ML number.

## Project structure

```
PROD/
  simulateScenario.m        Core: link budget, fading, node selection, allocation, energy, service state
  test_simulation.m         Entry point — single reference scenario
  monteCarloDriver.m        Batch driver over randomized topologies -> labeled CSV
  kpiRepeatedRuns.m         Repeats one topology to characterize channel-driven KPI variance
  temporalPassSimulation.m  Time-stepped LEO pass with handover cost accounting
  geometryValidation.m      Regression check on link geometry / antenna heights
  energyModelValidation.m   Check of how the energy metrics depend on user count
  saveRunVersion.m          Versioned run folders (parameters, commit, diff, outputs)
  array.m                   Prints the per-user results table
  visual.m                  3D plot of base stations, users, satellite, serving links
  istn.zip                  Archived snapshot of an earlier version
Βοηθητικά Έργαλεία/
  graph.m                   Standalone plot of measured RX power vs. distance (not part of the pipeline)
Results/runs/               Versioned outputs of every script (tracked in git)
Dataset/                    CSV output of monteCarloDriver.m (generated, gitignored)
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
