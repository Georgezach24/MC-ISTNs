# MC-ISTNs — Multi-Connectivity Integrated Satellite-Terrestrial Networks

MATLAB simulation of a joint terrestrial (5G NR) and non-terrestrial (LEO satellite) network, developed as Part 1 of a thesis on ML-aided multi-connectivity management for 6G. It models a set of users, terrestrial base stations, and a LEO satellite, computes the link budget from each user to each candidate node, and selects the best-serving node per user based on SNR.

This simulation is the reference generator for a dataset (Part 2, early proof-of-concept underway — see `Model/`) intended to train a machine-learning model that predicts/decides the best connectivity option per user, in place of the SNR-comparison rule used here.

## What it does

For a fixed set of geographic positions (base stations, users, one satellite):

1. Computes terrestrial path loss per user–BS pair using the 3GPP TR 38.901 channel model (`UMa`/`UMi` scenario, selectable), via MATLAB's 5G Toolbox — including a per-link LOS/NLOS draw from the TR 38.901 §7.4.2 distance-dependent LOS probability, and log-normal shadow fading from the model's own shadow-fading standard deviation (§7.4.1).
2. Computes the satellite path loss per user using free-space path loss, gated by a minimum elevation visibility mask.
3. Derives SNR for every candidate link (BS and satellite) from fixed transmit power / EIRP and thermal noise.
4. Selects each user's serving node as the candidate with the **highest SNR** (single-connectivity, per-user greedy selection).
5. Splits each node's bandwidth equally among the users it serves, and computes per-user Shannon capacity.
6. Computes an energy-per-bit proxy per user: node power draw (EARTH linear power model for the BS, a linear power-amplifier-efficiency model for the satellite) split equally across served users and divided by each user's bit rate.
7. Prints a results table and renders a 3D plot of the network topology and serving links.

Beyond the single static run above, two further modes reuse the same per-user link budget: repeated runs of the same topology to characterize KPI variance from channel stochasticity (`kpiRepeatedRuns.m`), and a time-stepped simulation of a satellite pass with spatially-correlated shadow fading across steps (`temporalPassSimulation.m`) — see below.

## Requirements

- MATLAB with:
  - **5G Toolbox** (`nrCarrierConfig`, `nrPathLossConfig`, `nrPathLoss`, `nrTDLChannel`, `nrCDLChannel`)
  - **Mapping Toolbox** or **Aerospace Toolbox** (`wgs84Ellipsoid`, `geodetic2enu`, `geodetic2aer`, `distance`)
  - **Phased Array System Toolbox** / **Communications Toolbox** (`fspl`, `physconst`)

## Running

Open MATLAB with `PROD/` on the path and run the entry script:

```matlab
run('PROD/test_simulation.m')
```

or from a shell with MATLAB on `PATH`:

```
matlab -batch "run('PROD/test_simulation.m')"
```

The script prints a per-user results table (serving node, distance, path loss, SNR, capacity, satellite elevation) and opens a 3D visualization figure.

### Generating a dataset (Monte-Carlo batch driver)

To generate a labeled dataset by running the same per-user link budget over many randomized topologies (BS/user counts and positions, satellite geometry, UMa/UMi mix), use `monteCarloDriver.m`:

```matlab
monteCarloDriver()          % 200 randomized scenarios -> Dataset/dataset.csv
monteCarloDriver(1000)      % 1000 scenarios -> Dataset/dataset.csv
```

Radio configuration (transmit power, bandwidth, power models) stays fixed — identical to `test_simulation.m` — only topology and satellite geometry are randomized, so results stay comparable to the single-run reference scenario. Output is one CSV row per user per scenario, with columns for scenario metadata (`ScenarioID`, `NumBS`, `NumUsers`, `ScenarioType`), geometry, the same per-user metrics as the single-run table, `NodeLoad` (how many users share the assigned node), and per-candidate diagnostics (`CandBS_*`/`CandSat_*`: best-BS and satellite SNR/distance/elevation/path-loss, independent of which one wins) so a model can learn the comparison rather than just read off the winner. `Dataset/` is gitignored since it's generated output, not source.

### Testing a model (Part 2, early proof-of-concept)

`Model/train_model.py` trains classifiers on the generated dataset to predict `ServingType` (Terrestrial vs Satellite) from the candidate-level features. See `Model/README.md` for setup, current metrics, and honest caveats about what this first pass does and doesn't prove.

### Characterizing KPI variance (repeated runs)

`kpiRepeatedRuns.m` runs the *same* static reference topology repeatedly with a different RNG seed each time, to separate KPI variance caused by channel stochasticity (LOS draw + shadow fading) from variance caused by topology:

```matlab
kpiRepeatedRuns()         % 500 runs -> Results/
kpiRepeatedRuns(1000)     % 1000 runs -> Results/
```

Writes per-run/per-user results and a mean/std/min/max summary grouped by `ServingType` to `Results/`, plus overlaid histograms per KPI.

### Simulating a satellite pass over time

`temporalPassSimulation.m` runs the same static BS/user topology across successive time steps while moving the satellite along a simplified ground track, so elevation/SNR/capacity evolve realistically over a LEO pass:

```matlab
temporalPassSimulation()      % dt = 5s -> Results/
temporalPassSimulation(2)     % dt = 2s
```

Shadow fading/LOS state carries spatial correlation across time steps (Gudmundson 1991 autocorrelation, TR 38.901 Table 7.5-6 correlation distances) instead of resampling independently at every step — this was added because independent-per-step sampling produced physically implausible handover "ping-pong" for otherwise-stationary users. Writes a per-step CSV, prints per-user handover counts, and writes time-series plots to `Results/`.

## Project structure

```
PROD/
  test_simulation.m       Entry point — scenario definition, calls simulateScenario, then reporting
  simulateScenario.m      Per-user link budget: LOS draw, shadow fading, SNR, node selection, capacity
  monteCarloDriver.m      Batch driver — runs simulateScenario over randomized topologies, writes a labeled CSV
  kpiRepeatedRuns.m       Repeats the static reference topology many times to characterize KPI variance
  temporalPassSimulation.m  Runs the static topology across time as the satellite moves through a pass
  array.m                 Formats and prints the per-user results table
  visual.m                3D plot of base stations, users, satellite, and serving links
  istn.zip                Archived snapshot of an earlier version
Βοηθητικά Έργαλεία/
  graph.m              Standalone plot of measured RX power vs. distance (not part of the simulation pipeline)
Plots/
  ...                  Saved figure exports from test_simulation.m/visual.m runs
Results/
  ...                  CSV/PNG output of kpiRepeatedRuns.m and temporalPassSimulation.m (tracked in git — feeds thesis figures)
Dataset/
  ...                  CSV output from monteCarloDriver.m (generated, gitignored)
Model/
  train_model.py       Part 2 proof-of-concept: trains/evaluates classifiers on Dataset/dataset.csv
  requirements.txt     Python dependencies (pandas, scikit-learn, numpy, matplotlib)
  results/             Saved metrics.json and plots from the last training run
Thesis/
  latex/               LaTeX thesis source (chapters, figures, main.tex/main.pdf)
```

## Current limitations / scope

- No user mobility — users are stationary within a run; only the satellite moves (in `temporalPassSimulation.m`), via a simplified constant-latitude ground track, not a full orbit propagator (e.g. SGP4).
- Best-node selection is a greedy, per-user SNR comparison with **no minimum-SNR/outage threshold** — a user with no usable candidate is still assigned to the least-bad one instead of being marked out-of-coverage (reproduced concretely in the temporal pass simulation). Not a joint network-wide optimization and not true dual/multi-connectivity (each user attaches to exactly one node), and node re-selection has no handover hysteresis/margin.
- No inter-cell interference — SNR is noise-limited only.
- The Monte-Carlo driver randomizes topology and satellite geometry per scenario, but not radio configuration (power, bandwidth) — each generated row is still a static single-snapshot link budget (no user mobility within a scenario).
- Energy-per-bit is a modeled proxy (EARTH power model / PA-efficiency model, not a hardware measurement), and only accounts for nodes actively serving at least one user — idle-node power is not yet tracked.
- Satellite channel is a deterministic free-space-path-loss model given fixed geometry — no atmospheric/scintillation loss or antenna-gain pattern yet, unlike the stochastic, TR 38.901-grounded terrestrial channel. This asymmetry shows up empirically as satellite KPIs having far lower run-to-run variance than terrestrial KPIs (`kpiRepeatedRuns.m`).

## Roadmap

In priority order, per the thesis discussion chapter:

1. Add an explicit outage state (minimum acceptable SNR/SINR threshold in node selection).
2. Add handover hysteresis/margin (e.g. akin to TS 38.331 Event A3), complementing the temporal shadow-fading correlation already implemented.
3. Replace the simplified ground-track model with a real orbit propagator (e.g. SGP4 on real two-line elements).
4. Add atmospheric/scintillation loss to the satellite link (TR 38.821 §6.1) and antenna gain modeling on the terrestrial side (TR 38.901 §7.3), so terrestrial and satellite links are on equal footing.
5. Real multi-connectivity — extend node selection beyond single-connectivity greedy argmax to simultaneous BS+satellite service.
6. Cap capacity to standard MCS/CQI spectral efficiency (TS 38.214) instead of unbounded Shannon capacity.
7. Part 2: a non-trivial ML target — the current `ServingType` classifier is a near-deterministic function of its own input features (see `Model/README.md`); next targets include prediction from noisy/partial observations, pre-handover prediction, or regression/ranking/fairness-aware formulations.
