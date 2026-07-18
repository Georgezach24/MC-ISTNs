# MC-ISTNs — Multi-Connectivity Integrated Satellite-Terrestrial Networks

MATLAB simulation of a joint terrestrial (5G NR) and non-terrestrial (LEO satellite) network, developed as Part 1 of a thesis on ML-aided multi-connectivity management for 6G. It models a set of users, terrestrial base stations, and a LEO satellite, computes the link budget from each user to each candidate node, and selects the best-serving node per user based on SNR.

This simulation is the reference generator for a future dataset (Part 2) intended to train a machine-learning model that predicts/decides the best connectivity option per user, in place of the SNR-comparison rule used here.

## What it does

For a fixed set of geographic positions (base stations, users, one satellite):

1. Computes terrestrial path loss per user–BS pair using the 3GPP TR 38.901 channel model (`UMa`/`UMi` scenario, selectable), via MATLAB's 5G Toolbox — including a per-link LOS/NLOS draw from the TR 38.901 §7.4.2 distance-dependent LOS probability, and log-normal shadow fading from the model's own shadow-fading standard deviation (§7.4.1).
2. Computes the satellite path loss per user using free-space path loss, gated by a minimum elevation visibility mask.
3. Derives SNR for every candidate link (BS and satellite) from fixed transmit power / EIRP and thermal noise.
4. Selects each user's serving node as the candidate with the **highest SNR** (single-connectivity, per-user greedy selection).
5. Splits each node's bandwidth equally among the users it serves, and computes per-user Shannon capacity.
6. Computes an energy-per-bit proxy per user: node power draw (EARTH linear power model for the BS, a linear power-amplifier-efficiency model for the satellite) split equally across served users and divided by each user's bit rate.
7. Prints a results table and renders a 3D plot of the network topology and serving links.

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

## Project structure

```
PROD/
  test_simulation.m   Entry point — scenario definition, calls simulateScenario, then reporting
  simulateScenario.m   Per-user link budget: LOS draw, shadow fading, SNR, node selection, capacity
  array.m              Formats and prints the per-user results table
  visual.m             3D plot of base stations, users, satellite, and serving links
  istn.zip             Archived snapshot of an earlier version
Βοηθητικά Έργαλεία/
  graph.m              Standalone plot of measured RX power vs. distance (not part of the simulation pipeline)
Plots/
  ...                  Saved figure exports from previous simulation runs
```

## Current limitations / scope

- Single static topology per run — no user mobility or satellite motion (elevation is fixed, not propagated over time).
- Best-node selection is a greedy, per-user SNR comparison, not a joint network-wide optimization and not true dual/multi-connectivity (each user attaches to exactly one node).
- No inter-cell interference — SNR is noise-limited only.
- No batch/Monte-Carlo driver yet — this repo does not currently generate a dataset, only a single labeled scenario per run.
- Energy-per-bit is a modeled proxy (EARTH power model / PA-efficiency model, not a hardware measurement), and only accounts for nodes actively serving at least one user — idle-node power is not yet tracked.

## Roadmap

- Add a Monte-Carlo batch driver on top of `simulateScenario.m` to generate a labeled dataset (varied topologies, load, and satellite geometry) as CSV/Parquet.
- Add atmospheric/scintillation loss to the satellite link (TR 38.821 §6.1) and antenna gain modeling on the terrestrial side (TR 38.901 §7.3), so terrestrial and satellite links are on equal footing.
- Cap capacity to standard MCS/CQI spectral efficiency (TS 38.214) instead of unbounded Shannon capacity.
- Add satellite pass dynamics (time-varying elevation).
- Part 2 (separate, future work): train an ML model on the generated dataset to predict the best connectivity option per user.
