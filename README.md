# MC-ISTNs — Multi-Connectivity Integrated Satellite-Terrestrial Networks

MATLAB simulation of a joint terrestrial (5G NR) and non-terrestrial (LEO satellite) network, developed as Part 1 of a thesis on ML-aided connectivity management for 6G. It models a set of users, terrestrial base stations, and LEO satellites; computes the full link budget from each user to every candidate node; selects the serving node; and classifies the resulting service state against a standards-derived requirement.

The simulation is the reference generator for a dataset (Part 2, proof-of-concept — see `Model/`) intended to train a machine-learning model that predicts the best connectivity option per user, in place of the rule used here.

> **Scope note.** This branch implements **single connectivity**: each user is served by exactly one node, or by none (outage). Simultaneous BS+satellite service and the fairness-aware load balancing of the abandoned `ml_v1` branch are **not** part of this code. Handover hysteresis *is* implemented here (since 2026-09-26), as a transposition of the 3GPP A3 criterion, and load balancing is implemented as the standard's own Mobility Load Balancing (a load-dependent A3 offset, since 2026-10-24), in two extra policies alongside the plain one.

## What it does

For a set of geographic positions (base stations, users, satellites):

**Scenario.** 19 base stations (centre + 2 rings) on a hexagonal grid with 500 m inter-site distance and 25 m height — the UMa layout of TR 38.901 Table 7.2-1 ("19 macro sites, ISD = 500 m"), **with wrap-around** as in the 3GPP system-level calibration: the 19-site cluster is repeated six times around itself (shifts of ISD·√19, axial (5,−2), (2,3), (−3,5) and their negatives), and every link — path loss, LOS state, shadowing, interference — is computed from the nearest copy of the station. Every site, including the outer ring, therefore sees the same surroundings. **Each site has 3 sectors** (TR 38.901 Table 7.8-1, "3 sectors per cell site: 30, 150 and 270 degrees"; ITU-R M.2412-0 §8.3, "19 sites, each of 3 TRxPs (cells)"), so the network has 57 cells; cell *c* = 3(site−1) + sector is node `BSc`. With neighbours at 0°, 60°, … every sector points between two neighbouring sites. Propagation (distance, LOS, shadowing, path loss) depends only on positions and is computed once per site; the three cells of a site differ only in antenna gain. `bsSectors = 1` gives the earlier omnidirectional model (one cell per site, constant 18 dBi), which reproduces earlier results byte for byte. At start-up the run checks that the copies tile the plane exactly (every site of the next ring is a copy of a station, no two copies overlap) and that a point at the same offset from any station sees the same set of distances; `wrapAround = false` gives the plain cluster. 10–30 users per pass start uniformly inside the 19 hexagonal cells (Table 7.2-1, "UT distribution: Uniform"), each with a uniformly random heading. The 35 m minimum BS–UT distance of the same table is enforced on the whole walk: a position/heading pair whose straight 750 m path would come closer than 35 m to any station is redrawn. A user starting near the edge can walk up to 750 m past the cluster; with wrap-around it is then simply in a copy of the network, and every link stays inside the model's 5 km validity range. (Earlier versions put 7 stations in the middle and all users 0.76–4.25 km away from them, which left the stations bunched in the centre of the user area.) 3 consecutive LEO-600 satellites in one orbital plane; their spacing (317 s, 19.7°, about 18 per plane) is chosen so that the three are exactly the satellites of the plane that can be seen inside the 900 s window — the next one rises just as the window closes, even on a zenith pass — so no satellite that would really be visible is missing. Consecutive satellites overlap on high passes (satellite↔satellite handover) and leave a gap on low ones. All of these are options (`bsRings`, `bsIsdM`, `numUsersRange`, `numSats`, `satFreqReuse`).

1. **Terrestrial link budget** per user–BS pair, using the 3GPP TR 38.901 `UMa`/`UMi` channel model via MATLAB's 5G Toolbox:
   - geodesic 2D distance and *physical* antenna heights (never ENU vertical coordinates, which drift with Earth curvature);
   - an applicability gate (TR 38.901 Table 7.4.1-1: 10 m ≤ d2D ≤ 5 km, 1.5 m ≤ hUT ≤ 22.5 m) — links outside the model's valid range are not computed at all, and are reported as such;
   - a spatially-consistent LOS/NLOS state (§7.4.2 probability, §7.6.3.3 consistency), used as the **soft LOS state** of §7.6.3.3 (eq. 7.6-18/19, option `softLosState`, default on): the TR itself notes that with mobility a LOS→NLOS change gives "a hard transition in the channel response" and defines `LOS_soft = 1/2 + (1/π)·arctan(√(d_LOS/λ)·(G + F(d)))`, `F(d) = √2·erf⁻¹(2P_LOS(d) − 1)`, with the channel computed once as LOS and once as NLOS and combined as `√LOS_soft·H_LOS + √(1−LOS_soft)·H_NLOS` (here in power, with independent phases). `d_LOS` is not defined in the clause; it is taken as the LOS-state correlation distance of Table 7.6.3.1-2 (50 m, UMa/UMi). With the hard state the latent variable crosses the threshold many times around each transition, so path-loss jumps of more than 10 dB grew from 2.8 to 7.5 per user-minute when the step went from 1 s to 0.2 s; with the soft state they are 0.18;
   - log-normal shadow fading with the model's own σ (§7.4.1), exponentially correlated across time steps (§7.4.4) — with soft LOS one process per state (LOS 37 m, NLOS 50 m);
   - small-scale fading: Rician (K = 9 dB, LOS) or Rayleigh (NLOS), flat.
2. **Satellite link budget**: free-space path loss + gaseous attenuation on the slant path (TR 38.811 §6.6.4 / ITU-R P.676), Shadowed-Rician fading (Abdi et al., 2003) parameterized by elevation, gated by a minimum-elevation visibility mask. The shadowing term (the Nakagami-m LOS amplitude) is **correlated along the user's path** with an exponential autocorrelation and the 37 m horizontal correlation distance of TR 38.811 Table 6.7.2-3a (Urban, S-band, LOS) — the same form as the terrestrial shadow fading. One latent Gaussian per user–satellite pair (each satellite is seen through different obstacles; TR 38.811 gives no cross-correlation) evolves as AR(1) with the user's displacement and is mapped to the model's Gamma distribution through its CDF, so the marginal distribution is exactly Abdi's and only the time correlation is new; as elevation changes the parameters change, but the shadowing percentile carries on smoothly. The scatter component stays independent per step. Decorrelation from the satellite's own motion is not modelled: the standard gives only the horizontal distance, and notes that real distances may be shorter.
3. **SINR** per candidate link. Terrestrial received power = conducted power (49 dBm) + **sector antenna gain towards the user** − path loss. The gain (`bsSectorGainDb.m`) is the element pattern of TR 36.873 Table 7.1-1 / TR 38.901 Table 7.3-1 (θ3dB = φ3dB = 65°, SLA_V = A_m = 30 dB, G_E,max = 8 dBi) plus the array factor of a 10-element vertical column, d_V = 0.5λ, all elements on one port with 102° electrical tilt (TR 38.901 Table 7.8-1, UMa); the peak is 8 + 10 = 18 dBi. The direction is the geometric (LOS) one from the site, or its nearest wrap-around copy, as in the large-scale calibration of Table 7.8-1 ("pathloss considering the LOS angle"); satellite EIRP is derived from the **EIRP density** of TR 38.821 Table 6.1.1.1-1 (LEO-600 S-band) and the channel bandwidth, with RF power kept distinct from amplifier electrical draw. Terrestrial candidates carry **inter-site interference** — every base station transmits continuously on the whole band (frequency reuse 1), so for candidate *b* the ratio is `P_b / (N + Σ_{j≠b} P_j)`. There is no interference between the two segments (2.0 vs 3.5 GHz). Between satellites, every other beam on the same band is interference (TR 38.821 §6.1.1.1, Fig. 6.1.1.1-2), again full buffer; only satellites above the elevation mask (those illuminating the area) count, so side-lobe interference is not modelled and the figure is a lower bound. The frequency reuse factor is an option (TR 38.821 Table 6.1.1.1-5: Option 1 = 1, default; Option 2 = 3, each beam then gets 30/FRF MHz). With one visible satellite per band, SINR equals SNR.
4. **Node selection** with **hysteresis**: the candidate with the highest SINR wins, but a *switch* away from the current serving node additionally requires the challenger to lead by a margin `Hys` and to keep that lead for a confirmation time `TTT` — the A3 criterion of TS 38.331 §5.5.4.4 transposed to this rule. The serving node is left only by such a handover or by a **radio link failure** (see *Radio link monitoring and RLF* below), never because one sample dipped. At first attachment, after outage and after an RLF, the best candidate by measurement is chosen among those at or above the minimum usable SINR (≈ −7.53 dB, the Shannon-equivalent of MCS 0, TS 38.214); if there is none the user is **out of coverage**.
5. **Resource allocation**: the serving node's bandwidth is split equally among its users; per-user capacity is Shannon, capped at the NR maximum spectral efficiency (5.5547 bit/s/Hz — MCS 28 of the ≤64QAM table, TS 38.214 Table 5.1.3.1-1).
6. **Service classification**: `Served` / `BelowTarget` / `Outage`, against the 5th-percentile user spectral-efficiency requirement (0.3 bit/s/Hz, TR 37.910 Table 5.4.1.1.1-1). Capacity and *delivered throughput* are reported separately, and the reason a candidate was unusable (below SINR floor / not visible / outside model range) is recorded per link.
7. **Energy**: per-user energy-per-bit (EARTH linear power model for the BS, Auer et al. 2011; linear PA-efficiency model for the satellite) plus a network-level power / bit-per-joule tally, including idle-BS power and a symmetric amplifier-only figure for fair BS-vs-satellite comparison.
8. **Paired policy comparison**: every step is additionally scored under two reference policies — terrestrial-only and satellite-only — and under the hybrid and terrestrial-only policies **with Mobility Load Balancing**, all on the *same* channel realization, with resource allocation recomputed for each. No extra random draws are taken, so the comparison is paired and the statistic is the per-pass difference.
9. **Reporting**: a per-user results table, a 3D plot of the topology and serving links, and a versioned run folder with the exact parameters used.

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
runSimulation(struct('numPasses', 10))            % quick check
runSimulation(struct('geometryGrid', false))      % random instead of stratified pass geometry
runSimulation(struct('hysteresisDb', 0, ...       % reproduce the no-hysteresis rule
                     'timeToTriggerMs', 0))
runSimulation(struct('label', 'v2'))              % tag the versioned run folder
```

**Run time.** Passes run in parallel (Parallel Computing Toolbox, option `numWorkers`, default 6; `0` runs them one after another). Each pass seeds its own generator (`rng(rngSeed + passIdx, 'twister')`) and resets channel and decision state, so it depends only on its index. Passes are computed in groups of `numWorkers`, then written to the dataset, the convergence curve and the checkpoint in pass order, so the output files are byte-identical to a serial run (checked on 2 and 12 passes: dataset, handover events and `simulation_results.mat`). Inside a pass, everything deterministic — distances to all station copies, path loss, LOS probability, sector gains, satellite geometry, SINR and L3 filtering — is computed for all users at once; only the random draws stay in the per-link loop, in the original order (user → site → satellite), which is what keeps results bit-identical. Measured on a 6-core laptop (i7-10750H): at Δt = 0.2 s, with 10–30 users, 6 passes take about 3.8 minutes, so a full run of 300 passes takes about 2.6 hours (it was half an hour at Δt = 1 s). With 570 users (10 per TRxP) the time grows about five-fold as well, from about 11 hours at Δt = 1 s. A checkpoint is written after every pass, so an interrupted run resumes from where it stopped rather than from the start.

### How the simulated time is structured

The run is a loop over **satellite passes**; inside each pass a 0.2 s time loop advances both the satellite and the users. Every timing choice is taken from a standard rather than picked:

| Quantity | Value | Source |
|---|---|---|
| Elementary window | one pass, 900 s, centred on the middle satellite's closest approach and containing the passes of all three | TR 38.821 Table 4.2-3 NOTE 1 — "a period of time corresponding to the visibility time of the satellite" |
| Time step | 0.2 s | the L1-to-L3 period of the TR 36.839 mobility calibration (Table 5.3.2.1, 200 ms) and the FR1 measurement period of TS 38.133 Table 9.2.5.2-1 (200 ms), so each step is one measurement and the L3 filter needs no rescaling; it also meets ITU-R M.2412-0 Annex 1 §5.3.2 (displacement below 1 m per step: 0.17 m) and TR 38.821 Table 7.3.2.1.4-1 (fastest LEO mobility timescale 6.61 s). At 1 s, TTT values below 1 s and the 1 s T310 could not be resolved |
| User speed / direction | 3 km/h, fixed per pass, uniformly random azimuth | ITU-R M.2412-0 §8.4, TABLE 5 b)/c) — "fixed and identical speed of all UEs of the same mobility class, randomly and uniformly distributed direction" |
| Total duration | fixed, 300 passes | ITU-R M.2412-0 §7.1 — "a sufficient number of drops … to ensure convergence" plus "the width of confidence intervals"; extended to satellite evaluation by ITU-R M.2514-0 §8.2.4 |

No standard fixes a total duration — TR 38.821's own NTN system-level calibration table (6.1.1.1-5) has no duration field at all. The run therefore uses a **fixed** number of passes, the same for every scenario, so that pass *n* carries the same randomness everywhere and scenarios can be compared pass-by-pass. The confidence interval no longer stops the run; it is **reported**, as §7.1 asks. The target is a 95% half-width, computed **across passes**, below 5% (relative, for rate and energy metrics) or 2 percentage points (absolute, for fractions). The previous convergence-driven run met it at pass 227, which is where the default of 300 comes from. At the end the run states whether the target was met and from which pass it held uninterrupted. The convergence trace is saved as a result in its own right (`convergence.png`, with the numbers in `simulation_results.mat`), which is what makes the sample size defensible instead of arbitrary.

Each pass is an independent *drop* in the ITU-R sense: users are re-dropped at new positions and the channel state is reset between passes, while being threaded step-to-step **within** a pass. `PassID` is therefore the unit for the ML train/test split, for confidence intervals, and for any per-repeat aggregation — users inside one pass are not independent observations.

Pass geometry varies: the RAAN is offset per pass so transits range from grazing (peak elevation at the 20° mask) to near-zenith. The largest useful offset is found numerically at startup rather than assumed. The offset is **stratified**: the range is split into `numPasses` equal strata and each pass takes the midpoint of one, in a shuffled order drawn from its own random stream. The target distribution is the same uniform one as before, but every part of the geometry range is guaranteed equal coverage, and the shuffle keeps the convergence curve from sweeping the geometry monotonically.

### Hysteresis and time-to-trigger

The decision rule carries a margin and a confirmation time, transposed from the A3 event of TS 38.331 §5.5.4.4 (`Mn − Hys > Mp` to enter, `Mn + Hys < Mp` to leave, with all per-cell offsets set to zero). Both parameter ranges come from §6.3.2 of the same specification: hysteresis is any multiple of 0.5 dB up to 15 dB, and time-to-trigger is one of sixteen enumerated values. The motivation is named by the standard itself — TR 38.821 §7.3.2.1.3 flags UE *ping-ponging* between cells as a high-priority NTN problem.

The defaults are **Set 1 of the TR 36.839 mobility calibration** (Table 5.3.2.1): TTT 480 ms, A3 offset 3 dB (here the hysteresis), L3 filter K = 4, L1-to-L3 period 200 ms, at the 3 km/h of the same table. TTT is counted in whole steps, rounded up, because the condition must *hold* for that long: 480 ms is 3 steps (0.6 s) at Δt = 0.2 s. (Earlier versions used 2,560 ms, picked from a sweep; with a 1 s ping-pong window that made a ping-pong between two A3 handovers impossible — every return within 1 s came from the old immediate switch on a weak link.)

### Radio link monitoring and RLF

A serving link is not dropped when one sample is weak. As in TS 38.331 §5.3.10, with the values of TR 36.839 Table 5.2.1.3.1: SINR below **Qout = −8 dB** is an out-of-sync indication, above **Qin = −6 dB** an in-sync one; **N310 = 1** out-of-sync starts **T310 = 1 s** (the 36.331 default), **N311 = 1** in-sync while it runs stops it, and if it expires a **radio link failure** is declared and the user re-establishes (§5.3.7) on the best usable candidate, or goes out of coverage. One indication per step; the measurement of a step is its 200 ms average, which is the Qout window of TR 36.839 (its Qin window is 100 ms). While T310 runs the user stays attached; below the MCS 0 SINR it gets no throughput. Handover failures follow Definition 3 of TR 36.839: the handover command arrives while T310 runs (the user stays and T310 keeps running, §5.4.1), an RLF occurs while the A3 entry condition holds, or the target is below Qout at execution. The recovery model is left open by the TR ("companies should have the flexibility to choose a realistic RLF recovery model", §5.4.1); here it is the same cell selection as first attachment, and its duration is not modelled. Only successful A3 handovers count as handovers (Definition 6: "successful handovers excl. handover failures"); RLFs (Definition 2, reported per user-minute) and the handover-failure share (Definition 4) are separate KPIs. After an RLF or failure the time-of-stay chain restarts (§5.2.2: "whenever there is a handover failure, the time of stay should not be logged").

### Mobility Load Balancing (the `*MLB` policies)

With a satellite visible almost all the time and 19 co-channel, full-buffer base stations interfering with each other, a rule that looks only at SINR sends about half the users to the satellite — even users inside the cells — whose bandwidth is then split among them. MLB is the standard's answer: TS 38.300 §15.5.1 defines it as distributing load "by means of optimization of cell reselection/handover parameters", with the change expressed as a delta on the handover trigger (§15.5.1.4) and kept within OAM limits. In the A3 formula that is the per-cell offset (Ocn/Ocp), zero in the plain policies. The decision stays SINR-based; only the threshold moves with the load.

The standard specifies the mechanism, not how much the offset should move. Two published pieces fill that in:

- **Offset rule** — Nasri & Altman, *Handover adaptation for dynamic load balancing in 3GPP LTE systems* (arXiv:1307.1212), eq. 4–5: `HM(e,k) = f(χe − χk)` with `f(x) = f(0) + (f(0) − HMmax)·x`, `f(0) = 6 dB`, `HMmax = 12 dB`. Only the load term is transposed (the planned margin here is the A3 hysteresis): the target gains an offset of `6 dB × (χ_serving − χ_target)`, the difference clipped to [−1, 1]. The offset is the Ocn that the serving cell gives the neighbour, so it is rounded to the values the field can actually take — `Q-OffsetRange` of TS 38.331: −24, −22, …, −8, −6…+6 in 1 dB steps, 8, …, 24 dB, i.e. whole dB within the ±6 dB used here (rounded on the magnitude, so the rule stays antisymmetric). The rule is antisymmetric (`f(x) + f(−x) = 2f(0)`), which is the condition the paper proves is needed to avoid ping-pong — the loaded cell pushes users out and the light cell delays sending them back.
- **Load** — the *virtual load* of Awada et al., *A game-theoretic approach to load balancing in cellular radio networks* (PIMRC 2010), eq. 1–3: each user needs the share `κ_u = D_u / (R(SINR_u)·N_tot)` of the node's resources to be satisfied, and `ρ_c = Σ κ_u` may exceed 1. With full buffer every PRB is always occupied, so PRB occupancy cannot tell a full node from an empty one; the virtual load can. The requirement `D_u` is the same service target used for classification, so `κ_u = 0.3 / SE_u` and `ρ_c ≤ 1` means the node can serve all its users.

Loads are reported periodically (load reporting, TS 38.300 §15.5.1) with a Reporting Periodicity of 1000 ms (option `mlbReportPeriodMs`; TS 38.423 §9.1.3.18 allows 500, 1000, 2000, 5000 and 10000 ms) and the last report holds in between, so a decision never sees the load it is about to create. At Δt = 1 s this was one report per step; at 0.2 s it is one every five. The offset applies only to the comparison against the serving node; first attachment and recovery from a lost link have no cell pair and use the plain measurement. The plain `Actual`, `TerrestrialOnly` and `SatelliteOnly` policies are unchanged — their results are bit-identical with and without the MLB code.

### Moving users, and what is checked every step

Users walk continuously during a pass. That is also what activates the correlated shadow fading and the spatially-consistent LOS state — with static users the per-step displacement was zero and both models were inert.

Because the users move, the run asserts on **every step of every pass** that their state is what it should be: antenna height still exactly 1.5 m, coordinates finite and in range, and every user still holding at least one terrestrial link inside the TR 38.901 Table 7.4.1-1 validity box (10 m ≤ d2D ≤ 5 km). It also checks the outputs for consistency — outage rows carry no capacity, zero throughput and zero load; attached rows carry finite, non-negative capacity — positive whenever the SINR is above the MCS 0 floor, zero only while T310 runs below it or when the serving link has left the model (SINR NaN) — and load ≥ 1. Any violation aborts the run naming the pass, step and user. It also asserts that no user is ever closer than 10 m to any base station. Initial positions and headings are sampled so that no walk can leave the validity box, so these assertions are a check rather than a correction.

### Outputs

| File | Content |
|---|---|
| `Dataset/dataset.csv` | one row per pass/step/user — geometry, per-candidate diagnostics, serving decision, capacity, delivered throughput, service state, link state, energy. Written every 5th step. Feeds both the results and the ML side. Gitignored. |
| `Dataset/handover_events.mat` | one record per successful A3 handover of the `Actual` and `ActualMLB` policies, at full 0.2 s resolution (the 5 s stride of the CSV can hide a quick return). Features are only what is known **at the decision**, and every column has a source. From Aggarwal, Shyamsundar & Saxena, *Predicting Unstable Handovers Using Machine Learning in LTE/5G Networks*, IEEE ICSigSys 2026, §IV-A/C: serving and candidate measurement (here the L3-filtered SINR the A3 rule uses), their gap, the 1 s change of each, their mean and least-squares slope over the last 3 s, and BS distance. From the UE History Information of TS 38.413 §9.3.1.97: time stayed in the source node ("Time UE Stayed in Cell") and whether the target is the previously visited node. From the Resource Status Update of TS 38.423 §9.1.3.21: load and number of active UEs (§9.2.2.62) of both nodes, as reported in the previous step — inputs listed by TR 37.817 §5.3.2.4. Also satellite elevation (the satellite's geometry, as distance is the BS's) and the MLB offset (an A3 parameter, like the hysteresis feature of Aggarwal et al.). The 1 s / 3 s features use the samples at t, t−1 s and t−2 s, as in the per-second dataset of that paper. What happened **after** is kept separately as outcome: time on the new node (`StayTo_s`), where the user went next (`NextNode`), how it left (`NextCause`: `HO`, or `RLF`/`HOF` for re-establishment after a radio link failure or a failed handover — an RLF that re-establishes on the same node also ends the stay), and `Censored` if the pass ended first. A ping-pong is `NextNode == FromNode && NextCause == "HO" && StayTo_s <= pingPongMaxStay_s` (in the file's metadata: 0.8 s, i.e. a stay *shorter* than the 1 s MTS); this reproduces the simulation's own count exactly. The unstable-handover label of Aggarwal et al. (stay ≤ 15 s or a return within 15 s) can be built from the same columns. One variable per column, `-v7`, so Python reads it with `scipy.io.loadmat`. Gitignored with the rest of `Dataset/`. |
| `Results/policy_gain.png` | how much better the hybrid policy is than terrestrial-only and satellite-only: mean paired per-pass difference with 95% CI. Rate and energy as relative change, user fractions in percentage points; positive always means the hybrid is better (outage sign flipped). Available spectrum is printed in the title as context, not plotted as a gain |
| `Results/mlb_gain.png` | what MLB adds: hybrid with vs without MLB, terrestrial-only with vs without MLB, and hybrid vs terrestrial-only when both have MLB — same paired statistic and sign convention as `policy_gain.png` (fewer handovers and a lower ping-pong share count as better) |
| `Results/handover_rates.png` | successful A3 handovers per user per minute for each policy, stacked by type (terrestrial↔terrestrial, terrestrial↔satellite, satellite↔satellite), with the 95% CI of the total; each bar is labelled with the ping-pong share, the handover-failure share (TR 36.839 Definition 4) and the RLF rate. A ping-pong is a return to the node before the previous handover after a stay shorter than `pingPongMtsS` (default 1 s, the MTS of TR 36.839). First attachment, reconnection after outage and re-establishment after RLF are not handovers |
| `Results/convergence.png` | per-pass CI half-widths of the tracked KPIs, against the precision target |
| `Results/simulation_results.mat` | the numbers behind every figure: reference pass at full 0.2 s resolution, convergence trace, per pass × policy KPIs, paired differences, summary |
| `Results/temporal_*.png` | time series over the reference pass |
| `Results/network_3d.png` | topology and serving links at closest approach |
| `Results/runs/<timestamp>_runSimulation[_label]/` | versioned copy of the `.mat`, the policy and 3D figures, the time series and the convergence plot, plus the exact parameters |

`Dataset/dataset.csv` is the **only** CSV the MATLAB side writes; it and `Dataset/handover_events.mat` are the inputs of the ML side.

### Verification scripts

```matlab
geometryValidation()                  % antenna-height & distance regression check
energyModelValidation()               % user-count dependence of the energy metrics
linkBudgetValidation()                % 75 checks over the link budget, resource allocation, soft LOS and RLF
```

All three print explicit pass/fail lines, write a PNG, store their result table in the run's `params.mat` (field `results`), and `error()` out on the first failure.

`geometryValidation` confirms the user height reaching the path-loss model never drifts with distance, that 3D range is consistent with 2D distance and heights, and that the 5 km validity gate fires — and quantifies the error the earlier ENU-based geometry introduced. `energyModelValidation` shows that per-user capacity scales as 1/L while energy-per-bit stays flat (the user count cancels algebraically), and that the network bit/J metric — unlike energy-per-bit — does respond to user composition.

`linkBudgetValidation` compares the code's value against an independently computed closed form, at double-precision tolerance, in ten groups: unit conversions (dBm/dBW/W and both EIRPs, including a round-trip that would fail if the satellite EIRP density were misread as amplifier power), `fspl` against `20·log10(4πd/λ)`, `kTB` noise, the identity `SINR = EIRP − PL − N` on the returned per-candidate values, known behaviour, bandwidth conservation (each user's share recovered as `C_u/SE_u` and summed per node), the network power balance with no double counting, four pre-specified edge cases — no candidate at all, a single user, a `−Inf` tie, and elevation exactly at the mask — and the sector antenna: pattern values at known directions (peak, −3 dB at φ = 32.5°, 30 dB back lobe), the array factor against a loop over the complex weights of TR 36.873 Table 7.1-1, and the SINR of a 3-sector site through `simulateScenario` against `max S_c / (N + other sectors)`; and soft LOS and RLF: `LOS_soft = 1` at P_LOS = 1, `LOS_soft` in [0, 1] and above 1/2 exactly when the hard state is LOS, and on a single link driven below Qout by a 60 dB power drop — the user stays with T310 started and no throughput, an RLF and outage follow when T310 expires without in-sync (counted as an RLF, not a handover), and T310 stops when the SINR is back above Qin in time.

Two of the known-behaviour checks are worth calling out, because they are what detects whether interference is really in the computation: with a **single** base station the identity `SINR = EIRP − PL − N` holds exactly and ±3 dB of transmit power moves SINR by exactly ±3 dB; with **two**, SINR is strictly below that identity for every user and ±3 dB moves it by strictly less — the signature of an interference-limited system. The same is checked between satellites: two visible satellites on one band give SINR strictly below SNR; with frequency reuse 3 they are on different bands and SINR equals SNR exactly; and a satellite below the mask does not interfere.

### Result versioning

Every run ends by calling `saveRunVersion`, producing:

```
Results/runs/<YYYYMMDD_HHMMSS>_<script>[_<label>]/
  params.mat      exact parameter restore
  params.txt      readable snapshot, incl. git commit, MATLAB version, dataset SHA-256
  changed.txt     parameter diff vs. the previous run of the same script
  <CSV/PNG output of that run>
```

The dataset itself is not copied in (it is tens of MB and this folder is tracked); its SHA-256 and byte size go into `params.txt` instead, which is what ties a set of results to one specific data file rather than to a filename.

### Part 2 (ML)

```
pip install -r Model/requirements.txt
python Model/train_model.py
```

Three passes on the generated dataset — exact SINR, noisy SINR, and geometry-only — predicting the serving type from candidate-level features, split **by pass** (`PassID`) so no pass contributes to both train and test. No decimation happens on the Python side: the CSV is already sampled every 5 s at source. Each run also writes `split.json`, per-sample `predictions_<Model>.csv` and the fitted `pipeline_<Model>.joblib`, so any published metric can be recomputed without retraining. See `Model/README.md` for the current numbers.

## Project structure

```
PROD/
  simulateScenario.m        Core: link budget, fading, node selection, allocation, energy, service state
  bsSectorGainDb.m          Sector antenna gain (TR 36.873 Table 7.1-1, TR 38.901 Table 7.8-1)
  runSimulation.m           The entry point — one run: passes, moving users, precision report, dataset
  geometryValidation.m      Regression check on link geometry / antenna heights
  energyModelValidation.m   Check of how the energy metrics depend on user count
  linkBudgetValidation.m    75 numeric checks on the link budget, resource allocation, soft LOS and RLF
  saveRunVersion.m          Versioned run folders (parameters, commit, diff, outputs)
  array.m                   Prints the per-user results table
  visual.m                  3D plot of base stations, users, satellite, serving links
Results/runs/               Versioned outputs of every run (tracked in git)
Dataset/                    CSV output of runSimulation.m (generated, gitignored)
Model/                      Part 2: Python training scripts, metrics and plots
```

