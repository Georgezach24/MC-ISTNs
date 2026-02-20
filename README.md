# MC-ISTNs
**Multi-Connectivity in Integrated Satellite–Terrestrial Networks (ISTNs)** — MATLAB time-slotted simulator for generating datasets and benchmarking multi-connectivity policies (URLLC + eMBB) over TN (terrestrial BS) and NTN (LEO satellite) links.

This repository is part of a thesis project on:

**ML-Aided Multi-Connectivity Management in Integrated Satellite–Terrestrial 6G Networks**

---

## Overview

This project implements a discrete-time simulator of an Integrated Satellite–Terrestrial Network (ISTN) environment with multi-connectivity support.

The objective is to:

- Model heterogeneous TN–NTN links
- Support URLLC and eMBB services simultaneously
- Evaluate multi-connectivity modes (TN, NTN, Split, Duplication)
- Generate datasets for ML / RL-based controllers

The simulator abstracts PHY-level behavior while preserving realistic performance trade-offs.

---

## System Model

The simulated topology includes:

- 1 Mobile User Equipment (UE)
- 2 Terrestrial Base Stations (TN domain)
- 2 LEO Satellites (NTN domain)
- Downlink (DL) and Uplink (UL) modeling
- Discrete-time simulation with slot duration Δt = 1 ms

### Channel Modeling

The simulator includes:

- TN fading: Rayleigh
- NTN fading: Nakagami-m (implemented without external toolboxes)
- Distance-based pathloss with configurable exponent
- Shannon-based achievable rate approximation
- Propagation delay differentiation between TN and NTN

The goal is not exact PHY accuracy, but realistic relative behavior between terrestrial and satellite links.

---

## Supported Services

### URLLC (Ultra-Reliable Low-Latency Communications)

- Packet size: 64 Bytes
- Periodic arrivals (100 packets per second)
- Strict latency deadline constraint
- Per-packet reliability evaluation

The simulator measures whether each URLLC packet:
- Is delivered
- Meets the latency deadline
- Contributes to overall reliability

---

### eMBB (Enhanced Mobile Broadband)

- Backlogged traffic model
- Throughput-oriented evaluation
- DL and UL accounting

eMBB performance is measured primarily in terms of throughput.

---

## Multi-Connectivity Modes

The simulator supports per-slot multi-connectivity decisions.

### URLLC Modes
- TN (Terrestrial only)
- NTN (Satellite only)
- SPLIT (Traffic distributed across TN and NTN)
- DUPLICATION (Same packet transmitted over both links)

### eMBB Modes
- TN
- NTN
- SPLIT (Aggregated throughput over both links)

These modes allow evaluation of:
- Reliability improvement through duplication
- Throughput enhancement through aggregation
- Trade-offs between latency, energy, and resource usage

---

## QoS Scenarios

Two QoS configurations are implemented to stress the system differently.

### Scenario A (Standard Profile)
- URLLC deadline: 20 ms
- Target reliability: 99.99%
- eMBB soft latency reference: 200 ms

This represents a realistic near-term 6G operating profile.

---

### Scenario C (Strict / 6G-Extreme Profile)
- URLLC deadline: 10 ms
- Target reliability: 99.999%
- eMBB soft latency reference: 150 ms

This configuration stresses the system under very tight URLLC constraints and highlights the value of intelligent multi-connectivity control.


