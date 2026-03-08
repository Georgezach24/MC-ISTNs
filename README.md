# MC-ISTNs

**Multi-Connectivity in Integrated Satellite–Terrestrial Networks (ISTNs)**

A MATLAB time-slotted simulator for studying traffic management, link selection, and handover dynamics in integrated terrestrial–non-terrestrial network environments.

This repository is developed in the context of the thesis project:

**ML-Aided Multi-Connectivity Management in Integrated Satellite–Terrestrial 6G Networks**

---

## Overview

This project implements a discrete-time MATLAB simulator for an Integrated Satellite–Terrestrial Network (ISTN) scenario, where a user can be served by either a terrestrial link (TN) or a non-terrestrial / satellite link (NTN).

The simulator is designed as a **system-level research tool** for analyzing how time-varying link quality affects:

- serving-link selection,
- queue evolution,
- traffic delivery,
- and handover behavior.

The current implementation provides a clean baseline for future extensions toward intelligent multi-connectivity control and machine learning-based decision mechanisms.

---

## Current Scope

At its current stage, the simulator focuses on the following core mechanisms:

- slot-based evolution of TN and NTN link conditions,
- achievable-rate computation for both links,
- traffic generation for heterogeneous services,
- queue-based packet servicing,
- persistent serving-link management,
- and handover execution with temporary transmission interruption.

In other words, the simulator has already moved beyond a simple “best-link-per-slot” policy and now includes a stateful handover process.

---

## System Model

The current implementation considers a **single-user, two-link access setting** consisting of:

- **1 User Equipment (UE)**
- **1 terrestrial access link (TN)**
- **1 non-terrestrial / satellite access link (NTN)**

The simulation operates in **discrete time** with slot duration:

- **Δt = 1 ms**

At each slot, the simulator updates the system state and records the resulting network behavior.

---

## Per-Slot Simulation Flow

For each simulation slot, the following steps are executed:

1. **Link update**  
   The quality of the TN and NTN links is updated according to time-varying traces.

2. **Rate computation**  
   The achievable transmission rate of each link is estimated.

3. **Traffic generation**  
   New traffic arrivals are generated and inserted into their corresponding queues.

4. **Handover update**  
   The handover controller evaluates whether the serving link should remain unchanged or be switched.

5. **Traffic servicing**  
   Queued traffic is served through the currently active serving link, subject to the available slot capacity.

6. **Logging**  
   Key performance variables are stored for later analysis and visualization.

This structure provides a modular baseline that can be extended gradually with more realistic or more advanced control logic.

---

## Link and Rate Modeling

The current simulator uses an abstraction-oriented link model rather than a full physical-layer implementation.

The implemented link modeling includes:

- time-varying TN and NTN quality evolution,
- distinct average conditions for terrestrial and satellite links,
- slot-level achievable-rate estimation,
- configurable bandwidth,
- and a configurable efficiency factor.

At this stage, link rates are computed through a lightweight Shannon-like approximation, which is sufficient to capture the relative differences between TN and NTN behavior while keeping the simulator simple, transparent, and computationally efficient.

The objective is not exact PHY-level replication, but a realistic-enough basis for studying control and scheduling behavior.

---

## Traffic Model

Two service classes are currently supported:

### URLLC
URLLC traffic represents delay-sensitive traffic and is assigned **higher scheduling priority**.  
New URLLC arrivals are inserted into a dedicated queue and are always served before eMBB traffic whenever transmission capacity is available.

### eMBB
eMBB traffic represents throughput-oriented traffic and is stored in a separate queue.  
It is served after URLLC, using any remaining capacity of the selected serving link.

This queue-based structure allows the simulator to capture basic service differentiation and congestion behavior.

---

## Handover Model

A basic TN ↔ NTN handover mechanism has already been incorporated into the simulator.

Instead of selecting the best link independently at every slot, the system now maintains a **persistent serving link**, and handover is triggered only when the alternative link remains sufficiently better for a required amount of time.

The current handover model includes:

- **filtered link metrics**,
- **candidate target-link tracking**,
- **Time-to-Trigger (TTT)**,
- and **handover interruption duration**.

This makes the simulation significantly more realistic than instantaneous slot-by-slot link switching, since it allows the study of:

- unnecessary switching reduction,
- temporary service interruption during handover,
- and queue buildup caused by handover-related transmission pauses.

---

## Logged Variables

The simulator currently records, at each time slot:

- TN achievable rate,
- NTN achievable rate,
- current serving link,
- handover active state,
- delivered bits per slot,
- URLLC queue occupancy,
- eMBB queue occupancy.

These variables are sufficient to analyze the basic dynamics of the system and validate the handover behavior.

---

## Current Outputs

The current implementation can produce plots such as:

- **TN and NTN rate evolution over time**
- **serving-link evolution over time**
- **handover interruption state**
- **URLLC and eMBB queue evolution**
- **running average delivered throughput**

These outputs are useful for verifying that the handover mechanism behaves as expected and for identifying the performance impact of link transitions.

---

## Implementation Status

### Implemented
- Discrete-time simulation framework
- Time-varying TN/NTN link model
- Achievable-rate calculation for both links
- URLLC and eMBB traffic queues
- Priority-based service mechanism
- Persistent serving-link state
- Basic handover controller
- Handover interruption modeling
- Per-slot logging and visualization

### Planned Extensions
- More realistic terrestrial and satellite channel models
- Explicit delay and deadline tracking for URLLC
- Reliability evaluation for URLLC traffic
- Split and duplication transmission modes
- Uplink modeling
- Multi-node topology (multiple BSs / satellites)
- Dataset generation for ML / RL training
- Learning-based handover and multi-connectivity control

---

## Research Direction

The long-term goal of this simulator is to serve as a flexible framework for evaluating intelligent multi-connectivity strategies in integrated 6G satellite–terrestrial systems.

More specifically, the project aims to study the trade-offs among:

- throughput,
- latency,
- reliability,
- handover cost,
- and resource utilization,

while progressively moving toward ML-assisted control policies.

---

## Notes

This repository is intentionally developed as a **research-oriented abstraction model**.

It does **not** currently aim at full 5G/6G PHY-level accuracy.  
Instead, it provides a modular and extensible platform for controlled experimentation, rapid prototyping, and system-level evaluation of TN/NTN management strategies.
