<!--
Section: Calibration to physical units
Update when: a calibration attempt is run, accepted or rejected, or the calibration tools change.
Code: chromatin_studio/calibration.py, chromatin_studio/kinetics.py, chromatin_studio/components/core.py (PROVISIONAL_CLOCK).
Sources: docs/PROJECT_STATUS.md "Physical calibration is unfinished", references/hansen/.
-->

# Calibration: from reduced units to nanometres and seconds

::: {.box}
**Short version.** Comparing the simulated motion of two points on a 2 Mb chain with live-cell tracking of the Fbn2 locus gives about **47 nm and 0.2 s per reduced unit**. The fit is good on the training data, but the held-out data set fails, the time estimate depends strongly on the assumed localisation noise, and the mapping has not been accepted. It is used only, and always labelled, as provisional.
:::

## The method

The calibration fixes length first, then time:

1. **Length.** Choose two beads separated by the same genomic distance as two labelled sites in the experiment (here 510–520 kb). Match the simulated mean-square separation to the measured one. That gives nanometres per reduced length unit.
2. **Time.** Compute the mean-square displacement (MSD) of that separation over lag time in the simulation and the experiment. Find the single time scale that makes the curves overlap. That gives seconds per reduced time unit.
3. **Hold out.** Test the fixed mapping, without refitting, on a separate set of movies.

The reference data are the public Hansen-lab Fbn2 measurements: 12 conditions and 2,894 trajectories, processed and kept separate. The fitting set has 424 trajectories from 74 movies; 113 trajectories from 18 other movies are held out.

## What the runs found

The baseline model (2,000 beads × 10 kb, cohesin depleted) was run with 8 replicas to 96,000 reduced time units on a GPU, discarding the first 48,000 as burn-in.

- **Stationarity.** Approximately reached after 48,000 units; results barely move for burn-ins of 24,000 to 64,000.
- **Length.** 47.2 nm per unit (95% interval 40.6–53.9). Pooling all interior bead pairs at the same separation gives 51.3 nm.
- **Time.** 0.205 s per unit (bootstrap 0.13–0.32). The MSD slope matches well (0.57 simulated against 0.54 measured).
- **Time step.** Converged: halving or quartering $\Delta t$ changes bead and pair MSD by less than 0.5%.
- **Held-out test fails** (log error 0.357 against a threshold of 0.25). But the two experimental acquisitions differ from **each other** by the same amount, so no single mapping could pass this test as written.

![**The clock depends on the noise model.** Fitted seconds per reduced time unit (left axis, bars) and fit error (right axis, line) when the assumed localisation-noise offset is scaled by $s$; $s = 1$ is the reference correction. A ±20% change in noise moves the clock by −21% to +48%, while length barely changes. Values from `docs/PROJECT_STATUS.md` (2026-11-08).](figures/cal-noise.svg){#fig:noise width=80%}

## What is needed to accept it

- A held-out criterion that allows for the measured difference between acquisitions.
- An independent measurement of localisation noise for this data set (this would constrain the clock more than more simulation).
- A decision on the length estimator (single pair or pooled interior pairs).

Even once accepted, one polymer clock gives consistent units; it does not calibrate protein binding, marks or transcription. Each of those needs its own data.
