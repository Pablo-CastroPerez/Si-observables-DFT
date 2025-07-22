# Silicon: DFT Simulation and Observable Extraction

This repository contains input files, scripts, and results for a density functional theory (DFT) study of **Silicon**, performed with **Quantum ESPRESSO**. The goal is to compute observables such as the bulk modulus, equilibrium volume, and equilibrium energy by fitting total energies to several equations of state (EOS).  

For details, see the results PDF [`Si_observables_results.pdf`](./Si_observables_results.pdf).

---

## Files

| File                      | Description                                                               |
|---------------------------|---------------------------------------------------------------------------|
| `si.scf.in`               | Self-consistent field (SCF) input file for QE                             |
| `QE-scf.slurm`            | SLURM job script for HPC submission ([Archer2](https://www.archer2.ac.uk))|
| `Si.pbe-n-van.UPF`        | Norm-conserving pseudopotential for Si (PBE GGA)                          |
| `Silicon_observables.py`  | Python script for fitting EOS and extracting observables from SCF data    |
| `si_energies`             | Raw energy vs. volume data file generated from QE outputs                 |
| `Si_observables_results.pdf` | PDF report summarising the extracted observables and fitting analysis  |

---

## Methods

- **DFT code:** [Quantum ESPRESSO](https://www.quantum-espresso.org/)
- **Exchange-correlation:** PBE functional (GGA)
- **Equations of State fitted:**
  - Murnaghan
  - Birch–Murnaghan
  - Vinet

---

## Results Summary

- Extracted observables include:
  - Equilibrium volume \( V_0 \)
  - Bulk modulus \( B_0 \)
  - Pressure derivative \( B_0' \)
  - Equilibrium energy \( E_0 \)

See [`Si_observables_results.pdf`](./Si_observables_results.pdf) for plots and fitted curves.

---

## Tools

Badges for tools used:

![Python](https://img.shields.io/badge/Python-3670A0?style=flat-square&logo=python&logoColor=white)
![Quantum ESPRESSO](https://img.shields.io/badge/Quantum--ESPRESSO-black?style=flat-square&logo=codeforces&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-11557C?style=flat-square&logo=matplotlib&logoColor=white)
![NumPy](https://img.shields.io/badge/NumPy-013243?style=flat-square&logo=numpy&logoColor=white)

---

## Note

This was done for educational purposes only and contain certain oversimplifications, such as calculating the sound speed directly from the bulk modulus instead of elastic coefficients. The results are not meant to be used as rigorous data.
