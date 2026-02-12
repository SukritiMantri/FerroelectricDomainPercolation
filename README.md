# Domain Percolation in Ferroelectric Polycrystals

This repository contains a C++ implementation for studying **domain-wall
percolation in ferroelectric polycrystalline microstructures**, with a focus on the
role of **grain boundaries, crystallographic misorientation, and special (Σ3) grain
boundary populations and Textured microstructures**.

The code generates synthetic polycrystalline microstructures, assigns crystallographic
orientations under orientation constraints, and evaluates
**domain percolation pathways across grains** based on geometric matching and
charge neutrality criteria.

---

## Scientific Background

- Ferroelectric domain walls propagate across grains through **grain boundary–mediated
  interactions**, which depend on:
  - crystallographic misorientation,
  - domain wall orientation,
  - polarization compatibility.
- In polycrystalline materials, long-range functional behavior (e.g. switching,
  conductivity, dielectric response) depends on the **connectivity of these domain-wall
  networks**, not just local domain structure.
- **Special grain boundaries**, such as Σ3 boundaries, are known to influence domain
  percolation and are explicitly modeled in this code.

This implementation enables:
- statistical analysis of **domain percolation length**,
- evaluation of **total domain percolating grain clusters**, and
- comparison between random, textured (MD), and Σ3-weighted microstructures.

---

## Reference and Citation

If you use this code, methodology, or any results derived from it, **please cite** the
following paper:

> **Mantri, S., & Daniels, J.**  
> *Percolation of domain walls in ferroelectric polycrystals.*  
> **Acta Materialia**, 234, 118615 (2022).  
> https://doi.org/10.1016/j.actamat.2022.118615

---

## Repository Structure

- `DomainPercolationFerroelectrics/src`  
  Core C++ source files implementing the domain percolation calculations.

- `DomainPercolationFerroelectrics.xcodeproj/`  
  Xcode project configuration.

- `python/`  
  Python scripts for post-processing and plotting results.

---
## License

This project is released under the **MIT License**.
See the `LICENSE` file for details.


--- 
## Notes

This codebase has been **carefully refactored from a long-standing research implementation**
to improve readability, reproducibility, and extensibility **without altering the
underlying physics or numerical behavior**.
