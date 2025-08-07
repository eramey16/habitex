<p align="center">
  <img src="/habitexlogo.png" width="200"/>
</p>

**HabitEx** is a Python-based tool designed to vet and characterize potentially habitable exoplanets using public data from the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/). The pipeline evaluates exoplanets based on their stellar and orbital properties to determine whether they reside in a **conservative** or **optimistic habitable zone**, and supports custom filtering for survey planning, target selection, and comparative exoplanetology.

---

## Objective

For each confirmed exoplanet, HabitEx:

- Retrieves key planetary and stellar parameters using `astroquery`
- Computes whether the planet lies within the conservative or optimistic habitable zone using models from Kopparapu et al. (2013, 2014)
- Estimates additional properties such as planet density (if mass and radius are available)
- Supports optional filtering by orbital, physical, or observational criteria
- Outputs structured data for follow-up analysis and survey planning


## Functional Overview

### Data Retrieval

- Queries the NASA Exoplanet Archive using `astroquery`
- Defaults to the most recent planet entry if multiple are available
- Allows optional filtering by paper or table (`pscomppars` (default), `cumulative`)

### Habitable Zone Assessment

- Calculates incident stellar flux using either semi-major axis or orbital period (via Kepler’s Third Law, if there are no values returned from astroquery)
- Evaluates whether the planet falls into the following categories, according to Kopparapu et al. 2013:
  - **Conservative Habitable Zone** (e.g. water loss to maximum greenhouse limits)
  - **Optimistic Habitable Zone** (e.g. recent Venus to early Mars)

### Planetary Density Estimation

- If both mass and radius are known, density is calculated
- Rocky planets, water worlds, and gas planets are identified based on density thresholds

### Custom Filtering

Users may apply custom filters on:

- Stellar effective temperature
- Orbital period
- Planet mass or minimum mass
- Radius
- Declination (for observability or site-based filtering)

### Custom Input Mode

- Users can supply their own stellar and planetary inputs (e.g., from simulations or mission concepts)
- Supports offline mode without querying the Exoplanet Archive

---

## Outputs

- Filtered and ranked CSV file of potentially habitable planets
- Flags for:
  - Conservative HZ inclusion
  - Optimistic HZ inclusion
  - Three potential populations (rocky planets, water worlds, and gas planets)
- Plots to visualize planet orbit compared to optimistic and conservative habitable zone
- Mass-Radius diagram for planets falling in the habitable zone

---

## Installation

HabitEx can be installed via pip:

```bash
pip install habitex
```

---

## Dependencies

- `astroquery`
- `pandas`
- `numpy`
- `matplotlib`

---

## References

- Kopparapu et al. (2013), ApJ, 765, 131  
- Kopparapu et al. (2014), ApJ, 787, L29
- Luque et al. (2022), Science, 377, 6611
- NASA Exoplanet Archive API: https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html

[![A rectangular badge, half black half purple containing the text made at Code Astro](https://img.shields.io/badge/Made%20at-Code/Astro-blueviolet.svg)](https://semaphorep.github.io/codeastro/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16756123.svg)](https://doi.org/10.5281/zenodo.16756123)
