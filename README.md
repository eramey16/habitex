<p align="center">
  <img src="https://github.com/user-attachments/assets/0a2d6c00-8ac1-4c0d-aa27-b71e496c3e9f" width="200"/>
</p>

**habitex** is a Python-based tool designed to vet and characterize potentially habitable exoplanets using public data from the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/). The pipeline evaluates exoplanets based on their stellar and orbital properties to determine whether they reside in a **conservative** or **optimistic habitable zone**, and supports custom filtering for survey planning, target selection, and comparative exoplanetology.

---

## Objective

For each confirmed exoplanet, habitex:

- Retrieves key planetary and stellar parameters using `astroquery`
- Computes whether the planet lies within the conservative or optimistic habitable zone using models from Kopparapu et al. (2013, 2014)
- Estimates additional properties such as planet density (if mass and radius are available)
- Supports optional filtering by orbital, physical, or observational criteria
- Outputs structured data for follow-up analysis and survey planning

---

## Parameters Retrieved

habitex queries the following fields from the NASA Exoplanet Archive:

| Parameter       | Description                          |
|----------------|--------------------------------------|
| `pl_orbper`     | Orbital Period [days]                |
| `pl_orbsmax`    | Semi-Major Axis [AU]                 |
| `pl_masse`      | Planet Mass [Earth masses]           |
| `pl_msinie`     | Minimum Mass [Earth masses]          |
| `pl_rade`       | Planet Radius [Earth radii]          |
| `pl_eqt`        | Planet Equilibrium Temperature [K]   |
| `pl_orbeccen`   | Orbital Eccentricity                 |
| `st_teff`       | Stellar Effective Temperature [K]    |
| `hostname`      | Stellar Host Name                    |
| `dec`           | Declination [deg]                    |

---

## Functional Overview

### Data Retrieval

- Queries the NASA Exoplanet Archive using `astroquery`
- Defaults to the most recent planet entry if multiple are available
- Allows optional filtering by paper or table (`pscomppars` (default), `cumulative`)

### Habitable Zone Assessment

- Calculates incident stellar flux using either semi-major axis or orbital period (via Kepler’s Third Law)
- Evaluates whether the planet falls into the following categories, according to Kopparapu et al. 2013:
  - **Conservative Habitable Zone** (e.g. water loss to maximum greenhouse limits)
  - **Optimistic Habitable Zone** (e.g. recent Venus to early Mars)

### Planetary Density Estimation

- If both mass and radius are known, density is calculated
- Rocky planets may be flagged based on radius or density thresholds

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
  - Rocky planet likelihood
- Plots such as orbital diagrams or HZ placement (future versions)

---

## Dependencies

- `astroquery`
- `astropy`
- `pandas`
- `numpy`
- `matplotlib`
- `scipy`

---

## References

- Kopparapu et al. (2013), ApJ, 765, 131  
- Kopparapu et al. (2014), ApJ, 787, L29
- Luque et al. (2022) arXiv:2209.03871v1
- NASA Exoplanet Archive API: https://exoplanetarchive.ipac.caltech.edu/docs/program_interfaces.html  
