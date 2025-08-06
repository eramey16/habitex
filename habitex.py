import numpy as np
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

class ArchiveExplorer:

    cols = ['hostname','dec', 'st_mass', 'pl_orbper','pl_orbsmax', 'pl_masse', 
            'pl_msinie','pl_rade','st_teff', 'pl_eqt','pl_orbeccen', 'pl_dens', 'pl_name']
    G = 6.743e-11 # m^3 kg^-1 s^-2
    m_earth = 5.9722e24 # mass of earth in kg
    m_sun = 1.989e30 # mass of sun in kg
    au = 1.496e11 # 1 AU in m
    day = 60 * 60 * 24 # day in seconds

    def __init__(self):
        pass

    def classify_planet_by_density(self, density_ratio):
        if density_ratio < 0.4:
            return "Gas planets"
        elif density_ratio < 0.7:
            return "Water worlds"
        else:
            return "Rocky planets"
    

    def query_exo(self, table='pscomppars', hostname=None, t_eff=None, dec=None, 
                 period=None, mandr=None, pl_dens=None, pl_name=None, pl_type=None):
        cuts = ["pl_orbeccen<0.3"]
        if mandr:
            cuts.append("pl_masse is not null and pl_rade is not null")
        if t_eff is not None:
            cuts.append(f"st_teff>{t_eff[0]} and st_teff<{t_eff[1]}")
        if hostname is not None:
            cuts.append(f"hostname like {hostname}")
        if dec is not None:
            cuts.append(f"dec>{dec[0]} and dec<{dec[1]}")
        if period is not None:
            cuts.append(f"pl_orbper>{period[0]} and pl_orbper<{period[1]}")
        if pl_dens is not None:
            cuts.append(f"pl_dens>{pl_dens[0]} and pl_dens<{pl_dens[1]}")
        if pl_name is not None:
            cuts.append(f"planetname like {pl_name}")
    
        tab = NasaExoplanetArchive.query_criteria(table="pscomppars", 
                                                  select=', '.join(self.cols),
                                                  where=' and '.join(cuts)
                                                  ).to_pandas()
        
        if 'pl_dens' in tab.columns:
            tab['pl_type'] = tab['pl_dens'].apply(lambda x: self.classify_planet_by_density(x) 
                                                  if pd.notnull(x) else None)
        else:
            tab['pl_type'] = None

        tab['pl_orbdist'] = self._orb_dist(tab)
        
        self.results = tab
        return tab
    
    def _orb_dist(self, data):
        """ Calculates orbital distance from orbital period """
        r = np.cbrt((self.G * data.st_mass * self.m_sun / (4 * np.pi**2)) 
                    * (data.pl_orbper * self.day)**2)
        return r / self.au # orbital distance in AU