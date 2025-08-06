import numpy as np
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

class ArchiveExplorer:

    cols = ['gaia_id', 'pl_pubdate', 'pl_name', 'hostname','dec', 'st_mass', 'pl_orbper','pl_orbsmax', 'pl_masse', 
            'pl_msinie','pl_rade','st_teff', 'pl_eqt','pl_orbeccen', 'pl_dens', 'st_lum']
    G = 6.743e-11 # m^3 kg^-1 s^-2
    m_earth = 5.9722e24 # mass of earth in kg
    m_sun = 1.989e30 # mass of sun in kg
    au = 1.496e11 # 1 AU in m
    day = 60 * 60 * 24 # day in seconds

    def __init__(self):
        pass

    def classify_planet_by_density(self, density_ratio):
#         if np.isnan(density_ration):
#           return None
        if density_ratio < 0.4:
            return "Gas planets"
        elif density_ratio < 0.7:
            return "Water worlds"
        else:
            return "Rocky planets"
    

    def query_exo(self, table='pscomppars', hostname=None, t_eff=None, dec=None, 
                 period=None, mandr=False):
        
        # Add default cuts, unless user specified
        _range = lambda param, minmax: f"{param}>{minmax[0]} and {param}<{minmax[1]}"

        # Cut on eccentricity (important for the equations)
        cuts = ["pl_orbeccen<0.3"]

        # Other cuts
        if mandr: cuts.append("pl_masse is not null and pl_rade is not null")
        if hostname is not None: cuts.append(f"hostname like '{hostname}'")
        if t_eff is not None: cuts.append(_range('st_teff', t_eff))
        if dec is not None: cuts.append(_range('dec', dec))
        if period is not None: cuts.append(_range('pl_orbper', period))

        # Query exoplanet archive
        tab = NasaExoplanetArchive.query_criteria(table=table, 
                                                  select=', '.join(self.cols),
                                                  where=' and '.join(cuts)
                                                  ).to_pandas()
        
        tab['pl_type'] = tab['pl_dens'].apply(lambda x: self.classify_planet_by_density(x) 
                                              if pd.notnull(x) else None)
        
        
        # Calculate orbital distance and add to table
        tab['pl_orbdist'] = self._orb_dist(tab)
        tab['pl_type'] = tab['pl_dens'].apply(lambda x: self.classify_planet_by_density(x) 
                                                  if pd.notnull(x) else None)

        # Drop duplicates (last first) if the table includes them
        if table!='pscomppars':
            tab.sort_values(by='pl_pubdate', ascending=False, ignore_index=True, inplace=True)
            tab.drop_duplicates(subset=['gaia_id', 'pl_name'], keep='first', inplace=True, ignore_index=True)

        self.results = tab
        return tab
    
    def _orb_dist(self, data):
        """ Calculates orbital distance from orbital period """
        r = np.cbrt((self.G * data.st_mass * self.m_sun / (4 * np.pi**2)) 
                    * (data.pl_orbper * self.day)**2)
        return r / self.au # orbital distance in AU