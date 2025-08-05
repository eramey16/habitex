import numpy as np
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

class ArchiveExplorer:

    cols = ['hostname','dec','pl_orbper','pl_orbsmax', 'pl_masse', 
            'pl_msinie','pl_rade','st_teff', 'pl_eqt','pl_orbeccen']
    G = 6.743e-11 # m^3 kg^-1 s^-2

    def __init__(self):
        pass

    def query_exo(self, table='pscomppars', hostname=None, t_eff=None, dec=None, 
                 period=None, mandr=None):
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
        tab = NasaExoplanetArchive.query_criteria(table="pscomppars", 
                                                  select=', '.join(self.cols),
                                                  where=' and '.join(cuts)
                                                  ).to_pandas()
        self.results = tab
        return tab
    
    def _orbital_distance(self):
        """ Calculates orbital distance from orbital period """

        r = np.cbrt((self.G * self.results.pl_masse / 4 / np.pi**2) 
                    * (self.results.pl_orbper**2))
        self.results['pl_orbdist'] = r
        return r