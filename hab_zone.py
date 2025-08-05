import numpy as np
import matplotlib.pyplot as plt
import astropy

class HabZoneEvaluator:

    def __init__(self):
        pass

    def conservative_habzone():
        pass
import numpy as np
import matplotlib.pyplot as plt
import astropy
import habitex

class HabZoneEvaluator:
    def __init__(self):
        pass

    def conservative_habzone(self, hostname=None, t_eff=None, dec=None, period=None, mandr=None):
        #Inner radius = runaway greenhouse (per Kopparapu et al. 2013)
        #Outer radius = maximum greenhouse
        sol_flux_in = 1.0512
        a_in = 1.3242e-4
        b_in = 1.5418e-8
        c_in = -7.9895e-12
        d_in = -1.8328e-15

        sol_flux_out = 0.3438
        a_out = 5.8942e-5
        b_out = 1.6558e-9
        c_out = -3.0045e-12
        d_out = -5.2983e-16
        data = habitex.ArchiveExplorer.query_exo(hostname=hostname, t_eff=t_eff, dec=dec, period=period, mandr=mandr)
        data['Conservative Inner Radius (AU)'] = np.full(data['hostname'], np.nan, dtype=float)
        data['Conservative Outer Radius (AU)'] = np.full(data['hostname'], np.nan, dtype=float)
        data['In Conservative Habitable Zone'] = np.full(data['hostname'], False, dtype=bool)
        for index, row in data.iterrows():
            semimajor = habitex._orb_dist(row)
            t_star = row['st_teff'] - 5780
            pl_stflux = (10**row['st_lum'])/semimajor**2 #Stellar luminosity is in units of log(L/L_sun)
            cons_inner_stflux = (sol_flux_in + a_in*t_star + b_in*(t_star**2) + c_in*(t_star**3) + d_in*(t_star**4))/np.sqrt(1 - row['pl_orbeccen']**2)
            cons_outer_stflux = (sol_flux_out + a_out*t_star + b_out*(t_star**2) + c_out*(t_star**3) + d_out*(t_star**4))/np.sqrt(1 - row['pl_orbeccen']**2)
            cons_inner_rad = np.sqrt((10**row['st_lum'])/cons_inner_stflux)
            cons_outer_rad = np.sqrt((10**row['st_lum'])/cons_outer_stflux)
            row['Conservative Inner Radius (AU)'] = cons_inner_rad
            row['Conservative Outer Radius (AU)'] = cons_outer_rad
            if pl_stflux > cons_outer_stflux and pl_stflux < cons_inner_stflux:
                row['In Conservative Habitable Zone'] = True
            else:
                row['In Conservative Habitable Zone'] = False

    def optimistic_habzone():
        pass