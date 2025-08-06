import numpy as np
import matplotlib.pyplot as plt
import astropy
import archive_explorer

class HabZoneEvaluator:
    """Class for Habitable Zone Calculations
    
    Contains functions to calculate whether the planet is in the habitable zone, per several different models of the habitable zone.

    """
    def __init__(self):
        pass

    def conservative_habzone(self, hostname=None, t_eff=None, dec=None, period=None, mandr=None):
        """Conservative Habitable Zone

        Evaluation of whether a given planet is in the conservative habitable zone, as given by the runaway greenhouse and maximum greenhouse limits in Kopparapu et al. 2013

        Args:
            hostname (tuple): Names of the host stars you are interested in (as strings with single quotations)
            t_eff (tuple): Tuple of minimum and maximum stellar effective temperatures you are interested in
            dec (tuple): Tuple of minimum and maximum declinations you are interested in
            period (tuple): Tuple of minimum and maximum planet orbital periods you are interested in (in days)
            mandr (Boolean, default=False): Require that there be measured mass and radius values

        Returns:
            DataFrame: Table of the planet systems and properties produced by ArchiveExplorer.query_exo(), with additional columns for the inner and outer radius of the conservative habitable zone and whether the planet is in the conservative habitable zone.
        """
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
        explorer = archive_explorer.ArchiveExplorer()
        data = explorer.query_exo(hostname=hostname, t_eff=t_eff, dec=dec, period=period, mandr=mandr)
        inner_rads = []
        outer_rads = []
        in_hz = []
        for index, row in data.iterrows():
            semimajor = explorer._orb_dist(row)
            t_star = float(row['st_teff']) - 5780
            pl_stflux = float((10**row['st_lum'])/semimajor**2) #Stellar luminosity is in units of log(L/L_sun)
            cons_inner_stflux = float((sol_flux_in + a_in*t_star + b_in*(t_star**2) + c_in*(t_star**3) + d_in*(t_star**4))/np.sqrt(1 - float(row['pl_orbeccen'])**2))
            cons_outer_stflux = float((sol_flux_out + a_out*t_star + b_out*(t_star**2) + c_out*(t_star**3) + d_out*(t_star**4))/np.sqrt(1 - float(row['pl_orbeccen'])**2))
            cons_inner_rad = np.sqrt(float((10**row['st_lum'])/cons_inner_stflux))
            cons_outer_rad = np.sqrt(float((10**row['st_lum'])/cons_outer_stflux))
            inner_rads.append(cons_inner_rad)
            outer_rads.append(cons_outer_rad)
            if pl_stflux > cons_outer_stflux and pl_stflux < cons_inner_stflux:
                in_hz.append(True)
            else:
                in_hz.append(False)
        data['Conservative Inner Radius (AU)'] = inner_rads
        data['Conservative Outer Radius (AU)'] = outer_rads
        data['In Conservative Habitable Zone'] = in_hz
        return data

    def optimistic_habzone(self, hostname=None, t_eff=None, dec=None, period=None, mandr=None):
        """Optimistic Habitable Zone

        Evaluation of whether a given planet is in the optimistic habitable zone, as given by the recent Venus and early Mars limits in Kopparapu et al. 2013

        Args:
            hostname (tuple): Names of the host stars you are interested in (as strings)
            t_eff (tuple): Tuple of minimum and maximum stellar effective temperatures you are interested in
            dec (tuple): Tuple of minimum and maximum declinations you are interested in
            period (tuple): Tuple of minimum and maximum planet orbital periods you are interested in (in days)
            mandr (Boolean, default=False): Require that there be measured mass and radius values

        Returns:
            DataFrame: Table of the planet systems and properties produced by ArchiveExplorer.query_exo(), with additional columns for the inner and outer radius of the conservative habitable zone and whether the planet is in the optimistic habitable zone.
        """
        #Inner radius = recent Venus (per Kopparapu et al. 2013)
        #Outer radius = early Mars
        sol_flux_in = 1.7753
        a_in = 1.4316e-4
        b_in = 2.9875e-9
        c_in = -7.5702e-12
        d_in = -1.1635e-15

        sol_flux_out = 0.3179
        a_out = 5.4513e-5
        b_out = 1.5313e-9
        c_out = -2.7786e-12
        d_out = -4.8997e-16
        explorer = archive_explorer.ArchiveExplorer()
        data = explorer.query_exo(hostname=hostname, t_eff=t_eff, dec=dec, period=period, mandr=mandr)
        inner_rads = []
        outer_rads = []
        in_hz = []
        for index, row in data.iterrows():
            semimajor = explorer._orb_dist(row)
            t_star = float(row['st_teff'])- 5780
            pl_stflux = float((10**row['st_lum'])/semimajor**2) #Stellar luminosity is in units of log(L/L_sun)
            opt_inner_stflux = float((sol_flux_in + a_in*t_star + b_in*(t_star**2) + c_in*(t_star**3) + d_in*(t_star**4))/np.sqrt(1 - float(row['pl_orbeccen'])**2))
            opt_outer_stflux = float((sol_flux_out + a_out*t_star + b_out*(t_star**2) + c_out*(t_star**3) + d_out*(t_star**4))/np.sqrt(1 - float(row['pl_orbeccen'])**2))
            opt_inner_rad = np.sqrt(float((10**row['st_lum'])/opt_inner_stflux))
            opt_outer_rad = np.sqrt(float((10**row['st_lum'])/opt_outer_stflux))
            inner_rads.append(opt_inner_rad)
            outer_rads.append(opt_outer_rad)
            if pl_stflux > opt_outer_stflux and pl_stflux < opt_inner_stflux:
                in_hz.append(True)
            else:
                in_hz.append(False)
        data['Optimistic Inner Radius (AU)'] = inner_rads
        data['Optimistic Outer Radius (AU)'] = outer_rads
        data['In Optimistic Habitable Zone'] = in_hz
        return data