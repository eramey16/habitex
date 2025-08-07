import numpy as np
import pandas as pd
from astroquery.ipac.nexsci.nasa_exoplanet_archive import NasaExoplanetArchive

class ArchiveExplorer:

    """ ArchiveExplorer class - Handles queries between the user and NASA's Exoplanet Archive
    """

    cols = ['gaia_id', 'ra', 'dec', 'pl_pubdate', 'pl_name', 'hostname', 'st_mass', 'st_teff', 'st_lum',
             'pl_orbper','pl_orbsmax', 'pl_masse', 'pl_msinie','pl_rade', 'pl_eqt','pl_orbeccen', 'pl_dens']
    G = 6.743e-11 # m^3 kg^-1 s^-2
    m_earth = 5.9722e24 # mass of earth in kg
    m_sun = 1.989e30 # mass of sun in kg
    au = 1.496e11 # 1 AU in m
    day = 60 * 60 * 24 # day in seconds

    cons_params = {
        'sol_flux_in': 1.0512, 'sol_flux_out': 0.3438,
        'a_in': 1.3242e-4, 'a_out': 5.8942e-5,
        'b_in': 1.5418e-8, 'b_out': 1.6558e-9,
        'c_in': -7.9895e-12, 'c_out': -3.0045e-12,
        'd_in': -1.8328e-15, 'd_out': -5.2983e-16,}
    
    opt_params = {
        'sol_flux_in': 1.7753, 'sol_flux_out': 0.3179,
        'a_in': 1.4316e-4, 'a_out': 5.4513e-5,
        'b_in': 2.9875e-9, 'b_out': 1.5313e-9,
        'c_in': -7.5702e-12, 'c_out': -2.7786e-12,
        'd_in': -1.1635e-15, 'd_out': -4.8997e-16,
    }

    def __init__(self, optimistic=False):
        pass

    def query_exo(self, table='pscomppars', hostname=None, t_eff=None, dec=None, 
                 period=None, mandr=False, paper=None, cols=None):
        """ Queries the NASA Exoplanet archive

        Calculates orbital distance and planet density and adds them to query results

        Args:
            table (str, optional): Table to pull from (typically 'ps' or 'pscomppars')
            hostname (str/list, optional): Specify hostname of a star or set of stars (e.g. 'Kepler-%')
            t_eff (tuple, optional): Range of effective temperatures [lo, hi]
            dec (tuple, optional): Declination range to search [lo, hi]
            period (tuple, optional): Period range to search [lo, hi]
            mandr (bool, optional): Specifies that both mass and radius must be non-null
            paper (str/list, optional): Reference name for discovery publication, only used with table 'ps' 
                (formatted as 'Author et al. Year')
            optimistic (bool, optional): Whether to use an optimistic habitable zone (False for conservative)
            cols (list, optional): List of additional column names as string. Default values:
                | *gaia_id*: Gaia ID of the star
                | *ra*: Right Ascension (star)
                | *dec*: Declination (star)
                | *pl_pubdate*: Initial date of publication of the planet's data
                | *pl_name*: Most commonly used planet name
                | *hostname*: Name of host star
                | *st_mass*: Mass of host star (in solar masses)
                | *st_teff*: Effective temperature of host star (in Kelvin)
                | *st_lum*: Log luminosity of host star (in log10(Solar))
                | *pl_orbper*: Orbital period of planet
                | *pl_orbsmax*: Orbital distance of planet
                | *pl_masse*: Mass of planet (in Earth masses)
                | *pl_msinie*: Minimum mass (in Earth masses)
                | *pl_rade*: Radius of planet (in Earth radii)
                | *pl_eqt*: Equilibrium temperature of the planet (in Kelvin)
                | *pl_orbeccen*: Orbital eccentricity of the planet
                | *pl_dens*: Density of planet (in g/cm^3)
        
        Returns:
            results (pd.DataFrame): Results of query as a pandas dataframe.
            Orbital distance will be a new column *pl_orbdist* (in AU), 
            planet density classification will be in column *pl_type* (as string)

            Habitable zone variables will be in 'hz_inner_<opt/cons>', 
            'hz_outer_<opt/cons>', and 'in_hz_<opt/cons>' columns
        """
        
        # Add default cuts, unless user specified
        _range = lambda param, minmax: f"{param}>{minmax[0]} and {param}<{minmax[1]}"

        # Cut on eccentricity (important for the equations)
        cuts = ["pl_orbeccen<0.3"]
        if cols is not None: [self.cols.append(col) for col in cols if col not in self.cols]

        # Other cuts
        if mandr: cuts.append("pl_masse is not null and pl_rade is not null")
        if hostname is not None:
                if isinstance(hostname, list):
                    host_cuts = " or ".join([f"hostname like '{h}'" for h in hostname])
                    cuts.append(f"({host_cuts})")
                else:
                    cuts.append(f"hostname like '{hostname}'")
        if t_eff is not None: cuts.append(_range('st_teff', t_eff))
        if dec is not None: cuts.append(_range('dec', dec))
        if period is not None: cuts.append(_range('pl_orbper', period))
        if paper is not None and table=='ps': cuts.append(f"disc_refname like '%{paper}%'")

        # Query exoplanet archive
        tab = NasaExoplanetArchive.query_criteria(table=table, 
                                                  select=', '.join(self.cols),
                                                  where=' and '.join(cuts)
                                                  ).to_pandas()
        
        
        # Drop duplicates (last first) if the table includes them
        if table!='pscomppars':
            tab.sort_values(by='pl_pubdate', ascending=False, ignore_index=True, inplace=True)
            tab.drop_duplicates(subset=['gaia_id', 'pl_name'], keep='first', inplace=True, ignore_index=True)
        
        new_data = self.calc_exo(tab)

        tab = tab.join(new_data)

        tab.reset_index(inplace=True, drop=True)
        self.results = tab
        return self.results
    
    def calc_exo(self, pl_data):
        """ Calculates exoplanet parameters based on user-input data

        Args:
            pl_data (pd.DataFrame): Table of planetary data (TODO - list required columns)
            optimistic (bool, optional): Whether to use an optimistic habitable zone (False for conservative)
        """
        # Calculate orbital distance and add to table
        new_data = pd.DataFrame()
        new_data['pl_orbdist'] = self._orb_dist(pl_data)
        new_data['pl_type'] = pl_data['pl_dens'].apply(lambda x: self._classify_planet_by_density(x) 
                                                  if pd.notnull(x) else None)

        # Calculate the habitable zone for conservative and optimistic and add to the table
        hz_data_opt = self._hab_zone(pl_data.join(new_data), optimistic=True)
        hz_data_cons = self._hab_zone(pl_data.join(new_data), optimistic=False)

        new_data = new_data.join(hz_data_opt).join(hz_data_cons)

        return new_data
    
    def _classify_planet_by_density(self, density_ratio):
        """ Classifies a planet's type based on its density

        Args:
            density_ratio (array-like): A pandas series or numpy array of density ratios (float)
        
        Returns:
            String classification of density (Gas planet, water, world, or Rocky planet)
        """
        if density_ratio < 0.4:
            return "Gas"
        elif density_ratio < 0.7:
            return "Water"
        else:
            return "Rocky"
    
    def _hab_zone(self, data, optimistic=False):
        """Calculates exoplanet habitable zone based on user-input dataframe

        The conservative habitable zone is given by the runaway greenhouse and 
        maximum greenhouse limits in Kopparapu et al. 2013.
        The optimistic habitable zone is given by the Venus and early Mars limits.

        Args:
            data (pd.DataFrame): A dataframe with planetary parameters
            optimistic (bool): Whether to use an optimistic habitable zone calculation (default False)
        
        Returns: pd.DataFrame
            'hz_inner' (AU), 'hz_outer' (AU), and 'in_hz' columns indicating whether 
            the planet is in the habitable zone
        """

        params = self.opt_params if optimistic else self.cons_params
        
        semimajor = data.pl_orbdist.astype(float)
        t_star = data.st_teff.astype(float) - 5780
        pl_stflux = (10**data.st_lum.astype(float))/(semimajor**2) #Stellar luminosity is in units of log(L/L_sun)
        inner_stflux = (params['sol_flux_in'] + 
                        params['a_in'] * t_star + 
                        params['b_in'] * (t_star**2) +
                        params['c_in'] * (t_star**3) +
                        params['d_in'] * (t_star**4)
                        ) / np.sqrt(1 - data.pl_orbeccen.astype(float)**2)
        outer_stflux = (params['sol_flux_out'] + 
                        params['a_out'] * t_star + 
                        params['b_out'] * (t_star**2) +
                        params['c_out'] * (t_star**3) +
                        params['d_out'] * (t_star**4)
                        ) / np.sqrt(1 - data.pl_orbeccen.astype(float)**2)
        inner_rad = np.sqrt((10**data.st_lum)/inner_stflux)
        outer_rad = np.sqrt((10**data.st_lum)/outer_stflux)

        new_data = pd.DataFrame()
        tag = '_opt' if optimistic else '_cons'
        new_data['hz_inner'+tag] = inner_rad
        new_data['hz_outer'+tag] = outer_rad
        new_data['in_hz'+tag] = (np.array(pl_stflux > outer_stflux) & 
                         np.array(pl_stflux < inner_stflux))
        return new_data
    
    def _orb_dist(self, data):
        """ Calculates orbital distance from orbital period 
        Args:
            data: A pandas dataframe obtained from query_expo
        Returns:
            A pandas series of orbital distance in AU
        """
        r = np.cbrt((self.G * data.st_mass * self.m_sun / (4 * np.pi**2)) 
                    * (data.pl_orbper * self.day)**2)
        return r / self.au # orbital distance in AU
