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

    def __init__(self):
        pass

    def _classify_planet_by_density(self, density_ratio):
        """ Classifies a planet's type based on its density

        Args:
            density_ratio: A pandas series or numpy array of density ratios (float)
        
        Returns:
            String classificatoin of density (Gas planet, water, world, or Rocky planet)
        """
#         if np.isnan(density_ration):
#           return None
        if density_ratio < 0.4:
            return "Gas"
        elif density_ratio < 0.7:
            return "Water"
        else:
            return "Rocky"
    

    def query_exo(self, table='pscomppars', hostname=None, t_eff=None, dec=None, 
                 period=None, mandr=False, cols=None):
        """ Queries the NASA Exoplanet archive

        Calculates orbital distance and planet density and adds them to query results

        Args:
            table (optional): Table to pull from (typically 'ps' or 'pscomppars')
            hostname (optional): Specify hostname of a star or set of stars (e.g. 'Kepler-%')
            t_eff (optional): Range of effective temperatures [lo, hi]
            dec (optional): Declination range to search [lo, hi]
            period (optional): Period range to search [lo, hi]
            mandr (optional): Specifies that both mass and radius must be non-null
            cols: (optional): List of additional column names as string
            Default columns:
                gaia_id: Gaia ID of the star
                ra: Right Ascension (star)
                dec: Declination (star)
                pl_pubdate: Initial date of publication of the planet's data
                pl_name: Most commonly used planet name
                hostname: Name of host star
                st_mass: Mass of host star (in solar masses)
                st_teff: Effective temperature of host star (in Kelvin)
                st_lum: Log luminosity of host star (in log10(Solar))
                pl_orbper: Orbital period of planet
                pl_orbsmax: Orbital distance of planet
                pl_masse: Mass of planet (in Earth masses)
                pl_msinie: Minimum mass (in Earth masses)
                pl_rade: Radius of planet (in Earth radii)
                pl_eqt: Equilibrium temperature of the planet (in Kelvin)
                pl_orbeccen: Orbital eccentricity of the planet
                pl_dens: Density of planet (in g/cm^3)
        
        Returns:
            results: Results of query as a pandas dataframe.
            Orbital distance will be a new column 'pl_orbdist' (in AU), 
            and planet density classification will be in 'pl_type' (string)
        """
        
        # Add default cuts, unless user specified
        _range = lambda param, minmax: f"{param}>{minmax[0]} and {param}<{minmax[1]}"

        # Cut on eccentricity (important for the equations)
        cuts = ["pl_orbeccen<0.3"]
        if cols is not None: [self.cols.append(col) for col in cols if col not in self.cols]

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
        
        
        # Calculate orbital distance and add to table
        tab['pl_orbdist'] = self._orb_dist(tab)
        tab['pl_type'] = tab['pl_dens'].apply(lambda x: self._classify_planet_by_density(x) 
                                                  if pd.notnull(x) else None)

        # Drop duplicates (last first) if the table includes them
        if table!='pscomppars':
            tab.sort_values(by='pl_pubdate', ascending=False, ignore_index=True, inplace=True)
            tab.drop_duplicates(subset=['gaia_id', 'pl_name'], keep='first', inplace=True, ignore_index=True)

        self.results = tab
        return tab
    
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
