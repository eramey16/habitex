import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import (AnchoredOffsetbox, AuxTransformBox,
                                  DrawingArea, TextArea, VPacker)
from matplotlib.patches import Circle, Ellipse, Annulus
from habitex import ArchiveExplorer


class PlotHZ:
    """Class: PlotHZ
    Useful Plots to visualize planets in the habitable zone

    """
    def __init__(self):
        pass
         
    
    def plot_hab(self, hostname=None, pl_name=None,sma=None, eccen=None, cons_in=None, cons_out=None, opt_in=None, opt_out=None):

        """Plot Habitable Zone
        Visual representation of the planet orbit and Habitable Zone around the star
        The conservative and optimistic habitable zone are plotted in concentric circles,
        while the planet orbit will be an ellipse depending on eccentricity

        Args:
            hostname: Name of host star
            pl_name: Most commonly used planet name
            sma (float): Semi-major axis in AU
            eccen (float): eccentricity
            cons_in: Inner bound of conservative habitable zone in AU
            cons_out: Outer bound of conservative habitable zone in AU
            opt_in: Inner bound of optimistic habitable zone in AU
            opt_out: Outer bound of optimistic habitable zone in AU
        
        Returns:
            matplotlib.pyplot
        """
        planets_data = []

        if hostname:
            exp = ArchiveExplorer()
            tab = exp.query_exo(hostname=hostname)

            if tab.empty:
                print("No matching exoplanet data found.")
                return

            if pl_name:
                tab = tab[tab["pl_name"] == pl_name]
                if tab.empty:
                    print(f"No planet named {pl_name} found for {hostname}.")
                    return

            for _, row in tab.reset_index(drop=True).iterrows():
                planets_data.append({
                    "name": row["pl_name"],
                    "sma": row["pl_orbsmax"],
                    "eccen": row["pl_orbeccen"],
                    "cons_in": row.get("hz_inner_cons"),
                    "cons_out": row.get("hz_outer_cons"),
                    "opt_in": row.get("hz_inner_opt"),
                    "opt_out": row.get("hz_outer_opt")
                })

        else:
            #for custom user data
            if None in (pl_name, sma, eccen):
                print("For custom data, at least pl_name, sma, and eccen must be provided.")
                return

            planets_data.append({
                "name": pl_name,
                "sma": sma,
                "eccen": eccen,
                "cons_in": cons_in,
                "cons_out": cons_out,
                "opt_in": opt_in,
                "opt_out": opt_out
            })

        #plot all planets
        for pdata in planets_data:
            #hz
            fig, ax = plt.subplots()

            if pdata["cons_in"] is not None and pdata["cons_out"] is not None:
                cons_zone = Annulus((0, 0), pdata["cons_out"],
                                    pdata["cons_out"] - pdata["cons_in"],
                                    color='green', alpha=0.8, label="Conservative HZ")
                ax.add_patch(cons_zone)

            if pdata["opt_in"] is not None and pdata["opt_out"] is not None:
                opt_zone = Annulus((0, 0), pdata["opt_out"],
                                pdata["opt_out"] - pdata["opt_in"],
                                color='green', alpha=0.4, label="Optimistic HZ")
                ax.add_patch(opt_zone)

            #orbit
            a = pdata["sma"]
            b = np.sqrt(1 - pdata["eccen"]**2) * a
            focus_offset = pdata["eccen"] * a
            orbit = Ellipse((-focus_offset, 0), 2 * a, 2 * b,
                            color='black', fill=False, label="Planet Orbit")
            ax.add_patch(orbit)

            #star
            ax.plot(0, 0, marker='*', markersize=10, color='gold', zorder=5, label="Host Star")

            ax.set_xlabel("Distance (AU)")
            ax.set_ylabel("Distance (AU)")
            ax.set_aspect('equal')

            #to see all content
            max_items = [a * (1 + pdata["eccen"])]
            if pdata["opt_out"] is not None:
                max_items.append(pdata["opt_out"])
            elif pdata["cons_out"] is not None:
                max_items.append(pdata["cons_out"])
            max_radius = max(max_items)

            ax.set_xlim(-1.2 * max_radius, 1.2 * max_radius)
            ax.set_ylim(-1.2 * max_radius, 1.2 * max_radius)
            ax.legend()

            plt.title(f"Habitable Zone and Orbit for {pdata['name']}")
            plt.show()
        return

    def plot_massradius_conservative(self):
        """Plot Mass-Radius Diagram 
        Plot the mass radius diagram for planets in the conservative habitable zone 

        Args:
            None
        
        Returns:
            matplotlib.pyplot
        """
        exp = ArchiveExplorer()
        tab = exp.query_exo()

        cons_table = tab[tab['in_hz_cons'] == True]
        cons_mass = cons_table['pl_msinie']
        cons_radii = cons_table['pl_rade']
        cons_temp = cons_table['st_teff'].values

        plt.scatter(cons_mass,cons_radii,c=cons_temp,cmap='inferno')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Minimum Mass (M$_{\oplus}$)')
        plt.ylabel('Planet Radius (R$_{\oplus}$)')
        plt.colorbar(label='Host Star T$_{eff}$')
        plt.title('Mass-Radius Relation for Planets in the Conservative HZ')
        plt.show()

        return
    
    def plot_massradius_optimistic(self):
        """Plot Mass-Radius Diagram 
        Plot the mass radius diagram for planets in the optimistic habitable zone 

        Args:
            None
        
        Returns:
            matplotlib.pyplot
        """
        exp = ArchiveExplorer()
        tab = exp.query_exo()
        
        opt_table = tab[tab['in_hz_opt'] == True]
        opt_mass = opt_table['pl_msinie']
        opt_radii = opt_table['pl_rade']
        opt_temp = opt_table['st_teff'].values

        plt.scatter(opt_mass,opt_radii,c=opt_temp,cmap='inferno')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Minimum Mass (M$_{\oplus}$)')
        plt.ylabel('Planet Radius (R$_{\oplus}$)')
        plt.colorbar(label='Host Star T$_{eff}$')
        plt.title('Mass-Radius Relation for Planets in the Optimistic HZ')
        plt.show()

        return
    
   