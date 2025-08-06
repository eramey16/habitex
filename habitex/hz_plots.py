import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import (AnchoredOffsetbox, AuxTransformBox,
                                  DrawingArea, TextArea, VPacker)
import astropy
from matplotlib.patches import Circle, Ellipse, Annulus
import habitex
import hab_zone


class PlotHZ:
    """Class: PlotHZ
    Useful Plots to visualize planets in the habitable zone

    """
    def __init__(self):
        pass
         
    
    def plot_hab(self, hostname=None,sma=None, eccen=None, cons_in=None, cons_out=None, opt_in=None, opt_out=None):

        """Plot HZ
        Visual representation of the planet orbit and Habitable Zone around the star
        The conservative and optimistic habitable zone are plotted in concentric circles,
        while the planet orbit will be an ellipse depending on eccentricity

        Args:
            hostname (string)
            sma (float): Semi-major axis in AU
            eccen (float): eccentricity
            cons_in: Inner bound of conservative habitable zone in AU
            cons_out: Outer bound of conservative habitable zone in AU
            opt_in: Inner bound of optimistic habitable zone in AU
            opt_out: Outer bound of optimistic habitable zone in AU
        
        Returns:
            pyplot
        """
        exp = habitex.ArchiveExplorer()
        if hostname:
            tab = exp.query_exo(hostname=hostname)
        else:
            tab = exp.query_exo()

        if tab.empty:
            print("No matching exoplanet data found.")
            return

        eccen = tab["pl_orbeccen"].iloc[0]
        sma = tab["pl_orbsmax"].iloc[0]

        eval = hab_zone.HabZoneEvaluator()

        cons_data = eval.conservative_habzone()
        cons_in = cons_data['Conservative Inner Radius (AU)'].min()
        cons_out = cons_data['Conservative Outer Radius (AU)'].max()

        opt_data = eval.optimistic_habzone()
        opt_in = opt_data['Optimistic Inner Radius (AU)'].min()
        opt_out = opt_data['Optimistic Outer Radius (AU)'].max()

        cons_zone = Annulus((0, 0), cons_out, cons_out - cons_in, color='green', alpha=0.8, label="Conservative HZ")
        opt_zone = Annulus((0, 0), opt_out, opt_out - opt_in, color='green', alpha=0.4, label="Optimistic HZ")

        smin = np.sqrt(1 - (eccen)**2)*sma
        orbit = Ellipse((0, 0), sma, smin, angle=30, color='black', fill=False, label="Planet Orbit")

        fig, ax = plt.subplots()

        ax.add_patch(cons_zone)
        ax.add_patch(opt_zone)

        ax.add_patch(orbit)

        ax.set_xlabel("Distance (AU)")
        ax.set_ylabel("Distance (AU)")

        ax.set_aspect('equal')

        ax.text(0, 0, '★', fontsize=20, ha='center', va='center', color='gold', zorder=5)

        max_radius = max(cons_out, opt_out, sma * (1 + eccen))  # add margin for orbit
        ax.set_xlim(-1.2*max_radius, 1.2*max_radius)
        ax.set_ylim(-1.2*max_radius, 1.2*max_radius)
        ax.legend()

        planet_name = tab["pl_name"].iloc[0]

        plt.title(f"Habitable Zone and Orbit for {planet_name}")
        plt.show()

        return

    def plot_massradius_conservative(self, hostname=None):
        """Plot Mass-Radius Diagram 
        Plot the mass radius diagram for the conservative habitable zone 

        Args:
            hostname (string)
        
        Returns:
            pyplot
        """
        exp = habitex.ArchiveExplorer()
        if hostname:
            tab = exp.query_exo(hostname=hostname)
        else:
            tab = exp.query_exo()

        if tab.empty:
            print("No matching exoplanet data found.")
            return

        eccen = tab["pl_orbeccen"].iloc[0]
        sma = tab["pl_orbsmax"].iloc[0]

        eval = hab_zone.HabZoneEvaluator()
        table = eval.conservative_habzone()
        cons_table = table[table['In Conservative Habitable Zone'] == True]
        cons_mass = cons_table['pl_msinie']
        cons_radii = cons_table['pl_rade']
        cons_temp = cons_table['st_teff'].values

        plt.scatter(cons_mass,cons_radii,c=cons_temp,cmap='inferno')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Minimum Mass ($M_{\earth}$)')
        plt.ylabel('Planet Radius ($R_{\earth}$)')
        plt.colorbar(label='Host Star $T_{eff}$')
        plt.title('Mass-Radius Relation for Planets in the Conservative HZ')
        plt.show()

        return
    
    def plot_massradius_optimistic(self, hostname=None):
        """Plot Mass-Radius Diagram 
        Plot the mass radius diagram for the optimistic habitable zone 

        Args:
            hostname (string)
        
        Returns:
            pyplot
        """
        exp = habitex.ArchiveExplorer()
        if hostname:
            tab = exp.query_exo(hostname=hostname)
        else:
            tab = exp.query_exo()

        if tab.empty:
            print("No matching exoplanet data found.")
            return

        eval = hab_zone.HabZoneEvaluator()
        table2 = eval.optimistic_habzone()
        opt_table = table2[table2['In Optimistic Habitable Zone'] == True]
        opt_mass = opt_table['pl_msinie']
        opt_radii = opt_table['pl_rade']
        opt_temp = opt_table['st_teff'].values

        plt.scatter(opt_mass,opt_radii,c=opt_temp,cmap='inferno')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Minimum Mass ($M_{\earth}$)')
        plt.ylabel('Planet Radius ($R_{\earth}$)')
        plt.colorbar(label='Host Star $T_{eff}$')
        plt.title('Mass-Radius Relation for Planets in the Optimistic HZ')
        plt.show()

        return