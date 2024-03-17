# import packages necessary for Task 3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs # for map projections
import cartopy.feature as cfeature # for map features like coastline, etc.
from cartopy.util import add_cyclic_point # get rid of the white line at 0 degree longitude for some models
import cmocean # a nice collection of colormap


# Define function to perform task 3
def season_avg_mld(mld, plot_year):

    '''
    Function for plotting seasonal mixed layer depth in the Southern Ocean
    and Antarctic region. The function provides a grid of subplots in which
    all three seasonally averaged datasets are plotted with the same color
    map, so that they can be compared with visual and analytic ease. Note,
    these subplots are not a return, but will automatically be displayed in
    the output when this function is called in a Jupyter notebook. The actual
    output of this function is a dictionary containing all of the seasonal
    average data, which is stored in xarray data array.
    ---
    Input:
    mld (xarray.core.dataarray.DataArray): An xarray data array of mixed-layer
                                           depth, created as the output of the
                                            cal_mld function.
    plot_year (int): The target year to be analyzed.
    ---
    Return:
    all_season_mld (dict): A dictionary with one column being the names of each
                           season and the other being the xarray data arrays for 
                           each seasonal mld averages.
    '''

    # Define each season as a set of the three months that make it up
    all_months = list(np.arange(12) + 1)
    JJA_months = all_months[5:8]
    SON_months = all_months[8:11]
    DJF_months = all_months[11:] + all_months[:2]
    MAM_months = all_months[2:5]
    all_seasons = [JJA_months, SON_months, DJF_months, MAM_months]
    season_names = ['JJA', 'SON', 'DJF', 'MAM']

    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), subplot_kw={'projection': ccrs.SouthPolarStereo()})
    axes = axes.flatten()

    # Define empty list to store all seasonal averages
    seasonal_avgs = []

    # Iterate over each season to calculate the seasonal average mld and plot
    for season, season_name, ax in zip(all_seasons, season_names, axes):

        # Define variable in which to store aggregate datasets for the season
        season_store = 0
        # Iterate over each month in the season to calculate mld
        for month in season:
            
            # Perform mld calculations for each individual month as defined by mld function
            plot_mld = mld.isel(time=((mld.time.dt.year == plot_year) & (mld.time.dt.month == month)))
            plot_mld = plot_mld.isel(time=0)
            
            # Add each month's mld data to the aggregate dataset
            season_store += plot_mld

        # Average the values over the three months in each season
        season_avg = season_store / 3
        season_avg.load()
        
        # Append seasonal average to list
        seasonal_avgs.append(season_avg)

        # Add geographic features
        ax.add_feature(cfeature.LAND, zorder=1)
        ax.add_feature(cfeature.COASTLINE, linewidth=0.5)

        # Set grid
        gl = ax.gridlines(
            xlocs = np.arange(-180,180,30),
            ylocs = np.arange(-90,90,20))

        # Set extent
        ax.set_extent([-180, 180, -90, -40], ccrs.PlateCarree())

        # Define color mesh
        im = ax.pcolormesh(season_avg.lon, season_avg.lat, season_avg,
                        transform=ccrs.PlateCarree(),
                        cmap=cmocean.cm.deep)

        ax.set_title('Mixed-Layer Depth [m] - Year: {}; Season: {}'.format(plot_year, season_name))

    # Combine all seasonal averages into a single array
    combined_seasonal_avgs = np.stack(seasonal_avgs)

    # Get the minimum and maximum values from all seasonal averages, excluding NaN values
    vmin = np.nanmin(combined_seasonal_avgs)
    vmax = np.nanmax(combined_seasonal_avgs)

    # Create a single colorbar for all subplots
    cbar_ax = fig.add_axes([0.1, 0.05, 0.8, 0.02])  # [left, bottom, width, height]
    cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal', extend='both')
    cbar.set_label('Mixed-Layer Depth [m]')
    cbar.set_ticks(np.linspace(vmin, vmax, num=10))

    # Adjust layout
    plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.show()

    # Create a dictionary in which to store the seasonal mld data and return it as the output of the function
    all_seasons_mld = {'Season' : season_names,
                       'Avg MLD' : seasonal_avgs}
    return all_seasons_mld