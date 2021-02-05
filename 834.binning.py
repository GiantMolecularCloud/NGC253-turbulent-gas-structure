##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data      scriptdir
dendrograms = load_dendrograms()


###################################################################################################
# binning function
###################################################################################################

def get_binned_percentiles(x_measure, y_measure, x_step, xdata, ydata):

    from astropy.table import Table

    def filter_window(x,y,window_min,window_max):
        """select values in the window"""
        x_window = x[(window_min<x) & (x<window_max)]
        y_window = y[(window_min<x) & (x<window_max)]
        return x_window,y_window


    # binned percentiles table
    bin_perc = Table(names = ['window min','window center','window max','x 16th','x median','x 84th','y 16th','y median','y 84th'],
                     meta = {'x': x_measure, 'y': y_measure, 'x_step': x_step}
                    )

    x_matched, y_matched = crossmatch(xdata,ydata)

    x_min = np.nanmin(x_matched)
    x_max = np.nanmax(x_matched)
    y_min = np.nanmin(y_matched)
    y_max = np.nanmax(y_matched)

    for log_window_cen in np.arange(np.log10(x_min), np.log10(x_max)+x_step, x_step):
        log_window_min = log_window_cen-x_step/2.
        log_window_max = log_window_cen+x_step/2.
        window_min = np.power(10,log_window_min)
        window_cen = np.power(10,log_window_cen)
        window_max = np.power(10,log_window_max)
        x_window,y_window = filter_window(x_matched,y_matched,window_min,window_max)

        # percentiles only make sense when at least two measurements are available
        if len(x_window)>1:
            x_perc = np.percentile(x_window, [16,50,84])
            y_perc = np.percentile(y_window, [16,50,84])
            bin_perc.add_row( flatten([window_min,window_cen,window_max,x_perc,y_perc]) )

    return bin_perc


###################################################################################################
# calculate binned size-linewidth relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            for lw_measure in lw_measures:

                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.1,
                                              xdata  = catalog[R_measure],
                                              ydata  = catalog[lw_measure]
                                             )

                fits[co][gal]['size-linewidth'][R_measure][lw_measure]['bins'] = bins


# leaves
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            for lw_measure in lw_measures:

                size      = catalog[R_measure]
                linewidth = catalog[lw_measure]
                filter    = (catalog['type']=='leaf')
                xdata = size[filter]
                ydata = linewidth[filter]

                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.1,
                                              xdata  = xdata,
                                              ydata  = ydata
                                             )

                fits[co][gal]['size-linewidth (leaves)'][R_measure][lw_measure]['bins'] = bins


###################################################################################################
# calculate binned density-virial relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['density-virial']
        except:
            fits[co][gal]['density-virial'] = {}

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            try:
                fits[co][gal]['density-virial'][R_measure]
            except:
                fits[co][gal]['density-virial'][R_measure] = {}

            for lw_measure in lw_measures:
                try:
                    fits[co][gal]['density-virial'][R_measure][lw_measure]
                except:
                    fits[co][gal]['density-virial'][R_measure][lw_measure] = {}

                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.25,
                                              xdata  = catalog['column density'],
                                              ydata  = catalog[lw_measure]**2/catalog[R_measure]
                                             )

                fits[co][gal]['density-virial'][R_measure][lw_measure]['bins'] = bins


# leaves
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['density-virial (leaves)']
        except:
            fits[co][gal]['density-virial (leaves)'] = {}

        catalog    = dendrograms[co][gal]['catalog']

        for R_measure in R_measures:
            try:
                fits[co][gal]['density-virial (leaves)'][R_measure]
            except:
                fits[co][gal]['density-virial (leaves)'][R_measure] = {}

            for lw_measure in lw_measures:
                try:
                    fits[co][gal]['density-virial (leaves)'][R_measure][lw_measure]
                except:
                    fits[co][gal]['density-virial (leaves)'][R_measure][lw_measure] = {}

                coldens   = catalog['column density']
                linewidth = catalog[lw_measure]
                size      = catalog[R_measure]
                filter = (catalog['type']=='leaf')

                xdata = coldens[filter]
                ydata = linewidth[filter]**2/size[filter]

                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.25,
                                              xdata  = xdata,
                                              ydata  = ydata
                                             )

                fits[co][gal]['density-virial (leaves)'][R_measure][lw_measure]['bins'] = bins


###################################################################################################
# calculate binned mass-linewidth relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        try:
            fits[co][gal]['mass-linewidth']
        except:
            fits[co][gal]['mass-linewidth'] = {}

        catalog    = dendrograms[co][gal]['catalog']
        bins = get_binned_percentiles('mass', 'linewidth (astrodendro)',
                                      x_step = 0.1,
                                      xdata  = catalog['mass'],
                                      ydata  = catalog['linewidth (astrodendro)']
                                     )
        fits[co][gal]['mass-linewidth']['bins'] = bins


###################################################################################################
# calculate binned size-mass relation
###################################################################################################

# all
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']

        bins = get_binned_percentiles('size (astrodendro)', 'mass',
                                      x_step = 0.1,
                                      xdata  = catalog['size (astrodendro)'],
                                      ydata  = catalog['mass']
                                     )
        fits[co][gal]['size-mass']['bins'] = bins



###################################################################################################
# save fits
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(fits[co][gal], join(compdir,gal+'.'+co+'.fits.pickle'))


###################################################################################################
#
###################################################################################################
