##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()
fits        = load_fits()


###################################################################################################
# plot fitted size-mass relation
###################################################################################################

# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-mass']['fit']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        size = catalog['size (astrodendro)']
        mass = catalog['mass']

        R_min = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max = np.nanmax(size)
        M_min = np.nanmin([x if x>1e-2 else np.nan for x in mass])
        m_max = np.nanmax(mass)

        ax.scatter(size, mass, marker='.', s=10, c=color, label=co+' '+gal, zorder=1)

        # fit
        x = np.logspace(np.log10(fits[co][gal]['size-mass']['fit range'][0]), np.log10(fits[co][gal]['size-mass']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['MCMC'][0],np.power(10,fit.c['MCMC'][0])), lw=0.5, color='k', zorder=3)

        # fit uncertainty
        samples = [powlaw(x,a,np.power(10,c)) for a,c in fit.samples]
        p16, p84 = np.percentile(samples, (16,84), axis=0)
        ax.fill_between(x, p16, p84, facecolor='grey', edgecolor=None, alpha=0.5, zorder=2)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('M [$_\odot$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*M_min, 2.0*M_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '07.size-mass', 'scatter_fit')
        fig.savefig(join(plotdir, '07.size-mass', 'scatter_fit', 'size-mass.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

R_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
R_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))
M_min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()])])
M_max = np.nanmax(flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-mass']['fit']

        mass = catalog['mass']
        size = catalog['size (astrodendro)']

        ax.scatter(size, mass, marker='.', s=10, c=color, label=co+' '+gal)
        x = np.logspace(np.log10(fits[co][gal]['size-mass']['fit range'][0]), np.log10(fits[co][gal]['size-mass']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['MCMC'][0],np.power(10,fit.c['MCMC'][0])), lw=0.5, color='k', zorder=3)

ax.set_xlabel('R [pc]')
ax.set_ylabel('M [$_\odot$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*M_min, 2.0*M_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '07.size-mass', 'scatter_fit')
fig.savefig(join(plotdir, '07.size-mass', 'scatter_fit', 'size-mass.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# plot binned size-mass relation
###################################################################################################

# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-mass']['binned fit']
        bins       = fits[co][gal]['size-mass']['bins']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        mass = catalog['mass']
        size = catalog['size (astrodendro)']

        R_min  = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max  = np.nanmax(size)
        M_min  = np.nanmin([x if x!=0. else np.nan for x in mass])
        M_max = np.nanmax(mass)

        i = 0
        for row in bins:
            window = [row['window min'], row['window max']]
            window_cen = row['window center']
            median = [row['y median'], row['y median']]
            perc16 = [row['y 16th'], row['y 16th']]
            perc84 = [row['y 84th'], row['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=2, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=1, label=co+' '+gal if i==0 else '')
            i += 1

        # fit
        x = np.logspace(np.log10(fits[co][gal]['size-mass']['fit range'][0]), np.log10(fits[co][gal]['size-mass']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

        # fit uncertainty
        low = powlaw(x,fit.a['ls'][0]-fit.a['ls'][1],np.power(10,fit.c['ls'][0]-fit.c['ls'][1]))
        up  = powlaw(x,fit.a['ls'][0]+fit.a['ls'][1],np.power(10,fit.c['ls'][0]+fit.c['ls'][1]))
        ax.fill_between(x, low, up, facecolor='grey', edgecolor=None, alpha=0.5, zorder=2)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('M [$_\odot$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*M_min, 2.0*M_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '07.size-mass', 'binned_fit')
        fig.savefig(join(plotdir, '07.size-mass', 'binned_fit', 'size-mass.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

R_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
R_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))
M_min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()])])
M_max = np.nanmax(flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-mass']['binned fit']
        bins       = fits[co][gal]['size-mass']['bins']

        mass = catalog['mass']
        size = catalog['size (astrodendro)']

        i = 0
        for row in bins:
            window = [row['window min'], row['window max']]
            window_cen = row['window center']
            median = [row['y median'], row['y median']]
            perc16 = [row['y 16th'], row['y 16th']]
            perc84 = [row['y 84th'], row['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=2, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=1, label=co+' '+gal if i==0 else '')
            i += 1

        x = np.logspace(np.log10(fits[co][gal]['size-mass']['fit range'][0]), np.log10(fits[co][gal]['size-mass']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

ax.set_xlabel('R [pc]')
ax.set_ylabel('M [$_\odot$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*M_min, 2.0*M_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '07.size-mass', 'binned_fit')
fig.savefig(join(plotdir, '07.size-mass', 'binned_fit', 'size-mass.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
