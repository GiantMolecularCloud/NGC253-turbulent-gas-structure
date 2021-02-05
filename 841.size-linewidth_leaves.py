scriptdir##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()


###################################################################################################
# plot fitted size-linewidth relation
###################################################################################################

# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['fit']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        leaves    = (catalog['type']=='leaf')
        size      = catalog['size (astrodendro)'][leaves]
        linewidth = catalog['linewidth (astrodendro)'][leaves]

        R_min  = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max  = np.nanmax(size)
        lw_min = np.nanmin([x if x>1e-2 else np.nan for x in linewidth])
        lw_max = np.nanmax(linewidth)

        ax.scatter(size, linewidth, marker='.', s=10, c=color, label=co+' '+gal, zorder=1)

        # fit
        x = np.logspace(np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][0]), np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['MCMC'][0],np.power(10,fit.c['MCMC'][0])), lw=0.5, color='k', zorder=3)

        # fit uncertainty
        samples = [powlaw(x,a,np.power(10,c)) for a,c in fit.samples]
        p16, p84 = np.percentile(samples, (16,84), axis=0)
        ax.fill_between(x, p16, p84, facecolor='grey', edgecolor=None, alpha=0.5, zorder=2)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*lw_min, 2.0*lw_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '04.size-linewidth', 'scatter_fit')
        fig.savefig(join(plotdir, '04.size-linewidth', 'scatter_fit', 'size-linewidth_leaves.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

R_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
R_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))
lw_min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
lw_max = np.nanmax(flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['fit']

        leaves    = (catalog['type']=='leaf')
        size      = catalog['size (astrodendro)'][leaves]
        linewidth = catalog['linewidth (astrodendro)'][leaves]

        ax.scatter(size, linewidth, marker='.', s=10, c=color, label=co+' '+gal)
        x = np.logspace(np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][0]), np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['MCMC'][0],np.power(10,fit.c['MCMC'][0])), lw=0.5, color='k', zorder=3)

ax.set_xlabel('R [pc]')
ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*lw_min, 2.0*lw_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '04.size-linewidth', 'scatter_fit')
fig.savefig(join(plotdir, '04.size-linewidth', 'scatter_fit', 'size-linewidth_leaves.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# plot binned size-linewidth relation
###################################################################################################

# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['binned fit']
        bins       = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['bins']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        leaves    = (catalog['type']=='leaf')
        size      = catalog['size (astrodendro)'][leaves]
        linewidth = catalog['linewidth (astrodendro)'][leaves]

        R_min  = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max  = np.nanmax(size)
        lw_min = np.nanmin([x if x>1e-2 else np.nan for x in linewidth])
        lw_max = np.nanmax(linewidth)

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
        x = np.logspace(np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][0]), np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

        # fit uncertainty
        low = powlaw(x,fit.a['ls'][0]-fit.a['ls'][1],np.power(10,fit.c['ls'][0]-fit.c['ls'][1]))
        up  = powlaw(x,fit.a['ls'][0]+fit.a['ls'][1],np.power(10,fit.c['ls'][0]+fit.c['ls'][1]))
        ax.fill_between(x, low, up, facecolor='grey', edgecolor=None, alpha=0.5, zorder=2)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*lw_min, 2.0*lw_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '04.size-linewidth', 'binned_fit')
        fig.savefig(join(plotdir, '04.size-linewidth', 'binned_fit', 'size-linewidth_leaves.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

R_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
R_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))
lw_min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
lw_max = np.nanmax(flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        fit        = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['binned fit']
        bins       = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['bins']

        leaves    = (catalog['type']=='leaf')
        size      = catalog['size (astrodendro)'][leaves]
        linewidth = catalog['linewidth (astrodendro)'][leaves]

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

        x = np.logspace(np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][0]), np.log10(fits[co][gal]['size-linewidth (leaves)']['fit range'][1]), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

ax.set_xlabel('R [pc]')
ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*lw_min, 2.0*lw_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '04.size-linewidth', 'binned_fit')
fig.savefig(join(plotdir, '04.size-linewidth', 'binned_fit', 'size-linewidth_leaves.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# size-linewidth fit histograms
###################################################################################################

# scatter fit
###################################################################################################

fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='none', sharey='none', figsize=(8,4))

# plot parameter histograms
for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        color = dendrograms[co][gal]['color']
        fit   = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['fit']
        # fit   = fits[co][gal]['size-linewidth (leaves)']['size (astrodendro)']['linewidth (astrodendro)']['binned fit']
        axes[0].hist(fit.samples[:,0], histtype='step',       bins=200, range=(0.5,1.5),   density=True, stacked=True, color=color, label=co+' '+gal)
        axes[0].hist(fit.samples[:,0], histtype='stepfilled', bins=200, range=(0.5,1.5),   density=True, stacked=True, color=color, alpha=0.5)
        axes[1].hist(fit.y10s,         histtype='step',       bins=200, range=(0,25), density=True, stacked=True, color=color, label=co+' '+gal)
        axes[1].hist(fit.y10s,         histtype='stepfilled', bins=200, range=(0,25), density=True, stacked=True, color=color, alpha=0.5)

axes[0].legend(bbox_to_anchor=(0.,1.02,2.2,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=4, fancybox=True, fontsize=10)
axes[0].set_xlabel('R -- $\sigma$ slope')
axes[1].set_xlabel(r'$\sigma_\mathrm{10\,pc}$ [km\,s$^{-1}$]')
axes[1].yaxis.set_label_position('right')
axes[1].yaxis.tick_right()
for ax in axes:
    ax.set_ylabel('frequency')
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
axes[0].set_xlim([0.6,1.4])
axes[1].set_xlim([0,25])
mkpath(plotdir, '04.size-linewidth', 'histogram')
fig.savefig(join(plotdir, '04.size-linewidth', 'histogram', 'histograms_leaves.slope_sigma10.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
