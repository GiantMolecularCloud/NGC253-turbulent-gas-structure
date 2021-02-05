##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data        = fnunpickle(join(datadir, 'data.pickle'))
dendrograms = load_dendrograms()


###################################################################################################
# comapre linewidth estimates
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        fig = plt.figure(figsize=(8,8))
        ax  = fig.subplots(1)

        colors = mpl.cm.inferno(np.linspace(0,1,len(lw_measures)))

        min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([catalog[lw_measure].data for lw_measure in lw_measures])])
        max = np.nanmax([catalog[lw_measure].data for lw_measure in lw_measures])

        for lw_measure,color in zip(lw_measures,colors):
            linewidth_ad = catalog['linewidth (astrodendro)'].data
            linewidth    = catalog[lw_measure].data
            ax.scatter(linewidth_ad, linewidth, marker='.', s=10, c=[color],
                       label=lw_measure.replace('linewidth (','').replace(')','').replace('_',' ').replace('%',r'\%'),
                       zorder=3)

        # maybe fit if there is a consistent over-/underprediction

        ax.plot([0.5*min,2.0*max],[0.5*min,2.0*max], ls='-', color='grey', lw=1, zorder=2)

        ax.set_xlabel(r'$\sigma_\mathrm{astrodendro}$ [pc]')
        ax.set_ylabel(r'$\sigma$ [pc]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*min, 2.0*max)
        ax.set_ylim(0.5*min, 2.0*max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=3, fancybox=True, fontsize=10)
        mkpath(plotdir, '03.compare', 'linewidth')
        fig.savefig(join(plotdir, '03.compare', 'linewidth', 'compare_linewidth.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
