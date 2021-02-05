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
# compare size estimates
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        fig = plt.figure(figsize=(8,8))
        ax  = fig.subplots(1)

        colors = mpl.cm.inferno(np.linspace(0,1,len(R_measures)))

        min = np.nanmin([x if x!=0. else np.nan for x in flatten([catalog[R_measure].data for R_measure in R_measures])])
        max = np.nanmax([catalog[R_measure].data for R_measure in R_measures])

        for R_measure,color in zip(R_measures,colors):
            size_ad = catalog['size (astrodendro)'].data
            size    = catalog[R_measure].data
            ax.scatter(size_ad, size, marker='.', s=10, c=[color],
                       label=R_measure.replace('size (','').replace(')','').replace('_',' '),
                       zorder=3)

        # maybe fit if there is a consistent over-/underprediction

        ax.plot([0.5*min,2.0*max],[0.5*min,2.0*max], ls='-', color='grey', lw=1, zorder=2)

        ax.set_xlabel(r'R$_\mathrm{astrodendro}$ [pc]')
        ax.set_ylabel(r'R [pc]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*min, 2.0*max)
        ax.set_ylim(0.5*min, 2.0*max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=3, fancybox=True, fontsize=10)
        mkpath(plotdir, '03.compare', 'size')
        fig.savefig(join(plotdir, '03.compare', 'size', 'compare_size.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
