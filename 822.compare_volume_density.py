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
# compare volume density estimates
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        fig = plt.figure(figsize=(8,8))
        ax  = fig.subplots(1)

        colors = mpl.cm.inferno(np.linspace(0,1,len(R_measures)))

        min = np.nanmin([x if x!=0. else np.nan for x in flatten([catalog['volume density ('+(R_measure.split('(')[1]).split(')')[0]+')'].data for R_measure in R_measures])])
        max = np.nanmax([catalog['volume density ('+(R_measure.split('(')[1]).split(')')[0]+')'].data for R_measure in R_measures])

        for R_measure,color in zip(R_measures,colors):
            density_ad = catalog['volume density (astrodendro)'].data
            density    = catalog['volume density ('+(R_measure.split('(')[1]).split(')')[0]+')'].data
            ax.scatter(density_ad, density, marker='.', s=10, c=[color],
                       label=R_measure.replace('size ','').replace(')','').replace('_',' '),
                       zorder=3)

        # maybe fit if there is a consistent over-/underprediction

        ax.plot([0.5*min,2.0*max],[0.5*min,2.0*max], ls='-', color='grey', lw=1, zorder=2)

        ax.set_xlabel(r'$\rho_\mathrm{astrodendro}$ [cm$^{-3}$]')
        ax.set_ylabel(r'$\rho$ ('+(R_measure.split('(')[1]).split(')')[0]+') [cm$^{-3}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*min, 2.0*max)
        ax.set_ylim(0.5*min, 2.0*max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
        mkpath(plotdir, '03.compare', 'volume_density')
        fig.savefig(join(plotdir, '03.compare', 'volume_density', 'compare_volume_density.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
