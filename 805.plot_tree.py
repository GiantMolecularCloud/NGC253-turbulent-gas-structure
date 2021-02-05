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
# plot dendrogram trees
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,5))

        # plot the tree
        plot = dendrograms[co][gal]['dendrogram'].plotter()
        plot.plot_tree(ax, color='black', lw=0.5)

        # labels
        ax.set_xlabel(r'structure \#', fontsize=12)
        ax.set_ylabel(r'T$_\mathrm{b}$ [K]', fontsize=12)
        ax.text(0.05,0.9, co+' '+gal, fontsize=12, va='top', ha='left', transform=ax.transAxes, bbox=props)
        fig.tight_layout()

        mkpath(plotdir, '02.trees')
        fig.savefig(join(plotdir, '02.trees', 'tree.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
