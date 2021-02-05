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
# assign colors
###################################################################################################

dendrograms['CO(1-0)']['NGC253']['color'] = 'gold'
dendrograms['CO(1-0)']['GC']['color']     = 'orange'
dendrograms['CO(3-2)']['NGC253']['color'] = 'aqua'
dendrograms['CO(3-2)']['GC']['color']     = 'royalblue'


###################################################################################################
# plot size-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        size      = catalog['size (astrodendro)']
        linewidth = catalog['linewidth (astrodendro)']

        R_min  = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max  = np.nanmax(size)
        lw_min = np.nanmin([x if x>1e-2 else np.nan for x in linewidth])
        lw_max = np.nanmax(linewidth)

        ax.scatter(size, linewidth, marker='.', s=10, c=color, label=co+' '+gal)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*lw_min, 2.0*lw_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '04.size-linewidth', 'scatter')
        fig.savefig(os.path.join(plotdir, '04.size-linewidth', 'scatter', 'size-linewidth.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
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

        size      = catalog['size (astrodendro)']
        linewidth = catalog['linewidth (astrodendro)']

        ax.scatter(size, linewidth, marker='.', s=10, c=color,
                   label=co+' '+gal)

ax.set_xlabel('R [pc]')
ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*lw_min, 2.0*lw_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '04.size-linewidth', 'scatter')
fig.savefig(os.path.join(plotdir, '04.size-linewidth', 'scatter', 'size-linewidth.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
