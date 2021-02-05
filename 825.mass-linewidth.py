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
# plot mass-linewidth relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        mass      = catalog['mass']
        linewidth = catalog['linewidth (astrodendro)']

        M_min  = np.nanmin([x if x!=0. else np.nan for x in mass])
        M_max  = np.nanmax(mass)
        lw_min = np.nanmin([x if x>1e-2 else np.nan for x in linewidth])
        lw_max = np.nanmax(linewidth)

        ax.scatter(mass, linewidth, marker='.', s=10, c=color, label=co+' '+gal)

        ax.set_xlabel('M [M$_\odot$]')
        ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*M_min, 2.0*M_max)
        ax.set_ylim(0.5*lw_min, 2.0*lw_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '06.mass-linewidth', 'scatter')
        fig.savefig(os.path.join(plotdir, '06.mass-linewidth', 'scatter', 'mass-linewidth.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

M_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()])])
M_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()]))
lw_min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
lw_max = np.nanmax(flatten([dendrograms[co][gal]['catalog']['linewidth (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        mass      = catalog['mass']
        linewidth = catalog['linewidth (astrodendro)']

        ax.scatter(mass, linewidth, marker='.', s=10, c=color,
                   label=co+' '+gal)

ax.set_xlabel('M [M$_\odot$]')
ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*M_min, 2.0*M_max)
ax.set_ylim(0.5*lw_min, 2.0*lw_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '06.mass-linewidth', 'scatter')
fig.savefig(os.path.join(plotdir, '06.mass-linewidth', 'scatter', 'mass-linewidth.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
