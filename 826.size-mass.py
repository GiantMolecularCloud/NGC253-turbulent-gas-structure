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
# plot size-mass relation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        size      = catalog['size (astrodendro)']
        mass      = catalog['mass']

        R_min  = np.nanmin([x if x!=0. else np.nan for x in size])
        R_max  = np.nanmax(size)
        M_min  = np.nanmin([x if x!=0. else np.nan for x in mass])
        M_max  = np.nanmax(mass)

        ax.scatter(size, mass, marker='.', s=10, c=color, label=co+' '+gal)

        ax.set_xlabel('R [pc]')
        ax.set_ylabel('M [M$_\odot$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*R_min, 2.0*R_max)
        ax.set_ylim(0.5*M_min, 2.0*M_max)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '07.size-mass', 'scatter')
        fig.savefig(os.path.join(plotdir, '07.size-mass', 'scatter', 'mass-linewidth.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

R_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()])])
R_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['size (astrodendro)'].data for co in lines.keys() for gal in galaxies.keys()]))
M_min  = np.nanmin([x if x!=0. else np.nan for x in flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()])])
M_max  = np.nanmax(flatten([dendrograms[co][gal]['catalog']['mass'].data for co in lines.keys() for gal in galaxies.keys()]))

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        size      = catalog['size (astrodendro)']
        mass      = catalog['mass']

        ax.scatter(size, mass, marker='.', s=10, c=color,
                   label=co+' '+gal)

ax.set_xlabel('R [pc]')
ax.set_ylabel('M [M$_\odot$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(0.5*R_min, 2.0*R_max)
ax.set_ylim(0.5*M_min, 2.0*M_max)
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '07.size-mass', 'scatter')
fig.savefig(os.path.join(plotdir, '07.size-mass', 'scatter', 'mass-linewidth.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
