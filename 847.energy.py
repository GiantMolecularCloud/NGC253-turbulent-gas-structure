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
# energy
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        mass      = catalog['mass'].quantity
        linewidth = catalog['linewidth (astrodendro)'].quantity

        Ekin = (0.5*mass*linewidth**2).to(u.erg)

        catalog.add_column( Column(name='kinetic energy', data=Ekin.value, dtype=np.float64, unit='erg') )


###################################################################################################
# energy histogram
###################################################################################################

# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        Ekin = catalog['kinetic energy']
        Ekin = Ekin[np.isfinite(Ekin)]

        bins  = np.logspace(np.log10(1e40),np.log10(1e60), 40+1)
        hist,bins = np.histogram(Ekin, bins=bins, density=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

ax.set_xlabel(r'E$_\mathrm{kin}$ [erg]')
ax.set_ylabel(r'PDF')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
mkpath(plotdir, '08.PDF', 'energy')
fig.savefig(os.path.join(plotdir, '08.PDF', 'energy', 'energy_PDF.all.pdf'), dpi=300, bbox_inches='tight')


# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        Ekin = catalog['kinetic energy']
        Ekin = Ekin[np.isfinite(Ekin)]

        bins  = np.logspace(np.log10(1e40),np.log10(1e60), 40+1)
        hist,bins = np.histogram(Ekin, bins=bins, density=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        ax.set_xlabel(r'E$_\mathrm{kin}$ [erg]')
        ax.set_ylabel(r'PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '08.PDF', 'energy')
        fig.savefig(os.path.join(plotdir, '08.PDF', 'energy', 'energy_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
