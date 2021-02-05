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
# plot cloud PDFs
###################################################################################################

# column density spectra
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        density = catalog['column density']
        density = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e19),np.log10(1e24), 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        ax.set_xlabel(r'$\Sigma$ [cm$^{-2}$]')
        ax.set_ylabel(r'cloud density PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlim([5e-2, 5e4])
        # ax.set_ylim(ylim)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '08.PDF', 'column_density')
        fig.savefig(os.path.join(plotdir, '08.PDF', 'column_density', 'column_density_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# volume density spectra
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        density = catalog['volume density (astrodendro)']
        density = density[np.isfinite(density)]

        bins  = np.logspace(np.log10(1e-2),np.log10(1e5), 50+1)
        hist,bins = np.histogram(density, bins=bins, density=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        ax.set_xlabel(r'$\rho_\mathrm{vol}$ [cm$^{-3}$]')
        ax.set_ylabel(r'cloud density PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlim([5e-2, 5e4])
        # ax.set_ylim(ylim)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '08.PDF', 'volume_density')
        fig.savefig(os.path.join(plotdir, '08.PDF', 'volume_density', 'volume_density_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# mass spectra
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        mass   = catalog['mass']
        mass = mass[np.isfinite(mass)]

        bins  = np.logspace(np.log10(1e1),np.log10(1e9), 50+1)
        hist,bins = np.histogram(mass, bins=bins, density=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        ax.set_xlabel('M [M$_\odot$]')
        ax.set_ylabel('cloud mass PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        # ax.set_xlim([1e1, 1e9])
        # ax.set_ylim(ylim)
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '08.PDF', 'mass')
        fig.savefig(os.path.join(plotdir, '08.PDF', 'mass', 'mass_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
