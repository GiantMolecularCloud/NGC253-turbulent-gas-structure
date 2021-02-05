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

        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        hist  = fits[co][gal]['column density PDF']['hist']
        bins  = fits[co][gal]['column density PDF']['bins']
        shape = fits[co][gal]['column density PDF']['shape']
        loc   = fits[co][gal]['column density PDF']['loc']
        scale = fits[co][gal]['column density PDF']['scale']

        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        x = np.logspace(np.log10(fits[co][gal]['column density PDF']['fit range'][0]),np.log10(fits[co][gal]['column density PDF']['fit range'][1]),100)
        pdf = lognorm.pdf(x, shape, loc, scale)
        ax.plot(x, pdf, color='k', ls='-')

        ax.set_xlabel(r'$\Sigma$ [cm$^{-3}$]')
        ax.set_ylabel(r'cloud density PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e19, 1e24])
        ax.set_ylim(0.5*np.min(hist[hist>0.0]), 2*np.max(hist[hist>0.0]))
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '06.PDF', 'column_density_fitted')
        fig.savefig(os.path.join(plotdir, '06.PDF', 'column_density_fitted', 'column_density_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')

# normal fit to logged data
for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        hist  = fits[co][gal]['column density PDF']['hist']
        bins  = fits[co][gal]['column density PDF']['bins']
        mu    = fits[co][gal]['column density PDF']['mu']
        sigma = fits[co][gal]['column density PDF']['sigma']

        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5)

        x = np.linspace(np.log10(fits[co][gal]['column density PDF']['fit range'][0]),np.log10(fits[co][gal]['column density PDF']['fit range'][1]),100)
        pdf = norm.pdf(x, mu, sigma)
        ax.plot(x, pdf, color='k', ls='-')

        ax.set_xlabel(r'log $\Sigma$ [cm$^{-3}$]')
        ax.set_ylabel(r'cloud density PDF')
        ax.set_xlim([19, 24])
        ax.set_ylim(0.5*np.min(hist[hist>0.0]), 2*np.max(hist[hist>0.0]))
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '06.PDF', 'column_density_fitted')
        fig.savefig(os.path.join(plotdir, '06.PDF', 'column_density_fitted', 'column_density_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# volume density spectra
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        density   = catalog['volume density (astrodendro)']
        density = density[np.isfinite(density)]

        bins  = fits[co][gal]['volume density PDF']['bins']
        hist  = fits[co][gal]['volume density PDF']['hist']
        shape = fits[co][gal]['volume density PDF']['shape']
        loc   = fits[co][gal]['volume density PDF']['loc']
        scale = fits[co][gal]['volume density PDF']['scale']

        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        x = np.logspace(np.log10(fits[co][gal]['volume density PDF']['fit range'][0]),np.log10(fits[co][gal]['volume density PDF']['fit range'][1]),100)
        pdf = lognorm.pdf(x, shape, loc, scale)
        ax.plot(x, pdf, color='k', ls='-')

        ax.set_xlabel(r'$\rho_\mathrm{vol}$ [cm$^{-3}$]')
        ax.set_ylabel(r'cloud density PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e-1, 1e5])
        ax.set_ylim(0.5*np.min(hist[hist>0.0]), 2*np.max(hist[hist>0.0]))
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '06.PDF', 'volume_density_fitted')
        fig.savefig(os.path.join(plotdir, '06.PDF', 'volume_density_fitted', 'volume_density_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


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

        bins      = fits[co][gal]['mass PDF']['bins']
        hist      = fits[co][gal]['mass PDF']['hist']
        slope     = fits[co][gal]['mass PDF']['slope']
        intercept = fits[co][gal]['mass PDF']['intercept']

        ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
        ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)

        fit_fn   = np.poly1d([slope['value'], intercept['value']])
        x  = np.log10(np.linspace(fits[co][gal]['mass PDF']['fit range'][0],fits[co][gal]['mass PDF']['fit range'][1], 10))
        ax.plot(np.power(10,x),  np.power(10,fit_fn(x)),  ls='-', color='k')

        ax.set_xlabel('M [M$_\odot$]')
        ax.set_ylabel('cloud mass PDF')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([1e1, 1e9])
        ax.set_ylim(0.5*np.min(hist[hist>0.0]), 2*np.max(hist[hist>0.0]))
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '06.PDF', 'mass_fitted')
        fig.savefig(os.path.join(plotdir, '06.PDF', 'mass_fitted', 'mass_PDF.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
