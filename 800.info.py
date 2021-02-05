###############################
# DENDROGRAM SUB-PROJECT INFO #
###############################

# settings for project 800: dendrogram cloud comparison in NGC253 and GC


####################################################################################################
# directories
####################################################################################################

subprojectdir  = 'NGC253/800.dendrogram/'
datadir        = join(subprojectdir, '01.data/')
compdir        = join(subprojectdir, '02.galaxy_comparison/')

plotdir        = 'plots/NGC253/800.dendrogram/'
compplotdir    = join(plotdir,'02.galaxy_comparison/')
paperplotdir   = 'plots/NGC253/paper_20b/'


####################################################################################################
# data
####################################################################################################

lines    = {'CO(1-0)': {'restfreq': 115.27120180*u.GHz,
                        'tex': 'CO(1--0)',
                        'excitation correction': 1.00,
                        'resolution': {'spatial': 32*u.pc, 'spectral': 5.0*u.km/u.s},
                        'fov': [1500*u.pc, 750*u.pc],
                        'rms': 38*u.mK,
                        'weight': 28*u.u},
            'CO(3-2)': {'restfreq': 345.79598990*u.GHz,
                        'tex': 'CO(3--2)',
                        'excitation correction': 0.67,
                        'resolution': {'spatial': 3.0*u.pc, 'spectral': 2.5*u.km/u.s},
                        'fov': [800*u.pc, 400*u.pc],
                        'rms': 115*u.mK,
                        'weight': 28*u.u}
           }
galaxies = {'NGC253': {'distance': 3.5*u.Mpc,
                       'vsys': 250*u.km/u.s,
                       'PA': 33*u.degree,
                       'center': SkyCoord('00h47m33.134s -25d17m19.68s'),
                       'inclination': 78*u.degree,
                       'Xco': 0.5e20},
            'GC':     {'distance': 8178*u.pc,
                       'vsys': 0*u.km/u.s,
                       'center': SkyCoord(0.*u.degree, 0.*u.degree, frame='galactic'),
                       'Xco': 1.0e20}
           }
data = {'CO(1-0)': {'NGC253': {'color': 'dodgerblue',
                               'original': {'rms': 0.080*u.K, 'file': datadir+'raw/NGC253.CO_1-0.orig.fits'},
                               'matched':  {'rms': None, 'file': None}},
                    'GC':     {'color': 'darkorange',
                               'original': {'rms': 0.092*u.K, 'file': datadir+'raw/COGAL_deep_interp.ordered.shift.2kpc.fits'},
                               'matched':  {'rms': None, 'file': None}}
                    },
        'CO(3-2)': {'NGC253': {'color': 'dodgerblue',
                               'original': {'rms': 0.330*u.K, 'file': datadir+'raw/NGC253.band7.TP+12m-mid+12m-high.CO_3-2.image.fits'},
                               'matched':  {'rms': None, 'file': None}},
                    'GC':     {'color': 'darkorange',
                               'original': {'rms': 0.000*u.K, 'file': datadir+'raw/CHIMPS2_CO_3-2.fits'},
                               'matched':  {'rms': None, 'file': None}}
                    }
       }

R_measures = ['size (astrodendro)',
              'size (area_exact)',
              'size (area_ellipse)',
              'size (manual)',
              'size (effective)']
lw_measures = ['linewidth (astrodendro)',
               'linewidth (mom2 mean)',
               'linewidth (mom2 median)',
               'linewidth (90% flux)',
               'linewidth (FWHM)',
               'linewidth (FW10%)']


###################################################################################################
# load dendrogram pipeline
###################################################################################################

if not 'casa' in globals():
    interactive = False
    from astrodendro import Dendrogram
    from spectral_cube import SpectralCube

    # increase the recursion limit or astrodendro with many structures will fail
    sys.setrecursionlimit(100000)


###################################################################################################
#
###################################################################################################
