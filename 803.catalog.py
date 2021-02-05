##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data = fnunpickle(join(datadir, 'data.pickle'))


###################################################################################################
# compute catalog
###################################################################################################

# catalog function
###################################################################################################

def get_catalog(cogal):
    from astrodendro import Dendrogram
    from astrodendro import ppv_catalog

    co  = cogal[0]
    gal = cogal[1]

    # get data and info
    file      = data[co][gal]['noise matched']['file']
    cube      = fits.open(file)[0]
    noise     = data[co][gal]['noise matched']['rms'].to(u.K)
    bmin      = cube.header['bmin']*u.degree
    bmaj      = cube.header['bmaj']*u.degree
    beam_area = 1.13309*bmin*bmaj
    pix1      = u.Quantity(str(np.abs(cube.header['cdelt1']))+cube.header['cunit1'])
    pix2      = u.Quantity(str(np.abs(cube.header['cdelt2']))+cube.header['cunit2'])
    pix3      = u.Quantity(str(np.abs(cube.header['cdelt3']))+cube.header['cunit3'])
    pix_beam  = (beam_area/pix1/pix2)
    d = Dendrogram.load_from(join(compdir,gal+'.'+co+'.dendrogram.fits'))

    # time execution
    print("\nComputing catalog: "+co+" "+gal+"\n")
    start = time.time()

    # calculate catalog
    metadata = {'data_unit': u.K,
                'spatial_scale': parsec_to_angle(pix1, source=gal),
                'velocity_scale': pix3,
                'beam_major': bmaj,
                'beam_minor': bmin,
                'vaxis': 0,
                'wavelength': lines[co]['restfreq'],
                'wcs': WCS(cube.header)
               }
    catalog = ppv_catalog(d, metadata)
    fnpickle(catalog, join(compdir,gal+'.'+co+'.catalog.pickle'))

    # time execution
    stop = time.time()
    exec_time = np.round(stop-start, 1)
    print("\nFinished catalog: "+co+" "+gal+"\nExecution took "+str(datetime.timedelta(seconds=exec_time))+ "hours for "+str(len(catalog))+" structures.\n")


# execute catalog function in parallel
###################################################################################################

pool = Pool(4)
pool.map(get_catalog, [['CO(1-0)','NGC253'], ['CO(1-0)','GC'], ['CO(3-2)','NGC253'], ['CO(3-2)','GC']])
pool.close()


###################################################################################################
#
###################################################################################################
