##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

execfile(join(scriptdir, '800.info.py'))
data = fnunpickle(join(datadir, 'data.pickle'))


###################################################################################################
# compute denddrograms
###################################################################################################

# dendrogram function
###################################################################################################

def get_dendrogram(cogal, min_snr=5., fac_delta=1., fac_npix=1.):
    import time,datetime
    from astrodendro import Dendrogram

    co  = cogal[0]
    gal = cogal[1]

    # get data and info
    file      = data[co][gal]['noise matched']['file']
    cube      = fits.open(file)[0]
    noise     = data[co][gal]['noise matched']['rms'].to(u.K)
    bmin      = angle_to_parsec(cube.header['bmin']*u.degree, source=gal)
    bmaj      = angle_to_parsec(cube.header['bmaj']*u.degree, source=gal)
    beam_area = 1.13309*bmin*bmaj
    pix1      = u.Quantity(str(np.abs(cube.header['cdelt1']))+cube.header['cunit1'])
    pix2      = u.Quantity(str(np.abs(cube.header['cdelt2']))+cube.header['cunit2'])
    pix3      = u.Quantity(str(np.abs(cube.header['cdelt3']))+cube.header['cunit3'])
    pix_beam  = (beam_area/pix1/pix2)

    # time execution
    print("\nComputing dendrogram: "+co+" "+gal+"\n")
    start = time.time()

    # calculate dendrogram
    d = Dendrogram.compute(cube.data,
                           #wcs = WCS(cube.header),           # causes a buffer overflow b/o wcslib errors: https://github.com/astropy/astropy/issues/7412
                           min_value = (min_snr*noise).value,
                           min_delta = (fac_delta*noise).value,
                           min_npix  = (fac_npix*pix_beam).value,
                           verbose   = True
                          )
    savepath = join(compdir,gal+'.'+co+'.dendrogram.fits')
    mkdir(os.path.dirname(escape_filename(savepath)))
    d.save_to(savepath)

    # time execution
    stop = time.time()
    exec_time = np.round(stop-start, 1)
    print("\nFinished dendrogram: "+co+" "+gal+"\nExecution took "+str(datetime.timedelta(seconds=exec_time))+ "hours.\n")


# execute dendrogram function in parallel
###################################################################################################

pool = Pool(4)
pool.map(get_dendrogram, [['CO(1-0)','NGC253'], ['CO(1-0)','GC'], ['CO(3-2)','NGC253'], ['CO(3-2)','GC']])
pool.close()


###################################################################################################
#
###################################################################################################
