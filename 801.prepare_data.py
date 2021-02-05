##############
# DENDROGRAM #
##############

###################################################################################################
# prepare data in CASA
###################################################################################################

execfile(join(scriptdir, '800.info.py'))


###################################################################################################
# execute data preparation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        cube = SpectralCube.read(data[co][gal]['original']['file'])
        cube.allow_huge_operations = True


# get image info
###################################################################################################

        bmaj = cube.header['bmaj']*u.degree
        bmin = cube.header['bmin']*u.degree
        bpa  = cube.header['bpa']*u.degree
        pix1 = u.Quantity(str(np.abs(cube.header['cdelt1']))+cube.header['cunit1'])
        pix2 = u.Quantity(str(np.abs(cube.header['cdelt2']))+cube.header['cunit2'])
        spec = u.Quantity(str(np.abs(cube.header['cdelt3']))+cube.header['cunit3'])

        res_spat = angle_to_parsec((bmaj+bmin)/2., source=gal)
        if spec.unit.is_equivalent(u.km/u.s):
            res_spec = spec.to(u.km/u.s)
        elif spec.unit.is_equivalent(u.Hz):
            res_spec = np.round(np.abs((CO['restfreq']+spec).to(u.km/u.s, equivalencies=u.doppler_radio(CO['restfreq']))), 1)


# smooth to circular common beam
###################################################################################################

        # target beam
        res = (np.sin((CO['resolution']['spatial']/GAL['distance']).to(u.dimensionless_unscaled).value)*u.rad).to(u.degree)

        # smooth cube to target beam
        import radio_beam
        beam = radio_beam.Beam(major=res, minor=res, pa=0*u.degree)
        cube1 = cube.convolve_to(beam)


# reproject onto common pixel grid
###################################################################################################

        # construct target image header
        reheader = copy.deepcopy(cube1.header)
        if cube.header['cdelt1']>0:
            pix_x    = (res/5.).to(u.degree).value
            origin_x = cube1.longitude_extrema[1].to(u.degree).value if gal=='GC' else cube1.longitude_extrema[0].to(u.degree).value
        else:
            pix_x    = -1.*(res/5.).to(u.degree).value
            origin_x = (cube1.longitude_extrema[0]).to(u.degree).value if gal=='GC' else (cube1.longitude_extrema[1]).to(u.degree).value
        if cube.header['cdelt2']>0:
            pix_y    = (res/5.).to(u.degree).value
            origin_y = cube1.latitude_extrema[0].to(u.degree).value
        else:
            pix_y    = -1.*(res/5.).to(u.degree).value
            origin_y = cube1.latitude_extrema[1].to(u.degree).value
        if gal=='NGC253':
            npix_x   = int(np.ceil(np.diff(cube1.longitude_extrema, n=1)[0]/np.abs(pix_x)).value)
            npix_y   = int(np.ceil(np.diff(cube1.latitude_extrema, n=1)[0]/np.abs(pix_y)).value)
        elif gal=='GC':
            npix_x   = int(np.ceil((cube1.longitude_extrema[0]-(cube1.longitude_extrema[1]-360*u.degree))/np.abs(pix_x)).value)
            npix_y   = int(np.ceil(np.diff(cube1.latitude_extrema, n=1)[0]/np.abs(pix_y)).value)

        reheader['cdelt1'] = pix_x
        reheader['cdelt2'] = pix_y
        reheader['naxis1'] = npix_x
        reheader['naxis2'] = npix_y
        reheader['crval1'] = origin_x
        reheader['crval2'] = origin_y
        reheader['crpix1'] = 0
        reheader['crpix2'] = 0
        try:
            del reheader['lonpole']
            del reheader['latpole']
            del reheader['wcsaxes']
        except:
            pass

        # regrid cube to target pixel size
        cube2 = cube1.reproject(reheader, order='bilinear', use_memmap=False, filled=True)


# smooth to common spectral resolution and interpolate onto common spectral grid
###################################################################################################

        from astropy.convolution import Gaussian1DKernel

        spectral_axis = GAL['vsys'] + np.arange(-250,251,CO['resolution']['spectral'].value)*CO['resolution']['spectral'].unit
        spectral_factor = (CO['resolution']['spectral']/res_spec).value
        fwhm_factor = np.sqrt(8*np.log(2))

        cube3 = cube2.spectral_smooth(Gaussian1DKernel(4/fwhm_factor))
        if cube3.spectral_axis.unit.is_equivalent(u.km/u.s):
            cube3_kms = cube3
        elif cube3.spectral_axis.unit.is_equivalent(u.Hz):
            cube3_kms = cube3.with_spectral_unit(u.km/u.s, velocity_convention='optical', rest_value=CO['restfreq'])
        cube4 = cube3_kms.spectral_interpolate(spectral_axis, suppress_smooth_warning=True)
        cube4 = cube4.hdu
        cube4.header['crval3'] = -250
        cube4 = SpectralCube.read(cube4)
        cube4.write(join(datadir,gal+'.'+co+'.regridded.fits'), overwrite=True)


###################################################################################################
# align orientation
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        cube4 = SpectralCube.read(join(datadir,gal+'.'+co+'.regridded.fits'))
        cube4.allow_huge_operations = True


# set reference pixels to galaxy center
###################################################################################################

        # Spectalcube does not allow to write header data! Need to use astropy instead.
        cube5 = cube4.hdu

        # get center coordinates
        try:
            center_x = GAL['center'].ra
        except:
            center_x = GAL['center'].l
        try:
            center_y = GAL['center'].dec
        except:
            center_y = GAL['center'].b

        # get image coordinates
        spine_x = cube4.world[0,0,:][2]
        spine_y = cube4.world[0,:,0][1]

        # find pixel index of center
        diff_x = [np.abs(x-center_x).value for x in spine_x]
        diff_y = [np.abs(y-center_y).value for y in spine_y]
        crpix1 = np.argmin(diff_x)
        crpix2 = np.argmin(diff_y)

        # update header
        cube5.header['crpix1'] = crpix1
        cube5.header['crpix2'] = crpix2
        cube5.header['crval1'] = 0.
        cube5.header['crval2'] = 0.


# rotate NGC253
###################################################################################################

        if gal=='NGC253':

            # rotation matrix for 33 degree
            rotation_matrix = {'PC1_1':  8.386705679454E-01,
                               'PC2_1': -5.446390350150E-01,
                               'PC3_1':  0.000000000000E+00,
                               'PC1_2': 5.446390350150E-01,
                               'PC2_2':  8.386705679454E-01,
                               'PC3_2':  0.000000000000E+00,
                               'PC1_3':  0.000000000000E+00,
                               'PC2_3':  0.000000000000E+00,
                               'PC3_3':  1.000000000000E+00}

            reheader = copy.deepcopy(cube5.header)
            for key,val in rotation_matrix.items():
                reheader[key] = val

            # rotate by regridding onto rotated header
            cube5 = SpectralCube.read(cube5)
            cube5.allow_huge_operations = True
            cube6 = cube5.reproject(reheader, order='bilinear', use_memmap=False, filled=True)
            cube6 = cube6.hdu

            # remove roation matrix
            for key in rotation_matrix.keys():
                try:
                    del cube6.header[key]
                except:
                    pass

        else:
            cube6 = cube5


# convert to pc grid
###################################################################################################

        cdelt1 = u.Quantity(str(cube6.header['cdelt1'])+cube6.header['cunit1'])
        cdelt2 = u.Quantity(str(cube6.header['cdelt2'])+cube6.header['cunit2'])

        cdelt1_pc = angle_to_parsec(cdelt1, source=gal)
        cdelt2_pc = angle_to_parsec(cdelt2, source=gal)

        cube6.header['cdelt1'] = cdelt1_pc.value
        cube6.header['cdelt2'] = cdelt2_pc.value
        cube6.header['cunit1'] = 'pc'
        cube6.header['cunit2'] = 'pc'
        cube6.header['ctype1'] = 'linear'
        cube6.header['ctype2'] = 'linear'


# extract matching subcube (take inclination into account)
###################################################################################################

        fov_x = CO['fov'][0]
        fov_y = CO['fov'][1]

        # get cube info
        center_x = int(cube6.header['crpix1'])
        center_y = int(cube6.header['crpix2'])
        cdelt_x = u.Quantity(str(np.abs(cube6.header['cdelt1']))+cube6.header['cunit1'])
        cdelt_y = u.Quantity(str(np.abs(cube6.header['cdelt2']))+cube6.header['cunit2'])

        # calculate pixel position of subcube
        npix_x = int(np.round(fov_x/2./cdelt_x))
        npix_y = int(np.round(fov_y/2./cdelt_y))
        cutout_x = [center_x-npix_x, center_x+npix_x]
        cutout_y = [center_y-npix_y, center_y+npix_y]

        # check image dimensions
        grow_x1 = 0-cutout_x[0] if cutout_x[0]<0 else 0
        grow_x2 = cutout_x[1]-cube6.header['naxis1'] if cutout_x[1]>cube6.header['naxis1'] else 0
        grow_x  = np.max([grow_x1, grow_x2])

        grow_y1 = 0-cutout_y[0] if cutout_y[0]<0 else 0
        grow_y2 = cutout_y[1]-cube6.header['naxis2'] if cutout_y[1]>cube6.header['naxis2'] else 0
        grow_y  = np.max([grow_y1, grow_y2])

        if grow_x>0 or grow_y>0:
            print("Growing image to allow cutting out symmetric subimage.")
            shape = cube6.shape
            grown = np.full((shape[0],shape[1]+2*grow_y,shape[2]+2*grow_x),np.nan)
            for v,V in tqdm(enumerate(cube6.data)):
                for y,Y in enumerate(V):
                    for x,X in enumerate(Y):
                        grown[v,y+grow_y,x+grow_x] = X
            cube6.data = grown
            cutout_x = [x+grow_x for x in cutout_x]
            cutout_y = [y+grow_y for y in cutout_y]

        # cut out subimage
        cube7   = cube6.data[:,cutout_y[0]:cutout_y[1],cutout_x[0]:cutout_x[1]]
        header7 = copy.deepcopy(cube6.header)

        # update reference pixels
        header7['crpix1'] = npix_x
        header7['crpix2'] = npix_y

        fits.writeto(join(datadir,gal+'.'+co+'.matched.fits'),
                     data      = cube7,
                     header    = header7,
                     overwrite = True
                    )
        data[co][gal]['matched']['file'] = join(datadir,gal+'.'+co+'.matched.fits')


fnpickle(data, join(datadir, 'data.pickle'))


###################################################################################################
# match noise
###################################################################################################

# measure noise by hand
###################################################################################################

data['CO(1-0)']['GC']['matched']['rms']     = 32*u.mK
data['CO(1-0)']['NGC253']['matched']['rms'] = 38*u.mK
data['CO(3-2)']['GC']['matched']['rms']     = 37*u.mK
data['CO(3-2)']['NGC253']['matched']['rms'] = 115*u.mK
fnpickle(data, join(datadir, 'data.pickle'))


# measure noise automatically
###################################################################################################

fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='col', sharey='row', figsize=(10,8))
for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        im   = fits.open(data[co][gal]['matched']['file'])[0]
        rms  = np.sqrt(np.nanmean(np.square(im.data), axis=(1,2)))*1000
        velo = [im.header['cdelt3']*(i-im.header['crpix3']+1)+im.header['crval3'] for i in np.arange(im.header['naxis3'])]
        ax.plot(velo, rms, ls='-' if co=='CO(3-2)' else '--', lw=2, color=data[co][gal]['color'], label=co+' '+gal)
        ax.axhline(y=data[co][gal]['matched']['rms'].to(u.mK).value, ls='-' if co=='CO(3-2)' else '--', lw=2, color=data[co][gal]['color'], alpha=0.5, zorder=2)

ax.set_xlim(-250, 250)
ax.set_yscale('log')
ax.xaxis.set_major_locator(MultipleLocator(50))
ax.xaxis.set_minor_locator(MultipleLocator(10))
ax.set_axisbelow(True)
ax.grid(axis='both')
ax.set_ylabel(r'rms [mK]', fontsize=12)
ax.set_xlabel(r'v [km\,s$^{-1}$]')
fig.legend(loc='best') #, bbox_to_anchor=(0.02, 0.975), bbox_transform=ax.transAxes)
fig.tight_layout()

savepath = join(plotdir, '01.data', 'noise.pdf')
os.system('mkdir -p '+os.path.dirname(savepath))
fig.savefig(savepath, dpi=300, bbox_inches='tight')


# add noise to get identical noise levels
###################################################################################################

# assume Gaussian noise and rms behave similar to Gaussian dispersion,
# i.e. dispersion (noise) adds in quadrature

# match GC data to NGC253 data
for co,CO in lines.items():

    im = fits.open(data[co]['GC']['matched']['file'])[0]
    npixels = np.product(im.data.shape)

    target_noise = data[co]['NGC253']['matched']['rms']
    actual_noise = data[co]['GC']['matched']['rms']
    additional_sigma = np.sqrt(target_noise**2 - actual_noise**2).to(u.K).value

    additional_noise = np.random.normal(0., additional_sigma, npixels)
    additional_noise = np.reshape(additional_noise, im.data.shape)

    new_data = im.data+additional_noise

    fits.writeto(join(datadir,'GC.'+co+'.noise_matched.fits'),
                 data   = new_data,
                 header = im.header,
                 overwrite = True
                )
    data[co]['GC']['noise matched'] = {'file': join(datadir,'GC.'+co+'.noise_matched.fits'), 'rms': target_noise}

# keep NGC253 data as is
    os.system('cp '+escape_filename(join(datadir,'NGC253.'+co+'.matched.fits'))+' '+escape_filename(join(datadir,'NGC253.'+co+'.noise_matched.fits')))
    data[co]['NGC253']['noise matched'] = {'file': join(datadir,'NGC253.'+co+'.noise_matched.fits'), 'rms': target_noise}


fnpickle(data, join(datadir, 'data.pickle'))

###################################################################################################
#
###################################################################################################
