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
# get linewidth
###################################################################################################

# astrodendro linewidth
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        lw     = catalog['v_rms']
        log_lw = [np.log10(l) for l in lw]

        catalog.add_column( Column(name='linewidth (astrodendro)',     data=lw,     dtype=np.float64, unit='km/s') )
        catalog.add_column( Column(name='log linewidth (astrodendro)', data=log_lw, dtype=np.float64, unit='log10 km/s') )


# calculate moment maps for all structures
###################################################################################################

def moment0(cube, chanwidth, mask_zeros=True):
    mom0 = chanwidth*np.nansum(cube.value, axis=0)*cube.unit
    if mask_zeros:
        mom0[mom0==0.0] = np.nan
    return mom0

def moment1(cube, moment0, chanwidth, pix_centers):
    mom1 = np.nansum(cube*pix_centers*chanwidth, axis=0)/moment0
    return mom1

def moment2(cube, moment0, moment1, chanwidth, pix_centers, mask_zeros=True):
    mom2 = np.full_like(moment1.value, 0.0) *cube.unit*pix_centers.unit**2*chanwidth.unit
    for chan in np.arange(cube.shape[0]):
        mom2 += cube[chan] *(pix_centers[chan]-moment1)**2 *chanwidth
    mom2 = np.sqrt(mom2/moment0)
    if mask_zeros:
        mom2[mom2==0.0] = np.nan
    return mom2

moms = {co:{gal: {'mom0': {}, 'mom1': {}, 'mom2': {}} for gal in galaxies.keys()} for co in lines.keys()}

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        cube = fits.open(data[co][gal]['noise matched']['file'])[0]
        chans,ny,nx = cube.data.shape
        bunit = u.Unit(cube.header['bunit'])
        cunit = u.Unit(cube.header['cunit3'])
        cdelt = cube.header['cdelt3']*cunit
        crval = cube.header['crval3']*cunit
        crpix = cube.header['crpix3']
        pix_cen = np.array([((i+crpix-1)*cdelt+crval).value for i in np.arange(chans)])*cunit
        pix_cen = np.array([np.full((ny,nx), cen) for cen in pix_cen])*pix_cen.unit

        # mask cube
        min_snr = 5.            # check dendrogram threshold!
        cube.data[cube.data*bunit<min_snr*data[co][gal]['noise matched']['rms']] = np.nan

        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask    = struct.get_mask()
            mcube   = np.nan_to_num(cube.data*mask)*bunit

            mom0 = moment0(mcube,             cdelt)
            mom1 = moment1(mcube, mom0,       cdelt, pix_cen)
            mom2 = moment2(mcube, mom0, mom1, cdelt, pix_cen)

            moms[co][gal]['mom0'][struct.idx] = mom0
            moms[co][gal]['mom1'][struct.idx] = mom1
            moms[co][gal]['mom2'][struct.idx] = mom2

        fnpickle(moms[co][gal], join(compdir,gal+'.'+co+'.moments.pickle'))

# Do NOT save as single file!
# Extremely large for this many structures!
#fnpickle(moms, join(compdir,'moments.pickle'))

# # load moms after restart
# moms = {co:{gal: {'mom0': {}, 'mom1': {}, 'mom2': {}} for gal in galaxies.keys()} for co in lines.keys()}
# for co,CO in lines.items():
#     for gal,GAL in galaxies.items():
#         moms[co][gal] = fnunpickle(join(compdir,gal+'.'+co+'.moments.pickle'))


# mean/median of per pixel moment 2
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']

        lw = {'mean': [], 'std': [], 'median': [], '16th': [], '84th': []}
        for idx in tqdm(moms[co][gal]['mom0'].keys()):
            mom2 = moms[co][gal]['mom2'][idx]
            mean = np.nanmean(mom2.value)
            std  = np.nanstd(mom2.value)
            p16,median,p84 = np.nanpercentile(mom2.value, (16,50,84))
            lw['mean'].append(mean)
            lw['std'].append(std)
            lw['median'].append(median)
            lw['16th'].append(p16)
            lw['84th'].append(p84)
        log_lw = {k:np.log10(v) for k,v in lw.items()}

        for k in lw.keys():
            catalog.add_column( Column(name='linewidth (mom2 '+k+')',   data=lw[k],   dtype=np.float64, unit='km/s') )

        for k in log_lw.keys():
            catalog.add_column( Column(name='log linewidth (mom2 '+k+')',   data=log_lw[k],   dtype=np.float64, unit='km/s') )


# linewidth containing 90% of all emission within structure contour
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        cube = fits.open(data[co][gal]['noise matched']['file'])[0]
        chans,ny,nx = cube.data.shape
        bunit = u.Unit(cube.header['bunit'])
        cunit = u.Unit(cube.header['cunit3'])
        cdelt = cube.header['cdelt3']*cunit
        crval = cube.header['crval3']*cunit
        crpix = cube.header['crpix3']
        pix_cen = np.array([((i+crpix-1)*cdelt+crval).value for i in np.arange(chans)])*cunit

        # mask cube
        min_snr = 5.            # check dendrogram threshold!
        cube.data[cube.data*bunit<min_snr*data[co][gal]['noise matched']['rms']] = np.nan

        def count_up(spectrum, level=0.9):
            allI = np.nansum(spec)
            peak = np.nanmax(spectrum)
            for i in np.linspace(10, peak, 100):
                sumI = np.nansum(spectrum[spectrum<i])
                if sumI/allI > 1.-level:
                    return i

        lw = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask    = struct.get_mask()
            mcube   = np.nan_to_num(cube.data*mask)*bunit

            # decrease threshold from top down in steps of 1K until 90% of the flux is captured
            spec = np.nansum(mcube, axis=(1,2)).value
            spec[spec==0.0] = np.nan
            trsh = count_up(spec, level=0.9)

            # get width at 90% level
            try:
                intersects = []
                for v,x1,x2 in zip(pix_cen[:-1],spec[:-1],spec[1:]):
                    if x1<trsh and x2>trsh:
                        intersects.append(v.value)
                    elif x1>trsh and x2<trsh:
                        intersects.append(v.value)
                width = np.max(intersects)-np.min(intersects)
            except:
                width = np.nan

            lw.append(width)

        log_lw = [np.log10(l) for l in lw]

        catalog.add_column( Column(name='linewidth (90% flux)',     data=lw,     dtype=np.float64, unit='km/s') )
        catalog.add_column( Column(name='log linewidth (90% flux)', data=log_lw, dtype=np.float64, unit='log10 km/s') )


# actual full width at half maximum
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        cube = fits.open(data[co][gal]['noise matched']['file'])[0]
        chans,ny,nx = cube.data.shape
        bunit = u.Unit(cube.header['bunit'])
        cunit = u.Unit(cube.header['cunit3'])
        cdelt = cube.header['cdelt3']*cunit
        crval = cube.header['crval3']*cunit
        crpix = cube.header['crpix3']
        pix_cen = np.array([((i+crpix-1)*cdelt+crval).value for i in np.arange(chans)])*cunit

        # mask cube
        min_snr = 5.            # check dendrogram threshold!
        cube.data[cube.data*bunit<min_snr*data[co][gal]['noise matched']['rms']] = np.nan

        lw = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask  = struct.get_mask()
            mcube = np.nan_to_num(cube.data*mask)*bunit
            spec  = np.nansum(mcube, axis=(1,2)).value
            hmax  = (spec.max()+spec.min())/2.

            # get min/max non-zero channel
            spec[spec<hmax] = 0.0
            nonzero = np.transpose(np.nonzero(spec))
            width = (pix_cen[np.nanmax(nonzero)]-pix_cen[np.nanmin(nonzero)]+cdelt).value       # add one channel because the methods cuts half a channel twice

            if width==0.:
                lw.append(np.nan)
            else:
                lw.append(width)

        log_lw = [np.log10(l) for l in lw]

        catalog.add_column( Column(name='linewidth (FWHM)',     data=lw,     dtype=np.float64, unit='km/s') )
        catalog.add_column( Column(name='log linewidth (FWHM)', data=log_lw, dtype=np.float64, unit='log10 km/s') )



# full width at 10% level
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        cube = fits.open(data[co][gal]['noise matched']['file'])[0]
        bunit = u.Unit(cube.header['bunit'])

        # mask cube
        min_snr = 5.            # check dendrogram threshold!
        cube.data[cube.data*bunit<min_snr*data[co][gal]['noise matched']['rms']] = np.nan

        lw = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask  = struct.get_mask()
            mcube = np.nan_to_num(cube.data*mask)*bunit
            spec  = np.nansum(mcube, axis=(1,2)).value
            p10   = 0.1*(spec.max()-spec.min())+spec.min()

            # get min/max non-zero channel
            spec[spec<p10] = 0.0
            nonzero = np.transpose(np.nonzero(spec))
            width = (pix_cen[np.nanmax(nonzero)]-pix_cen[np.nanmin(nonzero)]+cdelt).value       # add one channel because the methods cuts half a channel twice

            if width==0.:
                lw.append(np.nan)
            else:
                lw.append(width)

        log_lw = [np.log10(l) for l in lw]

        catalog.add_column( Column(name='linewidth (FW10%)',     data=lw,     dtype=np.float64, unit='km/s') )
        catalog.add_column( Column(name='log linewidth (FW10%)', data=log_lw, dtype=np.float64, unit='log10 km/s') )


###################################################################################################
# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
