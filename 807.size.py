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
# get size
###################################################################################################

# astrodendro size
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        size     = [angle_to_parsec(r, source=gal).value for r in catalog['radius'].quantity]
        log_size = [np.log10(s) for s in size]

        catalog.add_column( Column(name='size (astrodendro)',     data=size,     dtype=np.float64, unit='pc') )
        catalog.add_column( Column(name='log size (astrodendro)', data=log_size, dtype=np.float64, unit='log10 pc') )


# astrodendro area_exact size
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        area     = catalog['area_exact'].quantity
        size     = [angle_to_parsec(np.sqrt(a/np.pi), source=gal).value for a in area]
        log_size = [np.log10(s) for s in size]

        catalog.add_column( Column(name='size (area_exact)',     data=size,     dtype=np.float64, unit='pc') )
        catalog.add_column( Column(name='log size (area_exact)', data=log_size, dtype=np.float64, unit='log10 pc') )


# astrodendro area_ellipse size
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        area     = catalog['area_ellipse'].quantity
        size     = [angle_to_parsec(np.sqrt(a/np.pi), source=gal).value for a in area]
        log_size = [np.log10(s) for s in size]

        catalog.add_column( Column(name='size (area_ellipse)',     data=size,     dtype=np.float64, unit='pc') )
        catalog.add_column( Column(name='log size (area_ellipse)', data=log_size, dtype=np.float64, unit='log10 pc') )



# manual size
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        # cube = fits.open(data[co][gal]['noise matched']['file'])[0]
        # pix  = u.Quantity(str(np.abs(cube.header['cdelt1']))+cube.header['cunit1'])
        pix  = CO['resolution']['spatial']/5.

        size = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask = struct.get_mask()
            mom0 = np.sum(mask,axis=0)
            mom0[mom0>0] = 1
            A = np.sum(mom0)                                # identical to area_exact
            r = np.sqrt(A/np.pi)*pix.to(u.pc).value
            size.append(r)
        log_size = [np.log10(s) for s in size]

        catalog.add_column( Column(name='size (manual)',     data=size,     dtype=np.float64, unit='pc') )
        catalog.add_column( Column(name='log size (manual)', data=log_size, dtype=np.float64, unit='log10 pc') )


# effective size
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        major_sigma = angle_to_parsec(catalog['major_sigma'].quantity, source=gal)
        minor_sigma = angle_to_parsec(catalog['minor_sigma'].quantity, source=gal)
        size        = ((major_sigma+minor_sigma)/2.).value
        log_size    = [np.log10(s) for s in size]

        catalog.add_column( Column(name='size (effective)',     data=size,     dtype=np.float64, unit='pc') )
        catalog.add_column( Column(name='log size (effective)', data=log_size, dtype=np.float64, unit='log10 pc') )


# area
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        cat_area_exact = catalog['area_exact'].quantity
        cat_area_ellipse = catalog['area_ellipse'].quantity
        major_sigma = catalog['major_sigma'].quantity
        minor_sigma = catalog['minor_sigma'].quantity

        area_dendro    = np.pi*catalog['size (astrodendro)']**2
        area_exact     = angle_to_parsec(np.sqrt(cat_area_exact), source=gal)**2
        area_ellipse   = angle_to_parsec(np.sqrt(cat_area_ellipse), source=gal)**2
        area_effective = np.pi*angle_to_parsec(major_sigma, source=gal)*angle_to_parsec(minor_sigma, source=gal)

        try:
            del catalog['area (astrodendro)']
            del catalog['area (exact)']
            del catalog['area (ellipse)']
            del catalog['area (effective)']
        except:
            pass

        catalog.add_column( Column(name='area (astrodendro)', data=area_exact.value,     dtype=np.float64, unit='pc2') )
        catalog.add_column( Column(name='area (exact)',       data=area_exact.value,     dtype=np.float64, unit='pc2') )
        catalog.add_column( Column(name='area (ellipse)',     data=area_ellipse.value,   dtype=np.float64, unit='pc2') )
        catalog.add_column( Column(name='area (effective)',   data=area_effective.value, dtype=np.float64, unit='pc2') )


# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
