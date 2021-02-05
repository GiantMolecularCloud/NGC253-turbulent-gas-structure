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
# get column density
###################################################################################################

# integrated flux
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        try:
            del catalog['integrated Tb']
        except:
            pass

        chanwidth = CO['resolution']['spectral']

        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        flux = []
        for struct in tqdm(all_structs_ordered):
            flux.append( struct.values().sum()*chanwidth.value )

        # add column to catalog table
        catalog.add_column( Column(name='integrated Tb', data=flux, dtype=np.float64, unit='K km/s') )


# structure area
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        Apix = ((CO['resolution']['spatial']/5.)**2)

        npix = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]
        for struct in tqdm(all_structs_ordered):
            mask = struct.get_mask()
            n = np.sum(mask, axis=0)
            n[n>0] = 1.
            n = np.sum(n)

            npix.append(n)

        area_projected = npix*Apix

        catalog.add_column( Column(name='npix (projected)', data=npix,           dtype=np.int64) )
        catalog.add_column( Column(name='area (projected)', data=area_projected, dtype=np.float64, unit='pc^2') )


# effective pixel area
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        Apix = ((CO['resolution']['spatial']/5.)**2)

        npix_exact     = catalog['area (exact)']/Apix
        npix_ellipse   = catalog['area (ellipse)']/Apix
        npix_effective = catalog['area (effective)']/Apix

        catalog.add_column( Column(name='npix (exact)', data=npix_exact, dtype=np.float64) )
        catalog.add_column( Column(name='npix (ellipse)', data=npix_ellipse, dtype=np.float64) )
        catalog.add_column( Column(name='npix (effective)', data=npix_effective, dtype=np.float64) )

for i in ['projected','exact','ellipse','effective']:
    print(np.percentile(dendrograms['CO(1-0)']['NGC253']['catalog']['npix ('+i+')'], (1,50,99)))


# colunm density
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        Xco  = GAL['Xco'] *u.cm**-2 /(u.K*u.km/u.s)

        catalog['integrated Tb']*Xco/CO['excitation correction']

        # factor 1.36 due to Helium contribution
        coldens_projected = catalog['integrated Tb']*1.36*Xco/CO['excitation correction'] /catalog['npix (projected)']
        coldens_exact     = catalog['integrated Tb']*1.36*Xco/CO['excitation correction'] /catalog['npix (exact)']
        coldens_ellipse   = catalog['integrated Tb']*1.36*Xco/CO['excitation correction'] /catalog['npix (ellipse)']
        coldens_effective = catalog['integrated Tb']*1.36*Xco/CO['excitation correction'] /catalog['npix (effective)']

        for col in ['projected','exact','ellipse','effective']:
            try:
                del catalog['column density ('+col+')']
            except:
                pass

        catalog.add_column( Column(name='column density (projected)',   data=coldens_projected, dtype=np.float64, unit='cm^-2') )
        catalog.add_column( Column(name='column density (exact)',       data=coldens_exact,     dtype=np.float64, unit='cm^-2') )
        catalog.add_column( Column(name='column density (ellipse)',     data=coldens_ellipse,   dtype=np.float64, unit='cm^-2') )
        catalog.add_column( Column(name='column density (effective)',   data=coldens_effective, dtype=np.float64, unit='cm^-2') )


for i in ['projected','exact','ellipse','effective']:
    print(np.percentile(dendrograms['CO(1-0)']['NGC253']['catalog']['column density ('+i+')'], (1,50,99)))
for i in ['projected','exact','ellipse','effective']:
    print(np.percentile(dendrograms['CO(3-2)']['NGC253']['catalog']['column density ('+i+')'], (1,50,99)))



# mass
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        Apix = ((CO['resolution']['spatial']/5.)**2).to(u.cm**2)
        Xco  = GAL['Xco'] *u.cm**-2 /(u.K*u.km/u.s)

        # atomic weight of H2: 2; 1.36 to account for helium
        mass = ((catalog['integrated Tb']*Xco/CO['excitation correction'] *Apix *1.36*2.0*u.u).to(u.Msun)).value

        catalog.add_column( Column(name='mass', data=mass, dtype=np.float64, unit='Msun') )

        # log mass
        mass[(mass<1e0) & ~np.isfinite(mass)] = np.nan
        log_mass = np.log10(mass)

        catalog.add_column( Column(name='log mass', data=log_mass, dtype=np.float64, unit='log Msun') )


###################################################################################################
# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
