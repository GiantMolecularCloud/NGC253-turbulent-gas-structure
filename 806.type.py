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
# get structure type
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        from astropy.table.column import Column

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']

        stypes = []
        idx_list = [i.idx for i in dendrogram.all_structures]
        all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(dendrogram.all_structures)))]

        for struct in tqdm(all_structs_ordered):
            if struct.is_leaf:
                stypes.append('leaf')
            elif struct.is_branch:
                if struct.parent==None:
                    stypes.append('trunk')
                else:
                    stypes.append('branch')
            else:
                stypes.append('unknown')

        if 'type' in catalog.colnames:
            catalog.remove_column('type')
        catalog.add_column( Column(name='type', data=stypes, dtype=str) )


# save catalog
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():
        fnpickle(dendrograms[co][gal]['catalog'], join(compdir,gal+'.'+co+'.catalog.pickle'))


###################################################################################################
#
###################################################################################################
