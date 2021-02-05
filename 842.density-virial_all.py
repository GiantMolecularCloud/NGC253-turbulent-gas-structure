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
# plot density-virial relation
###################################################################################################

# separate figures
###################################################################################################

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        bins       = fits[co][gal]['density-virial']['size (astrodendro)']['linewidth (astrodendro)']['bins']

        fig = plt.figure(figsize=(4,4))
        ax  = fig.subplots(1)

        size      = catalog['size (astrodendro)']
        linewidth = catalog['linewidth (astrodendro)']
        s2R       = linewidth**2/size
        coldens   = catalog['column density']

        s2R_min = np.nanmin([x if x>1e-5 else np.nan for x in s2R])
        s2R_max = np.nanmax(s2R)
        cd_min  = np.nanmin([x if x>1e15 else np.nan for x in coldens])
        cd_max  = np.nanmax(coldens)

        i = 0
        for row in bins:
            window = [row['window min'], row['window max']]
            window_cen = row['window center']
            median = [row['y median'], row['y median']]
            perc16 = [row['y 16th'], row['y 16th']]
            perc84 = [row['y 84th'], row['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=3, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=co+' '+gal if i==0 else '')
            i += 1

        ax.set_xlabel(r'$\Sigma$ [cm$^{-2}$]')
        ax.set_ylabel(r'$\frac{\sigma^2}{R}$ [km$^2$\,s$^{-2}$\,pc$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(0.5*cd_min, 2.0*cd_max)
        ax.set_ylim(0.5*s2R_min, 2.0*s2R_max)
        for alpha in [1,10,100]:
            ax.plot([1e19, 1e24],
                    [(10/9*alpha*c.G*np.pi/5 *(2*u.u *r/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value for r in [1e19, 1e24]],
                    ls=':', lw=1, color='darkgrey', zorder=1
                   )
            ax.text(1e23,
                    1.7*(10/9*alpha*c.G*np.pi/5 *(2*u.u *1e23/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value,
                    r'$\alpha_\mathrm{vir} = '+str(alpha)+'$',
                    rotation = np.arctan((np.log10(ax.get_ylim()[1])-np.log10(ax.get_ylim()[0]))/(np.log10(ax.get_xlim()[1])-np.log10(ax.get_xlim()[0])))/2/np.pi*360,
                    color = 'darkgrey'
                   )
        ax.text(0.95,0.05, 'all sizes', color='darkgrey', transform=ax.transAxes, ha='right', va='bottom')
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')
        mkpath(plotdir, '05.coldens-sigma2size', 'binned_fit')
        fig.savefig(os.path.join(plotdir, '05.coldens-sigma2size', 'binned_fit', 'coldens-sigma2size_binned.'+co+'.'+gal+'.pdf'), dpi=300, bbox_inches='tight')


# single figure
###################################################################################################

fig = plt.figure(figsize=(4,4))
ax  = fig.subplots(1)

for co,CO in lines.items():
    for gal,GAL in galaxies.items():

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        color      = dendrograms[co][gal]['color']
        bins       = fits[co][gal]['density-virial']['size (astrodendro)']['linewidth (astrodendro)']['bins']

        i = 0
        for row in bins:
            window = [row['window min'], row['window max']]
            window_cen = row['window center']
            median = [row['y median'], row['y median']]
            perc16 = [row['y 16th'], row['y 16th']]
            perc84 = [row['y 84th'], row['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=3, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=co+' '+gal if i==0 else '')
            i += 1

ax.set_xlabel(r'$\Sigma$ [cm$^{-2}$]')
ax.set_ylabel(r'$\frac{\sigma^2}{R}$ [km$^2$\,s$^{-2}$\,pc$^{-1}$]')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(1e19, 1e24)
ax.set_ylim(1e-2, 1e3)
for alpha in [1,10,100]:
    ax.plot([1e19, 1e24],
            [(10/9*alpha*c.G*np.pi/5 *(2*u.u *r/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value for r in [1e19, 1e24]],
            ls=':', lw=1, color='darkgrey', zorder=1
           )
    ax.text(1e23,
            1.7*(10/9*alpha*c.G*np.pi/5 *(2*u.u *1e23/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value,
            r'$\alpha_\mathrm{vir} = '+str(alpha)+'$',
            rotation = np.arctan((np.log10(ax.get_ylim()[1])-np.log10(ax.get_ylim()[0]))/(np.log10(ax.get_xlim()[1])-np.log10(ax.get_xlim()[0])))/2/np.pi*360,
            color = 'darkgrey'
           )
ax.text(0.95,0.05, 'all sizes', color='darkgrey', transform=ax.transAxes, ha='right', va='bottom')
ax.set_axisbelow(True)
ax.grid(ls=':', c='lightgrey')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
mkpath(plotdir, '05.coldens-sigma2size', 'binned_fit')
fig.savefig(os.path.join(plotdir, '05.coldens-sigma2size', 'binned_fit', 'coldens-sigma2size_binned.all.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# density-virial relation as same size scale
###################################################################################################

def select_range(size, larger_than, smaller_than, *lists):
    selection = tuple([np.logical_and(size>larger_than, size<smaller_than)])
    selected_lists = [size[selection]]
    for l in lists:
        selected_lists.append(l[selection])
    return selected_lists

# size_bin = [8,10]
bounds = 2**np.arange(0,8)
for size_bin in tqdm(zip(bounds[:-1],bounds[1:])):

    fig = plt.figure(figsize=(4,4))
    ax  = fig.subplots(1)

    for co,CO in lines.items():
        for gal,GAL in galaxies.items():

            catalog    = dendrograms[co][gal]['catalog']
            dendrogram = dendrograms[co][gal]['dendrogram']
            color      = dendrograms[co][gal]['color']

            size      = catalog['size (astrodendro)'].data
            linewidth = catalog['linewidth (astrodendro)'].data
            s2R       = linewidth**2/size
            coldens   = catalog['column density'].data

            size, linewidth, s2R, coldens = select_range(size, size_bin[0], size_bin[1], linewidth, s2R, coldens)

            if not len(size)==0:
                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.25,
                                              xdata  = coldens,
                                              ydata  = s2R
                                             )

                i = 0
                for row in bins:
                    window = [row['window min'], row['window max']]
                    window_cen = row['window center']
                    median = [row['y median'], row['y median']]
                    perc16 = [row['y 16th'], row['y 16th']]
                    perc84 = [row['y 84th'], row['y 84th']]

                    ax.plot(window, median,                 color=color,                                zorder=3, label='')
                    ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=co+' '+gal if i==0 else '')
                    i += 1

    ax.set_xlabel(r'$\Sigma$ [cm$^{-2}$]')
    ax.set_ylabel(r'$\frac{\sigma^2}{R}$ [km$^2$\,s$^{-2}$\,pc$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(1e19, 1e24)
    ax.set_ylim(1e-2, 1e3)
    for alpha in [1,10,100]:
        ax.plot([1e19, 1e24],
                [(10/9*alpha*c.G*np.pi/5 *(2*u.u *r/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value for r in [1e19, 1e24]],
                ls=':', lw=1, color='darkgrey', zorder=1
               )
        ax.text(1e23,
                1.7*(10/9*alpha*c.G*np.pi/5 *(2*u.u *1e23/u.cm**2) ).to(u.km**2/u.s**2/u.pc).value,
                r'$\alpha_\mathrm{vir} = '+str(alpha)+'$',
                rotation = np.arctan((np.log10(ax.get_ylim()[1])-np.log10(ax.get_ylim()[0]))/(np.log10(ax.get_xlim()[1])-np.log10(ax.get_xlim()[0])))/2/np.pi*360,
                color = 'darkgrey'
               )
    ax.text(0.95,0.05, 'size bin '+str(size_bin[0])+'--'+str(size_bin[1])+r'\,pc', color='darkgrey', transform=ax.transAxes, ha='right', va='bottom')
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
    ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=2, fancybox=True, fontsize=10)
    mkpath(plotdir, '05.coldens-sigma2size', 'binned_fit')
    fig.savefig(os.path.join(plotdir, '05.coldens-sigma2size', 'binned_fit', 'coldens-sigma2size_binned.all.'+str(size_bin[0])+'-'+str(size_bin[1])+'pc.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
