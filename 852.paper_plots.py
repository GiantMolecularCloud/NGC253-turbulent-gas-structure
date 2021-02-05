##############
# DENDROGRAM #
##############

###################################################################################################
# load data
###################################################################################################

# execfile(join(scriptdir, '800.info.py'))
# data        = fnunpickle(join(datadir, 'data.pickle'))
# dendrograms = load_dendrograms()

gal_colors = ['dodgerblue','darkorange']


###################################################################################################
# load paperfits
###################################################################################################

# execfile(join(scriptdir, '830.MCMC.py'))
# paperfits = {co: {gal: {} for gal in galaxies.keys()} for co in lines.keys()}
# for co,CO in lines.items():
#     for gal,GAL in galaxies.items():
#         paperfits[co][gal] = fnunpickle(join(compdir,gal+'.'+co+'.paperfits.pickle'))


###################################################################################################
# table 1: dataset info
###################################################################################################

# done manually


###################################################################################################
# table 2: size-linewidth fits
###################################################################################################

table = r"""\begin{deluxetable*}{lcrrrcccc}
	\tablewidth{\textwidth}
	\tablecaption{Dendrogram statistics and fit results for power law fits to the binned size--linewidth and size--mass relations shown in Figure~\ref{fig1} and \ref{fig3}.}
	\label{table2}
	\tablehead{\colhead{source} & \colhead{line} & \multicolumn{3}{c}{dendrogram structures} & \multicolumn{2}{c}{size--linewidth} & \multicolumn{2}{c}{size--mass}\\
	&& total & branches & leaves & $b$ & \sigmaten & $d$ & log $\mten$\\\\
    &&&&& (1) & (2) & (3) & (4)
	}
	\startdata
"""
for co in lines:
     for gal in galaxies:
         catalog    = dendrograms[co][gal]['catalog']
         dendrogram = dendrograms[co][gal]['dendrogram']
         fit  = paperfits[co][gal]['size-linewidth']['binned fit']
         # fit2 = paperfits[co][gal]['mass-linewidth']['binned fit']
         fit2 = paperfits[co][gal]['size-mass']['binned fit']

         total    = len(catalog)
         branches = len(catalog[catalog['type']=='branch'])
         leaves   = len(catalog[catalog['type']=='leaf'])

         b     = fit.a['ls'][0]
         b_err = fit.a['ls'][1]
         s     = fit.y10['ls'][0]
         s_err = fit.y10['ls'][1]

         b2     = fit2.a['ls'][0]
         b2_err = fit2.a['ls'][1]
         # a2     = np.power(10,fit2.c['ls'][0])
         # a2_err = np.power(10,fit2.c['ls'][1])*np.log(10)
         a2     = fit2.y10['ls'][0]
         a2_err = fit2.y10['ls'][1]
         loga2  = np.log10(a2)
         loga2_err = np.log10(1+a2_err/a2)

         table += ('{:<6}'.format(gal)+' & '+
                   '{:<7}'.format(co)+' & '+
                   '{:>5d}'.format(total)+' & '+
                   '{:>4d}'.format(branches)+' & '+
                   '{:>4d}'.format(leaves)+' & '+
                   r'${:4.2f}\pm{:.2f}$'.format(b,b_err)+' & '+
                   r'${:4.1f}\pm{:.1f}$'.format(s,s_err)+' & '+
                   r'${:4.2f}\pm{:.2f}$'.format(b2,b2_err)+' & '+
                   # r'${:4.1f}\pm{:.1f}$'.format(a2,a2_err)+r' \\'+'\n'
                   r'${:4.2f}\pm{:.2f}$'.format(loga2,loga2_err)+r' \\'+'\n'
                  )
table += r"""    \enddata
    \tablecomments{(1) Exponent $b$ of the power law fit to the size--linewidth relation according to Equation~\ref{equation: size-linewidth}.
                   (2) Characteristic linewidth at 10\,pc according to the power law fit to the size--linewidth relation (Equation~\ref{equation: size-linewidth}) in \kms.
                   (3) Exponent $d$ of the power law fit to the size--mass relation according to Equation~\ref{equation: size-mass}.
                   (4) Characteristic mass at 10\,pc according to Equation~\ref{equation: size-mass} in $\log \mathrm{M}_\odot$.}
\end{deluxetable*}
"""
print(table)


###################################################################################################
# fig 1: size-linewidth binned
###################################################################################################

vol_thresh = {'CO(3-2)': [[2.4e0,1e-1],[2e-1,9e0]],
              'CO(1-0)': [[3e1,1e-1],[6e-1,3e2]]}

fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='none', sharey='none', figsize=(8,5))
fig.subplots_adjust(hspace=0., wspace=0.)

for (co,CO),ax in zip(reversed(list(lines.items())),axes):
    for idx,(gal,color) in enumerate(zip(reversed(list(galaxies.keys())),gal_colors)):

        fit        = paperfits[co][gal]['size-linewidth']['binned fit']
        bins       = paperfits[co][gal]['size-linewidth']['bins']

        i = 0
        for bin in bins:
            window = [bin['window min'], bin['window max']]
            window_cen = bin['window center']
            median = [bin['y median'], bin['y median']]
            perc16 = [bin['y 16th'], bin['y 16th']]
            perc84 = [bin['y 84th'], bin['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=3, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=co+' '+gal if i==0 else '')
            i += 1

        x = np.logspace(np.log10(np.min(bins['window min'])), np.log10(np.max(bins['window max'])), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=4)

        ax.text(0.05,0.85+idx*0.05, CO['tex']+' '+gal, color=color, transform=ax.transAxes, ha='left', va='bottom')

        vth = vol_thresh[co]
        ax.fill_between([vth[0][0],vth[1][0]],[vth[0][1],vth[1][1]], [vth[0][1],vth[0][1]], facecolor='lightgrey', edgecolor=None, alpha=0.5, zorder=1)
        ax.plot([vth[0][0],vth[1][0]],[vth[0][1],vth[1][1]], ls='-', color='darkgrey', alpha=0.5, zorder=1)

    ax.set_xlabel('R [pc]')
    ax.set_ylabel('$\sigma$ [km\,s$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(6e-1, 8e1)
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
axes[0].set_xlim(3e-1, 5e1)
axes[1].set_xlim(3e0, 5e2)
axes[1].yaxis.set_label_position('right')
axes[1].tick_params(axis='y', which='both', labelleft='off', labelright='on')
fig.savefig(join(plotdir, 'paper', 'fig2.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fig 2: size-luminosity binned
###################################################################################################

fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='none', sharey='none', figsize=(8,5))
fig.subplots_adjust(hspace=0., wspace=0.)

for (co,CO),ax in zip(reversed(list(lines.items())),axes):
    for idx,(gal,color) in enumerate(zip(reversed(list(galaxies.keys())),gal_colors)):

        fit        = paperfits[co][gal]['size-luminosity']['binned fit']
        bins       = paperfits[co][gal]['size-luminosity']['bins']

        i = 0
        for bin in bins:
            window = [bin['window min'], bin['window max']]
            window_cen = bin['window center']
            median = [bin['y median'], bin['y median']]
            perc16 = [bin['y 16th'], bin['y 16th']]
            perc84 = [bin['y 84th'], bin['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=2, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=1, label=co+' '+gal if i==0 else '')
            i += 1

        x = np.logspace(np.log10(np.min(bins['window min'])), np.log10(np.max(bins['window max'])), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

        ax.text(0.05,0.85+idx*0.05, CO['tex']+' '+gal, color=color, transform=ax.transAxes, ha='left', va='bottom')

    x = np.logspace(-1,3,100)
    for v,l,a in [[2,'N',38],[3,r'\rho',47]]:
        for i in [-1,1,3,5]:
            y = 10**i *x**v
            line = ax.plot(x, y,
                           ls=':', lw=1, color='darkgrey', zorder=1
                          )
            if ax==axes[1]:
                label_line(ax, line, 4e2, 1.1* 10**i *4e2**v,
                           r'$'+l+' = const$',
                           angle=a,
                           **{'color': 'darkgrey', 'ha':'right'}
                          )

    ax.set_xlabel(r'R [pc]')
    ax.set_ylabel(r'L [K\,km\,s$^{-1}$\,pc$^2$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(8e0, 1.5e8)
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
axes[0].set_xlim(3e-1, 5e1)
axes[1].set_xlim(3e0, 5e2)
axes[1].yaxis.set_label_position('right')
axes[1].tick_params(axis='y', which='both', labelleft='off', labelright='on')
fig.savefig(join(plotdir, 'paper', 'fig3.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# figure 3: density-virial binned
###################################################################################################

def alpha_x_to_y(x):
        return (5/3*alpha*c.G*np.pi/5 *(1.36*2.*u.u *x) ).to(u.km**2/u.s**2/u.pc)

def Pek_x_to_y(x, Pek, Gamma=3./5.):
    return ( 1./3. *( np.pi*Gamma*c.G*x*1.36*2.*u.u + 4*c.k_B/(x*1.36*2.*u.u)*Pek) ).to(u.km**2 /u.s**2 /u.pc)


fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(8,4))
fig.subplots_adjust(hspace=0., wspace=0.)

# format axes
for ax in axes:
    ax.set_xlim(1.5e20, 2.0e24)
    ax.set_ylim(8e-2, 3e2)
    ax.set_xlabel(r'$N$ [cm$^{-2}$]')
    ax.set_ylabel(r'$\frac{\sigma^2}{R}$ [km$^2$\,s$^{-2}$\,pc$^{-1}$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')

# alpha lines
for ax in axes:
    for alpha in [1,10,100]:
        # plot alpha lines
        x = np.array([1e19, 1e25])*u.cm**-2
        y = alpha_x_to_y(x)
        line = ax.plot(x.value, y.value,
                       ls=':', lw=1, color='darkgrey', zorder=1
                      )
        # label only second panel
        if ax==axes[1]:
            x = 1.2e24*u.cm**-2
            y = 0.5*alpha_x_to_y(x)
            label_line(ax, line, x.value, y.value,
                       r'$\alpha_\mathrm{vir} = '+str(alpha)+'$',
                       angle=0.,
                       **{'color': 'darkgrey', 'ha':'right'}
                      )

# Pe/k lines
for ax in axes:
    for Pek in np.array([1e4,1e5,1e6,1e7])*u.K/u.cm**3:
        # plot Pe/k lines
        x = np.logspace(np.log10(1e20), np.log10(5e24),100)*u.cm**-2
        y = Pek_x_to_y(x, Pek)
        line = ax.plot(x.value, y.value,
                       ls='--', lw=1, color='darkgrey', zorder=1
                      )
        # label only second panel
        if ax==axes[1]:
            x = 1e21*u.cm**-2
            y = Pek_x_to_y(x, Pek)
            ax.text(x.value, y.value,
                    'P/k\,=\,$10^{'+str(int(np.log10(Pek.value)))+r'}$\,K\,cm$^{\textnormal{-} 3}$',
                    color = 'grey',
                    ha    = 'left',
                    fontsize = 8
                   )

# plot data
for col,co in enumerate(reversed(list(lines.keys()))):
    ax = axes[col]
    for gal,color in zip(reversed(list(galaxies.keys())),gal_colors):

        catalog    = dendrograms[co][gal]['catalog']
        dendrogram = dendrograms[co][gal]['dendrogram']
        size      = catalog['size (astrodendro)'].data
        linewidth = catalog['linewidth (astrodendro)'].data
        s2R       = linewidth**2/size
        coldens   = catalog['column density (effective)'].data

        if not len(size)==0:
            bins = get_binned_percentiles(R_measure, lw_measure,
                                          x_step = 0.25,
                                          xdata  = coldens,
                                          ydata  = s2R
                                         )

            i = 0
            for bin in bins:
                window = [bin['window min'], bin['window max']]
                window_cen = bin['window center']
                median = [bin['y median'], bin['y median']]
                perc16 = [bin['y 16th'], bin['y 16th']]
                perc84 = [bin['y 84th'], bin['y 84th']]

                ax.plot(window, median,                 color=color,                                zorder=3, label='')
                ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=gal if i==0 else '')
                i += 1

# finish formatting
axes[0].tick_params(axis='y', left='on', top='off', right='off', bottom='off', labelleft='on', labeltop='off', labelright='off', labelbottom='off')
axes[1].tick_params(axis='y', left='off', top='off', right='on', bottom='off', labelleft='off', labeltop='off', labelright='on', labelbottom='off')
axes[1].yaxis.set_label_position('right')
axes[0].text(0.5,1.02, 'CO(3--2)', color='k', transform=axes[0].transAxes, ha='center', va='bottom', fontsize=14)
axes[1].text(0.5,1.02, 'CO(1--0)', color='k', transform=axes[1].transAxes, ha='center', va='bottom', fontsize=14)
axes[0].legend(bbox_to_anchor=(0.52,0.08,0.4,0.4), loc='lower right', borderaxespad=0., ncol=1, fancybox=True, fontsize=10)
fig.savefig(join(plotdir, 'paper', 'fig4.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fig A.1: example structure (show size definition)
###################################################################################################

def moment0(cube, mask_zeros=True, return_HDU=True):

    # get info
    bunit = u.Unit(cube.header['bunit'])
    cunit = u.Unit(cube.header['cunit3'])
    cdelt = cube.header['cdelt3']*cunit

    # calculate moment
    mom0 = cdelt*np.nansum(cube.data, axis=0)*bunit
    if mask_zeros:
        mom0[mom0==0.0] = np.nan

    # construct header
    mom_head = copy.deepcopy(cube.header)
    mom_head['naxis'] = 2
    for k in cube.header.keys():
        if '3' in k:
            del mom_head[k]

    if return_HDU:
        mom0_hdu = fits.PrimaryHDU(mom0.astype('short').value, mom_head)
        return mom0_hdu
    else:
        return mom0, mom_head


# dataset to use:
###################################################################################################

co  = 'CO(3-2)'
gal = 'NGC253'
idx = {'idx': 4822, 'vmin': 5000, 'vmax': 12000}


# get moment 0 map of matched image
###################################################################################################

cube = fits.open(data[co][gal]['noise matched']['file'])[0]

if not os.path.exists(join(subprojectdir,'paper',co+'.'+gal+'.moment_0.fits')):

    # mask cube as in dendrogram calculation
    min_snr = 5.
    cube.data[cube.data*bunit<min_snr*data[co][gal]['noise matched']['rms']] = np.nan

    # get moment 0
    mom0 = moment0(cube)
    mom0.writeto(join(subprojectdir,'paper',co+'.'+gal+'.moment_0.fits'))


# select structures and get masks
###################################################################################################

if not os.path.exists(join(subprojectdir,'paper',co+'.'+gal+'.struct_'+str(idx['idx'])+'.moment_0.fits')):
    structure   = dendrograms[co][gal]['dendrogram'][idx['idx']]
    struct_mask = structure_08.get_mask()
    struct_mom0 = moment0( fits.PrimaryHDU(struct_mask.astype('short'), cube.header) )
    struct_mom0.writeto(join(subprojectdir,'paper',co+'.'+gal+'.struct_'+str(idx['idx'])+'.moment_0.fits'))
else:
    struct_mom0 = fits.open(join(subprojectdir,'paper',co+'.'+gal+'.struct_'+str(idx['idx'])+'.moment_0.fits'))[0]


# plot map and structures
###################################################################################################

# figure
fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(5,3))

# labels
ax.set_xlabel('galactic longitude [pc]')
ax.set_ylabel('galactic latitude [pc]')

xaxis = get_pixel_locations(mom0, axis=1).value
yaxis = get_pixel_locations(mom0, axis=2).value
xlabels = [int(i) for i in np.arange(np.max(np.round(xaxis,-1)), np.min(np.round(xaxis,-1)), -2)]
ylabels = [int(i) for i in np.arange(np.min(np.round(yaxis,-1)), np.max(np.round(yaxis,-1)), 2)]
xticks  = coordinate_to_pixel(mom0, 1, xlabels*u.pc).value
yticks  = coordinate_to_pixel(mom0, 2, ylabels*u.pc).value

ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_xticklabels(xlabels);
ax.set_yticklabels(ylabels);

# show image
im = ax.imshow(mom0.data,
               origin        = 'lower',
               interpolation = 'nearest',
               cmap          = plt.cm.Blues,
               aspect        = 'equal',
               vmin          = idx['vmin'],
               vmax          = idx['vmax']
              )
cbar = fig.colorbar(im)
cbar.set_label('intensity [K\,km\,s$^{-1}$]')

# show contour
d = dendrograms[co][gal]['dendrogram']
p = d.plotter()
p.plot_contour(ax, structure=idx['idx'], lw=1, colors='orange')

# show ellipse
idx_list = [i.idx for i in d.all_structures]
all_structs_ordered = [x for _,x in sorted(zip(idx_list,list(d.all_structures)))]
s = PPVStatistic(all_structs_ordered[idx['idx']])
ellipse = s.to_mpl_ellipse(edgecolor='red', facecolor='none')
ax.add_patch(ellipse)

# zoom in on structure
center   = ellipse.get_center()
vertices = ellipse.get_patch_transform().transform( ellipse.get_path().vertices.copy() )
bounds   = (np.min(vertices[:,0]), np.max(vertices[:,0]), np.min(vertices[:,1]), np.max(vertices[:,1]))
extent   = (bounds[1]-bounds[0], bounds[3]-bounds[2])

xmin = center[0]-2.5*extent[0]
xmax = center[0]+2.5*extent[0]
ymin = center[1]-2.5*extent[1]
ymax = center[1]+2.5*extent[1]
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)

# show half axes
dist = np.sqrt( (vertices[:,0]-center[0])**2 + (vertices[:,1]-center[1])**2 )
p_maj = vertices[np.argmax(dist)+1]
p_min = vertices[np.argmin(dist)]

a_maj = ax.plot((center[0],p_maj[0]), (center[1],p_maj[1]), c='red', ls='--', lw=0.5)
a_min = ax.plot((center[0],p_min[0]), (center[1],p_min[1]), c='red', ls='--', lw=0.5)

# plot beam
from matplotlib.patches import Ellipse
bmaj_pix = ( angle_to_parsec(cube.header['bmaj']*u.deg, source=gal) / cdelt ).value
bmin_pix = ( angle_to_parsec(cube.header['bmin']*u.deg, source=gal) / cdelt ).value
bpa      = cube.header['bpa']
beam = Ellipse(xy     = (xmin+0.15*(xmax-xmin), ymin+0.15*(xmax-xmin)),
               width  = bmaj_pix,
               height = bmin_pix,
               angle  = bpa,
               ec     = None,
               fc     = 'black',
               alpha  = 1.
              )
ax.add_artist(beam)

# save figure
fig.savefig(join(plotdir, 'paper', 'fig1.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fig A.2: structure definition comparison
###################################################################################################

catalog    = dendrograms['CO(3-2)']['NGC253']['catalog']
dendrogram = dendrograms['CO(3-2)']['NGC253']['dendrogram']

fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='none', sharey='none', figsize=(8,4))
fig.subplots_adjust(hspace=0., wspace=0.1)

R_measures = ['size (astrodendro)',
              'size (area_ellipse)',
              'size (manual)']
lw_measures = ['linewidth (astrodendro)',
               'linewidth (mom2 mean)',
               'linewidth (mom2 median)',
               'linewidth (90% flux)',
               'linewidth (FWHM)',
               'linewidth (FW10%)']


# A1: size definition comparison
###################################################################################################

ax = axes[0]
A1colors = mpl.cm.YlOrRd(np.linspace(0.1,0.9,len(R_measures)))

min = np.nanmin([x if x!=0. else np.nan for x in flatten([catalog[R_measure].data for R_measure in R_measures])])
max = np.nanmax([catalog[R_measure].data for R_measure in R_measures])

for R_measure,color in zip(R_measures,A1colors):
    size_ad = catalog['size (astrodendro)'].data
    size    = catalog[R_measure].data
    ax.scatter(size_ad, size, marker='.', s=8, c=[color],
               label=R_measure.replace('size (','').replace(')','').replace('_',' ').replace('astrodendro',r'R$_\mathrm{astrodendro}$').replace('manual',r'R$_\mathrm{circular}$').replace('area ellipse',r'R$_\mathrm{ellipse}$'),
               zorder=3,
               rasterized=True)
ax.plot([0.5*min,2.0*max],[0.5*min,2.0*max], ls='-', color='grey', lw=1, zorder=2)

ax.set_xlabel(r'R$_\mathrm{astrodendro}$ [pc]')
ax.set_ylabel(r'R [pc]')
ax.set_xlim(0.5*min, 2.0*max)
ax.set_ylim(0.5*min, 2.0*max)
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=3, scatterpoints=1, handletextpad=0., fancybox=True, fontsize=10)


# A2: linewidth definition comparison
###################################################################################################

ax = axes[1]
A2colors = mpl.cm.YlGnBu(np.linspace(0.1,1,len(lw_measures)))

min = np.nanmin([x if x>1e-2 else np.nan for x in flatten([catalog[lw_measure].data for lw_measure in lw_measures])])
max = np.nanmax([catalog[lw_measure].data for lw_measure in lw_measures])

for lw_measure,color in zip(lw_measures,A2colors):
    linewidth_ad = catalog['linewidth (astrodendro)'].data
    linewidth    = catalog[lw_measure].data
    ax.scatter(linewidth_ad, linewidth, marker='.', s=8, c=[color],
               label=lw_measure.replace('linewidth (','').replace(')','').replace('_',' ').replace('astrodendro',r'$\sigma_\mathrm{astrodendro}$').replace('mom2 mean',r'$\sigma_\mathrm{mom2\ mean}$').replace('mom2 median',r'$\sigma_\mathrm{mom2\ median}$').replace('90% flux',r'$\sigma_\mathrm{90\%}$').replace('FWHM',r'$\sigma_\mathrm{FWHM}$').replace('FW10%',r'$\sigma_\mathrm{FW10\%}$'),
               zorder=3,
               rasterized=True)
ax.plot([0.5*min,2.0*max],[0.5*min,2.0*max], ls='-', color='grey', lw=1, zorder=2)

ax.set_xlabel(r'$\sigma_\mathrm{astrodendro}$ [km\,s$^{-1}$]')
ax.set_ylabel(r'$\sigma$ [km\,s$^{-1}$]')
ax.set_xlim(0.5*min, 2.0*max)
ax.set_ylim(0.5*min, 2.0*max)
ax.yaxis.set_label_position('right')
ax.tick_params(axis='y', which='both', labelleft='off', labelright='on')
ax.legend(bbox_to_anchor=(0.,1.02,1.,0.05), loc='lower left', mode='expand', borderaxespad=0., ncol=3, scatterpoints=1, handletextpad=0., fancybox=True, fontsize=10)


# save figure
for ax in axes:
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
fig.savefig(join(plotdir, 'paper', 'figA.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fig B: size mass relation
###################################################################################################

fig,axes = plt.subplots(nrows=1, ncols=2, squeeze=True, sharex='none', sharey='none', figsize=(8,5))
fig.subplots_adjust(hspace=0., wspace=0.)

for (co,CO),ax in zip(reversed(list(lines.items())),axes):
    for idx,(gal,color) in enumerate(zip(reversed(list(galaxies.keys())),gal_colors)):

        fit        = paperfits[co][gal]['size-mass']['binned fit']
        bins       = paperfits[co][gal]['size-mass']['bins']

        i = 0
        for bin in bins:
            window = [bin['window min'], bin['window max']]
            window_cen = bin['window center']
            median = [bin['y median'], bin['y median']]
            perc16 = [bin['y 16th'], bin['y 16th']]
            perc84 = [bin['y 84th'], bin['y 84th']]

            ax.plot(window, median,                 color=color,                                zorder=2, label='')
            ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=1, label=co+' '+gal if i==0 else '')
            i += 1

        x = np.logspace(np.log10(np.min(bins['window min'])), np.log10(np.max(bins['window max'])), 100)
        ax.plot(x, powlaw(x,fit.a['ls'][0],np.power(10,fit.c['ls'][0])), color=color, zorder=3)

        ax.text(0.05,0.85+idx*0.05, CO['tex']+' '+gal, color=color, transform=ax.transAxes, ha='left', va='bottom')

    x = np.logspace(-1,3,100)
    for v,l,a in [[2,'N',38],[3,r'\rho',47]]:
        for i in [-1,1,3,5]:
            y = 10**i *x**v
            line = ax.plot(x, y,
                           ls=':', lw=1, color='darkgrey', zorder=1
                          )
            if ax==axes[1]:
                label_line(ax, line, 4e2, 1.1* 10**i *4e2**v,
                           r'$'+l+' = const$',
                           angle=a,
                           **{'color': 'darkgrey', 'ha':'right'}
                          )

    ax.set_xlabel(r'R [pc]')
    ax.set_ylabel(r'M$_\mathrm{lum}$ [M$_\odot$]')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim(8e0, 1.5e8)
    ax.set_axisbelow(True)
    ax.grid(ls=':', c='lightgrey')
axes[0].set_xlim(3e-1, 5e1)
axes[1].set_xlim(3e0, 5e2)
axes[1].yaxis.set_label_position('right')
axes[1].tick_params(axis='y', which='both', labelleft='off', labelright='on')
fig.savefig(join(plotdir, 'paper', 'figB.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# fig C: size splitted virial state
###################################################################################################

def select_range(size, larger_than, smaller_than, *lists):
    selection = tuple([np.logical_and(size>larger_than, size<smaller_than)])
    selected_lists = [size[selection]]
    for l in lists:
        selected_lists.append(l[selection])
    return selected_lists


fig,axes = plt.subplots(nrows=6, ncols=2, squeeze=True, sharex='col', sharey='row', figsize=(6,10))
fig.subplots_adjust(hspace=0., wspace=0.)

bounds = 2.**np.arange(-1,6)
for row,size_bin in tqdm(enumerate(zip(bounds[:-1],bounds[1:]))):
    for col,co in enumerate(lines.keys()):
        ax = axes[row,col]

        for gal,color in zip(reversed(list(galaxies.keys())),gal_colors):

            catalog    = dendrograms[co][gal]['catalog']
            dendrogram = dendrograms[co][gal]['dendrogram']
            size      = catalog['size (astrodendro)'].data
            linewidth = catalog['linewidth (astrodendro)'].data
            s2R       = linewidth**2/size
            coldens   = catalog['column density (effective)'].data
            size, linewidth, s2R, coldens = select_range(size, size_bin[0], size_bin[1], linewidth, s2R, coldens)

            if not len(size)==0:
                bins = get_binned_percentiles(R_measure, lw_measure,
                                              x_step = 0.25,
                                              xdata  = coldens,
                                              ydata  = s2R
                                             )

                i = 0
                for bin in bins:
                    window = [bin['window min'], bin['window max']]
                    window_cen = bin['window center']
                    median = [bin['y median'], bin['y median']]
                    perc16 = [bin['y 16th'], bin['y 16th']]
                    perc84 = [bin['y 84th'], bin['y 84th']]

                    ax.plot(window, median,                 color=color,                                zorder=3, label='')
                    ax.fill_between(window, perc16, perc84, facecolor=color, edgecolor=None, alpha=0.5, zorder=2, label=gal if i==0 else '')
                    i += 1

        ax.set_xlim(1.5e20, 2.0e24)
        ax.set_ylim(8e-2, 3e2)
        ax.set_xlabel(r'$N$ [cm$^{-2}$]')
        ax.set_ylabel(r'$\frac{\sigma^2}{R}$ [km$^2$\,s$^{-2}$\,pc$^{-1}$]')
        ax.set_xscale('log')
        ax.set_yscale('log')

        # alpha lines
        for alpha in [1,10,100]:
            # plot alpha lines
            x = np.array([1e19, 1e25])*u.cm**-2
            y = alpha_x_to_y(x)
            line = ax.plot(x.value, y.value,
                           ls=':', lw=1, color='darkgrey', zorder=1
                          )
            # label only second panel
            if ax==axes[1][0]:
                x = 8e21*u.cm**-2
                y = 0.3*alpha_x_to_y(x)
                label_line(ax, line, x.value, y.value,
                           r'$\alpha_\mathrm{vir} = '+str(alpha)+'$',
                           angle=0.,
                           **{'color': 'darkgrey', 'ha':'left', 'fontsize': 8}
                          )

        # Pe/k lines
        for Pek in np.array([1e4,1e5,1e6,1e7])*u.K/u.cm**3:
            # plot Pe/k lines
            x = np.logspace(np.log10(1e20), np.log10(5e24),100)*u.cm**-2
            y = Pek_x_to_y(x, Pek)
            line = ax.plot(x.value, y.value,
                           ls='--', lw=1, color='darkgrey', zorder=1
                          )
            # label only second panel
            if ax==axes[0][0]:
                ax.text(1e23, (Pek.value/1e5)**0.5,
                        'P/k\,=\,$10^{'+str(int(np.log10(Pek.value)))+r'}$\,K\,cm$^{\textnormal{-} 3}$',
                        color = 'grey',
                        ha    = 'left',
                        fontsize = 6
                       )

        # label bin size
        size_bin = [int(x) if x>=1. else x for x in size_bin]
        ax.text(0.95,0.05, 'size bin '+str(size_bin[0])+'--'+str(size_bin[1])+r'\,pc', color='darkgrey', transform=ax.transAxes, ha='right', va='bottom')
        ax.set_axisbelow(True)
        ax.grid(ls=':', c='lightgrey')

for row in axes:
    row[0].tick_params(axis='y', left='on', top='off', right='off', bottom='off', labelleft='on', labeltop='off', labelright='off', labelbottom='off')
    row[1].tick_params(axis='y', left='off', top='off', right='on', bottom='off', labelleft='off', labeltop='off', labelright='on', labelbottom='off')
    row[1].yaxis.set_label_position('right')
    for ax in row:
        if not np.all(row==axes[-1]):
            ax.tick_params(axis='x', left='off', top='on', right='off', bottom='on', labelleft='off', labeltop='off', labelright='off', labelbottom='off')
            ax.set_xlabel('')
        else:
            ax.tick_params(axis='x', left='off', top='on', right='off', bottom='on', labelleft='off', labeltop='off', labelright='off', labelbottom='on')

axes[0][0].text(0.5,1.02, 'CO(1--0)', color='k', transform=axes[0][0].transAxes, ha='center', va='bottom', fontsize=14)
axes[0][1].text(0.5,1.02, 'CO(3--2)', color='k', transform=axes[0][1].transAxes, ha='center', va='bottom', fontsize=14)
axes[0][1].legend(bbox_to_anchor=(-0.95,0.5,0.4,0.4), loc='upper left', borderaxespad=0., ncol=1, fancybox=True, fontsize=10)
fig.savefig(join(plotdir, 'paper', 'figD.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# table X: machine-readable binned data
###################################################################################################

datfile = 1
for relation,symx,symy,unitx,unity in [['size-linewidth', 'R', 'sigma',   'pc',   'km/s'],
                                       ['size-luminosity','R', 'L',       'pc',   'K.km/s.pc2'],
                                       ['density-virial', 'N', 'sigma2/R','1/cm2','km2/s2/pc']]:

    # set precision
    if symx=='R' or symx=='sigma':
        precx = '6.2f'
        fmatx = 'F3.2'
        lengthx = 6
    else:
        precx = '8.2e'
        fmaty = 'E1.2'
        lengthx = 8

    if symy=='R' or symy=='sigma' or symy=='sigma2/R':
        precy = '6.2f'
        fmaty = 'F3.2'
        lengthy = 6
    else:
        precy = '8.2e'
        fmaty = 'E1.2'
        lengthy = 8

    # set up table
    tab = open(join(plotdir,'paper','MRF.'+relation+'.txt'),'w+')
    tab.write('Title: The turbulent gas structure in the centers of NGC253 and the Milky Way\n')
    tab.write('Authors: Krieger N., Bolatto A.D., Koch E., Leroy A.K., Rosolowsky E., Walter F., Weiss A., Eden D.J., Levy R.C., Meier D.S., Mills E.A.C., Moore T., Ott J., Su Y., Veilleux S.\n')
    tab.write('Table: '+relation+' relation\n')

    # byte header
    tab.write('================================================================================\n')
    tab.write('Byte-by-byte Description of file: datafile'+str(datfile)+'.txt\n')
    tab.write('--------------------------------------------------------------------------------\n')
    tab.write('{:<5} {:<6} {:<10} {:<12} {}\n'.format('Bytes','Format','Units','Label','Explanations'))
    tab.write('--------------------------------------------------------------------------------\n')

    # count up bytes
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(1, 6,  'A6', '---', 'galaxy', 'source'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(9, 16, 'A7', '---', 'CO',     'CO transition'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(19,                        18+lengthx,                 fmatx, unitx, symx+'-min', 'Lower edge of '+symx+' bin'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(19+lengthx+3,              18+2*lengthx+3,             fmatx, unitx, symx+'-max', 'Upper edge of '+symx+' bin'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(19+2*lengthx+6,            18+2*lengthx+6+lengthy,     fmaty, unity, symy+'-min', 'Lower bound of '+symy+' distribution (16th percentile)'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(19+2*lengthx+6+lengthy+3,  18+2*lengthx+6+2*lengthy+3, fmaty, unity, symy+'-med', 'Median of '+symy+' distribution (50th percentile)'))
    tab.write('{:>2}-{:2} {:<4}   {:<10} {:<12} {}\n'.format(19+2*lengthx+6+2*lengthy+6,18+2*lengthx+6+3*lengthy+6, fmaty, unity, symy+'-max', 'Upper bound of '+symy+' distribution (84th percentile)'))

    tab.write('--------------------------------------------------------------------------------\n')
    tab.write('Note (1): \n')
    tab.write('--------------------------------------------------------------------------------\n')

    for co,CO in lines.items():
        for gal,GAL in galaxies.items():

            # data to be written to table
            bin_min = paperfits[co][gal][relation]['bins']['window min'].data
            bin_max = paperfits[co][gal][relation]['bins']['window max'].data
            y_min   = paperfits[co][gal][relation]['bins']['y 16th'].data
            y_med   = paperfits[co][gal][relation]['bins']['y median'].data
            y_max   = paperfits[co][gal][relation]['bins']['y 84th'].data

            # data
            for a,b,c,d,e in zip(bin_min,bin_max,y_min,y_med,y_max):
                tab.write('{:<6}   {:<7}   {:{precx}}   {:{precx}}   {:{precy}}   {:{precy}}   {:{precy}}\n'.format(gal, co, a,b,c,d,e, precx=precx, precy=precy))

    # save file
    tab.close()

    datfile += 1


###################################################################################################
# table X: example of machine-readable tables
###################################################################################################

table = r"""\begin{deluxetable*}{llccccc}
	\tablewidth{\linewidth}
	\tablecaption{Sample of the binned size -- line width data for \co32 in NGC253. All data of the size -- line width, size -- luminosity and column density -- $\sigma^2$/R relations for both tracers and both sources is available in the machine-readable format in the online journal. The table shown here provides guidance regarding the form and content.}
	\label{tableX}
	\tablehead{\colhead{galaxy} & \colhead{CO} & \colhead{R$_\mathrm{min}$} & \colhead{R$_\mathrm{max}$} & \colhead{$\sigma_\mathrm{16th}$} & \colhead{$\sigma_\mathrm{median}$} & \colhead{$\sigma_\mathrm{84th}$} \\
	&& [pc] & [pc] & [km\,s$^{-1}$] & [km\,s$^{-1}$] & [km\,s$^{-1}$]\\ \relax
    && (1) & (2) & (3) & (4) & (5)
	}
	\startdata
"""
co = 'CO(3-2)'
gal = 'NGC253'
relation = 'size-linewidth'

for bin_min,bin_max,y_min,y_med,y_max in zip(paperfits[co][gal][relation]['bins']['window min'].data,
                                             paperfits[co][gal][relation]['bins']['window max'].data,
                                             paperfits[co][gal][relation]['bins']['y 16th'].data,
                                             paperfits[co][gal][relation]['bins']['y median'].data,
                                             paperfits[co][gal][relation]['bins']['y 84th'].data):
    table += '{:<6} & {:<7} & {:5.2f} & {:5.2f} & {:5.2f} & {:5.2f} & {:5.2f}\\\\\n'.format(gal,co,bin_min,bin_max,y_min,y_med,y_max)

table += r"""    \enddata
    \tablecomments{(1) Lower edge of the size bin.
                   (2) Upper edge of the size bin.
                   (3) Lower bound of the line width distribution (16th percentile).
                   (4) Median of the line width distribution (50th percentile).
                   (5) Upper bound of the line width distribution (84th percentile).
                  }
\end{deluxetable*}
"""
print(table)



###################################################################################################
# UNUSED fig 5: size-linewidth slope comparison
###################################################################################################

# import matplotlib.patches as patches
# from matplotlib.ticker import MultipleLocator
#
# slw_gal   = [{'label': 'Galactic plane CO(1--0) S+87',   'value': 0.50, 'error': 0.05,        'color': 'dodgerblue'},
#              {'label': 'Galactic plane CO(1--0) MD+17',  'value': 0.63, 'error': 0.30,        'color': 'dodgerblue'},
#              {'label': 'Galactic plane $^{13}$CO (3--2) R+19', 'value': 0.63, 'error': 0.03,  'color': 'dodgerblue'},
#              {'label': 'outer galaxy CO(1-0) R+16',      'value': 0.49, 'error': 0.04,        'color': 'dodgerblue'},
#              {'label': 'inner galaxy CO(1-0) R+16',      'value': 0.52, 'error': 0.03,        'color': 'dodgerblue'},
#              {'label': 'Perseus HCO$^+$(1--0) S+16',     'value': 0.50, 'error': 0.07,        'color': 'lightskyblue'},
#              {'label': 'Perseus HCN(1--0) S+16',         'value': 0.60, 'error': 0.10,        'color': 'deepskyblue'},
#              {'label': 'Perseus N$_2$H$^+$(1--0) S+16',  'value': 0.77, 'error': 0.21,        'color': 'cornflowerblue'},
#              {'label': 'GC N$_2$H$^+$(3--2) K+17',       'value': 0.66, 'error': 0.18,        'color': 'cornflowerblue'},
#              {'label': 'GC N$_2$H$^+$(1--0) S+12',       'value': 0.78, 'range': [0.41,1.13], 'color': 'cornflowerblue'},
#              {'label': 'GC HCN(1--0) S+12',              'value': 0.62, 'range': [0.27,1.02], 'color': 'deepskyblue'},
#              {'label': 'GC H$^{13}$CN(1--0) S+12',       'value': 0.66, 'range': [0.21,1.66], 'color': 'deepskyblue'},
#              {'label': 'GC HCO$^+$(1--0) S+12',          'value': 0.79, 'range': [0.45,1.30], 'color': 'lightskyblue'}]
# slw_own   = [{'label': 'GC CO(1--0)',                    'value': 0.74, 'error': 0.04,        'color': 'dodgerblue'},
#              {'label': 'GC CO(3--2)',                    'value': 0.72, 'error': 0.03,        'color': 'dodgerblue'},
#              {'label': 'NGC253 CO(1--0)',                'value': 0.82, 'error': 0.02,        'color': 'orange'},
#              {'label': 'NGC253 CO(3--2)',                'value': 0.62, 'error': 0.01,        'color': 'orange'}]
# slw_extra = [{'label': 'LMC $^{12}$CO(1--0) W+19',       'value': 0.64, 'range': [0.32,0.97], 'color': 'darkorange'},
#              {'label': 'LMC $^{13}$CO(1--0) W+19',       'value': 0.59, 'ramge': [0.36,1.00], 'color': 'darkorange'},
#              # {'label': 'LMC Planck cold cloud $^{12}$CO(2-1) W+19', 'value': 0.51, 'error': 0.02, 'color': 'darkorange'},
#              # {'label': 'LMC Planck cold cloud $^{13}$CO(2-1) W+19', 'value': 0.91, 'error': 0.14, 'color': 'darkorange'},
#              {'label': '30Doradus $^{12}$CO(1--0) W+19', 'value': 0.60, 'error': 0.03,        'color': 'darkorange'},
#              {'label': '30Doradus $^{13}$CO(1--0) W+19', 'value': 0.43, 'error': 0.08,        'color': 'darkorange'},
#              {'label': '30Doradus CO(1--0) N+16',        'value': 0.65, 'error': 0.04,        'color': 'darkorange'},
#              {'label': 'nearby galaxies CO(1--0 or 2--1) B+08', 'value': 0.60, 'error': 0.10, 'color': 'goldenrod'},
#              {'label': 'NGC300 CO(2--1) F+18',           'value': 0.48, 'error': 0.05,        'color': 'goldenrod'}]
#
#
# fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(4,4))
#
# for idx,slw in enumerate(slw_gal+slw_own+slw_extra):
#     if 'range' in slw.keys():
#         value = slw['value']
#         low   = slw['range'][0]
#         high  = slw['range'][1]
#     elif 'error' in slw.keys():
#         value = slw['value']
#         low   = value-slw['error']
#         high  = value+slw['error']
#     ax.add_patch(patches.Rectangle((idx-0.4,low), width=0.8, height=high-low, lw=0, edgecolor=slw['color'], facecolor=slw['color'], alpha=0.6, zorder=2))
#     ax.plot(idx, value, marker='_', color=slw['color'], lw=10, zorder=3)
#
# ax.axvspan(len(slw_gal)-0.5, len(slw_gal)+len(slw_own)-0.5, color='lightgrey', lw=0, alpha=1., zorder=0)
# ax.text(len(slw_gal)-0.5+len(slw_own)/2,                 1.5, 'this work',                   color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
# ax.text(len(slw_gal)/2-0.5,                              1.5, 'literature:\nGalactic',       color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
# ax.text(len(slw_gal)+len(slw_own)-0.5+len(slw_extra)/2., 1.5, 'literature:\nextra-galactic', color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
#
# ax.set_xlabel('')
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.set_xticklabels(['','']+[slw['label'] for slw in slw_gal+slw_own+slw_extra]+[''], rotation=90, fontsize=6)
# ax.set_ylabel('size -- line width exponent $b$')
# ax.set_xlim(-1, len(slw_gal+slw_own+slw_extra))
# ax.set_ylim(0.15, 1.70)
# ax.set_axisbelow(True)
# ax.grid(ls=':', c='darkgrey', axis='y', zorder=1)
# fig.savefig(join(plotdir, 'paper', 'fig5.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# UNUSED fig 6: size-mass slope comparison
###################################################################################################

# import matplotlib.patches as patches
# from matplotlib.ticker import MultipleLocator
#
# sm_gal   = [{'label': 'Galactic plane CO(1--0) MD+17',         'value': 2.20, 'error': 0.20, 'color': 'dodgerblue'},
#             {'label': 'Galactic plane CO(3-2) full C+19',      'value': 2.26, 'error': 0.06, 'color': 'deepskyblue'},
#             {'label': 'Galactic plane CO(3-2) fiducial C+19',  'value': 2.17, 'error': 0.38, 'color': 'deepskyblue'},
#             {'label': 'Galactic plane CO(3-2) full R+19',      'value': 2.26, 'error': 0.02, 'color': 'deepskyblue'},
#             {'label': 'Galactic plane CO(3-2) 8-12\,kpc R+19', 'value': 2.42, 'error': 0.05, 'color': 'deepskyblue'}]
# sm_own   = [{'label': 'GC CO(1--0)',                     'value': 3.25, 'error': 0.13, 'color': 'dodgerblue'},
#             {'label': 'GC CO(3--2)',                     'value': 2.69, 'error': 0.02, 'color': 'deepskyblue'},
#             {'label': 'NGC253 CO(1--0)',                 'value': 2.92, 'error': 0.07, 'color': 'darkorange'},
#             {'label': 'NGC253 CO(3--2)',                 'value': 2.89, 'error': 0.02, 'color': 'orange'}]
# sm_extra = [{'label': 'nearby galaxies B+08',            'value': 2.54, 'error': 0.20, 'color': 'orange'},
#             {'label': 'LMC Planck cold cloud CO(2-1) W+17', 'value': 2.52, 'error': 0.05, 'color': 'goldenrod'},
#             {'label': '30Doradus CO(2-1) W+17',          'value': 3.10, 'error': 0.11, 'color': 'goldenrod'},
#             {'label': 'M51 molecular ring CO(1-0) C+14', 'value': 2.20, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 material arms CO(1-0) C+14',  'value': 1.80, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 nuclear bar CO(1-0) C+14',    'value': 1.70, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 spiral arm DWI CO(1-0) C+14', 'value': 2.00, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 spiral arm DWO CO(1-0) C+14', 'value': 2.00, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 arm downstream CO(1-0) C+14', 'value': 1.50, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'M51 arm upstream CO(1-0) C+14',   'value': 1.50, 'error': 0.00, 'color': 'darkorange'},
#             {'label': 'NGC300 CO(2-1) F+18',             'value': 2.00, 'error': 0.12, 'color': 'goldenrod'},
#             {'label': 'CenA northern filament CO(1-0) S+17', 'value': 3.5, 'error': 0.5, 'color': 'darkorange'}]
#
#
# fig,ax = plt.subplots(nrows=1, ncols=1, squeeze=True, sharex='none', sharey='none', figsize=(4,4))
#
# for idx,sm in enumerate(sm_gal+sm_own+sm_extra):
#     if 'range' in sm.keys():
#         value = sm['value']
#         low   = sm['range'][0]
#         high  = sm['range'][1]
#     elif 'error' in sm.keys():
#         value = sm['value']
#         low   = value-sm['error']
#         high  = value+sm['error']
#     ax.add_patch(patches.Rectangle((idx-0.4,low), width=0.8, height=high-low, lw=0, edgecolor=sm['color'], facecolor=sm['color'], alpha=0.6, zorder=2))
#     ax.plot(idx, value, marker='_', color=sm['color'], lw=10, zorder=3)
#
# ax.axvspan(len(sm_gal)-0.5, len(sm_gal)+len(sm_own)-0.5, color='lightgrey', lw=0, alpha=1., zorder=0)
# ax.text(len(sm_gal)-0.5+len(sm_own)/2,                 3.75, 'this work',                   color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
# ax.text(len(sm_gal)/2-0.5,                             3.75, 'literature:\nGalactic',       color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
# ax.text(len(sm_gal)+len(sm_own)-0.5+len(sm_extra)/2.,  3.75, 'literature:\nextra-galactic', color='dimgrey', ha='center', va='center', weight='bold', fontsize=10, zorder=4)
#
# ax.set_xlabel('')
# ax.xaxis.set_major_locator(MultipleLocator(1))
# ax.set_xticklabels(['','']+[sm['label'] for sm in sm_gal+sm_own+sm_extra]+[''], rotation=90, fontsize=6)
# ax.set_ylabel('size--luminosity (mass) exponent $d$')
# ax.set_xlim(-1, len(sm_gal+sm_own+sm_extra))
# ax.set_ylim(1.3, 4.2)
# ax.set_axisbelow(True)
# ax.grid(ls=':', c='darkgrey', axis='y', zorder=1)
# fig.savefig(join(plotdir, 'paper', 'fig6.pdf'), dpi=300, bbox_inches='tight')

###################################################################################################
# UNUSED fig C: density mass PDFS
###################################################################################################

# plt.clf();plt.cla()
# fig  = plt.figure(figsize=(8,5))
# width  = (1.-2*0.05-0.025-0.075)/3
# height = (1.-2*0.05-0.025)/2
# left   = [0.05, 0.05+width+0.025, 0.05+width+0.05+width+0.05]
# bottom = [0.05+height+0.025, 0.05]
# axes = np.array([[plt.axes([left[0],bottom[0],width,height]), plt.axes([left[1],bottom[0],width,height]), plt.axes([left[2],bottom[0],width,height])],
#                  [plt.axes([left[0],bottom[1],width,height]), plt.axes([left[1],bottom[1],width,height]), plt.axes([left[2],bottom[1],width,height])]])
#
# for row,co in enumerate(lines.keys()):
#     for idx,(gal,color) in enumerate(zip(reversed(list(galaxies.keys())),gal_colors)):
#
#         catalog    = dendrograms[co][gal]['catalog']
#         dendrogram = dendrograms[co][gal]['dendrogram']
#
# # column density
# ###################################################################################################
#
#         ax = axes[row][0]
#         density   = catalog['column density']
#         density = density[np.isfinite(density)]
#
#         fitrange = paperfits[co][gal]['column density PDF']['fit range']
#         hist  = paperfits[co][gal]['column density PDF']['hist']
#         bins  = paperfits[co][gal]['column density PDF']['bins']
#         slope = paperfits[co][gal]['column density PDF']['slope']
#         intercept = paperfits[co][gal]['column density PDF']['intercept']
#
#         ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
#         ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)
#
#         fit_fn   = np.poly1d([slope['value'], intercept['value']])
#         x  = np.log10(np.linspace(fitrange[0],fitrange[1], 10))
#         ax.plot(np.power(10,x),  np.power(10,fit_fn(x)),  ls='-', color=color)
#
#         # x = np.array([3e21,3e22])
#         # for b in [-1.5, -2.0, -2.5]:
#         #     a = 10**(-np.log10(3e22)*b+np.log10(5e-26))                          # to make the lines cross (10^22, 10^-24.75)
#         #     y = a*x**b
#         #     ax.plot(x, y, 'grey', ls='--', lw=0.5)
#         #     ax.text(2e21, 10**(np.log10(5e-26)-b), str(b), color='grey', va='center', ha='right', fontsize=6)
#
#         ax.set_xlabel(r'$N$ [cm$^{-2}$]')
#         ax.set_ylabel(r'PDF')
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#         ax.set_axisbelow(True)
#         ax.grid(ls=':', c='lightgrey')
#         ax.text(0.95,0.88+idx*0.08, gal+' '+co, color=color, transform=ax.transAxes, ha='right', va='top')
#
#
# # volume density
# ###################################################################################################
#
#         ax = axes[row][1]
#         density   = catalog['volume density (astrodendro)']
#         density = density[np.isfinite(density)]
#
#         fitrange = paperfits[co][gal]['volume density PDF']['fit range']
#         hist  = paperfits[co][gal]['volume density PDF']['hist']
#         bins  = paperfits[co][gal]['volume density PDF']['bins']
#         shape = paperfits[co][gal]['volume density PDF']['shape']
#         loc   = paperfits[co][gal]['volume density PDF']['loc']
#         scale = paperfits[co][gal]['volume density PDF']['scale']
#
#         ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
#         ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)
#
#         x = np.logspace(np.log10(fitrange[0]),np.log10(fitrange[1]),100)
#         pdf = lognorm.pdf(x, shape['value'], loc['value'], scale['value'])
#         ax.plot(x, pdf, color=color, ls='-')
#
#         ax.set_xlabel(r'$\rho_\mathrm{vol}$ [cm$^{-3}$]')
#         ax.set_ylabel(r'PDF')
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#         ax.set_axisbelow(True)
#         ax.grid(ls=':', c='lightgrey')
#         ax.text(0.95,0.88+idx*0.08, gal+' '+co, color=color, transform=ax.transAxes, ha='right', va='top')
#
#
# # mass spectra
# ###################################################################################################
#
#         ax = axes[row][2]
#         mass   = catalog['mass']
#         mass = mass[np.isfinite(mass)]
#
#         fitrange = paperfits[co][gal]['mass PDF']['fit range']
#         hist  = paperfits[co][gal]['mass PDF']['hist']
#         bins  = paperfits[co][gal]['mass PDF']['bins']
#         slope = paperfits[co][gal]['mass PDF']['slope']
#         intercept = paperfits[co][gal]['mass PDF']['intercept']
#
#         ax.hist(bins[:-1], bins, weights=hist, histtype='step',       color=color,            log=True)
#         ax.hist(bins[:-1], bins, weights=hist, histtype='stepfilled', color=color, alpha=0.5, log=True)
#
#         fit_fn   = np.poly1d([slope['value'], intercept['value']])
#         x  = np.log10(np.linspace(fitrange[0],fitrange[1], 10))
#         ax.plot(np.power(10,x),  np.power(10,fit_fn(x)),  ls='-', color=color)
#
#         ax.set_xlabel('M [M$_\odot$]')
#         ax.set_ylabel('PDF')
#         ax.set_xscale('log')
#         ax.set_yscale('log')
#         ax.set_xlim(1e1, 1e9)
#         ax.set_axisbelow(True)
#         ax.grid(ls=':', c='lightgrey')
#         ax.text(0.95,0.88+idx*0.08, gal+' '+co, color=color, transform=ax.transAxes, ha='right', va='top')
#
#
# for row in axes:
#     row[1].tick_params(axis='y', left='on', top='off', right='on', bottom='off', labelleft='off', labeltop='off', labelright='on', labelbottom='off')
#     # row[1].yaxis.set_label_position('right')
#     row[1].yaxis.label.set_visible(False)
#     row[2].tick_params(axis='y', left='on', top='off', right='on', bottom='off', labelleft='off', labeltop='off', labelright='on', labelbottom='off')
#     row[2].yaxis.set_label_position('right')
#     row[0].set_xlim(2e20, 2e24)
#     row[1].set_xlim(2.5e-1, 5e4)
#     row[2].set_xlim(1e1, 1e9)
# for col in axes.transpose():
#     col[0].tick_params(axis='x', left='off', top='on', right='off', bottom='on', labelleft='off', labeltop='on', labelright='off', labelbottom='off')
#     col[0].xaxis.set_label_position('top')
#
# axes[0][0].set_ylim([4e-26, 5e-22])
# axes[1][0].set_ylim([2e-26, 5e-22])
# axes[0][1].set_ylim([8e-5, 2e-1])
# axes[1][1].set_ylim([2e-7, 3e-2]) #[8e-9, 3e-2])
# axes[0][2].set_ylim([3e-10, 1e-4])
# axes[1][2].set_ylim([3e-10, 1e-2])
#
# fig.savefig(join(plotdir, 'paper', 'figC.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
# UNUSED table 3: PDF fits
###################################################################################################

# table = r"""\begin{deluxetable}{lcccccccc}
#     \tablewidth{\textwidth}
#     \tablecaption{Results for the lognormal fits to the volume density PDFs and power law fits to the mass PDFs shown in figure~\ref{figC}.}
#     \label{table3}
#     \tablehead{\colhead{source} & \colhead{line} & \multicolumn{2}{c}{column density} & \multicolumn{2}{c}{volume density} & \multicolumn{2}{c}{mass}\\
#     && $\gamma$ & $a$\ & $\mu$ & $\sigma$ & $\gamma$ & $a$\\
#     && (1) & (2) & (3) & (4) & (5) & (6)
#     }
#     \startdata
# """
# for co in lines:
#      for gal in galaxies:
#
#          n_slope     = paperfits[co][gal]['column density PDF']['slope']['value']
#          n_intercept = paperfits[co][gal]['column density PDF']['intercept']['value']
#          d_mu        = paperfits[co][gal]['volume density PDF']['mu']['value']
#          d_sigma     = paperfits[co][gal]['volume density PDF']['sigma']['value']
#          m_slope     = paperfits[co][gal]['mass PDF']['slope']['value']
#          m_intercept = paperfits[co][gal]['mass PDF']['intercept']['value']
#
#          table += ('{:<6}'.format(gal)+' & '+
#                    '{:<7}'.format(co)+' & '+
#                    r'${:4.2f}$'.format(n_slope)+' & '+
#                    r'${:4.2f}$'.format(n_intercept)+' & '+
#                    r'${:4.2f}$'.format(d_mu)+' & '+
#                    r'${:4.2f}$'.format(d_sigma)+r' & '+
#                    r'${:4.2f}$'.format(m_slope)+' & '+
#                    r'${:4.2f}$'.format(m_intercept)+r' \\'+'\n'
#                   )
# table += r"""    \enddata
#     \tablecomments{(1) Exponent and (2) normalisation of the power law fit to the mass PDF. (3) Mean and (4) standard deviation of the lognormal fit to the volume density PDF. (5) Exponent and (6) normalisation of the power law fit to the mass PDF.}
# \end{deluxetable}
# """
# print(table)


###################################################################################################
#
###################################################################################################
