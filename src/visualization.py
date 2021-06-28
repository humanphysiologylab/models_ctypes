import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def plot_waveforms(phenotype_model, phenotype_control, *,
                   fig=None, axes=None,
                   color_model='C0', color_control='0.3',
                   ylabels=None, xlabel=None, titles=None, suptitle=None,
                   ylim_list=None, yticks_list=None,
                   xlim=None, xticks=None,
                   inset_row_indices=None, xlim_inset=None, xticks_inset=None,
                   grid=True, grid_inset=True,
                   points_per_ms=1):

    ncols = len(phenotype_model)
    nrows = len(phenotype_model[0])

    if fig is None or axes is None:
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols,
                                 figsize=plt.figaspect(nrows * 1.1 / ncols),
                                 sharex='col', sharey='row',
                                 )

    axes = axes.reshape((nrows, ncols))

    if inset_row_indices is None:
        inset_row_indices = []

    if ylabels is None:
        ylabels = [None] * nrows

    if titles is None:
        titles = [None] * ncols

    if ylim_list is None:
        ylim_list = [None] * nrows

    if yticks_list is None:
        yticks_list = [None] * nrows

    for i_col in range(ncols):
        for i_row in range(nrows):

            waveform_model   = phenotype_model[i_col][i_row]
            waveform_control = phenotype_control[i_col][i_row]

            ax = axes[i_row, i_col]

            if i_row in inset_row_indices:
                for child in ax.get_children():
                    if isinstance(child, matplotlib.axes._axes.Axes):
                        ax_inset = child
                        break
                else:
                    ax_inset = ax.inset_axes([0.5, 0.5, 0.4, 0.4])
                axes_current = [ax, ax_inset]
            else:
                axes_current = [ax]

            kw_list  = [{'marker': None, 'ls': '-', 'lw': 1.5},
                        {'marker': '.', 'ls': '-', 'lw': 1}]

            def plot_fancy(x, ax, **kw):
                t = np.arange(len(x)) / points_per_ms
                ax.plot(t, x, color='w', lw=kw['lw']+0.5, zorder=kw.get('zorder', -1))
                ax.plot(t, x, **kw)

            for ax_, kw in zip(axes_current, kw_list):

                kw['color'] = color_control
                plot_fancy(waveform_control, ax_, **kw, zorder=-1)

                kw['color'] = color_model
                plot_fancy(waveform_model, ax_, **kw)

            # # # #
            #

            ylim_current = ax.get_ylim() if ylim_list[i_row] is None else ylim_list[i_row]
            xlim_current = ax.get_xlim() if xlim is None else xlim
            ax.set_ylim(ylim_current)  # needed for ticks
            ax.set_xlim(xlim_current)

            yticks_current = ax.get_yticks() if yticks_list[i_row] is None else yticks_list[i_row]
            xticks_current = ax.get_xticks() if xticks is None else xticks

            ax.set_yticks(yticks_current)
            ax.set_xticks(xticks_current)
            ax.set_ylim(ylim_current)
            ax.set_xlim(xlim_current)

            ax.tick_params(axis='both', labelsize='x-small')

            ax.grid(grid)

            if i_col == 0:
                ax.set_ylabel(ylabels[i_row])

            if i_row == nrows-1:
                ax.set_xlabel(xlabel)

            if i_row == 0:
                ax.set_title(titles[i_col], size='medium')

            if i_row in inset_row_indices:

                xlim_inset_current = ax_inset.get_xlim() if xlim_inset is None else xlim_inset
                ax_inset.set_xlim(xlim_inset_current)
                ax_inset.set_ylim(ylim_current)

                xticks_inset_current = ax_inset.get_xticks() if xticks_inset is None else xticks_inset

                ax_inset.set_yticks(yticks_current) #  ax_inset.set_yticks(yticks_current[1:-1])
                ax_inset.set_xticks(xticks_inset_current)
                ax_inset.set_xlim(xlim_inset_current)
                ax_inset.set_ylim(ylim_current)

                ax_inset.tick_params(axis='both', labelsize='xx-small')

                ax_inset.grid(grid_inset)

                  #
            # # # #

    fig.align_labels()
    plt.tight_layout()

    plt.subplots_adjust(top=0.85)
    fig.suptitle(suptitle)

    return fig
