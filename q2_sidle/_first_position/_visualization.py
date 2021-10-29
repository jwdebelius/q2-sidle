import os
import warnings


import matplotlib.pyplot as plt
import numpy as np
import pkg_resources
import pandas as pd
import q2templates
import seaborn as sn
import skbio

from qiime2 import Metadata

from q2_sidle._primerless import (_generate_align_mask,
                          )


TEMPLATES = pkg_resources.resource_filename('q2_sidle._first_position',
                                            'assets')


def summarize_alignment_positions(output_dir,
                                  alignment: pd.Series,
                                  position_summary: Metadata,
                                  sort_cols: str='starting-position',
                                  weight_by_abundance: bool=True,
                                  colormap: str=None,
                                  heatmap_maskcolor: str=None,
                                  heatmap_grid: bool=True,
                                  tick_interval: int=100,
                                  ):
    """
    Builds a visualization describing alignment between a first position 
    and reference
    """
    sort_cols = sort_cols.split(',')
    summary = position_summary.to_dataframe()
    summary['starting-position'] = summary['starting-position'].astype(int)
    dir_ = summary['direction'].unique()[0]

    # Builds the heatmap
    weighted_coverage = _build_cover_matrix(alignment=alignment,
                                            summary=summary,
                                            sort_cols=sort_cols,
                                            weight=weight_by_abundance,
                                            )
    heatmap = _make_alignment_heatmap(weighted_coverage=weighted_coverage,
                                      cmap=colormap,
                                      maskcolor=heatmap_maskcolor,
                                      grid=heatmap_grid,
                                      tick_interval=tick_interval,
                                      direction=dir_
                                      )
    # Saves the heatmap 
    for ext in ['png', 'svg']:
        img_fp = os.path.join(output_dir, 'alignment-start-heatmap.%s' % ext)
        heatmap.savefig(img_fp)


    # Builds the regional summary and formats to be a pretty(?) table
    regional_summary = _regional_alignment_results(summary)

    context = {'regional_summary': q2templates.df_to_html(regional_summary)}
    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context=context)


def _build_cover_matrix(alignment,
                        summary,
                        sort_cols=['starting-position'],
                        weight=True,
                        ):
    """
    Builds heatmap of sequences for plotting
    """
    # We'll start by checking for an overlap between the alignment ids and
    # summarized ids

    seq_ids = set(summary.index)
    align_ids = set(alignment.index)
    if ((seq_ids == align_ids) | (align_ids.issuperset(seq_ids))) == False:
        raise ValueError('The alignmnent does not contain all summarized '
                         'sequences. Please check that you have the correct '
                         'alignment and summary.')

    summary.sort_values(sort_cols, ascending=True, inplace=True)
    coverage = _generate_align_mask(alignment, summary.index)
    if weight:
        weighted_coverage = \
            (coverage.T * np.log10(summary['sequence-counts'].values)).T
    else:
        weighted_coverage = coverage

    return weighted_coverage


def _make_alignment_heatmap(weighted_coverage, 
                           cmap=None, 
                           maskcolor=None,
                           grid=True,
                           tick_interval=100,
                           test=False,
                           direction='fwd',
                           ):
    """
    Generates a heatmap showing alignment positions
    """
    fig = plt.figure(constrained_layout=True,
                     dpi=150,
                     )
    hm_ax = fig.add_subplot(
        1,1,1, facecolor={None: 'white'}.get(maskcolor, maskcolor)
        )
    hm_mask = hm_mask = ((weighted_coverage < 1) & (maskcolor != None))

    sn.heatmap(weighted_coverage,#.replace({0: np.nan}), 
               mask=hm_mask,
               cmap=None,
               yticklabels=[],
               ax=hm_ax,
               cbar_kws={'label': '$\\log_{10}$(Frequency)'},
               )
    hm_ax.set_ylabel("ASV", size=12)
    hm_ax.set_xlabel("Alignment Position", size=12)

    xticks = np.arange(0, len(weighted_coverage.T), tick_interval)
    if direction == 'fwd':
        hm_ax.set_xticks(xticks)
    else:
        hm_ax.set_xlim(hm_ax.get_xlim()[::-1])
        hm_ax.set_xticks(hm_ax.get_xlim()[0] - xticks)
    hm_ax.set_xticklabels(xticks)
    hm_ax.xaxis.set_tick_params(bottom=True, 
                                top=False, 
                                labelbottom=True, 
                                labeltop=False, 
                                labelsize=10,
                                rotation=90,
                                )
    hm_ax.grid(grid)

    if test:
        return fig, hm_mask
    return fig


def _regional_alignment_results(summary): 
    """
    Identifies information about positions
    """
    per_region = \
        summary.groupby(['starting-position', 'direction']).describe()
    per_region = per_region['sequence-counts']
    per_region.rename(columns={"count": 'number of mapped ASVs',
                               'min': 'Frequency of ASV with fewest counts',
                               'max': 'Frequency of ASV with most counts',
                               '50%': 'Median Frequency of mapped ASVs',
                               },
                      inplace=True)
    per_region.index.set_names(['Starting Position', 'Sequence Direction'],
                               inplace=True)
    per_region.reset_index('Sequence Direction', inplace=True, drop=False)

    keep_cols = ['Sequence Direction', 
                 'number of mapped ASVs',
                 'Frequency of ASV with fewest counts',
                 'Median Frequency of mapped ASVs',
                 'Frequency of ASV with most counts',
                 ]

    return per_region[keep_cols]

