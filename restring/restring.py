# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com
# 051120

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd

from .gears import (
    manzlog,
    get_dirs,
    aggregate_results,
    tableize_aggregated,
    summary,
    write_all_aggregated,
    write_all_summarized
    )


def draw_clustermap(data, figsize=None, sort_values=None,
                    log_transform=True, log_base=10, log_na=0,
                    pval_min=None, custom_index=None, custom_cols=None,
                    unwanted_terms=None, title=None, title_size=24,
                    savefig=False, outfile_name="aggregated results.png",
                    dpi=300, readable=False, return_table=False,
                    **kwargs):
    
    """

    Draws a clustermap of an 'aggregated'-type table.
    
    This functions expects this table layout (example):
    
    terms │ exp. cond 1 | exp. cond 2 | exp cond n ..
    ──────┼─────────────┼─────────────┼────────────..
    term 1│     0.01    |      1      |  0.00023
    ------┼-------------┼-------------┼------------..
    term 2│      1      |    0.05     |     1    
    ------┼-------------┼-------------┼------------..
      ..  │     ..      |      ..     |     ..
      
    *terms* must be indices of the table. If present, 'common' column
    will be ignored.
    
    Params
    ------
    
    data      A <pandas.DataFrame object>, or a filename. If a filename is
              given, this will try to load the table as if it were produced
              by tableize_aggregated() and saved with pandas.DataFrame.to_csv()
              as a .csv file, with ',' as separator.
              
    sort_values   A <str> or <list> of <str>. Sorts the table accordingly. Please
              *also* set row_cluster=False (see below, seaborn.clustermap() 
              additional parameters), otherwise columns will still be clustered.              
    
    log_transform   If True, values will be log-transformed. Defaults to True.
              note: values are turned into -log values, thus p-value of 0.05
              gets transformed into 1.3 (with default parameters), as
              10^-1.3 ~ 0.05
    
    log_base  If log_transform, this base will be used as the logarithm base.
              Defaults to 10
    
    log_na    When unable to compute the logarithm, this value will be used
              instead. Defaults to 0 (10^0 == 1; actually, to -0 so as to
              have 0 in the final table)
    
    pval_min  Trims values to the ones matching *at least* this p value. If
              log_transform with default values, this needs to be set accordingly
              For example, if at least p=0.01 is desired, then aim for 
              pval_min=2 (10^-2 == 0.01)
              
    custom_cols   It is possible to pass a <list> of <str> to draw only specified 
              columns of input table
    
    custom_index  It is possible to pass a <list> of <str> to draw only specified
              rows of input table
    
    unwanted_terms   If a <list> of <str> is supplied, 
              
    readable  If set to False (default), the generated heatmap will be of reasonable
              size, possibly not showing all term descriptors. Setting readable=True
              will result into a heatmap where all descriptors are visible (this might
              result into a very tall heatmap)
              
    savefig   If True, a picture of the heatmap will be saved in the current directory.
    
    outfile_name  The file name of the picture saved.
    
    dpi       dpi resolution of saved picture. Defaults to 300.

    return_table   If set to True, it also returns the table that was manipulated
              internally to draw the heatmap (with all modifications applied).
              Returns: <seaborn.matrix.ClusterGrid>, <pandas.core.frame.DataFrame>

    **kwargs    The drawing is performed by seaborn.clustermap(). All additional
              keyword arguments are passed directly to it, so that the final picture
              can be precisely tuned. More at:
              https://seaborn.pydata.org/generated/seaborn.clustermap.html
    
    """
    
    if isinstance(data, str):
        table = pd.read_csv(data, index_col=0)
    else:
        table = data.copy()
    
    if "common" in table.columns:
        del table["common"]
    
    if sort_values is not None:        
        table = table.sort_values(by=sort_values)
    
    try:    
        if log_transform:
            table = table.applymap(lambda x: -manzlog(x, log_base, log_na))
            # as log(1) == 0, but -log(1) == -0, we need to switch back
            # the sign of zeroes
            table = table.replace({-0: 0})
    except TypeError:
        print("*error* Check the table layout, something's not right.")
        return table
    
    # custom column content --------
    
    if custom_cols is not None:
        table = table[custom_cols]
    
    if custom_index is not None:
        try:
            table = table.reindex(custom_index)
        except NameError:
            print(f"*error* Something went wrong with custom index: {custom_index}")
            print("      Table index was NOT modified.")
            pass
    
    if unwanted_terms is not None:
        try:
            if isinstance(unwanted_terms, str):
                table = table.drop(unwanted_terms, axis=0)
            else:
                for x in unwanted_terms:
                    try:
                        table = table.drop(x, axis=0)
                    except:
                        print(f"*warning* while processing unwanted_terms")
                        print(f"*warning* Something went wrong with {x}")
        except:
            print(f"*error* Something went wrong with unwanted terms: {unwanted_terms}")
            print("         Table index was NOT modified as intended.")
    
    
    # we do this last, major modification of table content
    if pval_min is not None:
        bser = [any(table.loc[x, :] > pval_min) for x in table.index]
        table = table.loc[bser,:]

    # end of custom column content --------
    
    # draw --------------------------------
    
    if not readable:
        clus = sns.clustermap(
            table,
            **kwargs
        )
    else:
        height = len(table.index) * 0.289
        FIGSIZE = (10, height)
        
        clus = sns.clustermap(
            table,
            figsize=FIGSIZE,
            **kwargs
        )
    
    # end of draw --------------------------
    
    if title is not None:
        plt.title(title, size=title_size)
    
    if savefig:
        clus.savefig(
            outfile_name,
            dpi=dpi,
            bbox_inches="tight"
        )
    
    if return_table:
        return clus, table
    else:
        return clus