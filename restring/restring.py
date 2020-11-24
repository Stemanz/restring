# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com
# 171120

from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import os
from time import time, sleep

from .gears import (
    manzlog,
    get_dirs,
    aggregate_results,
    tableize_aggregated,
    summary,
    write_all_aggregated,
    write_all_summarized,
    keep_start,
    keep_end,
    keep_inside,
    prune_start,
    prune_end,
    prune_inside,

    # String API
    session_ID,
    get_functional_enrichment,
    write_functional_enrichment_tables
    )

from .settings import sep


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


class Aggregation:

    """This takes care of running the whole analysis.

    TODO: more doc here
    """

    def __init__(self, working_directory=None, overwrite=False, verbose=True):

        if working_directory is None:
            self.working_directory = os.path.abspath(".")
            self._original_path = os.path.abspath(".")
        else:
            if isinstance(working_directory, str):
                self.working_directory = working_directory
                os.chdir(self.working_directory)
            else:
                raise TypeError(f"*Error* : unrecognized path {working_directory}")

        self.verbose = verbose
        self.overwrite = overwrite

        self.say(f"Current working directory: {self.working_directory}")
        #self.show_params()


    def say(self, *args, **kwargs):
        if self.verbose:
            print(*args, **kwargs)


    def setwd(self, stringlike):
    
        if os.path.exists(stringlike):
            os.chdir(stringlike)
            self.working_directory = stringlike
            self.say(f"Working directory changed to:\n{self.working_directory}")
        else:
            self.say(f"Invalid path: {stringlike}")        

        
    def show_params(self):
        print(f"Current working directory: {self.working_directory}")


    def get_files(self, startswith=None, endswith=None, inside=None):
        """retrieves all files in current path, then keeps or discards
        them based on input parameters. Folders starting with "__" and "."
        are discarded

        To better define the file names to be analyzed, one or more strings
        can be specified to be contained in different parts of the filename.

        Params:
        =======
        startswith:  a <str> or <list> of <str> containing the identifiers.
               file names starting with any one of the identifiers will be kept,
               the others discarded.
        
        endswith:  a <str> or <list> of <str> containing the identifiers
               file names ending with any one of the identifiers will be kept,
               the others discarded.

        inside:  a <str> or <list> of <str> containing the identifiers
               file names containing any one of the identifiers will be kept,
               the others discarded.


        Returns:
        ========

        <list>
        """

        files = sorted(os.listdir("."))
        files = [x for x in files if not os.path.isdir(x)]
        files = [x for x in files if not x.startswith("__") and not x.startswith(".")]

        if startswith is not None:
            files = keep_start(files, startswith)

        if endswith is not None:
            files = keep_end(files, endswith)

        if inside is not None:
            files = keep_inside(files, inside)

        self._files = files
        return files


    def get_dirs(self):
        """retrieves all folders in current path, then keeps or discards
        them based on input parameters. Folders starting with "__" and "."
        are discarded
        """

        dirs = sorted(os.listdir("."))
        dirs = [x for x in dirs if os.path.isdir(x)]
        dirs = [x for x in dirs if not x.startswith("__") and not x.startswith(".")]
        
        self._dirs = dirs
        return dirs


    def file_analysis(self, kind="csv", sep=sep, species="mouse",
                      reverse_direction=False, query_wait_time=1,
        ):

        if not hasattr(self, "_files"):
            self.say(f"Retrieving the files in the current directory..")
            self._files = self.get_files()
            self._files = prune_end(self._files, "ipynb")
        else:
            self.say(f"{len(self._files)} files were set up to be analyzed.")

        self.say(f"These {len(self._files)} files will be processed:")
        for f in self._files:
            self.say(f)
        self.say(f"{'='*80}\n")

        # here we go!
        t0 = time()
        for f in self._files:
            if any([f.endswith(x) for x in (".txt", ".csv", ".xls", ".tsv", ".doc", ".tdt")]):
                extension = -4
            elif f.endswith(".xlsx"):
                extension = -5
            else:
                extension = None

            folder_name = f[:extension]

            # can't have both filename and dir named the same way
            # if we don't strip the extension, we need to change the dir name
            if extension is None:
                temppath = self.working_directory + f"/{folder_name}.restring"
            else:
                temppath = self.working_directory + f"/{folder_name}"

            if os.path.exists(temppath):
                if self.overwrite:
                    self.say(f"*Notice*: Path {f} exists, but we're going to overwrite it.\n")
                    os.chdir(temppath)
                else:
                    print(f"*Error* : path '{temppath}' already exists, and overwrite='False'.")
                    return None
            else:
                os.mkdir(folder_name)
                os.chdir(temppath)
                self.say(f"Now in: {temppath}")

            # we assume that the input files are tabular in nature.
            # first column: gene identifiers
            # second columns: fold change values

            if kind == "csv":
                tempdf = pd.read_csv(f"../{f}", index_col=0, sep=sep)
                self.say(f"Reading from: {f}")
            elif kind == "xls":
                tempdf = pd.read_excel(f"../{f}", index_col=0)
                self.say(f"Reading from: {f}")
            else:
                raise NotImplementedError(f"*Error*: can't handle: {kind}.")

            # just a bunch of quick checks
            rownumber, colnumber = tempdf.shape
            if colnumber != 1:
                err_message = f"*Error*: Wrong format for '{f}'."
                err_message += f"\nColumns retrieved from table: {tempdf.columns}"
                if rownumber > 0:
                    err_message += f"First index element: '{tempdf.index[0]}'."
                raise NotImplementedError(err_message)
            if rownumber == 0:
                print(f"*Notice*: no genes to process in '{f}'.")
                continue

            col = tempdf.columns[0]
            all_gene_list  = list(tempdf.index)

            # a brief note on UP and DOWN genes. the convention is this:
            #
            # cond1  | cond2  | cond1_vs_cond2 |  log2FC |  
            # -------------------------------------------|
            #   143  |  748   |      0.191     |  -2.38  |
            # -------------------------------------------|
            #   50   |   4    |      12.5      |   3.64  |
            #
            # UP means: UP in cond1 vs cond2 (logFC > 0)
            # DOWN means: DOWN in cond1 vs cond2 (logFC < 0)
            #
            # To reverse this, call with reverse_direction=True
            
            if not reverse_direction:
                up_gene_list   = list(tempdf[tempdf[col] > 0].index)
                down_gene_list = list(tempdf[tempdf[col] < 0].index)
            else:
                up_gene_list   = list(tempdf[tempdf[col] < 0].index)
                down_gene_list = list(tempdf[tempdf[col] > 0].index)

            string_params = {
                "species": species,
                "caller_ID": session_ID,
                "allow_pubmed": 0
            }

            up_df = get_functional_enrichment(up_gene_list, **string_params)
            sleep(query_wait_time)
            write_functional_enrichment_tables(up_df, prefix="UP_")

            down_df = get_functional_enrichment(down_gene_list, **string_params) 
            sleep(query_wait_time)
            write_functional_enrichment_tables(down_df, prefix="DOWN_")

            # TODO: implement this analysis
            #all_df = get_functional_enrichment(all_gene_list, **string_params)
            #sleep(query_wait_time)
            #write_functional_enrichment_tables(all_df, prefix="ALL_")

            self.say()
            os.chdir(self.working_directory)

        t1 = time()
        self.say(f"{'='*80}\nFinished making functional enrichment tables.")
        self.say(f"{round(t1-t0, 2)} seconds elapsed.")
        
