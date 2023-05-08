import os
import pandas as pd
from os.path import isdir
from math import log
from io import StringIO
from glob import glob
import requests
from io import StringIO
from random import choice
from time import time, sleep
from matplotlib import pyplot as plt
import seaborn as sns


from .settings import (
    file_types,
    API_file_types,
    header_table,
    sep,
    PATH
)


session_ID = "".join(
    [choice(("abcdefghijklmnopqrstuvwxyz0123456789_")) for x in range(8)]
)
session_ID = f"restring-{session_ID}"


# basic gears ===============================================================

def manzlog(numlike, base=2, err=0):
    """This is a ZeroDivisionError safe version of math.log().
    Differently from numpy.log(), we are not looking for -inf, but
    looking to intercept the error and provide a custom value instead.
    """
    if numlike != 0:
        return log(numlike, base)
    else:
        return err


def get_dirs():
    """retrieves all folders in current path, then keeps or discards
    them based on input parameters. Folders starting with "__" and "."
    are discarded
    """

    dirs = sorted(os.listdir("."))
    dirs = [x for x in dirs if isdir(x)]
    dirs = [x for x in dirs if not x.startswith("__") and not x.startswith(".")]
    return dirs


def clean_dirs(listlike, startswith=["Ctrl", "Treated"], contains=None):
    """Clears a <list> of strings on unwanted terms. Used to obtain a subset
    of directories to process

    Params:
    =======
    listlike: <list> of strings to process


    startswith: <list> of textual identifiers for experimental condition 1; for example
                ["Ctrl", "Treated"]
                Directories names * must start * with these textual identifiers
    
    contains:    <list> of textual identifiers for experimental condition 2; for example
                ["24h", "48h"]
                Directories names * must contain * these textual identifiers
    """
    if not isinstance(startswith, list):
        say(f"Problems with param cond1: {startswith} of type {type(startswith)}")
        say("cond1 must be a <list>. Supplied a single <str> intead?")
        raise TypeError

        
    if not isinstance(contains, list):
        say(f"Problems with param cond2: {contains} of type {type(contains)}")
        say("cond2 must be a <list>. Supplied a single <str> intead?")
        raise TypeError

    wanted_dirs = []
    for elem in listlike:
        check = [elem.startswith(x) for x in startswith]
        if any(check):
            wanted_dir.append(elem)

    finally_wanted_dirs = []
    if contains is not None:
        for elem in wanted_dirs:
            check = [x in elem for x in contains]
            if any(check):
                finally_wanted_dirs.append(elem)

    return finally_wanted_dirs


def prune_start(listlike, unwanted):
    """
    Discards any of the elements of <listlike> if they start with
    any of the identifiers contained in <terms>
    """

    if isinstance(unwanted, str):
        unwanted = [unwanted]

    returnlist = []
    for x in listlike:
        for y in unwanted:
            err = 0
            if x.startswith(y):
                err += 1
                break
        if err == 0:    
            returnlist.append(x)

    return returnlist


def keep_start(listlike, wanted):
    """
    Keeps any of the elements of <listlike> if they start with
    any of the identifiers contained in <terms>
    """

    if isinstance(wanted, str):
        wanted = [wanted]

    returnlist = []
    for x in listlike:
        for y in wanted:
            if x.startswith(y):
                returnlist.append(x)
                break

    return returnlist


def prune_end(listlike, unwanted):
    """
    Discards any of the elements of <listlike> if they end with
    any of the identifiers contained in <terms>
    """

    if isinstance(unwanted, str):
        unwanted = [unwanted]

    returnlist = []
    for x in listlike:
        for y in unwanted:
            err = 0
            if x.endswith(y):
                err += 1
                break
        if err == 0:    
            returnlist.append(x)

    return returnlist


def keep_end(listlike, wanted):
    """
    Keeps any of the elements of <listlike> if they end with
    any of the identifiers contained in <terms>
    """

    if isinstance(wanted, str):
        wanted = [wanted]

    returnlist = []
    for x in listlike:
        for y in wanted:
            if x.endswith(y):
                returnlist.append(x)
                break

    return returnlist


def prune_inside( listlike, unwanted):
    """
    Discards any of the elements of <listlike> if they contain
    any of the identifiers contained in <terms>
    """

    if isinstance(unwanted, str):
        unwanted = [unwanted]

    returnlist = []
    for x in listlike:
        for y in unwanted:
            err = 0
            if y in x:
                err += 1
                break
        if err == 0:    
            returnlist.append(x)

    return returnlist


def keep_inside(listlike, wanted):
    """
    Keeps any of the elements of <listlike> if they contain
    any of the identifiers contained in <terms>
    """

    if isinstance(wanted, str):
        wanted = [wanted]

    returnlist = []
    for x in listlike:
        for y in wanted:
            if y in x:
                returnlist.append(x)
                break

    return returnlist


# ===============================================================

# ===============================================================

def aggregate_results(
    directories,
    kind="KEGG",
    directions=["UP", "DOWN"],
    verbose=True,

    # -- settings.py --

    file_types=file_types, 
    PATH=PATH,
    header_table=header_table,
):
    
    """
    Walks the given <directories> list, and reads the String .tsv files of
    defined <kind>.

    Params:
    =======
    directories: <list> of directories where to look for String files

    kind:    <str> Defines the String filetype to process. Kinds defined in
             settings.file_types

    directions: <list> contaning the up- or down-regulated genes in a comparison.
             Info is retrieved from either UP and/or DOWN lists.
             * Prerequisite *: generating files form String with UP and/or DOWN regulated
             genes separately.

    verbose: <bool>; turns verbose mode on or off


    Returns: <dict>
    ========
    
    Returned dict structure:
    ========================
    
    dictlike   keys      keys(2)
    {bestof} - {term1} - "exp condition" : <float>
                         "exp condition2": <float>
                         "hightes pval"  : <float>
             - {term2} - ...
             
    Call tableize_aggregated() on this dict to build a table
    """
    
    os.chdir(PATH)
    
    def say(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)
    
    say("Start walking the directory structure.\n")
    say(f"Parameters\n{'-'*10}\nfolders: {len(directories)}\nkind={kind}\ndirections={directions}\n")
    
    PROCESSED_DIRS = 0
    PROCESSED_FILES = 0
    
    if not isinstance(directions, list):
        try:
            # supplied a string? check if that's a valid direction
            temp = directions
            directions = []
            directions.append(temp)
            if directions[0] not in ["UP", "DOWN"]:
                raise TypeError
        except:
            say(f"Problems with param directions: {directions} of type {type(directions)}")
            say("directions must either be ['UP, DOWN'] 'UP' or 'DOWN'.")
            raise
    
    if kind not in file_types:
        raise TypeError(f"STRING analysis type must be one of these:\n{file_types}")
    
    bestof = {}
    
    # ID : the column name of the wanted name (process, KEGG, ..)
    # score: the column name or the wanted score (likely the false discovery rate or pval)
    # TERM_ID: the actual term to be retrieved
    # SCORE: the actual float score
    # gene_names: the column where the genes per term are stored
    
    # picking the header handles form settings.header_table
    # as of now, this is superfluous as all tables share the same layout. Should
    # this change in the future, just update settings.header_table
    ID, score, gene_names = header_table[kind].values()
                    
    # start walking the directories
    # =============================
    for d in directories:
        PROCESSED_DIRS += 1
        say(f"Processing directory: {d}")
        os.chdir(d)
        
        files = glob("*enrichment."+kind+".tsv")
        if len(files) > 0:
            for file in files: # either UP and/or DOWN
                
                if file[:2] in [x[:2] for x in directions]:
                    PROCESSED_FILES += 1
                    say(f"\tProcessing file {file}")

                    df = pd.read_csv(file, sep=sep, index_col=ID)
                    for TERM_ID in df.index:
                        # we are adding the first key to <bestof>.
                        # this key is the retrieved term.
                        # in this sub-dict, the first keys to be added
                        # are the following, then followed by all
                        # <d> (directory names, == experimental groups)
                        #
                        #
                        # bestof -- TERM -- "highest score" <float>
                        #               --- "genes" <set>
                        #               --- "common_temp" <dict>
                        #               --- "directory_1" <float>
                        #               --- "directory 2" <float>
                        #               --- ...
                        #        -- TERM2 - ...
                        #
                        bestof.setdefault(TERM_ID, {})
                        bestof[TERM_ID].setdefault("highest score", 1)
                        bestof[TERM_ID].setdefault("genes", set())
                        bestof[TERM_ID].setdefault("common_temp", dict())
                        
                        # adding a key (d - the directory) to the sub dict
                        # storing the pvalue of that term, for future heatmap
                        SCORE = df.loc[TERM_ID, score]
                        bestof[TERM_ID][d] = SCORE
                        if bestof[TERM_ID]["highest score"] > SCORE:
                            bestof[TERM_ID]["highest score"] = SCORE
                        
                        GENES = df.loc[TERM_ID, gene_names]
                        GENES = GENES.split(",")
                        bestof[TERM_ID]["genes"].update(set(GENES))
                        
                        # things get tricky. I'm looping over the directories,
                        # not terms. Now I have to store all genes for every
                        # directory, then only at the end I have to find
                        # the common elements
                        GENES = df.loc[TERM_ID, gene_names]
                        GENES = GENES.split(",")
                        bestof[TERM_ID]["common_temp"].setdefault(d, set(GENES))
        
        os.chdir("..")
    
    # final thing to do: loop over TERM_ID (bestof keys) and find the common genes,
    # then deleting individual sets
    
    for k in bestof.keys():
        curr = bestof[k]["common_temp"] # <dict>
        conditions = list(curr.keys())
        
        first_key = conditions.pop() # <set>
        first_set = curr[first_key]
        
        if len(conditions) > 0:
            for other_keys in conditions:
                first_set = first_set.intersection(curr[other_keys])

            # now first_set contains all genes that are commonly picked up in all conditions
            if len(first_set) == 0:
                first_set = set(list(["No common gene"]))
            
            
            del bestof[k]["common_temp"]
            bestof[k]["common"] = first_set
        else: # there is just one condition
            del bestof[k]["common_temp"]
            bestof[k]["common"] = set(list(["n/a (just one condition)"]))
        
        
    say(f"\nProcessed {PROCESSED_DIRS} directories and {PROCESSED_FILES} files.")
    say(f"Found a total of {len(bestof)} {kind} elements.")
    
    return bestof


def summary(dictlike):
    
    returnstring="ID\tscore\toccurrence\tall_genes\tcommon_genes\n"
    keys = sorted(dictlike.keys())
    for k in keys:
        occurrences = len(dictlike[k]) -3 # there are three hard-coded keys

        pval = dictlike[k]["highest score"]
        all_genes = ",".join(sorted(dictlike[k]["genes"]))
        common_genes = ",".join(sorted(dictlike[k]["common"]))
        
        returnstring += k+"\t"+str(pval)+"\t"+str(occurrences)+"\t"+all_genes+"\t"+common_genes+"\n"
    
    df = pd.read_csv(StringIO(returnstring), index_col=0, sep="\t").sort_values\
    (by="score", ascending=False)
    
    return df


def tableize_aggregated(dictlike, terms=None, not_found=1):
    
    """The structure of dictlike is as follows:
    
     dictlike   keys      keys(2)
    {bestof} - {term1} - "exp condition": float
                         "exp condition2": float
                         "hightes pval": float
             - {term2} - ...
    """
    if not isinstance(dictlike, dict):
        print(f"'dictlike' must be a dictionary, as the name suggests.")
        return -1
    
    from io import StringIO
    
    returnstring=""
    keys = sorted(dictlike.keys())
    
    # pulling all possible conditions that have been picked up
    exp_conditions = set()
    for retrieved_term in keys: # keys are the retrieved terms
        exp_conditions = exp_conditions.union(set(dictlike[retrieved_term].keys()))
    exp_conditions = exp_conditions.difference(set(["highest score"]))
    exp_conditions = exp_conditions.difference(set(["genes"]))
    
    
    exp_conditions = sorted(list(exp_conditions))
    returnstring += "term\t"
    returnstring += "\t".join(exp_conditions)
    returnstring += "\n"
    
    if terms is None:
        for retrieved_term in keys:
            returnstring += retrieved_term + "\t"
            for exp_condition in exp_conditions:
                returnstring += str(
                    dictlike[retrieved_term].get(exp_condition, not_found)
                )
                returnstring += "\t"
            returnstring = returnstring[:-1] + "\n"
    else:
        for retrieved_term in keys:
            if retrieved_term in terms:
                returnstring += retrieved_term + "\t"
                for exp_condition in exp_conditions:
                    returnstring += str(
                        dictlike[retrieved_term].get(exp_condition, not_found)
                    )
                    returnstring += "\t"
                returnstring = returnstring[:-1] + "\n"
            else:
                continue
    
    df = pd.read_csv(StringIO(returnstring), index_col=0, sep="\t").sort_values\
    (by="term", ascending=True)
    
    return df


def write_all_aggregated(
    directories,
    ways=[["UP"], ["DOWN"], ["UP", "DOWN"]],
    kind="all", prefix="aggregated",
    # wanted = "all", # TODO: enable custom table slicing
    ):

    """
    """

    print(f"Start invoking aggregate_results() with following parameters:\n{'-'*60}\n")
    print(f"Batch processing directories: {directories}")
    print(f"Enrichment tables following those patterns: {ways}")
    print(f"Enrichment tables for the following types: {kind}")
    print(f"Tables will be written prepended with this prefix: '{prefix}'\n")
    
    if kind == "all":
        kinds = file_types
    else:
        if isinstance(kind, str):
            kinds = [kinds]
        elif isinstance(kind, list):
            kinds = kind
        else:
            raise TypeError(f"kind parameter must be one of: {file_types}")

    TABLES = 0
    for way in ways:
        for k in kinds:
            db = aggregate_results(
                directions=way,
                kind=k,
                directories=directories
            )

            if len(db) == 0:
                continue

            table = tableize_aggregated(db)
            TABLES += 1
            
            # TODO: enable custom table slicing (not based on type)
            #if wanted != "all":
            #    table = table.loc[wanted:]
                
            del table["common"] # detailed in summary()
            
            outfile_name = f"{prefix}_{'_'.join(way)}_{k}.csv"
            table.to_csv(outfile_name, sep=sep)

    print(f"\n{'-'*60}")
    print(f"Finished. A total of {TABLES} tables were produced.")


def write_all_summarized(
    directories,
    ways=[["UP"], ["DOWN"], ["UP", "DOWN"]],
    kind="all", prefix="summary",
    # wanted = "all", # TODO: enable custom table slicing
    ):

    """
    """

    print(f"Start invoking aggregate_results() with following parameters:\n{'-'*60}\n")
    print(f"Batch processing directories: {directories}")
    print(f"Enrichment tables following those patterns: {ways}")
    print(f"Enrichment tables for the following types: {kind}")
    print(f"Tables will be written prepended with this prefix: '{prefix}'\n")

    if kind == "all":
        kinds = file_types
    else:
        if isinstance(kind, str):
            kinds = [kinds]
        elif isinstance(kind, list):
            kinds = kind
        else:
            raise TypeError(f"kind parameter must be one of: {file_types}")

    TABLES = 0
    for way in ways:
        for k in kinds:
            db = aggregate_results(
                directions=way,
                kind=k,
                directories=directories
            )

            if len(db) == 0:
                continue

            table = summary(db)
            TABLES += 1
            
            # TODO: enable custom table slicing (not based on type)
            #if wanted != "all":
            #    table = table.loc[wanted:]
                            
            outfile_name = f"{prefix}_{'_'.join(way)}_{k}.csv"
            table.to_csv(outfile_name, sep=sep)

    print(f"\n{'-'*60}")
    print(f"Finished. A total of {TABLES} tables were produced.")


def get_functional_enrichment(
    genes=None, species=None, caller_ID=session_ID,
    allow_pubmed=0, statistical_background=None, verbose=True,
    string_api_url = "https://string-db.org/api" # defaults to latest STRING release
):

    """
    Requests String functional enrichment via STRING API.
    Please see: https://string-db.org/help//api/

    Returns:
    ========
    pandas.core.frame.DataFrame: retrieved results

    """

    def say(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    species_book = {
        "mouse": 10090,
        "mus musculus": 10090,
        "10090": 10090,
        "human": 9606,
        "homo sapiens": 9606,
        "9606": 9606,
    }

    if genes is None:
        raise TypeError("A list of gene/protein identifiers is required.")

    if isinstance(species, str):
        species = species_book.get(species, None)

    if species is None:
        raise TypeError("Organism species must be provided. Mouse (10090)? Human (9606)?")

    say(f"Querying STRING. Session ID: {caller_ID}, TaxID: {species}, {len(genes)} genes/proteins.")

    output_format = "tsv"
    method = "enrichment"
    request_url = "/".join([string_api_url, output_format, method])

    if statistical_background is None:
        say("Running the analysis against a statistical",
            "background of the entire genome (default)."
            )

        params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : species,               # species NCBI identifier 
        "caller_identity" : caller_ID,     # your app name
        "allow_pubmed": 0,                 # this just seems to be ignored
        }
    else:
        say("Running the analysis against a statistical",
            "background of user-supplied terms."
            )

        params = {
        "identifiers" : "%0d".join(genes), # your proteins
        "species" : species,               # species NCBI identifier 
        "caller_identity" : caller_ID,     # your app name
        "allow_pubmed": 0,                 # this just seems to be ignored
        "background_string_identifiers": "%0d".join(statistical_background)
        }

    t0 = time()
    response = requests.post(request_url, data=params)
    t1 = time()

    say(f"STRING replied in {round((t1-t0)*1000, 2)} milliseconds.")

    df = pd.read_csv(StringIO(response.text.strip()), sep="\t", index_col=0)

    return df


def write_functional_enrichment_tables(df, databases="defaults", skip_empty=True,
                                       prefix=None, verbose=True):
    """
    For each type of functional enrichment, this **writes** a table.
    
    Params:
    =======

    databases   A <list> of <str> of wanted functional enrichments, as
         defined in settings.API_file_types.

         If "defaults", then only settings.file_types databases are produced
         (Component, Function, KEGG, Process, RCTM).

         If "all", then all possibile types of tables are produced.


    Returns:
    =======

    None
    """

    def say(*args, **kwargs):
        if verbose:
            print(*args, **kwargs)

    if databases != "all":
        if databases == "defaults":
            wanted = file_types

        elif isinstance(databases, (list, tuple)) and len(databases) != 0:
            wanted = databases

            for x in wanted:
                if x not in API_file_types:
                    print(f"*warning*: unknown database {x}")
                    wanted.pop(x)
            if len(x) == 0:
                raise TypeError("No valid database provided.")

        else:
            raise TypeError("A list of wanted databases is required")
    else:
        wanted = API_file_types

    for term in wanted:  # term like "KEGG", "Function", ...
        tempindex = df.index == term
        tempdf = df.loc[tempindex]
        tempname = f"enrichment.{term}.tsv"

        if prefix is not None:
            tempname = prefix + tempname

        # now we need to maquillage this table into the same layout of
        # tables that users retrieve via the web interface (it's not the same)

        new_col_names = {
            # old : new
            # category: None,
            "term": "#term ID",
            "description": "term description",
            "number_of_genes": "observed gene count",
            "number_of_genes_in_background": "background gene count",
            # "ncbiTaxonId": None,
            "inputGenes": "matching proteins in your network (labels)",  # guesswork
            "preferredNames": "matching proteins in your network (IDs)",  # guesswork
            # "p_value": None,
            "fdr": "false discovery rate",
        }

        new_col_order = [
            #"#term ID", #this is the index now
            "term description",
            "observed gene count",
            "background gene count",
            "false discovery rate",
            "matching proteins in your network (IDs)",
            "matching proteins in your network (labels)",
        ]

        tempdf = tempdf.rename(columns=new_col_names)
        tempdf = tempdf.set_index("#term ID")
        tempdf = tempdf[new_col_order]

        if skip_empty:
            if len(tempdf.index) == 0:
                say(f"*Notice*: skipping {tempname}: it's empty.")
                continue

        tempdf.to_csv(tempname, sep="\t")
        say(f"Table written: {tempname}")


def draw_clustermap(data, figsize=None, sort_values=None,
                    log_transform=True, log_base=10, log_na=0,
                    pval_min=None, custom_index=None, custom_cols=None,
                    unwanted_terms=None, title=None, title_size=24,
                    savefig=False, outfile_name="aggregated results.png",
                    dpi=300, readable=False, return_table=False,
                    savefigGUI=False, **kwargs
    ):
    
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
    # TODO FIX: if NOT log transforming, this should be < pval_min
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

    # just for the GUI
    if savefigGUI:
        clus.savefig(savefigGUI, dpi=dpi, bbox_inches="tight")

    if return_table:
        return clus, table
    else:
        return clus


def draw_bubbleplot(data, figsize=None, title=None, title_size=24, dpi=300,
                    savefig=False, outfile_name="aggregated summary.png",
                    return_table=False, savefigGUI=False, term_height=.25, 
                    sizes=(5,200), sort_values=None, **kwargs
    ):
    
    """
    
    Params
    ------
    
    data      A <pandas.DataFrame object>, or a filename. If a filename is
              given, this will try to load the table as if it were produced
              by tableize_aggregated() and saved with pandas.DataFrame.to_csv()
              as a .csv file, with ',' as separator.
              
    sort_values   A <str> or <list> of <str>. Sorts the table accordingly. Please
              *also* set row_cluster=False (see below, seaborn.clustermap() 
              additional parameters), otherwise columns will still be clustered.
              
    savefig   If True, a picture of the heatmap will be saved in the current directory.
    
    outfile_name  <str> The file name of the picture saved.
    
    return_table   If set to True, it also returns the table that was manipulated
              internally to draw the heatmap (with all modifications applied).
              Returns: <seaborn.matrix.ClusterGrid>, <pandas.core.frame.DataFrame>

    sizes    <tuple>. Defines minimum and maximum bubble sizes for seaborn's
             scatterplot.

    term_height  <float>. This picture is not managed, as all other figures,
              by plt.figure(figsize=(<float>, <float>)), but some parameter
              measured in inches instead. Defaults to 0.25.

    **kwargs    The drawing is performed by seaborn.clustermap(). All additional
              keyword arguments are passed directly to it, so that the final picture
              can be precisely tuned. More at:
              https://seaborn.pydata.org/generated/seaborn.clustermap.html

    """

    if isinstance(data, str):
        table = pd.read_csv(data, index_col=0)
    else:
        table = data.copy()

    if sort_values is not None:
        table = table.sort_values(by=sort_values)

    # draw --------------------------------
    
    # there is no "readable" here. The whole thing is drawn
    # by thinking in inches, and not by specifying picture size.
    # we therefore adapt the height in funcion of how many
    # terms we need to accommodate, to have all of them readable.
    # <term_height> can be tuned to adjust this (defaults to .25)

    pair = sns.PairGrid(
        data=data.sort_values(ascending=True, by="score"),
        height=len(data.index) * term_height, #this nightmare is in inches
        hue="score",
        x_vars=["all genes"],
        y_vars=["ID"],
    )

    pair.map(sns.scatterplot, data=data, size="common genes", sizes=sizes)
    plt.legend()

    # end of draw --------------------------

    if title is not None:
        plt.title(title, size=title_size)

    if savefig:
        pair.savefig(outfile_name, dpi=dpi, bbox_inches="tight")

    # just for the GUI
    # if it is true, it gets assigned a <str> value, from the GUI
    if savefigGUI:
        pair.savefig(savefigGUI, dpi=dpi, bbox_inches="tight")

    if return_table:
        return pair, table
    else:
        return pair


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
        self.say(f"{'=' * 80}\n")

        #here we go!
        t0 = time()
        print("*DEBUG* Here we go!!!")
        print("=" * 80, "\n")
        print("=" * 80, "\n")
        print("=" * 80, "\n")
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
            all_gene_list = list(tempdf.index)

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

        #now automatically running the aggregation of functional enrichment

        self.say(f"Getting directories to process.")
        self.say(f"*Python*: dirs = get_dirs()")
        dirs = get_dirs()
        produced_tables = []

        for term in file_types: # these are legit and recognized by the other functions

            self.say(f"\nAggregating data for: {term}")
            self.say(f"{'='*35}")

            self.say(f"*Python*: db = aggregate_results(dirs, kind='{term}')")
            db = aggregate_results(dirs, kind=term)

            self.say(f"*Python*: tableize_aggregated(db)")
            df = tableize_aggregated(db)
            outfile_results_name = f"{term}_results.csv"
            self.say(f"*Python*: df.to_csv('{outfile_results_name}')")
            df.to_csv(outfile_results_name)
            produced_tables.append(outfile_results_name)

            self.say(f"*Python*: res = summary(db)")
            res = summary(db)
            outfile_summary_name = f"{term}_summary.csv"
            self.say(f"*Python*: res.to_csv('{outfile_summary_name}')")
            res.to_csv(outfile_summary_name)
            produced_tables.append(outfile_summary_name)

        self.say(f"\n{'='*80}\nFinished aggregating all terms. Tables produced:")
        for name in sorted(produced_tables):
            self.say(name)


# ======================= advanced API features =========================

# def remap_identifiers(listlike, url="https://string-db.org/api",
#                       SPECIES=10090, session_ID="dummy session",
#     ):
#     """
#     From STRING doc:
    
#     You can call our STRING API with common gene names,
#     various synonyms or even UniProt identifiers and accession
#     numbers. However, STRING may not always understand them
#     which may lead to errors or inconsistencies. Before
#     using other API methods it is always advantageous to map
#     your identifiers to the ones STRING uses. In addition,
#     STRING will resolve its own identifiers faster, therefore
#     your tool/website will see a speed benefit if you use them.
#     For each input protein STRING places the best matching
#     identifier in the first row, so the first line will usually
#     be the correct one.
    
#     This function translates IDs according to the instructions
#     provided in:
#     https://string-db.org/cgi/help.pl?subpage=api%23mapping-identifiers
#     """
    
#     try:
#         assert isinstance(listlike, list)
#         assert len(listlike) != 0
#     except AssertionError:
#         return None
    
#     string_api_url = url
#     output_format = "tsv-no-header"
#     method = "get_string_ids"
    
#     params = {
#     "identifiers" : "\r".join(listlike), # your protein list
#     "species" : SPECIES,
#     "limit" : 1, # only one (best) identifier per input protein
#     "echo_query" : 1, # see your input identifiers in the output
#     "caller_identity" : session_ID,
#     }
        
#     request_url = "/".join([string_api_url, output_format, method])
#     results = requests.post(request_url, data=params)
    
#     translated = []
#     for line in results.text.strip().split("\n"):
#         l = line.split("\t")
#         #input_identifier, string_identifier = l[0], l[2]
#         #print("Input:", input_identifier, "STRING:", string_identifier, sep="\t")
#         translated.append(l[2])
    
#     return translated

def remap_identifiers(listlike, url="https://string-db.org/api",
                      SPECIES=10090, session_ID="dummy session",
    ):
    return listlike
