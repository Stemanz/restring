import os
import pandas as pd
from os.path import isdir
from math import log
from io import StringIO
from glob import glob
import requests
from io import StringIO
from random import choice
from time import time

from .settings import (
    file_types,
    API_file_types,
    header_table, 
    sep,
    PATH
)


session_ID = "".join([choice(("abcdefghijklmnopqrstuvwxyz0123456789_")) for x in range(8)])
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
    
    def say(stringlike, **kwargs):
        if verbose:
            print(stringlike, **kwargs)
    
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
    (by="score", ascending="False")
    
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
    (by="term", ascending="True")
    
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


def get_functional_enrichment(genes=None, species=None, caller_ID=session_ID,
                              allow_pubmed=0, verbose=True):

    """
    Requests String functional enrichment via STRING API.
    Please see: https://string-db.org/help//api/

    Returns:
    ========
    pandas.core.frame.DataFrame: retrieved results

    """

    def say(stringlike, **kwargs):
        if verbose:
            print(stringlike, **kwargs)

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

    string_api_url = "https://string-db.org/api"
    output_format = "tsv"
    method = "enrichment"
    request_url = "/".join([string_api_url, output_format, method])

    params = {
    "identifiers" : "%0d".join(genes), # your proteins
    "species" : species,               # species NCBI identifier 
    "caller_identity" : caller_ID,     # your app name
    "allow_pubmed": 0,                 # this just seems to be ignored
    }

    t0 = time()
    response = requests.post(request_url, data=params)
    t1 = time()

    say(f"STRING replied in {round((t1-t0)*1000, 2)} milliseconds.")

    df = pd.read_csv(StringIO(response.text.strip()), sep="\t", index_col=0)

    return df


def write_functional_enrichment_tables(df, databases="defaults",
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

    def say(stringlike, **kwargs):
        if verbose:
            print(stringlike, **kwargs)

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

    for term in wanted: # term like "KEGG", "Function", ...
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
            "preferredNames": "matching proteins in your network (IDs)", # guesswork
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

        tempdf.to_csv(tempname, sep="\t")
        say(f"Table written: {tempname}")

