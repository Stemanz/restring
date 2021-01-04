# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com

__version__ = "0.2.040121 beta"  # dalla buonanima di __init__.py

from gears import (
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

    draw_clustermap,
    Aggregation,

    # String API
    session_ID,
    get_functional_enrichment,
    write_functional_enrichment_tables,
    )

from guigears import ReorderableListbox

from settings import (
    file_types,
    API_file_types,
    header_table, 
    sep,
    PATH
)

import tkinter as tk
from tkinter import messagebox
from tkinter.scrolledtext import ScrolledText
import webbrowser
import re
import os
from os.path import isdir
from io import StringIO
import pandas as pd
from math import log
from glob import glob
import requests
from random import choice
from time import time, sleep
from matplotlib import pyplot as plt
import seaborn as sns

#TODO: unify error messages window popups

# super-globals
files = None
working_directory = None
SPECIES = 10090  #defaults to mouse
ANALYSIS_TYPE = ["UP", "DOWN"]

about_text = """\

         ╔════════════════════════════════════════════╗
         ║                 reString                   ║
         ║      functional enrichment aggregator      ║
         ╚════════════════════════════════════════════╝

                     Stefano Manzini, PhD
                      Marco Busnelli, PhD
                         Alice Colombo

                    Head of the laboratory:
                    Prof. Giulia Chiesa, PhD

               Università degli Studi di Milano
Laboratorio di farmacologia delle dislipidemie e dell’aterosclerosi

"""

instructions = """\

Restring is a tool for aggregating functional enrichment results from
multiple comparisons at once.

For a quick and comprehensive tour of reString:  Help > Online doc
To get sample data to familiarize with reString: File > Download sample data

Quick start:
============
1 - Select the files that contain your DE genes [Open files..]
2 - Choose a folder to save your results [Set folder]
3 - Automatically retrieve functional enrichment data [New analysis]

reString will retrieve automatically functional enrichment data and
will produce aggregate summaries of all findings.

After the analysis, you can draw and personalize clustermaps choosing
Analysis > Draw clustermap

"""

def say(*args, sep=" ", end="\n"):
    output_window_text.configure(state='normal')

    if len(args) > 1:
        for arg in args[:-1]:
            output_window_text.insert(tk.END, arg + sep)
        output_window_text.insert(tk.END, args[-1] + end)
    elif len(args) == 1:
        output_window_text.insert(tk.END, args[0] + end)
    elif len(args) == 0:
        output_window_text.insert(tk.END, end)
    
    output_window_text.configure(state='disabled')
    output_window_text.see(tk.END)
    root.update()


def display_image(imagepath, *args, **kwargs):
    webbrowser.open("file:///"+imagepath)


def dummy_command():
    say("*debug*: A dummy command was issued.")


# def display_about_screen():
#     messagewin = tk.Toplevel(root) # defines a new window, on top
#     messagewin.resizable(True, True)
#     messagewin_text = tk.Text(
#         messagewin,
#         width=68,
#         height=16,
#         font=("consolas", 18)
#     )
#     messagewin_text.pack()
#     messagewin_text.insert(tk.END, about_text)
#     messagewin_text.configure(state='disabled')
#     button = tk.Button(messagewin, text="OK", command=messagewin.destroy)
#     button.pack()
#     messagewin.title("about reString")


# TODO: set button background!!! it doesn't seem to work
def display_about_screen():
    background_color = "#212836"
    displayed_image_window = tk.Toplevel()
    displayed_image_window.resizable(False, False)
    #displayed_image_window.minsize(width=600, height=800)
    displayed_image_window.title("reString: about")
    
    #tk.Toplevel does not have .create_image method, tk.Canvas does
    canvas = tk.Canvas(
        displayed_image_window,
        width=600, height=800,
        bg=background_color,
        )
    #canvas.pack(expand=tk.YES, fill=tk.BOTH)
    image = tk.PhotoImage(file="credits.png")
    canvas.create_image(0, 0, image=image, anchor=tk.NW)
    canvas.image = image # to avoid garbage collection!
    canvas.pack()

    ok_button = tk.Button(
        displayed_image_window,
        text="OK",
        command=displayed_image_window.destroy,
        anchor="w",
        bg=background_color, # not working: MacOS
        bd=0, #Border width in pixels. Default is 2. (this parameter works)

    )
    ok_button_window = canvas.create_window(
        275, 750, 
        anchor="nw", window=ok_button,
    )
    

# this is basically a whole new application -.-'
def window_draw_heatmap():
    data = False
    data_vanilla = False
    CUSTOMINDEX = False
    GLOBALCUSTOMINDEX = False
    COLORDER = False

    window_drawmap = tk.Toplevel(root) # defines a new window, on top
    window_drawmap.resizable(True, True)
    window_drawmap.geometry("800x600+50+50")
    
    # TODO: is it necesseray to super-wrap the text into a tk.Frame?
    window_drawmap_text_frame = tk.Frame(
        window_drawmap,
        width=682, # 682 px, 68 tk.Text points?
        height=338, # 338 px, 16 tk.Text points?
        )
    window_drawmap_text_frame.pack(side=tk.TOP)

    window_drawmap_text = ScrolledText(
        window_drawmap_text_frame,
        width=68,
        height=16,
        font=("consolas", 18)
    )
    window_drawmap_text.pack()

    def initialize_dataframe():
        nonlocal data
        nonlocal data_vanilla
        nonlocal CUSTOMINDEX
        nonlocal GLOBALCUSTOMINDEX

        if input_clustermap_filename.get() == "No clustermap input file chosen yet":
            write_to_window_drawmap("Choose an input 'results'-type table first.")
            return

        data = pd.read_csv(
            input_clustermap_filename.get(), index_col=0, sep="\t"
        )

        data_vanilla = data.copy()

        # TODO: make a nicer version of this just placing the following
        # CUSTOMINDEX and GLOBALCUSTOMINDEX in the DataFrame object
        
        CUSTOMINDEX = [True for _ in range(len(data.index))]
        GLOBALCUSTOMINDEX = {k:v for k, v in zip(data.index, CUSTOMINDEX)}

        show_dataframe_stats(data)


    def show_dataframe_stats(dataframe):
        # assumes a results-type table
        # only considers numeric columns (to the exception of table index)

        cols = [isinstance(x, float) for x in dataframe.min()]
        MIN = min(dataframe.loc[:,cols].min())
        MAX = max(dataframe.loc[:,cols].max())
        MIN_log = manzlog(MIN, base=log_base.get())
        if MIN_log != 0:
            MIN_log = -MIN_log
        MAX_log = manzlog(MAX, base=log_base.get())
        if MAX_log != 0:
            MAX_log = -MAX_log

        write_to_window_drawmap(
            f"\nInput table name: {input_clustermap_filename.get()}"
        )
        write_to_window_drawmap(
            f"{len(dataframe.index)} total terms."
        )
        write_to_window_drawmap(
            f"{len(dataframe.loc[:,cols].columns)} total comparisons."
        )

        write_to_window_drawmap(
            f"Min pval: {round(MIN, 2)}\tMax pval:{round(MAX, 2)}\n",
            "If log-transformed with current settings:\n",
            f"Min pval: {round(MIN_log, 2)}\tMax pval:{round(MAX_log, 2)}",
            sep=""
        )


    def choose_terms():
        nonlocal CUSTOMINDEX
        nonlocal GLOBALCUSTOMINDEX
        nonlocal data

        if data is False:
            write_to_window_drawmap("Load a 'results'-type table first.")
            return

        HORIZONTAL_POINT_PER_ENTRY = 10
        VERTICAL_POINT_PER_ENTRY = 25
        SCROLLX = max([len(x) for x in data.index]) * HORIZONTAL_POINT_PER_ENTRY +10
        SCROLLY = len(data.index)*VERTICAL_POINT_PER_ENTRY + 2*VERTICAL_POINT_PER_ENTRY
        if SCROLLX < 800:
            WIDTH = SCROLLX
        else:
            WIDTH = 300 #pixel
        HEIGHT = 600 #pixel

        choose_terms_window = tk.Toplevel()
        choose_terms_window.title("Choose terms to show")
        
        frame=tk.Frame(choose_terms_window, width=WIDTH, height=HEIGHT)
        frame.pack(expand=True, fill=tk.BOTH)
        
        canvas = tk.Canvas(
            frame, width=WIDTH, height=HEIGHT,
            scrollregion=(0, 0, SCROLLX, SCROLLY),
            confine=True # stay within scrollregion bounds
        )
        
        hbar=tk.Scrollbar(frame, orient=tk.HORIZONTAL)
        hbar.config(command=canvas.xview)
        hbar.pack(side=tk.BOTTOM, fill=tk.X)
        
        vbar=tk.Scrollbar(frame, orient=tk.VERTICAL)
        vbar.config(command=canvas.yview)
        vbar.pack(side=tk.RIGHT, fill=tk.Y)
        
        canvas.config(width=WIDTH, height=HEIGHT)
        canvas.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        
        # read the DataFrame content
        var_drawer = []
        COORD1 = VERTICAL_POINT_PER_ENTRY *.8
        COORD2 = VERTICAL_POINT_PER_ENTRY *.8
        for i, k in enumerate(data.index):
            tempvar = tk.BooleanVar()
            tempvar.set(CUSTOMINDEX[i])
            var_drawer.append(tempvar)
            
            temp_checkbutton = tk.Checkbutton(canvas, variable=var_drawer[i], text=k)
            canvas.create_window(COORD1, COORD2, window=temp_checkbutton, anchor="w")
            COORD2 += VERTICAL_POINT_PER_ENTRY
        
        canvas.pack(side=tk.LEFT, expand=True, fill=tk.BOTH)
       
               
        def apply_customindex():
            nonlocal CUSTOMINDEX
            nonlocal GLOBALCUSTOMINDEX
            CUSTOMINDEX = [var.get() for var in var_drawer]

            # now we need to merge newly chosen/unchosen
            # terms over to GLOBALCUSTOMINDEX
            tempindex = {k:v for k, v in zip(data.index, CUSTOMINDEX)}
            GLOBALCUSTOMINDEX = {**GLOBALCUSTOMINDEX, **tempindex}
            choose_terms_window.destroy()
            
        apply_button = tk.Button(
            choose_terms_window, text="Apply & OK",
            command=apply_customindex
        )
        apply_button.pack(side=tk.BOTTOM)


    def choose_cols_order():
        nonlocal COLORDER
        nonlocal data

        if data is False:
            write_to_window_drawmap("Choose an input 'results'-type table first.")
            return

        choose_cols_order_window = tk.Toplevel()
        choose_cols_order_window.title("Drag and drop columns to reorder")
        
        frame=tk.Frame(
            choose_cols_order_window,
        )
        frame.pack(expand=True, fill=tk.BOTH)

        if COLORDER is False: # never attempted to edit col order
            COLORDER = tk.StringVar()
            listbox = ReorderableListbox(
                frame, listvariable=COLORDER, height=30, width=50,
                font=("consolas", "16"),
            )
            for col in data.columns:
                if col != "common": # TODO: just add numeric columns only
                    listbox.insert(tk.END, col)
        else:
            temp_col_list = eval(COLORDER.get())
            COLORDER = tk.StringVar() # reinitialize
            listbox = ReorderableListbox(
                frame, listvariable=COLORDER, height=30, width=50,
                font=("consolas", "16"),
            )
            for col in temp_col_list:
                listbox.insert(tk.END, col)
        listbox.pack(fill=tk.BOTH, expand=True)
            
        OK_button = tk.Button(
            choose_cols_order_window, text="OK",
            command=choose_cols_order_window.destroy
        )
        OK_button.pack(side=tk.BOTTOM)


    def write_to_window_drawmap(*args, sep=" ", end="\n"):
        window_drawmap_text.configure(state='normal')

        if len(args) > 1:
            for arg in args[:-1]:
                window_drawmap_text.insert(tk.END, arg + sep)
            window_drawmap_text.insert(tk.END, args[-1] + end)
        elif len(args) == 1:
            window_drawmap_text.insert(tk.END, args[0] + end)
        elif len(args) == 0:
            window_drawmap_text.insert(tk.END, end)
        
        window_drawmap_text.configure(state='disabled')
        window_drawmap_text.see(tk.END)
        window_drawmap.update()


    window_drawmap_text_description=(
        #"==================      limit, # included     =====================",
        "Load 'results'-type .tsv files to draw a clustermap.",
        "\n\n",
        "After aggregating functional enrichment info from different\n",
        "comparisons, reString produces two kind of tables: 'results' and \n",
        "'summary'.",
        "\n\n",
        "'results'-type tables contain the FDR values for each term in every\n",
        "comparison, and clustermaps are great to get an overall idea of the\n",
        "whole analysis.",
        "\n\n",
        "To draw a heatmap:\n==================\n",
        "1: select a 'results'-type table\n",
        "2: choose the output image filename\n",
        "3: hit 'Draw clustermap'",
    )

    write_to_window_drawmap("".join(window_drawmap_text_description))

    # == [input file Button] & actions ==
    input_clustermap_filename = tk.StringVar()
    input_clustermap_filename.set("No clustermap input file chosen yet")

    def get_clustermap_filename():
        returnfile = tk.filedialog.askopenfilename(
            title='Choose a file to draw a clustermap from',
            filetypes=[
            ('reString tab separated values', '.tsv'),
            ],
        )
        if not len(returnfile) == 0:
            input_clustermap_filename.set(returnfile)
            # we also initialize the table!
            initialize_dataframe()

        if len(returnfile) == 0: # if one hits "cancel"
            return
    choose_file_button = tk.Button(
        window_drawmap,
        text="Choose input table file",
        command=get_clustermap_filename
    )
    choose_file_button.pack()

    input_filename_frame =  tk.Frame(
        window_drawmap,
        width=580,
        #height=40,
        bd=0,
        highlightbackground="black", highlightcolor="black",
        highlightthickness=0,
    )
    input_filename_frame.pack()
    input_filename_frame_entry = tk.Entry(
        input_filename_frame,
        width=70,  # that's in characters!
        textvariable=input_clustermap_filename, justify="left"
    )
    input_filename_frame_entry.grid(row=0, column=0)

    # == [output file Button] & actions =
    output_clustermap_filename = tk.StringVar()
    output_clustermap_filename.set("No output file defined yet")

    def get_clustermap_outfilename():
        returnfile = tk.filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[
                ("Portable Network Graphics", "*.png"),
                ("Joint Photographic experts Group", "*.jpg"),
                ("Portable Document Format", "*.pdf")
            ],
            #initialdir="this can be set!",
            title="Choose output filename",
        )
        if not len(returnfile) == 0:
            output_clustermap_filename.set(returnfile)
        else:
            return
    save_as_button = tk.Button(
        window_drawmap,
        text="Choose output filename",
        command=get_clustermap_outfilename
    )
    save_as_button.pack()

    output_filename_frame =  tk.Frame(
        window_drawmap,
        width=580,
        #height=40,
        bd=0,
        highlightbackground="black", highlightcolor="black",
        highlightthickness=0,
    )
    output_filename_frame.pack()
    output_filename_frame_entry = tk.Entry(
        output_filename_frame,
        width=70, # that's in characters!
        textvariable=output_clustermap_filename, justify="left"
    )
    output_filename_frame_entry.grid(row=0, column=0)


    # -- -- -- flaggable stuff -- -- --

    flags_area = tk.Frame(
        window_drawmap,
        height=75,
        pady=10
    )
    flags_area.pack()

    readable=tk.BooleanVar()
    readable_checkbutton = tk.Checkbutton(
        flags_area,
        text="readable",
        variable=readable,
    )
    readable_checkbutton.grid(row=0, column=0)
    #readable_checkbutton.pack()

    log_transform=tk.BooleanVar()
    logtransform_checkbutton = tk.Checkbutton(
        flags_area,
        text="log transform",
        variable=log_transform,
    )
    logtransform_checkbutton.grid(row=0, column=1)
    #logtransform_checkbutton.pack()

    row_cluster=tk.BooleanVar()
    row_cluster_checkbutton = tk.Checkbutton(
        flags_area,
        text="Cluster rows",
        variable=row_cluster,
    )
    row_cluster_checkbutton.grid(row=0, column=2)
    #logtransform_checkbutton.pack()

    col_cluster=tk.BooleanVar()
    col_cluster_checkbutton = tk.Checkbutton(
        flags_area,
        text="Cluster columns",
        variable=col_cluster,
    )
    col_cluster_checkbutton.grid(row=0, column=3)
    #logtransform_checkbutton.pack()

    # -- -- -- /flaggable stuff -- -- --

    # -- -- --    Input area    -- -- --
    input_area = tk.Frame(window_drawmap)
    input_area.pack()

    # P-value cutoff
    pval_min = tk.DoubleVar()
    pval_min.set(1)
    pval_min_label = tk.Label(
        input_area, text="P-value cutoff ", bd=0,
    )
    pval_min_label.grid(row=0, column=0)
    pval_min_entry = tk.Entry(
        input_area, textvariable=pval_min, bd=1, width=5
    )
    pval_min_entry.grid(row=0, column=1)


    log_base = tk.IntVar()
    log_base.set(10)
    log_base_label = tk.Label(
        input_area, text="  Log base ", bd=0,
    )
    log_base_label.grid(row=0, column=2)
    log_base_entry = tk.Entry(
        input_area, textvariable=log_base, bd=1, width=3
    )
    log_base_entry.grid(row=0, column=3)


    output_dpi = tk.IntVar()
    output_dpi.set(300)
    output_dpi_label = tk.Label(
        input_area, text="  DPI ", bd=0,
    )
    output_dpi_label.grid(row=0, column=4)
    output_dpi_entry = tk.Entry(
        input_area, textvariable=output_dpi, bd=1, width=4
    )
    output_dpi_entry.grid(row=0, column=5)

    # -- -- --   /Input area    -- -- --

    # -- -- --  Action Buttons  -- -- -- 

    action_buttons_area = tk.Frame(window_drawmap)
    action_buttons_area.pack()

    # [Apply] Button
    def apply_params():
        nonlocal data
        nonlocal CUSTOMINDEX

        if data is False:
            write_to_window_drawmap(f"Load a 'results'-type table first.\n")
            return

        # 0: <data> is made from scratch, from data_vanilla
        # 1: <data> is modified
        # 2: <CUSTOMINDEX> is built on the (eventually) remaining
        #    terms on the basis of the info contained in GLOBALCUSTOMINDEX
        # 3: New stats are printed

        data = data_vanilla.copy()
        if "common" in data.columns:
            del data["common"]

        # applying log-transform, if necessary
        if log_transform.get() is True:
            write_to_window_drawmap("\nThe data table is log transformed.")
            data = data.applymap(lambda x: -manzlog(x, log_base.get()))
            data = data.replace({-0: 0}) # -0 is ugly

        # scissoring away values that fall outside allowed pvalue
        if log_transform.get() is False:
            bser = [any(data.loc[x, :] <= pval_min.get()) for x in data.index]
            data = data.loc[bser,:]
        else:
            bser = [any(data.loc[x, :] >= pval_min.get()) for x in data.index]
            data = data.loc[bser,:]

        # now, this processed <data> table must be reminiscent of the choices
        # that might have been made by choose_terms()
        CUSTOMINDEX = [GLOBALCUSTOMINDEX.get(x) for x in data.index]

        # now showing some stats

        cols = [isinstance(x, float) for x in data.min()] # numeric cols
        MIN = min(data.loc[:,cols].min())
        MAX = max(data.loc[:,cols].max())
        if log_transform.get() is False:
            write_to_window_drawmap(
                f"\nMin pvalue: {MIN}\tMax pvalue: {MAX}",
                f"{sum(CUSTOMINDEX)} of {len(data.index)} total terms selected.",
                sep="\n"
            )
        else:
            write_to_window_drawmap(
                f"\nMin pvalue: {MAX}\tMax pvalue: {MIN}",
                f"{sum(CUSTOMINDEX)} of {len(data.index)} total terms selected.",
                sep="\n"
            )            

    apply_params_button = tk.Button(
        action_buttons_area, text="Apply", command=apply_params
    )
    apply_params_button.grid(row=0, column=0)


    # [Choose terms..] Button.
    choose_terms_button = tk.Button(
        action_buttons_area, text="Choose terms..", command=choose_terms
    )
    choose_terms_button.grid(row=0, column=1)

    # [Choose col order] Button.
    choose_terms_button = tk.Button(
        action_buttons_area, text="Choose col order", command=choose_cols_order
    )
    choose_terms_button.grid(row=0, column=2)


    # [Draw clustermap] Button
    def draw():
        nonlocal data
        nonlocal CUSTOMINDEX
        nonlocal COLORDER

        if COLORDER is False:
            cols_to_draw = [x for x in data.columns]
        else:
            cols_to_draw = list(eval(COLORDER.get()))
        
        if data is False: # 'the truth value of a dataframe is ambiguous..'
            write_to_window_drawmap(f"Load a 'results'-type table first.\n")
            return

        # ultra-basic check
        if output_clustermap_filename.get() == "No output file defined yet":
            write_to_window_drawmap(f"Choose an output filename first.\n")
            return

        apply_params() # applying latest settings, if unapplied

        if len(data.index) == 0:
            write_to_window_drawmap(
                f"The data table appears to be empty with current params.\n"
            )
            return

        write_to_window_drawmap(f"Drawing heatmap")

        cmap = "rocket"
        if log_transform.get() is True:
            cmap += "_r"

        draw_clustermap(
            data.loc[CUSTOMINDEX,cols_to_draw],
            row_cluster=row_cluster.get(),
            col_cluster=col_cluster.get(),
            #pval_min=pval_min.get(), #already managed by apply_params()
            #log_base=log_base.get(), #already managed by apply_params()
            readable=readable.get(),

            #log_transform=log_transform.get(), #already managed by apply_params()
            log_transform=False, # True is default, we *don't* want <data> modified

            savefigGUI=output_clustermap_filename.get(),
            dpi=output_dpi.get(),
            cmap=cmap,
            )

        write_to_window_drawmap(
            f"Heatmap in: {output_clustermap_filename.get()}"
        )
        display_image(output_clustermap_filename.get())

    draw_button = tk.Button(
        action_buttons_area, text="Draw clustermap", command=draw
    )
    draw_button.grid(row=0, column=3
        #, columnspan=2
    )


    #TODO: da implementare!
    # [Reset] Button. reloads <data> & reinitializes vars
    reset_button = tk.Button(
        action_buttons_area, text="Reset", command=initialize_dataframe
    )
    reset_button.grid(row=1, column=0)


    #TODO: da implementare!
    # [Help] Button. Triggers a Messagebox
    show_variables_button = tk.Button(
        action_buttons_area, text="Help", command=dummy_command
    )
    show_variables_button.grid(row=1, column=1)


    #TODO: da implementare!
    # [Online manual] Button. Triggers a Messagebox
    experimental_feature_button = tk.Button(
        action_buttons_area, text="Online manual", command=dummy_command,
    )
    experimental_feature_button.grid(row=1, column=2)


    # [Close] window button
    close_window_button = tk.Button(
        action_buttons_area, text="Close", command=window_drawmap.destroy
    )
    close_window_button.grid(row=1, column=3)
    
    window_drawmap.title("Draw clustermap")
    # TODO: when finished, set DRAWING_HEATMAP back to False


# TODO: make generic. this just works for <output_window_text>
# TODO: this looks a lot like say()
# TODO: replace this with say??? TODO TODO TODO and see if that works
def write_to_textwall(message="write_to_textwall(): unspecified message\n"):
    output_window_text.configure(state='normal')
    output_window_text.insert(tk.END, message)
    output_window_text.configure(state='disabled')
    output_window_text.see(tk.END)
    output_window.update()
    

def go_to_my_github():
    webbrowser.open("https://github.com/Stemanz/restring")

    
def download_sample_data():
    data_url="https://github.com/Stemanz/restring/raw/main/data/restring_sample_tables.zip"
    webbrowser.open(data_url)


filename_pattern = r"name='(.*?)'"
filename_regex = re.compile(filename_pattern)
def decode_filepath_from_askopenfiles(stringlike):
    match = re.search(filename_regex, stringlike)
    
    try:
        return match.group(1)
    except:
        return None


# modified version to write stuff to the GUI
def get_functional_enrichmentGUI(genes=None, species=None, caller_ID=session_ID,
                              allow_pubmed=0, verbose=True):

    """
    Requests String functional enrichment via STRING API.
    Please see: https://string-db.org/help//api/

    Returns:
    ========
    pandas.core.frame.DataFrame: retrieved results

    """

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

    say(f"STRING replied in {round((t1-t0) * 1000, 2)} milliseconds.")

    df = pd.read_csv(StringIO(response.text.strip()), sep="\t", index_col=0)

    return df


# modified version to write to the GUI
def write_functional_enrichment_tablesGUI(df, databases="defaults", skip_empty=True,
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

    if databases != "all":
        if databases == "defaults":
            wanted = file_types

        elif isinstance(databases, (list, tuple)) and len(databases) != 0:
            wanted = databases

            for x in wanted:
                if x not in API_file_types:
                    say(f"*warning*: unknown database {x}")
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

        if skip_empty:
            if len(tempdf.index) == 0:
                say(f"*Notice*: skipping {tempname}: it's empty.")
                continue

        tempdf.to_csv(tempname, sep="\t")
        say(f"Table written: {tempname}")


def aggregate_resultsGUI(
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
             Info is retrieved from either UP and/or DOWN lists, *or*
             from the list of ALL genes together.
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
            say("directions must either be ['UP, DOWN'] 'UP', 'DOWN' or 'ALL'.")
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
            for file in files: # either UP and/or DOWN; or ALL
                
                if file[:2] in [x[:2] for x in directions]:
                    PROCESSED_FILES += 1
                    say(f"\tProcessing file {file}")

                    df = pd.read_csv(file, sep="\t", index_col=ID) # sep HAS TO BE "\t"
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


# modified version to fit in the GUI
class AggregationGUI:

    """This takes care of running the whole analysis.

    TODO: more doc here
    """

    def __init__(self, files, working_directory, overwrite=False, verbose=True):

        self.working_directory = working_directory
        self.files = files
        self._files = files # compatibility, TODO: eventually make one
        self.verbose = verbose
        self.overwrite = overwrite


    #TODO: replace with safer *args version
    def say(self, message):
        message = message + "\n"
        output_window_text.configure(state='normal')
        output_window_text.insert(tk.END, message)
        output_window_text.configure(state='disabled')
        output_window_text.see(tk.END)
        root.update()


    def file_analysis(
        self,
        #species="mouse",  # now there's a global SPECIES
        reverse_direction=False,
        query_wait_time=1,
    ):


        #here we go!
        t0 = time()

        os.chdir(self.working_directory)
        # TODO: now this stuff of getting the extension may be unnecessary
        # as we have the full path for all files
        for f in self._files:
            if f.endswith(".xlsx"):
                kind = "xls"
                sep = None
                extension = -5
            elif f.endswith(".xlsx"):
                kind = "xls"
                sep = None
                extension = -4
            elif f.endswith(".csv"):
                kind = "flat"
                sep = ","
                extension = -4
            elif f.endswith(".tsv"):
                kind = "flat"
                sep = "\t"
                extension = -4
            # the following should not be openable: not liked by Libre Office
            elif f.endswith(".tdt"):
                kind = "flat"
                sep = "\t"
                extension = -4

            # the folders we create stem from the names of the files
            file_basename = os.path.basename(f)
            folder_name = file_basename[:extension]

            # can't have both filename and dir named the same way
            # if we don't strip the extension, we need to change the dir name
            # update: with the GUI, this should not happen as files must
            # have an extension. TODO: safe to remove?
            if extension is None:
                temppath = self.working_directory + f"/{folder_name}.restring"
            else:
                temppath = self.working_directory + f"/{folder_name}"

            if os.path.exists(temppath):
                if self.overwrite:
                    self.say(f"*Notice*: Path {f} exists, but we're going to overwrite it.\n")
                    os.chdir(temppath)
                else:
                    self.say(f"*Error* : path '{temppath}' already exists, and overwrite='False'.")
                    return None
            else:
                os.mkdir(folder_name)
                os.chdir(temppath)
                self.say(f"Now in: {temppath}")

            # we assume that the input files are tabular in nature.
            # first column: gene identifiers
            # second columns: fold change values

            # note for self: dropping the ../ relative positions
            # as now we have the full filepath for each file
            if kind == "flat":
                tempdf = pd.read_csv(f"{f}", index_col=0, sep=sep)
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
                self.say(f"*Notice*: no genes to process in '{f}'.")
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
                "species": SPECIES,
                "caller_ID": session_ID,
                "allow_pubmed": 0
            }

            if "UP" in ANALYSIS_TYPE:
                up_df = get_functional_enrichmentGUI(up_gene_list, **string_params)
                sleep(query_wait_time)
                write_functional_enrichment_tablesGUI(up_df, prefix="UP_")

            if "DOWN" in ANALYSIS_TYPE:
                down_df = get_functional_enrichmentGUI(down_gene_list, **string_params) 
                sleep(query_wait_time)
                write_functional_enrichment_tablesGUI(down_df, prefix="DOWN_")

            if "ALL" in ANALYSIS_TYPE:  #not checking if unique: managed by GUI
                all_df = get_functional_enrichmentGUI(all_gene_list, **string_params)
                sleep(query_wait_time)
                write_functional_enrichment_tablesGUI(all_df, prefix="ALL_")

            self.say("")
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
            #db = aggregate_resultsGUI(dirs, kind=term, PATH=self.working_directory)
            db = aggregate_resultsGUI(
                dirs,
                kind=term,
                PATH=self.working_directory,
                directions=ANALYSIS_TYPE,
            )

            self.say(f"*Python*: tableize_aggregated(db)")
            df = tableize_aggregated(db)
            outfile_results_name = f"{term}_results.tsv"
            self.say(f"*Python*: df.to_csv('{outfile_results_name}')")
            df.to_csv(outfile_results_name, sep="\t")
            produced_tables.append(outfile_results_name)

            self.say(f"*Python*: res = summary(db)")
            res = summary(db)
            outfile_summary_name = f"{term}_summary.tsv"
            self.say(f"*Python*: res.to_csv('{outfile_summary_name}')")
            res.to_csv(outfile_summary_name, sep="\t")
            produced_tables.append(outfile_summary_name)

        self.say(f"\n{'='*80}\nFinished aggregating all terms. Tables produced:")
        for name in sorted(produced_tables):
            self.say(name)


def get_files():
    global files
    files_to_open = tk.filedialog.askopenfiles(
        filetypes=[
            ('comma separated values', '.csv'),
            ('tab separated text', '.tsv'),
            ('old Excel format', '.xls'),
            ('new Excel format', '.xlsx'),
            ],
        title='Open files',
        #mode="r",
        multiple=True,
    )
    
    files_to_open_decoded = [
        decode_filepath_from_askopenfiles(f.__repr__()) for f in files_to_open
    ]
    
    if len(files_to_open_decoded) > 0:
        files = files_to_open_decoded
        for i, f in enumerate(files):
            write_to_textwall(f"Adding file {i+1}: '{f}'.\n")
    else:
        files = None
    # TODO: add some other check if one hits "cancel":is this needed?


def set_working_dir():
    global working_directory
    working_directory = tk.filedialog.askdirectory(
        mustexist=True,
        title='Set destination directory'
    )
    if not len(working_directory) == 0:
        write_to_textwall(f"'{working_directory}' set as working directory.\n")
    if len(working_directory) == 0: # if one hits "cancel"
        working_directory = None


def de_settings():
    global ANALYSIS_TYPE
    # in here we set the analysis type with respect to UP, DOWN, UP+DOWN,
    # and UP,DOWN genes. More info in the doc
    window_de_settings = tk.Toplevel(root)
    window_de_settings.resizable(True, True)
    window_de_settings.geometry("400x200+200+200")
    window_de_settings.title("DE genes analysis settings")

    window_de_settings_frame = tk.Frame(window_de_settings)
    window_de_settings_frame.pack(side=tk.TOP)

    analysis_type = tk.StringVar()
    analysis_type.set(ANALYSIS_TYPE.__repr__())
    r1 = tk.Radiobutton(
        window_de_settings_frame,
        text="Upregulated genes only (logFC > 0)",
        variable=analysis_type,
        value="['UP']",
    )
    r1.pack(anchor=tk.N)

    r2 = tk.Radiobutton(
        window_de_settings_frame,
        text="Downregulated genes only (logFC < 0)",
        variable=analysis_type,
        value="['DOWN']",
    )
    r2.pack(anchor=tk.N)

    r3 = tk.Radiobutton(
        window_de_settings_frame,
        text="Upregulated and Downregulated, separately",
        variable=analysis_type,
        value="['UP', 'DOWN']",
    )
    r3.pack(anchor=tk.N)

    # TODO: not implemented!!!
    r4 = tk.Radiobutton(
        window_de_settings_frame,
        text="All genes together",
        variable=analysis_type,
        value="['ALL']",
    )
    r4.pack(anchor=tk.N)

    def set_analysis_type():
        global ANALYSIS_TYPE
        nonlocal analysis_type
        ANALYSIS_TYPE = eval(analysis_type.get())
        say(f"Analysis type set to: {analysis_type.get()}")
        window_de_settings.destroy()
    set_analysis_type_button = tk.Button(
        window_de_settings_frame, text="Set",
        command=set_analysis_type
    )
    set_analysis_type_button.pack()


def set_species():
    global SPECIES
    window_species = tk.Toplevel(root)
    window_species.resizable(True, True)
    window_species.geometry("400x200+200+200")
    window_species.title("Set species TAXID")

    window_species_frame = tk.Frame(window_species)
    window_species_frame.pack(side=tk.TOP)


    def set_fruitfly():
        global SPECIES
        SPECIES = 7227
        say(f"Species set to Drosophila melanogaster, TAXID 7227.")
        window_species.destroy()
    fruitfly_button = tk.Button(
        window_species_frame, text="fruit fly",
        command=set_fruitfly
    )
    fruitfly_button.pack(side=tk.BOTTOM)

    def set_human():
        global SPECIES
        SPECIES = 9606
        say(f"Species set to Homo sapiens, TAXID 9606.")
        window_species.destroy()
    human_button = tk.Button(
        window_species_frame, text="human",
        command=set_human
    )
    human_button.pack(side=tk.BOTTOM)

    def set_mouse():
        global SPECIES
        SPECIES = 10090
        say(f"Species set to Mus musculus, TAXID 10090.")
        window_species.destroy()
    mouse_button = tk.Button(
        window_species_frame, text="mouse",
        command=set_mouse
    )
    mouse_button.pack(side=tk.BOTTOM)

    def set_rat():
        global SPECIES
        SPECIES = 10116
        say(f"Species set to Rattus norvegicus, TAXID 10116.")
        window_species.destroy()
    rat_button = tk.Button(
        window_species_frame, text="rat",
        command=set_rat
    )
    rat_button.pack(side=tk.BOTTOM)

    def set_zebrafish():
        global SPECIES
        SPECIES = 7955
        say(f"Species set to Danio rerio, TAXID 7955.")
        window_species.destroy()
    zebrafish_button = tk.Button(
        window_species_frame, text="zebrafish",
        command=set_zebrafish
    )
    zebrafish_button.pack(side=tk.BOTTOM)

# -- custom species selection --

    window_custom_species_frame = tk.Frame(window_species)
    window_custom_species_frame.pack(side=tk.TOP)
    input_species = tk.IntVar()
    input_species.set(SPECIES)
    input_species_label = tk.Label(
        window_custom_species_frame, text="Set manually: ", bd=0,
    )
    input_species_label.grid(row=0, column=0)
    input_species_entry = tk.Entry(
        window_custom_species_frame, textvariable=input_species, bd=1, width=5
    )
    input_species_entry.grid(row=0, column=1)

    def set_custom():
        global SPECIES
        nonlocal input_species
        SPECIES = input_species.get()
        say(f"Species set to custom: {SPECIES}")
        window_species.destroy()
    custom_button = tk.Button(
        window_custom_species_frame, text="Set",
        command=set_custom
    )
    custom_button.grid(row=0, column=2)


# core function!
def new_analysis():
    # getting files
    #files = get_files()

    if files is None:
        message = "Error: files containing DE genes must be provided.\
\n\n To download sample input files to inspect how they need to be structured, \
choose:\n\nFile > Download sample data"
        tk.messagebox.showinfo("Alert Message", message)
        return  

    # setting working dir
    #working_directory = set_working_dir()
    if working_directory is None:
        message = "Error: please choose an output folder."
        tk.messagebox.showinfo("Alert Message", message)
        return


    write_to_textwall(f"\n{'='*60}\n")
    write_to_textwall("Aggregation started.\n")
    run = AggregationGUI(files, working_directory)
    run.file_analysis()


# == GUI ==
# TODO: put every GUI-specific fuction/class into guigears.py,
# and import them ONLY if needed (if the program is directly called)
if __name__ == "__main__":

    root = tk.Tk()
    root.geometry("1024x768")
    root.resizable(False, False)
    root.title("reString")
    root.after(1, lambda: root.focus_force())

    # -- file menu --

    menubar = tk.Menu(root)

    # File
    filemenu = tk.Menu(menubar, tearoff=0)
    #filemenu.add_command(label="Open...", command=donothing)
    filemenu.add_command(label="Open...", command=get_files)
    filemenu.add_command(label="Set output folder", command=set_working_dir)
    filemenu.add_separator()
    filemenu.add_command(label="Download sample data", command=download_sample_data)
    filemenu.add_separator()
    filemenu.add_command(label="Exit", command=root.destroy)
    menubar.add_cascade(label="File", menu=filemenu)

    # Analysis menu
    analysismenu = tk.Menu(menubar, tearoff=0)
    analysismenu.add_command(label="New analysis", command=new_analysis)
    analysismenu.add_separator()
    analysismenu.add_command(label="Set species", command=set_species)
    analysismenu.add_command(label="DE genes settings", command=de_settings)
    analysismenu.add_separator()
    analysismenu.add_command(label="Draw clustermap", command=window_draw_heatmap)
    menubar.add_cascade(label="Analysis", menu=analysismenu)

    # Edit
    #editmenu = tk.Menu(menubar, tearoff=0)
    #editmenu.add_command(label = "Set species", command = donothing)
    #editmenu.add_command(label = "Overwrite", command = donothing)
    #menubar.add_cascade(label = "Edit", menu = editmenu)

    # Help menu
    helpmenu = tk.Menu(menubar, tearoff=0)
    helpmenu.add_command(label="Online doc", command=go_to_my_github)
    helpmenu.add_separator()
    helpmenu.add_command(label="About...", command=display_about_screen)
    menubar.add_cascade(label="Help", menu=helpmenu)


    # -- app window --
    # ============================================

    # top banner -- -- -- -- -- -- 
    BANNER_FRAME_HEIGHT = 150
    #banner_frame = tk.Frame(
    banner_frame = tk.Canvas(
        root,
        width=1024,
        height=BANNER_FRAME_HEIGHT,
        bd=0,      # borderwidth
        bg="blue", # background
        highlightbackground="black",
        highlightcolor="black",
        highlightthickness=0,
    )

    banner_frame.pack(side=tk.TOP)

    banner_image = tk.PhotoImage(file="restring_banner_writings.png")
    banner_frame.create_image(0, 0, image=banner_image, anchor=tk.NW)

    # menu frame -- -- -- -- -- -- 
    menu_frame = tk.Frame(
        root,
        width=256,
        height=768 - BANNER_FRAME_HEIGHT,
        bd=2,      # borderwidth
        bg="white", # background
        highlightbackground="black",
        highlightcolor="black",
        highlightthickness=0,
    )

    menu_frame.pack(side=tk.LEFT)

    BUTTON_WIDTH = 8
    button_open_files = tk.Button(
        menu_frame,
        text="Open files..",
        pady=6,
        width=BUTTON_WIDTH,
        command=get_files) ####
    button_open_files.pack()

    button_set_workdir = tk.Button(
        menu_frame,
        text="Set folder",
        pady=6,
        width=BUTTON_WIDTH,
        command=set_working_dir) ####
    button_set_workdir.pack()

    button_new_analysis = tk.Button(
        menu_frame,
        text="New analysis",
        pady=6,
        width=BUTTON_WIDTH,
        command=new_analysis)
    button_new_analysis.pack()



    # button_2 = tk.Button(
    #     menu_frame,
    #     text="TEST_write",
    #     pady=6,
    #     width=BUTTON_WIDTH,
    #     command=write_to_textwall
    # )
    # button_2.pack()

    button_go_to_my_github = tk.Button(
        menu_frame,
        text="Online doc",
        pady=6,
        width=BUTTON_WIDTH,
        command=go_to_my_github
    )
    button_go_to_my_github.pack()

    button_about = tk.Button(
        menu_frame,
        text="About",
        pady=6,
        width=BUTTON_WIDTH,
        command=display_about_screen
    )
    button_about.pack()


    #output window -- -- -- -- -- --  (768x618)
    output_window = tk.Frame(
    #output_window = tk.Canvas(
        root,
        width=768,
        height=618,
        bg="red", # background
    )
    output_window.pack(
        side=tk.RIGHT,
        #fill=tk.BOTH
    )

    output_window_text = ScrolledText(
        output_window,
        height=618,
        width=400,
        state=tk.DISABLED,
    )
    output_window_text["font"] = ('consolas', '16')
    output_window_text.pack(
        side=tk.LEFT,
        fill=tk.BOTH
        #expand=True,
    )

    # GUI is complete here --- --- --- --- --- ---
    # ============================================

    write_to_textwall(f"Welcome to reString version {__version__}\n")
    write_to_textwall(f"{'='*48}\n")

    write_to_textwall(instructions)

    root.config(menu = menubar)

    # from tkinter source:
    #root.iconify()
    root.update()
    root.deiconify()
    root.mainloop()
