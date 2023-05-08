# TODO: installation via Pip needs to be fixed
__version__ = "0.1.20"

# the following tests if and how I can import stuff
#from .restring import files # this works
#from restring.restring import working_directory #this also works
#from .gears import manzlog # this works
#from restring.gears import get_dirs # this also works

from restring.restring import restring_gui

from restring.gears import (
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
    draw_bubbleplot,
    Aggregation,

    # String API
    session_ID,
    get_functional_enrichment,
    write_functional_enrichment_tables,
    remap_identifiers,
)

from restring.settings import(
    file_types,
    API_file_types,
    header_table, 
    sep,
    PATH
)

__all__ = (
    # gears
    "restring_gui",
    "manzlog",
    "get_dirs",
    "aggregate_results",
    "tableize_aggregated",
    "summary",
    "write_all_aggregated",
    "write_all_summarized",
    "keep_start",
    "keep_end",
    "keep_inside",
    "prune_start",
    "prune_end",
    "prune_inside",
  
    "draw_clustermap",
    "draw_bubbleplot",
    "Aggregation",

    # String API
    "session_ID",
    "get_functional_enrichment",
    "write_functional_enrichment_tables",
    "remap_identifiers",
    
    # settings
    "file_types",
    "API_file_types",
    "header_table", 
    "sep",
    "PATH"
)
