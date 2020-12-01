# -*- coding: utf-8 -*-
# Author: Manzini Stefano; stefano.manzini@gmail.com
# 171120



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

    draw_clustermap,
    Aggregation,

    # String API
    session_ID,
    get_functional_enrichment,
    write_functional_enrichment_tables
    )

from .settings import sep, file_types



