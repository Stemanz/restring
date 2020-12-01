import os

# separator for text tables
sep = "\t"

# working directory
PATH = os.getcwd()

# String file types
file_types = (
    "Component",
    "Function",
    "KEGG",
    "Process",
    "RCTM"
)

# Functional Enrichment from API
API_file_types = (
    "Process",
    "Component",
    "Function",
    "Keyword",
    "KEGG",
    "SMART",
    "InterPro",
    "Pfam",
    "PMID",
    "RCTM",
    "NetworkNeighborAL"
    )

# String files headers
header_table = {
    "KEGG": {
        "ID": "term description",
        "score": "false discovery rate",
        "gene name": "matching proteins in your network (labels)"
    },
    "Component" : {
        "ID": "term description",
        "score": "false discovery rate",
        "gene name": "matching proteins in your network (labels)"
    },
    "Function": {
        "ID": "term description",
        "score": "false discovery rate",
        "gene name": "matching proteins in your network (labels)"
    },
    "Process": {
        "ID": "term description",
        "score": "false discovery rate",
        "gene name": "matching proteins in your network (labels)"
    },
    "RCTM": {
        "ID": "term description",
        "score": "false discovery rate",
        "gene name" : "matching proteins in your network (labels)"
    }
}