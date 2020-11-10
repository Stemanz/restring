# restring
Functional enrichment terms aggregator

```restring``` is designed to aggregate functional enrichment data, preferably from String-compliant table delimited text (```.tdt```) tables. ```restring``` is currently under **active** development, including this doc, and it's still **not ready for use**.

doc: coming soon

examples: coming soon

## Overview
```restring``` aggregates functional enrichment data from the table separated values (```.tsv```) tables produced by [String](https://string-db.org/). It returns, in table-friendly format, aggregated results from analyses of multiple comparisons.

What KEGG pathway was found in which comparisons? What pvalues? What DE genes annotated in that pathway were shared in those comparisons? How can I simultaneously show results for all my experimental groups, for all terms, all at once? This can all be managed by ```restring```.

---

Modern high-throughput -omic approaches generate huge lists of differentially expressed (DE) genes/proteins, which can in turn be used for functional enrichment studies. Manualy reviewing a large number of such analyses is time consuming, especially for experimental designs with more than a few groups. Let's consider this experimental setup:

![](https://github.com/Stemanz/restring/raw/main/images/Figure%201.jpg)

This represents a fairly common experimental design, but manually inspecting functional enrichment results for such all possible combinations would require substantial effort. Let's take a look at how we can tackle this issue with ```restring```.

Our sample experimental setup has **two treatments**, given at **two time points** to **three different sample types**. Let's assume those samples are cells of different genotypes, and we'd like to mainly investigate genotype comparisons. After quantifying gene expression by RNAseq, we have DE genes for every comparison. As in many experimental pipelines, each list of DE genes is investigated with functional enrichment tools, such as [String](https://string-db.org/). But every comparison generates one or more tables. ```restring``` makes it easy to generate summary reports from all of them, automatically.

![](https://github.com/Stemanz/restring/raw/main/images/Figure%202.jpg)

---
## Procedure
### Prepping the files

Head over to [String](https://string-db.org/), and analyze your gene/protein list. Please refer to the [String documentation](https://string-db.org/cgi/help) for help. After running the analysis, hit the ![](https://github.com/Stemanz/restring/raw/main/images/analysis.png) button at the bottom of the page. This allows to download the results as tab delimited text files.

![](https://github.com/Stemanz/restring/raw/main/images/export_supported.png)

```restring``` is designed to work with the results types highlighted in green. For each one of your experimental settings, create a folder with a name that will serve as a label for it. Here's how our [sample data](https://github.com/Stemanz/restring/tree/main/sample_data) is arranged:

![](https://github.com/Stemanz/restring/raw/main/images/files.png)

Not all comparisons resulted in a DE gene list that long enough to generate functional enrichment results (see image above), thus a few comparisons _(folders)_ are missing. When the DE gene list was sufficiently long to generate results for all analyses, this is what the folder content looks like _(example of one folder)_:

![](https://github.com/Stemanz/restring/raw/main/images/folder%20content.png)

For each enrichment (```KEGG```, ```Component```, ```Function```, ```Process``` and ```RCTM```), we fed String with DE genes that were either up- or downregulated with respect of one of the genotypes of the analysis. ```restring``` **makes use** of the ```UP``` and ```DOWN``` labels in the filenames to know what direction the analysis went _(it's possible to aggregate_ ```UP``` _and_ ```DOWN``` _DE genes together)_.

It's OK to have folders that don't contain all files (if there were insufficient DE genes to produce some), like in the folder ```ctrl_t1_KO_vs_DKO```:

![](https://github.com/Stemanz/restring/raw/main/images/missing%20analyses%20folder.png)

### Aggregating the results

Once everything is set up, we can run ```restring``` to aggregate info from all the sparse results. The following example makes use of the String results that can be found in [sample data](https://github.com/Stemanz/restring/tree/main/sample_data).

```python
import restring

dirs = restring.get_dirs()
print(dirs)
```

```python
['ctrl_t0_wt_vs_DKO', 'ctrl_t0_wt_vs_KO', 'ctrl_t1_KO_vs_DKO', 'ctrl_t1_wt_vs_DKO', 'ctrl_t1_wt_vs_KO', 'treatment_t0_wt_vs_DKO', 'treatment_t0_wt_vs_KO', 'treatment_t1_KO_vs_DKO', 'treatment_t1_wt_vs_DKO', 'treatment_t1_wt_vs_KO']
```

```get_dirs()``` returns a ```list``` of all folders within the current directory, to the excepion of folders beginning with ```__``` or ```.```. We can start aggregating results with default parameters (```KEGG``` pathways for both ```UP``` and ```DOWN``` regulated genes).

```python
db = restring.aggregate_results(dirs)
```

```python
Start walking the directory structure.

Parameters
----------
folders: 10
kind=KEGG
directions=['UP', 'DOWN']

Processing directory: ctrl_t0_wt_vs_DKO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: ctrl_t0_wt_vs_KO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: ctrl_t1_KO_vs_DKO
Processing directory: ctrl_t1_wt_vs_DKO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: ctrl_t1_wt_vs_KO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: treatment_t0_wt_vs_DKO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: treatment_t0_wt_vs_KO
	Processing file DOWN_wt_enrichment.KEGG.tsv
Processing directory: treatment_t1_KO_vs_DKO
	Processing file UP_KO_enrichment.KEGG.tsv
	Processing file DOWN_KO_enrichment.KEGG.tsv
Processing directory: treatment_t1_wt_vs_DKO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv
Processing directory: treatment_t1_wt_vs_KO
	Processing file DOWN_wt_enrichment.KEGG.tsv
	Processing file UP_wt_enrichment.KEGG.tsv

Processed 10 directories and 17 files.
Found a total of 165 KEGG elements.
```

Running ```aggregate_results()``` with other parameters is possible:

```python
help(restring.aggregate_results)
```

```
# truncated output
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
```

The ```kind``` parameter is picked from the 5 supported String result tables:

```python
print(restring.settings.file_types)
```

```python
('Component', 'Function', 'KEGG', 'Process', 'RCTM')
```

To manipulate the aggregated results, it's convenient to put them into a table:

```
df = restring.tableize_aggregated(db)
```

This functions wraps the results into a handy ```pandas.DataFrame``` object, that can be saved as a table for further inspection:

```
df.to_csv("results.csv")
```

![](https://github.com/Stemanz/restring/raw/main/images/Figure%203.png)

The table contains all terms cumulatively retrieved from all comparisons _(directories/columns)_. For every term, common genes (if any) are listed. These common genes only include comparisons where the term actually shows up. If the term just appears in exactly one comparison, this is explicitly stated: ```n/a (just one condition)```. P-values are the ones retrieved from the String tables _(the lower, the better)_. Missing p-values are represented with ```1``` (this can be set to anything ```str()``` accepts when calling ```tableize_aggregated()``` with the ```not_found``` parameter).

These 'aggregated' tables are useful for charting the results (see later). There's other info that can be extracted form aggregated results, in the form of a 'summary' table:

```python
res = restring.summary(db)
res.to_csv("summary.csv")
```

![](https://github.com/Stemanz/restring/raw/main/images/Figure%204.png)

These can be useful to find the most interesting terms across all comparisons: better p-value, presence in most/selected comparisons), as well as finding the most recurring DE genes for each term.
