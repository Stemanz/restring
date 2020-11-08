# restring
Functional enrichment terms aggregator

```restring``` is designed to aggregate functional enrichment data, preferably from String-compliant table delimited text (```.tdt```) tables. ```restring``` is currently under **active** development, including this doc, and it's still **not ready for use**.

doc: coming soon

examples: coming soon

## Overview
```restring``` aggregates functional enrichment data from the table delimited text (```.tdt```) tables produced by [String](https://string-db.org/). It returns, in table-friendly format, aggregated results from analyses of multiple comparisons.

What KEGG pathway was found in which comparisons? What pvalues? What DE genes annotated in that pathway were shared in those? How can I simultaneously show results for all my experimental groups, for all terms, all at once? This can all be managed by ```restring```.

---

Modern high-throughput -omic approaches generate huge lists of differentially expressed (DE) genes/proteins, which can in turn be used for functional enrichment studies. Manualy reviewing a large number of such analyses is time consuming, especially for experimental designs with more than a few groups. Let's consider this experimental setup:

![](https://github.com/Stemanz/restring/raw/main/images/exp_design.jpg)

This represents a fairly common experimental design, but manually inspecting functional enrichment results for such all possible combinations would require substantial effort. Let's take a look at how we can tackle this issue with ```restring```.
