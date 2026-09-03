# The hidden links between animal weapons, fighting style, and their effect on contest success: a meta-analysis

Authors: Alexandre V. Palaoro, Paulo Enrique Cardoso Peixoto <br>
Journal: Biological Reviews, DOI: 10.1111/brv.12877 (<a href="https://doi.org/10.1111/brv.12877" target="_blank">link</a>) <br>
Contact about code, data, and analyses: alexandre.palaoro@gmail.com; alexandre.palaoro@ufpr.br

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22286622.svg)](https://doi.org/10.5281/zenodo.22286622)
[![paper](https://img.shields.io/badge/paper-10.1111%2Fbrv.12877-blue)](https://doi.org/10.1111/brv.12877)
[![code license](https://img.shields.io/badge/code%20license-MIT-green)](LICENSE)
[![data license](https://img.shields.io/badge/data%20license-CC0%201.0-brightgreen)](LICENSE-DATA)

---

### This readme has been divided in four parts. First, we will talk about file structure, then the code, then the dataset, then how to reproduce the analyses.

The paper has two halves that are independent of each other. The first is an individual-based
simulation that asks how the mean difference between winners and losers changes as a trait
becomes more important in determining who wins. The second is a phylogenetic meta-analysis
that uses that mean difference as a proxy for how much weapons influence contest success.
The code in this repository is split the same way.

##### File structure:

Everything sits at the root of the repository. There are three code files and four data files.

| File | What it is |
| --- | --- |
| `Simulation model - biol rev.R` | The individual-based simulation (first half of the paper) |
| `meta_analysis-weapons.R` | The full meta-analysis: tree building, all models, heterogeneity, publication bias |
| `PLOT-FORESTPLOT-FORPAPER.R` | The forest plot (Figure 3) and the counts behind Table 3 |
| `full_dataset_weapons2.csv` | The extracted effect sizes — the main dataset |
| `correl_matrix_weap_LIN-PERF.csv` | Correlation matrix among effect sizes, weapon linear + performance subset |
| `correl_matrix_display.csv` | Correlation matrix among effect sizes, display subset |
| `correl_matrix_weap_LIN.csv` | Correlation matrix among effect sizes, fighting-style subset |

The figures are written to the working directory when you run the code; we are not uploading
them because they are all in the paper.

##### Code:

**`Simulation model - biol rev.R`** — the individual-based model. It creates two groups of 100
individuals with the same mean trait value, pairs them, and lets the pair with the larger trait
difference win more often. The maximum probability that the *weaker* rival wins is controlled by
`div1` and `div2` at the top of the file: they are set to 4 and 2, which correspond to maximum
inversion probabilities of 0.25 and 0.5. The whole thing is repeated 1000 times. The output is
the mean difference and the mean standard deviation between winners and losers under each
scenario, plotted at the end of the file. This is what justifies using the winner–loser mean
difference as a proxy for weapon importance in the meta-analysis.

**`meta_analysis-weapons.R`** — the main analysis file, and the one to start with. In order, it:

1. Builds the phylogeny. Species names are matched against the Open Tree of Life with
   `rotl::tnrs_match_names()`, three beetle entries are replaced by a higher taxon
   (`Librodor` by its subfamily, *Onthophagus taurus* by its family, *acuminatus* by
   Dynastinae), the induced subtree is pulled down, branch lengths are added with Grafen's
   method, and the tip labels are renamed. This produces Figure S2.
2. Turns the tree into a correlation matrix with `ape::vcv()`, to be used as the phylogenetic
   random effect in `metafor`.
3. Runs the three main models, each an `rma.mv()` with random effects for phylogeny, study,
   and taxon nested in pairing method, and with the study-level correlation matrix supplied
   through `R = list(study = ...)`:
   - **linear versus performance measures** of the weapon (Figure 4),
   - **display type** — whether the animal displays the weapon, the body, or nothing (Figure 5),
   - **fighting style** — how the weapon is actually used in the fight (Figure 6).
4. Calculates I² (total, phylogenetic, species, and residual) with confidence intervals for each
   model, and runs Egger's test for publication bias.
5. Repeats the linear-versus-performance comparison within Decapoda and lizards only, which is
   the only group where both measure types are well represented (Figure S3).

**`PLOT-FORESTPLOT-FORPAPER.R`** — builds the forest plot in Figure 3 and the species-by-style
counts that became Table 3. It reads the same dataset and subsets it to weapon traits.

##### Data:

###### `full_dataset_weapons2.csv`

**This file is semicolon-separated.** Read it with `read.csv("full_dataset_weapons2.csv", h=T, sep=';')`.

Each row is one effect size. There are 173 rows: 107 for weapon traits and 66 for body traits,
drawn from 49 studies and 52 species. The 107 weapon effect sizes are the ones reported in the
paper; the body rows are kept for reference and are filtered out at the top of every analysis.

Columns A–S describe the study and the species. Columns T–AO hold the numbers as they were
reported in the original papers. Columns AP–AR hold the effect size we calculated from them.

COLUMN A: art - internal code we gave each source article during screening. <br>
COLUMN B: ref - short reference for the source study (first author, or first author plus "et al"). <br>
COLUMN C: year - publication year of the source study. <br>
COLUMN D: type.measure - the trait as it was named in the original paper, kept verbatim so the extraction can be checked. <br>
COLUMN E: category - broad trait category. Levels: massa_corpo (body mass), comprimento_corpo (body length), comprimento_arma (weapon length), massa_arma (weapon mass), forca_arma (weapon force), assimetria_arma (weapon asymmetry). <br>
COLUMN F: measure2 - finer measurement type. Levels: area, asym, height, index, length, mass, strength, width. <br>
COLUMN G: measure - the measurement class used in the models. Levels: linear, performance, area, mass, asymmetry, index. The linear-versus-performance model uses only the first two. <br>
COLUMN H: between_ef - unique identifier for each effect size (1 to 173). <br>
COLUMN I: study - study identifier. Effect sizes taken from the same study share a value. This is the grouping variable for the study random effect and the key that the correlation matrices are aligned to. <br>
COLUMN J: environ - where the contests were staged. Levels: field, lab. <br>
COLUMN K: dyad.type - how contestants were paired. Levels: random, paired (i.e., size-matched). <br>
COLUMN L: order - broad taxonomic group. Levels: Acari, Arachnida, Aves, Coleoptera, Decapoda, Dermaptera, Hemiptera, Heteroptera, Hoplocarida, Lizard, Mammalia, Orthoptera. <br>
COLUMN M: genus - the genus of the study species. <br>
COLUMN N: sp - the study species. These names are what gets matched against the Open Tree of Life. <br>
COLUMN O: trait - whether the effect size refers to the weapon or to the body. Levels: weapon, body. <br>
COLUMN P: disp - what the animal displays before contact, as coded during extraction. Levels: no, tact, visual, weapon. <br>
COLUMN Q: display - display category used in the models, collapsed from column P. Levels: weapon, body, no. <br>
COLUMN R: fight.cat - fighting-function category used in the models. Levels: weapon (the weapon does the work), body (the body does the work), interm (intermediate). <br>
COLUMN S: fighting.style - how the weapon is used during the fight. Levels are single functions (push, squeeze, impact) and combinations of them (e.g. lift.push, squeeze.pierce, push.squeeze.pull.lift). <br>
COLUMN T: RHP_losers - mean trait value of losers, as reported. <br>
COLUMN U: S_losers - standard deviation of the trait in losers, as reported. <br>
COLUMN V: N_losers - number of losers. <br>
COLUMN W: RHP_winners - mean trait value of winners, as reported. <br>
COLUMN X: S_winners - standard deviation of the trait in winners, as reported. <br>
COLUMN Y: N_winners - number of winners. <br>
COLUMN Z: Chi_sq - chi-square statistic, for studies that reported one instead of means and standard deviations. <br>
COLUMN AA: N_chi - sample size for the chi-square. <br>
COLUMN AB: z-score - z statistic, same purpose. <br>
COLUMN AC: N_z - sample size for the z. <br>
COLUMN AD: p_z - p-value for the z. <br>
COLUMN AE: rs - Spearman correlation, same purpose. <br>
COLUMN AF: N_rs - sample size for the correlation. <br>
COLUMN AG: paired_t - paired t statistic. <br>
COLUMN AH: pairedt_N - sample size for the paired t. <br>
COLUMN AI: pairedt_p - p-value for the paired t. <br>
COLUMN AJ: wilcox_z - Wilcoxon z statistic. <br>
COLUMN AK: wilcox_n - sample size for the Wilcoxon. <br>
COLUMN AL: wilcox_p - p-value for the Wilcoxon. <br>
COLUMN AM: F_value - F statistic. <br>
COLUMN AN: F_N - sample size for the F. <br>
COLUMN AO: F_p - p-value for the F. <br>
COLUMN AP: yi - the effect size: the standardised mean difference between winners and losers. This is the response variable in every model. <br>
COLUMN AQ: vi - the sampling variance of yi. <br>
COLUMN AR: N - total sample size behind the effect size. <br>

"NA" cells mean the original paper did not report that quantity. Most studies reported means
and standard deviations (columns T–Y), so most of columns Z–AO are empty; those columns exist
for the minority of studies where we had to back-calculate `yi` from a test statistic instead.

> [!NOTE]
> Columns T–Y were exported from a spreadsheet that used a different decimal convention, so
> some values appear with more than one dot (for example `1.037.466.667` for 1037.466667).
> This affects the raw reported values only. The analyses never touch these columns — they use
> `yi` and `vi`, which are clean — so results are unaffected. If you want to reuse the raw
> means, parse those columns carefully.

###### The three correlation matrices

`correl_matrix_weap_LIN-PERF.csv`, `correl_matrix_display.csv` and `correl_matrix_weap_LIN.csv`
hold the correlations between effect sizes that came from the same study — the non-independence
created when one paper reports several weapon measurements on the same animals. They are passed
to `metafor` through the `R = list(study = ...)` argument of `rma.mv()`. Off-diagonal entries are
zero between studies and non-zero within a study.

There is one matrix per analysis because each analysis uses a different subset of the data, and
the matrix has to have exactly as many rows as the subset:

| Matrix | Used by | Size | Separator |
| --- | --- | --- | --- |
| `correl_matrix_weap_LIN-PERF.csv` | linear vs. performance model | 97 × 97 | comma |
| `correl_matrix_display.csv` | display model | 74 × 74 | comma |
| `correl_matrix_weap_LIN.csv` | fighting-style model | 74 × 74 | **semicolon** |

> [!WARNING]
> **The matrices are aligned by position, not by name.** They carry no row names; the code
> assigns them with `colnames(study.mat) <- weap$study` after the data have been filtered. If you
> reorder or filter the rows of `full_dataset_weapons2.csv`, the matrices will silently
> mis-align with the data and the models will be wrong. Keep the row order as it is, or rebuild
> the matrices alongside any change.

##### How to reproduce:

Run `meta_analysis-weapons.R` from the repository root, with the working directory set there, so
the `read.csv()` calls find the data files. Then run `PLOT-FORESTPLOT-FORPAPER.R` for Figure 3.
`Simulation model - biol rev.R` is standalone and can be run at any point.

The code was run in R. <br>
Packages used: <br>
metafor <br>
ape <br>
rotl <br>
tidyverse <br>
scales <br>
faux <br>
sciplot <br>
boot <br>

> [!NOTE]
> The tree-building step in `meta_analysis-weapons.R` queries the Open Tree of Life live through
> `rotl`, so it needs an internet connection, and the taxonomy it returns can change as the Open
> Tree is updated. The three manual corrections in the script are made by row index
> (`spp[11,4]`, `spp[48,4]`, `spp[51,4]`), so if the matching ever returns a different ordering,
> check `spp` before running the rest. If you want a fully frozen pipeline, save the matched
> `spp` object and the resulting tree to an `.RDS` file and load that instead.

---

## Citation

If you use anything in this repository, please cite the paper:

> Palaoro, A.V. & Peixoto, P.E.C. (2022) The hidden links between animal weapons, fighting style, and their effect on contest success: a meta-analysis. *Biological Reviews* 97(5): 1948–1966. [https://doi.org/10.1111/brv.12877](https://doi.org/10.1111/brv.12877)

If you reuse the code or the archived files directly, please also cite the archive:

> Palaoro, A.V. & Peixoto, P.E.C. (2022) *Data and code for: The hidden links between animal weapons, fighting style, and their effect on contest success — a meta-analysis* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.22286622

`CITATION.cff` in this repository holds both citations in machine-readable form — GitHub's
**Cite this repository** button (top right of the repository page) will generate APA or
BibTeX from it for you.

## License

This repository is released under two licenses, because code and data are different things.

| Content | License | File |
| --- | --- | --- |
| Analysis code — the three `.R` files | [MIT](https://opensource.org/licenses/MIT) | [`LICENSE`](LICENSE) |
| Data — `full_dataset_weapons2.csv` and the three correlation matrices | [CC0 1.0 Universal](https://creativecommons.org/publicdomain/zero/1.0/) | [`LICENSE-DATA`](LICENSE-DATA) |

In short: do whatever you like with the data, no permission needed and no attribution
legally required; reuse the code freely as long as you keep the copyright notice. Academic
norms still apply — if the data or code are useful to you, cite the paper.

## Updating the archive

This repository is linked to Zenodo. Every new GitHub release is archived automatically and
gets its own version DOI, while the concept DOI above always resolves to the newest version —
so the DOI in this README and in `CITATION.cff` never needs changing again.

To publish an update: **Releases → Draft a new release**, bump the tag (`v1.1.0`), publish.
Zenodo picks it up within a minute or two. `.zenodo.json` supplies the title, authors, ORCIDs,
keywords and the link to the article, so there is nothing to retype in the Zenodo form.

Cite the concept DOI in papers, never a version DOI.

## Reproducibility

- Everything needed to reproduce the published analyses is in this repository.
- An earlier version of this work is on bioRxiv: [10.1101/2020.08.26.268185](https://doi.org/10.1101/2020.08.26.268185). Note that the preprint title uses "contest resolution"; the published title uses "contest success".
- The archived release on Zenodo is the version of record. GitHub history may move on; the
  DOI will not.
- Two things will break a re-run if you are not careful: the row order of the dataset relative
  to the correlation matrices, and the live Open Tree of Life query. Both are flagged above.
