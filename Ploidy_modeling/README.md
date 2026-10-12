# Viral ploidy and allele frequency modeling

This repository contains an R Markdown analysis of how the number of viral genomes in an infectious unit affects selection on a mutant genome. The model accompanies the allele frequency modeling section of an influenza A virus aggregation study. It compares fixed ploidy, variation in ploidy among infectious units, and ploidy that changes between generations.

The repository contains a theoretical model and its figures. It does not contain the experimental measurements, sequencing data, or analysis pipelines for the other experiments described in the study's Methods.

## Model

At generation `t`, `q` is the frequency of mutant genome `A2` and `p = 1 - q` is the frequency of reference genome `A1`. Genomes are randomly assembled into infectious units. Selection acts on each unit, and its mutant contribution is proportional to the fraction of `A2` genomes it contains. The model assumes an effectively infinite population and equal baseline progeny output per infectious unit before the unit's fitness is applied. It does not explicitly simulate genetic recombination or stochastic drift.

For a unit of fixed ploidy `c`, an `A2` mutation is **recessive** when mixed units have reference fitness `W1`; it is **dominant** when mixed units have mutant fitness `W2`. The mutant frequency in the next generation is:

| Mode | Frequency update |
| --- | --- |
| Recessive `A2` | `q_next = [W1*q + (W2 - W1)*q^c] / [W1 + (W2 - W1)*q^c]` |
| Dominant `A2` | `q_next = W2*q / [W1*p^c + W2*(1 - p^c)]` |

`W1` is set to 1. The main fixed-ploidy comparisons use `W2 = 1.5` for beneficial `A2` and `W2 = 0.2` for deleterious `A2`. They begin at `q = 0.5`, cover integer ploidies 1 through 10, and run for 10 generations. Separate recessive lethal examples use `W2 = 0`.

### Variation among infectious units

The fixed-versus-Poisson comparison uses fixed ploidy `c = 2` and a zero-truncated Poisson distribution, which includes only units with at least one genome. Its rate `lambda` is chosen so that the mean *nonempty* unit has ploidy 2: `lambda / (1 - exp(-lambda)) = 2`, giving `lambda ≈ 1.593624`. This is a comparison at matched mean ploidy. For the recessive lethal mutant, the probability `q^c` of an all-mutant fixed unit is replaced by `S(q) = [exp(-lambda*(1 - q)) - exp(-lambda)] / [1 - exp(-lambda)]`.

### Variation between generations

In a separate recessive lethal simulation, one ploidy is drawn independently from the integers 1 through 10 for each of 10 generations (`set.seed(123)`). Every unit within a given generation uses that generation's ploidy. The resulting frequency trajectory is compared with trajectories at constant arithmetic-mean ploidy, constant harmonic-mean ploidy, and *effective ploidy*. Effective ploidy is found numerically as the constant value that reproduces the variable-ploidy trajectory's **final** allele frequency; it need not match intermediate generations. The fitted and mean ploidies may be nonintegers because they are used as continuous parameters in the frequency equations.

## Run the analysis

The analysis was developed with R 4.3.3. Install the packages required to render the notebook and export its plots:

```r
install.packages(c("rmarkdown", "ggplot2", "tidyr", "dplyr", "svglite"))
```

From the repository root, run:

```r
rmarkdown::render("combined_ploidy_models_A4_8pt_ploidy2_comparison_ZTP_HM_updated.Rmd")
```

Rendering writes the notebook output and exports individual figures as SVG, PDF, and 600 dpi TIFF files to `ploidy_individual_plots_A4_8pt_final_layout/`. PDF export requires an R graphics installation with Cairo support. No external input dataset is required for these simulations.

## Current notebook scope

The supplied combined notebook also generates exploratory overdominance, underdominance, dominant-negative, and diploid comparison plots. Those analyses are outside the dominant/recessive manuscript model summarized above. 
