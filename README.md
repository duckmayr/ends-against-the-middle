# ends-against-the-middle

This archive contains all code and data necessary to reproduce the results in
"Ends Against the Middle: Measuring Latent Traits When Opposites Respond the
Same  Way for Antithetical Reasons".

## Software environment

The code was last run by the authors using the following software environment:

- Manjaro Linux 21
- R 4.1

We also used RStudio version 1.4.1717. If you are using RStudio, there is an
R Project file in the root of this directory that should be of use. However,
using RStudio should not be required to reproduce our results; you may instead
run the included R scripts (for example) via the command line `Rscript` tool.

## Instruction

To reproduce our results:

 1. Install the `GGUM` package included in the `code/` directory. This is a
    lightly patched version of the v0.4-2 previously available on CRAN that
    handles some changes to R beginning with R 4.1 which broke that version of
    `GGUM`. The version of the package currently available on CRAN, v0.4-3,
    should be equivalent to our lightly patched version of version 0.4-2.
    All other required R extension packages should be available on CRAN.
 2. Run the R script `code/01-theory-plots.R`. This reproduces Figures 1, 2, 3,
    A1, A2, A3, and B1.
 3. Run the R script `code/02-house-application.R`. This reproduces Figures 9,
    10, E1, E2, E3, E4, E5, I1, I2, I3, and I4 and Table E1.
 4. Run the R script `code/03-court-application.R`. This reproduces Figures 7
    and 8 and Table H1.
 5. Run the R script `code/04-immigration-application.R`. This reproduces
    Figure 6.
 6. Run the R script `code/05-likelihood-comparison.R`. This reproduces Table 2
 7. Run the R script `code/06-IFE-application.R`. This reproduces Figure K1 and
    Table K1.
 8. Run the R script `code/07-temperature-demonstration.R`. This reproduces
    Figure 4.
 9. Run the R script `code/08-second-dimension-concerns.R`. This reproduces
    Figures 5, F1, F2, F3, and F4 and Tables F1 and F2.
10. Run the R script `code/09-monotonic-simulation.R`. This reproduces Tables 1
    and D1.
11. Run the R script `code/10-software-comparison.R`. This reproduces Table C2.
12. Run the R script `code/11-non-monotonicity.R`. This reproduces Figure J1.

