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

We used a machine with a CPU with 12 logical cores (6 physical cores with
hyperthreading) and 64 GB of RAM (as well as a 64 GB swap file). On such a
machine, reproducing all results takes approximately 21 days, or 3 weeks;
approximately 7 days can be shaved off this time if you avoid running the
software comparison for Appendix C.3 which occurs in the file
`code/10-software-comparison.R`.

## Instructions for Reproducing Results

There are two ways you can reproduce our results.
First, we have included a [`bash`](https://www.gnu.org/software/bash/) script that you can run, which will run all R scripts required to reproduce our results without further input or work on your part.
From the terminal, simply run

```sh
./reproduce-results.sh
```

This script has the following optional arguments:

- `--help` or `-h`: Produce a message in the terminal explaining how to use the script
- `--install` or `-i`: Install the R package versions we used; this will help ensure results reproduce exactly, but please note that you may need to update your R packages after reproducing our results to restore your R addon packages to the latest versions. Please also note that this requires tools to compile R packages from source. If you are using a Mac, this may fail for the `oc` package depending on what version of `gfortran` you have installed. If so, we have included Mac x86 binaries for the appropriate version of the `oc` package.
- `--no-comparison` or `-n`: This option means that the R script to reproduce Table C2 will **not** be run; if you are uninterested in the software comparison, this will save you about a week of computation time or so.

Otherwise, you can do things piecemeal:

 1. Install the `GGUM` package included in the `code/` directory. This is a
    lightly patched version of the v0.4-2 previously available on CRAN that
    handles some changes to R beginning with R 4.1 which broke that version of
    `GGUM`. The version of the package currently available on CRAN, v0.4-3,
    should be equivalent to our lightly patched version of version 0.4-2.
    Install version 1.2 of the `oc` package; at the time of this writing it
    can't be installed using `install.packages()`, but you should be able to do
    so via the R command `devtools::install_version("oc", version = "1.2", repos = "http://cloud.r-project.org")`.
    Again, if you're using a Mac, you may need to instead install from the
    binary we have included (`oc_1.2.tgz` in this repository's root directory).
    All other required R extension packages should be available on CRAN;
    if you would like to use the same versions we used, they are as follows:
    - `bggum`: version 1.0.2
    - `coda`: version 0.19-4
    - `devtools`: version 2.4.2
    - `dplyr`: 1.0.7
    - `ggplot2`: 3.3.5
    - `ltm`: 1.1-1
    - `MCMCpack`: 1.5-0
    - `pROC`: 1.17.0.1
    - `Rcpp`: 1.0.7
    - `tibble`: 3.1.4
    - `tidyr`: 1.1.3
    - `wnominate`: 1.2.5
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
13. Run the R script `code/12-software-comparison2.R`.
    This reproduces Figures C1 and C2.
    This R script relies on intermediate output from a closed-source,
    Windows only program: `data/sim01-item-draws.txt`,
    `data/sim02-item-draws.txt`, and `data/sim03-item-draws.txt`.
    If you have access to a Windows machine and wish to recreate this
    intermediate output:
    1. Run the installer (`MCMCGGUM-1.2.0-Installer.exe`),
    2. Copy the files `data/sim-study-data.txt`, `data/sim01.mcs`,
       `data/sim02.mcs`, and `data/sim03.mcs` to the directory the program was
       installed to (by default `C:/Program Files (x86)/UCF IST/MCMC GGUM/`)
    3. Open the program
    4. Click "Import Syntax..." on the right-hand side and select `sim01.mcs`
    5. Click "Run" on the right-hand side
    6. When the run is complete, save the output via "File" > "Save As"
    7. Repeat steps iii-vi for `sim02.mcs` and `sim03.mcs`
    8. Then you can verify equality with the files `data/sim01-item-draws.txt`,
       `data/sim02-item-draws.txt`, and `data/sim03-item-draws.txt`.
