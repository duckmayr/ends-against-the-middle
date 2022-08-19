#!/bin/bash

function usage {
    echo "Usage: ./reproduce-results.sh [--options]"
    echo ""
    echo "--options accepted are"
    echo "    --help             Print usage and exit"
    echo "    --install          Install R package versions used by the authors"
    echo "    --no-comparison    Do not reproduce Table C2 comparing to 'bggum'"
    echo "                       to other software packages"
    echo ""
    echo "These options can be shortened to -h, -i, and -n respectively"
    exit 1
}

for arg in "$@"
do
    case "$arg" in
        "-h" | "--help")
            usage
            ;;
        "-i" | "--install")
            echo "Installing required packages"
            Rscript --no-echo "code/zz-install.R"
            ;;
        "-n" | "--no-comparison")
            COMPARISON="DON'T"
            ;;
        *)
            usage
            ;;
    esac
done

echo "Reproducing Figures 1-3, A1-A3, and B1"
Rscript "code/01-theory-plots.R"

echo "Reproducing Figures 9-10, E1-E5, and I1-I4 and Table E1"
# Rscript "code/02-house-application.R"

echo "Reproducing Figures 7-8 and Table H1"
# Rscript "code/03-court-application.R"

echo "Reproducing Figure 6"
# Rscript "code/04-immigration-application.R"

echo "Reproducing Table 2"
# Rscript "code/05-log-likelihood-comparison.R"

echo "Reproducing Figure K1 and Table K1"
# Rscript "code/06-IFE-application.R"

echo "Reproducing Figure 4"
# Rscript "code/07-temperature-demonstration.R"

echo "Reproducing Figures 5 and F1-F4 and Tables F1 and F2"
# Rscript "code/08-second-dimension-concerns.R"

echo "Reproducing Tables 1 and D1"
# Rscript "code/09-monotonic-simulation.R"

if [ "$COMPARISON" == "DON'T" ]
then
    echo "Skipping the software comparison due to option -n or --no-comparison"
else
    echo "Reproducing Table C2"
    # Rscript "code/10-software-comparison.R"
fi


echo "Reproducing Figure J1"
# Rscript "code/11-non-monotonicty.R"

