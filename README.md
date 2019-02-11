# inhibitionPaperCode
This repository contains all code needed to reproduce the results and analysis presented in Keenan et al. 2019 (K19), when combined with the REddyProc repository and a copy of the FLUXNET data archive (details below).

To reproduce the analysis in K19, the FLUXNET data repository (download from http://fluxnet.fluxdata.org) is first processed with partitioning algorithms. 
This processing is performed by calling:
runPartitioning_master_batch.sh
and
runPartitioning_master_batch_wInhib.sh

Each script runs batch partitioning using a modified version of REddyProc (https://github.com/trevorkeenan/REddyProc), either with or without allowing for the inhibition of daytime respiration in the light.

Post-processing code, including code to generate all figures from Keenan et al. 2019 Nature Ecology and Evolution, can be found in:
./code_batch/code_postPartitioning/

Note that due to data policy restrictions, the original FLUXNET data used in K19 is not redistributed here, but is freely available for download at www.fluxnet.fluxdata.org.
Similarly, the Isotope data from Wehr et al., used in Figure 4, should be obtained from Rick Wehr and is published in Wehr et al., 2016: https://www.nature.com/articles/nature17966



