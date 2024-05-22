TissueSpecificTscKnockouts
==========================

<!-- Open Knowledge Link -->
 <a href="http://opendefinition.org/">
 <img alt="This material is Open Knowledge" border="0"
  src="http://assets.okfn.org/images/ok_buttons/ok_80x15_blue.png" /></a>
<!-- /Open Knowledge Link -->

# Code and Raw Data for Tissue Specific Studies on mTORC1 Function

This repository contains raw data for studies done by the [Bridges Lab](http://bridgeslab.sph.umich.edu) and our collaborators on the tissue specific roles of mTORC1.  This repository contains the data for the manuscripts detailed below.  The tag column refers the state of the dataset at the indicated time of the manuscript's publishing:

| Publication | Dataset | Tag |
|-------------|---------|-----|
| Stephenson EJ, Redd JR, Snyder DS, Tran QT, Lu B, Peloquin MJ, Mulcahy MC, Harvey I, Fisher K, Han JC, Qi N, Saltiel AR, and Bridges D. Skeletal Muscle mTORC1 Activation Increases Energy Expenditure and Reduces Longevity in Mice. bioRxiv. [doi:10.1101/720540](https://doi.org/10.1101/720540). | [![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.22193.svg)](http://dx.doi.org/10.5281/zenodo.3359627) | [muscleTSC-v1.0](https://github.com/BridgesLab/TissueSpecificTscKnockouts/releases/tag/muscleTSC-v1.0 ) |
| El Habbal N, Meyer AC, Hafner H, Redd JR, Carlson Z, Mulcahy MC, Gregg B, and Bridges D. A mouse model of adipocyte mTORC1 activation increases milk fat. bioRxiv. [doi:10.1101/2021.07.01.450596](https://www.biorxiv.org/content/10.1101/2021.07.01.450596v1). | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5106602.svg)](https://doi.org/10.5281/zenodo.5106602) | [mammaryTSC-v0.1](https://github.com/BridgesLab/TissueSpecificTscKnockouts/releases/tag/mammaryTSC-v0.1) |
| Cousineau CM, Snyder D, Redd JR, Turner S, Carr T, and Bridges D.  Reduced beta-hydroxybutyrate disposal after ketogenic diet feeding in mice. | | |

# Struture of this Repository

In general data are provided in raw form, with a R script (as a Rmd file) that analyses the data, generating html and md outputs along with figures.  There are several subfolders to this repository:

Manuscripts
: Contains all text and image files for the manuscript.  See the README file in each to identify location of underlying files and scripts

Mouse Data
: Contains mouse data including analysis and figure generation scripts.

RNASeq
: Contains RNAseq data and analyses, including figure generation scripts.

Other Published Data
: Contains analyses of data from other sources, including GEO and re-analyses of other published data.

# Licence

This TissueSpecificTscKnockouts data is made available under the Open Data Commons Attribution License: http://opendatacommons.org/licenses/by/1.0.  For more information see [LICENSE](https://github.com/BridgesLab/TissueSpecificTscKnockouts/blob/master/LICENSE)

