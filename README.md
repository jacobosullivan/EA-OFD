# EA-OFD

## Jacob O'Sullivan
## j.osullivan@qmul.ac.uk | j.osullivan@zoho.com

Temporally robust spatial structure in ecosystems explained by local biodiversity regulation:

Environment Agency data on Macroinvertebrate, Macrophyte and Diatom observations in England and R scripts to reproduce analyses of the catchment-scale temporal occupancy frequency distribution

<p align="center">
<img width="450" height="310" src="https://github.com/jacobosullivan/EA-OFD/blob/master/EA-OFD_icon.png?raw=true">
</p>

Raw data are publicly available at https://environment.data.gov.uk/ecology-fish/downloads/ while catchment identifiers were taken from https://environment.data.gov.uk/catchment-planning/. Here we include both raw data and cleaned databases following various merging steps.

## Directory structure:

```
EA-OFD
| EA-OFD_icon.png
| LICENCE
| README.md
| README_database.md
└─── Diatoms
|    | DIAT_DATA_CLEANDED.zip # cleaned data generated in accompanying R script
|    | DIAT_OPEN_DATA.zip # raw data downloaded from environment.data.gov.uk
|    └─── README.md # description of raw data structure
|
└─── Macroinverts
|    | INV_DATA_CLEANDED.zip
|    | INV_OPEN_DATA.zip
|    └─── README.md
|
└─── Macrophytes
|    | MACP_DATA_CLEANDED.zip
|    | MACP_OPEN_DATA.zip
|    └─── README.md
|
└─── TaxonInfo
|    | OPEN_TAXON_INFO.zip # database of taxon identifiers and available taxonomic information
|    └─── README.md
|
└─── scripts
|    | EA_diatoms.R # clean and analyse diatom data
|    | EA_macroinverts.R # clean and analyse macroinvertebrate data
|    | EA_macrophytes.R # clean and analyse macrophyte data
|    | LSPD_realisation.R # single realisation of LSPDM
|    | analysis_functions.R # various functions used to analyse the data, applicable to all datasets
|    | fitLSPD.R # procedure for fitting LSPD to binary species-by-site table
|    | fitTimeScales.R # procedure for estimating m/alpha from binary species-by-site-by-time array
|    | tokeshi_method.R # tokeshi's (1992) test of biomodality
|    └─── README.md

```
