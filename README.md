# Notes on Immunology-related tools and databases

- `Haemosphere` - mostly murine immune cell signatures, downloadable data. http://haemosphere.org/datasets/show
 
## `data`

- `EINAV_INTERFERON_SIGNATURE_IN_CANCER.txt` - A gene expression signature found in a subset of cancer patients suggestive of a deregulated immune or inflammatory response. http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=EINAV_INTERFERON_SIGNATURE_IN_CANCER

## `ESTIMATE`

Yoshihara, Kosuke, Maria Shahmoradgoli, Emmanuel Martínez, Rahulsimham Vegesna, Hoon Kim, Wandaliz Torres-Garcia, Victor Treviño, et al. “Inferring Tumour Purity and Stromal and Immune Cell Admixture from Expression Data.” Nature Communications 4 (2013): 2612. doi:10.1038/ncomms3612. https://www.nature.com/articles/ncomms3612#supplementary-information

- `ncomms3612-s2.xlsx` - A gene list of stromal and immune signatures
- `ncomms3612-s3.xlsx` - A list of stromal, immune, and ESTIMATE scores in TCGA data sets. All cancers, all gene expression plaforms.

### `ImmQuant` data

- Frishberg, Amit, Avital Brodt, Yael Steuerman, and Irit Gat-Viks. “ImmQuant: A User-Friendly Tool for Inferring Immune Cell-Type Composition from Gene-Expression Data.” Bioinformatics 32, no. 24 (December 15, 2016): 3842–43. https://doi.org/10.1093/bioinformatics/btw535. - Deconvolution of immune cell lineages. http://csgi.tau.ac.il/ImmQuant/downloads.html

The log2-scaled reference data files. http://csgi.tau.ac.il/ImmQuant/download.html  

- `ImmGen.txt` - mouse reference data (Heng and Painter, 2008)
- `DMAP.txt` - human reference data (Novershtern et al., 2011)
- `IRIS.txt` - human reference data (Abbas et al., 2005)

### `quanTIseq` data

- Finotello, Francesca, Clemens Mayer, Christina Plattner, Gerhard Laschober, Dietmar Rieder, Hubert Hackl, Anne Krogsdam, et al. “QuanTIseq: Quantifying Immune Contexture of Human Tumors.” BioRxiv, January 1, 2017. https://doi.org/10.1101/223180. https://www.biorxiv.org/content/early/2017/11/22/223180

- `223180-4.xlsx` - Immune cell signatures, 170 genes x 10 immune cell types. [Source](https://www.biorxiv.org/highwire/filestream/68173/field_highwire_adjunct_files/3/223180-4.xlsx)


## `tools`

### `ImmQuant`

Download from http://csgi.tau.ac.il/ImmQuant/download.html, run as `java -jar ImmQuant.jar`


## `web`

- `TCIA` - The Cancer Immunome Atlas, https://tcia.at/home. Immunophenograms, cell type fraction table of TCGA samples. Survival analysis based on immune cell signatures. All analyses are on TCGA data.

