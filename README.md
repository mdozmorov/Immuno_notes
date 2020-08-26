# Notes on Immunology-related tools and databases

[![MIT License](https://img.shields.io/apm/l/atomic-design-ui.svg?)](https://github.com/tterb/atomic-design-ui/blob/master/LICENSEs) [![PR's Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat)](http://makeapullrequest.com) 

Mostly cancer-related. Please, [contribute and get in touch](CONTRIBUTING.md)! See [MDmisc notes](https://github.com/mdozmorov/MDmisc_notes) for other programming and genomics-related notes.

# Table of content

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Deconvolution](#deconvolution)
  - [Stand-alone tools](#stand-alone-tools)
  - [Web-based tools](#web-based-tools)
- [Signatures](#signatures)
- [Purity](#purity)
- [Data folder](#data-folder)
- [Misc Notes](#misc-notes)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Deconvolution

- Review and benchmarking of 7 immune deconvolution tools. Table 1 - overview of deconvolution methods, Table 2 - recommendations for different immune cell types. Using mix of gene expression from known proportions and types of single cells, mean expression. Access to the methods is implemented in  https://grst.github.io/immunedeconv/
    - Sturm, Gregor, Francesca Finotello, Florent Petitprez, Jitao David Zhang, Jan Baumbach, Wolf H Fridman, Markus List, and Tatsiana Aneichyk. “Comprehensive Evaluation of Computational Cell-Type Quantification Methods for Immuno-Oncology.” BioRxiv, January 27, 2019. https://doi.org/10.1101/463828. 

- Deconvolution using methylation profiles. Reference-free and reference-based. Table 1 overview of published methods. Description of Houseman algorithm, CIBERSORT. Problems with reference signatures. SVA performs well. Methylation is highly tissue-specific.
    - Teschendorff, Andrew E, and Shijie C Zheng. “Cell-Type Deconvolution in Epigenome-Wide Association Studies: A Review and Recommendations.” Epigenomics, March 14, 2017. https://doi.org/10.2217/epi-2016-0153. 

- Classical deconvolution paper, estimating cell proportions using cell signatures. IRIS cell signatures. Condition number, selection of biomarkers to minimize (improve) it. 
    - Abbas, Alexander R., Kristen Wolslegel, Dhaya Seshasayee, Zora Modrusan, and Hilary F. Clark. “Deconvolution of Blood Microarray Data Identifies Cellular Activation Patterns in Systemic Lupus Erythematosus.” Edited by Patrick Tan. PLoS ONE 4, no. 7 (July 1, 2009): e6098. https://doi.org/10.1371/journal.pone.0006098.

- Statistical derivation of deconvolution.
    - Venet, D., F. Pecasse, C. Maenhaut, and H. Bersini. “Separation of Samples into Their Constituents Using Gene Expression Data.” Bioinformatics (Oxford, England) 17 Suppl 1 (2001): S279-287.

### Stand-alone tools

- `immunedeconv` - A unified interface to immune deconvolution methods (CIBERSORT, EPIC, quanTIseq, TIMER, xCell, MCPcounter) https://icbi-lab.github.io/immunedeconv, https://github.com/icbi-lab/immunedeconv

- `MIXTURE` - an R package and a Shiny app to deconvolve cell proportions, Based on Support Vector Regression with a noise-constrained recursive feature selection. Compared with ABBAS, ABIS, CIBERSORT, DTANGLE, it detects less cell types, avoids overestimation of noisy detections. https://github.com/elmerfer/MIXTURE.App
    - Fernandez, Elmer, Yamil D Mahmoud, Florencia Veigas, Dario Rocha, Monica Balzarini, Hugo D Lujan, Gabriel A Rabinovich, and Maria Romina Girotti. “MIXTURE: An Improved Algorithm for Immune Tumor Microenvironment Estimation Based on Gene Expression Data.” BioRxiv, August 7, 2019. https://doi.org/10.1101/726562.

- `Linseed` - Mathematically elegant approach for cell type deconvolution when the number of cell types and cell signatures are unknown (Methods). Based on mutual linearity of cell type-specific genes. Mixed gene expression profiles form a simplex structure in the normalized sample-specific expression space, with corners representing normalized cell proportions. Implemented in R, Linseed (LINear Subspace identification for gene Expresion Deconvolution) is a package that provides tools and interface to explore gene expression datasets in linear space and perform complete gene expression deconvolution. https://github.com/ctlab/linseed
    - Zaitsev, Konstantin, Monika Bambouskova, Amanda Swain, and Maxim N. Artyomov. “Complete Deconvolution of Cellular Mixtures Based on Linearity of Transcriptional Signatures.” Nature Communications 10, no. 1 (May 17, 2019): 2209. https://doi.org/10.1038/s41467-019-09990-5.

- `MuSiC` - cell type deconvolution method to use scRNA-seq data (pre-determined cell types) to deconvolve bulk RNA-seq data. Gene weighting to prioritize stable and reliably expressed genes to build cell signatures. Compared with CIBERSORT, Nonnegative least squares (NNLS), BSEQ-sc. https://github.com/xuranw/MuSiC
    - Wang, Xuran, Jihwan Park, Katalin Susztak, Nancy R. Zhang, and Mingyao Li. “Bulk Tissue Cell Type Deconvolution with Multi-Subject Single-Cell Expression Reference.” Nature Communications 10, no. 1 (December 2019). https://doi.org/10.1038/s41467-018-08023-x.

- `MethylCIBERSORT` - methylation-based cell type deconvolution, using CIBERSORT model, reformatting methylation matrices for it. https://zenodo.org/record/1298968#.W9-iaHpKj-Y
    - Chakravarthy, Ankur, Andrew Furness, Kroopa Joshi, Ehsan Ghorani, Kirsty Ford, Matthew J. Ward, Emma V. King, et al. “Pan-Cancer Deconvolution of Tumour Composition Using DNA Methylation.” Nature Communications 9, no. 1 (December 2018). https://doi.org/10.1038/s41467-018-05570-1.

- [CDSeq](https://github.com/kkang7/CDSeq) - Complete Deconvolution of cell proportions using bulk RNA-seq only, without prior knowledge of cell type-specific profiles, only need the number of cells. Uses Latent Dirichlet Annotation, extends it to account for gene length and differences of RNA per cell if cell sizes differ. Similar to deriving abstract but meaningful topics from text. Compared with CIBERSORT, csSAM using synthetic and experimental datasets with known proportions. MATLAB implementation, https://github.com/kkang7/CDSeq
    - Kang, Kai, Qian Meng, Igor Shats, David M. Umbach, Melissa Li, Yuanyuan Li, Xiaoling Li, and Leping Li. “[A Novel Computational Complete Deconvolution Method Using RNA-Seq Data](https://doi.org/10.1101/496596).” BioRxiv, January 1, 2018

- `ImmQuant` - Deconvolution of human/mouse gene expression, output - immune cell proportions. Download from http://csgi.tau.ac.il/ImmQuant/download.html, run as `java -jar ImmQuant.jar`
    - Frishberg, Amit, Avital Brodt, Yael Steuerman, and Irit Gat-Viks. “ImmQuant: A User-Friendly Tool for Inferring Immune Cell-Type Composition from Gene-Expression Data.” Bioinformatics 32, no. 24 (December 15, 2016): 3842–43. https://doi.org/10.1093/bioinformatics/btw535.

- `DeconRNAseq` - deconvolution of RNA-seq datasets into cell proportions using cell signatures. Non-negative decomposition algorithm (X = AS) solved using quadratic programming.  https://bioconductor.org/packages/release/bioc/html/DeconRNASeq.html
    - Gong, Ting, and Joseph D. Szustakowski. “DeconRNASeq: A Statistical Framework for Deconvolution of Heterogeneous Tissue Samples Based on MRNA-Seq Data.” Bioinformatics (Oxford, England) 29, no. 8 (April 15, 2013): 1083–85. https://doi.org/10.1093/bioinformatics/btt090.

### Web-based tools

- `ABIS-seq` - ABsolute Immune Signal (ABIS) deconvolution, Shiny app https://giannimonaco.shinyapps.io/ABIS/ and local installationm https://github.com/giannimonaco/ABIS.
    - Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang, Christophe Carré, Nicolas Burdin, et al. “RNA-Seq Signatures Normalized by MRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types.” Cell Reports 26, no. 6 (February 2019): 1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041.

- `TCIA` - The Cancer Immunome Atlas, https://tcia.at/home. Immunophenograms, cell type fraction table of TCGA samples. Survival analysis based on immune cell signatures. All analyses are on TCGA data.

- `TIMER` - immune cell-oriented exploration of TCGA cancers. prehensive resource for systematical analysis of immune infiltrates across diverse cancer types. Exploring the abundances of six immune infiltrates (B cells, CD4+ T cells, CD8+ T cells, Neutrphils, Macrophages and Dendritic cells) with gene expression, survival, mutations, copy number variants and more. Six analysis modules: Gene correlation with immune cell proportions, immune proportions and survival, and mutations, and somatic copy number alterations, simple boxplot expression of a gene across all cancer/normal samples, correlation between two genes adjusted for tumor purity or age, deconvolution of user-provided gene expression, estimation of immune proportions in all TCGA samples. https://cistrome.shinyapps.io/timer/. Video tutorial at https://youtu.be/94v8XboCrXU
    - Li, Taiwen, Jingyu Fan, Binbin Wang, Nicole Traugh, Qianming Chen, Jun S. Liu, Bo Li, and X. Shirley Liu. “TIMER: A Web Server for Comprehensive Analysis of Tumor-Infiltrating Immune Cells.” Cancer Research 77, no. 21 (November 1, 2017): e108–10. https://doi.org/10.1158/0008-5472.CAN-17-0307. 

- `CIBERSORT` - cell type identification. Support Vector Regression. Excellent methods description. Non-log-linear space. p-value for the overall goodness of deconvolution (H0 - no cell types are present in a given gene expression profile), also Pearson and RMSE for estimating goodness of fit. References to datasets for benchmarking. https://cibersort.stanford.edu/index.php
    - Newman, Aaron M., Chih Long Liu, Michael R. Green, Andrew J. Gentles, Weiguo Feng, Yue Xu, Chuong D. Hoang, Maximilian Diehn, and Ash A. Alizadeh. “Robust Enumeration of Cell Subsets from Tissue Expression Profiles.” Nature Methods 12, no. 5 (May 2015): 453–57. https://doi.org/10.1038/nmeth.3337.

## Signatures

- Assessment of predictive biomarkers for immuno-checkpoint inhibitor response. ImmunoPhenoScore, IFNgamma signature don't perform well. Table 3 - gene lists for immune signatures. Data (FPKM values) https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220, code https://gitlab.com/bioinfo-ire-release/ici-biomarker-review/-/tree/master/
    - Pallocca, Matteo, Davide Angeli, Fabio Palombo, Francesca Sperati, Michele Milella, Frauke Goeman, Francesca De Nicola, et al. “Combinations of Immuno-Checkpoint Inhibitors Predictive Biomarkers Only Marginally Improve Their Individual Accuracy.” Journal of Translational Medicine 17, no. 1 (December 2019): 131. https://doi.org/10.1186/s12967-019-1865-8.

- `CellMaker` - cell type-specific markers, human, mouse, single-cell data. Browse and download at http://bio-bigdata.hrbmu.edu.cn/CellMarker/download.jsp, enrichment analysis at  https://yulab-smu.github.io/clusterProfiler-book/chapter3.html#cell-marker
    - Zhang, Xinxin, Yujia Lan, Jinyuan Xu, Fei Quan, Erjie Zhao, Chunyu Deng, Tao Luo, et al. “CellMarker: A Manually Curated Resource of Cell Markers in Human and Mouse.” Nucleic Acids Research 47, no. D1 (January 8, 2019): D721–28. https://doi.org/10.1093/nar/gky900.

- `10K immunomes project` - immunology reference dataset from 83 studies, 10 data types (CyTOF, proteomics, gene expression, others). Formatted (standard units of measurement) and normalized (batch-corrected, ComBat) data for visualization and download. http://10kimmunomes.org/
    - Kelly A. Zalocusky et al., “The 10,000 Immunomes Project: Building a Resource for Human Immunology,” Cell Reports 25, no. 2 (October 2018): 513-522.e3, https://doi.org/10.1016/j.celrep.2018.09.021.

- Brain immune atlas scRNA-seq resource. Border-associated macrophages from discrete mouse brain compartments, tissue-specific transcriptional signatures. http://www.brainimmuneatlas.org/index.php, https://github.com/saeyslab/brainimmuneatlas/
    - Van Hove, Hannah, Liesbet Martens, Isabelle Scheyltjens, Karen De Vlaminck, Ana Rita Pombo Antunes, Sofie De Prijck, Niels Vandamme, et al. “A Single-Cell Atlas of Mouse Brain Macrophages Reveals Unique Transcriptional Identities Shaped by Ontogeny and Tissue Environment.” Nature Neuroscience, May 6, 2019. https://doi.org/10.1038/s41593-019-0393-4.

- `DICE` - Database of Immune Cell eQTLs, Expression, Epigenomics. https://dice-database.org/

- `Haemosphere` - mostly murine immune cell signatures, downloadable data. http://haemosphere.org/datasets/show

- `63_immune_cells` - Gene expression profiles of 63 immune cell types. https://github.com/mdozmorov/63_immune_cells
 
## Purity

- `ESTIMATE` - tumor-stroma purity detection. 141 immune and stromal genes. single-sample GSEA analysis. ESTIMATE score as a combination of immune and stromal scores. - R package http://bioinformatics.mdanderson.org/estimate/rpackage.html
    - Yoshihara, Kosuke, Maria Shahmoradgoli, Emmanuel Martínez, Rahulsimham Vegesna, Hoon Kim, Wandaliz Torres-Garcia, Victor Treviño, et al. “Inferring Tumour Purity and Stromal and Immune Cell Admixture from Expression Data.” Nature Communications 4 (2013): 2612. https://doi.org/10.1038/ncomms3612.

- `ISOpureR` - Deconvolution of Tumour Profiles to purify tumor samples. Regression-based, uses purified tumor profile to estimate the proportion of tumor samples. Discussion of overfitting due to overparametrization. https://cran.r-project.org/web/packages/ISOpureR/index.html
    - Quon, Gerald, Syed Haider, Amit G Deshwar, Ang Cui, Paul C Boutros, and Quaid Morris. “Computational Purification of Individual Tumor Gene Expression Profiles Leads to Significant Improvements in Prognostic Prediction.” Genome Medicine 5, no. 3 (2013): 29. https://doi.org/10.1186/gm433.

- `PurityEst` - proportion of somatic mutations, averaged across all chromosomes, to estimate tumor purity. https://odin.mdacc.tmc.edu/~xsu1/PurityEst.html
    - Su, Xiaoping, Li Zhang, Jianping Zhang, Funda Meric-Bernstam, and John N. Weinstein. “PurityEst: Estimating Purity of Human Tumor Samples Using next-Generation Sequencing Data.” Bioinformatics (Oxford, England) 28, no. 17 (September 1, 2012): 2265–66. https://doi.org/10.1093/bioinformatics/bts365.

- `ABSOLUTE` - infers tumor purity, ploidy from SNPs, CNVs. Also detects subclonal heterogeneity. http://archive.broadinstitute.org/cancer/cga/ABSOLUTE
    - Carter, Scott L., Kristian Cibulskis, Elena Helman, Aaron McKenna, Hui Shen, Travis Zack, Peter W. Laird, et al. “Absolute Quantification of Somatic DNA Alterations in Human Cancer.” Nature Biotechnology 30, no. 5 (May 2012): 413–21. https://doi.org/10.1038/nbt.2203.

## [Data folder](data)

- [BRCA](data/BRCA) - Data from Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. “Single-Cell Immune Map of Breast Carcinoma Reveals Diverse Phenotypic States Driven by the Tumor Microenvironment.” BioRxiv, January 1, 2017. https://doi.org/10.1101/221994. - scRNA-seq of immune cells in BRCA. inDrop single-cell technology. SEQC processing pipeline, Bisquit Bayesian clustering and normalization that removes confounding technical effects. Heterogeneity of immune cell composition, clusters of immune cell subpopulations, covariance among them. Supplementary Material at https://www.biorxiv.org/content/early/2017/11/25/221994.figures-only
    - `221994-2.xlsx` - Table S2. Annotations of clusters inferred in full breast immune atlas (across all patients and tissues) and their proportions across tissues and patients.
    - `221994-3.xlsx` - Table S3. List of differentially expressed genes in clusters listed in Table S2 (sheet 1); the subset of differentially expressed immune-related genes (sheet 2).
    - `221994-2.xlsx` - Table S4. List of gene signatures (sources listed in STAR Methods)

- [Cibersort](data/Cibersort) - Data from Newman, Aaron M., Chih Long Liu, Michael R. Green, Andrew J. Gentles, Weiguo Feng, Yue Xu, Chuong D. Hoang, Maximilian Diehn, and Ash A. Alizadeh. “Robust Enumeration of Cell Subsets from Tissue Expression Profiles.” Nature Methods 12, no. 5 (May 2015): 453–57. https://doi.org/10.1038/nmeth.3337. - CIBERSORT - cell type identification. Support Vector Regression. Methods description. Non-log-linear space. p-value for the overall goodness of deconvolution (H0 - no cell types are present in a given gene expression profile), also Pearson and RMSE for estimating goodness of fit. https://cibersort.stanford.edu/index.php
    - `LM22.txt` - 547 genes X 22 immune cell types matrix of cell type specific gene signatures

- [ESTIMATE](data/ESTIMATE) - Yoshihara, Kosuke, Maria Shahmoradgoli, Emmanuel Martínez, Rahulsimham Vegesna, Hoon Kim, Wandaliz Torres-Garcia, Victor Treviño, et al. “Inferring Tumour Purity and Stromal and Immune Cell Admixture from Expression Data.” Nature Communications 4 (2013): 2612. doi:10.1038/ncomms3612. https://www.nature.com/articles/ncomms3612#supplementary-information
    - `ncomms3612-s2.xlsx` - A gene list of stromal and immune signatures
    - `ncomms3612-s3.xlsx` - A list of stromal, immune, and ESTIMATE scores in TCGA data sets. All cancers, all gene expression plaforms.

- [ImmQuant](data/ImmQuant) - Frishberg, Amit, Avital Brodt, Yael Steuerman, and Irit Gat-Viks. “ImmQuant: A User-Friendly Tool for Inferring Immune Cell-Type Composition from Gene-Expression Data.” Bioinformatics 32, no. 24 (December 15, 2016): 3842–43. https://doi.org/10.1093/bioinformatics/btw535. - Deconvolution of immune cell lineages. http://csgi.tau.ac.il/ImmQuant/downloads.html. The log2-scaled reference data files. http://csgi.tau.ac.il/ImmQuant/download.html  
    - `ImmGen.txt` - mouse reference data (Heng and Painter, 2008)
    - `DMAP.txt` - human reference data (Novershtern et al., 2011)
    - `IRIS.txt` - human reference data (Abbas et al., 2005)

- [immune_cell_signature_genes](https://github.com/caleblareau/immune_cell_signature_genes) - Repository for signature genes from Immune Cell Atlas. Also, R code for [Exploratory data analysis of SCSig collection: Signatures of Single Cell Identities](https://gist.github.com/mdozmorov/50aabcb6f5ada0d9068de3f0904b45f7#file-scsig-collection-eda)

- [quanTIseq](data/quanTIseq) - Finotello, Francesca, Clemens Mayer, Christina Plattner, Gerhard Laschober, Dietmar Rieder, Hubert Hackl, Anne Krogsdam, et al. “QuanTIseq: Quantifying Immune Contexture of Human Tumors.” BioRxiv, January 1, 2017. https://doi.org/10.1101/223180. https://www.biorxiv.org/content/early/2017/11/22/223180
    - `223180-4.xlsx` - Immune cell signatures, 170 genes x 10 immune cell types. [Source](https://www.biorxiv.org/highwire/filestream/68173/field_highwire_adjunct_files/3/223180-4.xlsx)

- [TIMER](data/TIMER) - Li, Taiwen, Jingyu Fan, Binbin Wang, Nicole Traugh, Qianming Chen, Jun S. Liu, Bo Li, and X. Shirley Liu. “TIMER: A Web Server for Comprehensive Analysis of Tumor-Infiltrating Immune Cells.” Cancer Research 77, no. 21 (November 1, 2017): e108–10. https://doi.org/10.1158/0008-5472.CAN-17-0307.
    - `ImmuneEstimation.xlsx` - Proportions of six immune cell types in all TCGA samples, downloaded from "Estimates" tab at TIMER web site https://cistrome.shinyapps.io/timer/

- `29_signatures.xlsx` - Well-Conditioned Signature Matrices for RNA-Seq (ABIS-Seq) and Microarray (ABIS-Microarray) Deconvolution. [Table S5](https://www.cell.com/cell-reports/fulltext/S2211-1247(19)30059-2#secsectitle0260). 
    - Monaco, Gianni, Bernett Lee, Weili Xu, Seri Mustafah, You Yi Hwang, Christophe Carré, Nicolas Burdin, et al. “RNA-Seq Signatures Normalized by MRNA Abundance Allow Absolute Deconvolution of Human Immune Cell Types.” Cell Reports 26, no. 6 (February 2019): 1627-1640.e7. https://doi.org/10.1016/j.celrep.2019.01.041. - Expression signatures of 29 immune subsets (FACS sorted). Modules of co-expressed, housekeeping genes (Table S3). Their robust normalization method (RLM) better suited for normalizing heterogeneous cell populations. Deconvolution for PBMC transcriptomic data. RNA-seq (ABIS-seq, 1296 genes) and microarray (ABIS-microarray, 819 genes) deconvolution panels. Outperforms five other methods (LM, non-negative LM, RLM, QP, CIBERSORT).TPM download at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011. 

- `EINAV_INTERFERON_SIGNATURE_IN_CANCER.txt` - A gene expression signature found in a subset of cancer patients suggestive of a deregulated immune or inflammatory response. [Source](http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=EINAV_INTERFERON_SIGNATURE_IN_CANCER)

- `Immune_signatures.xlsx` - List of gene signatures for "Treg", "CD8 T Cell Activation", "Anti-inflammatory", "Anergy", "Pro inflammatory", "Lipid mediators", "Glycolysis", "TCA cycle", "Pentose Phosphate Pathway", "Glycogen Metabolism", "Glucose Deprivation", "M1 Macrophage Polarization", "M2 Macrophage Polarization", "Cytolytics effector pathway", "Type I Interferon response", "Type II Interferon Response", "Hypoxia/HIF regulated", "TCell Terminal Differentiation", "G1/S", "G2/M". Sheets 2 and 3 - macrophage M1 and M2 (suppressive) signatures. Table S4 from [https://doi.org/10.1016/j.cell.2018.05.060](https://www.sciencedirect.com/science/article/pii/S0092867418307232). 

- `IRIS.xlsx` - Gene signatures of six immune cell types (T-cells, NK cells, B cells, monocytes and macrophages, Dendritic cells, Neutrophils). Microarray data). Latest gene expression data for twelve cell types, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22886. [Table S1](https://www.nature.com/articles/6364173#supplementary-information). 
    - Abbas, A. R., D. Baldwin, Y. Ma, W. Ouyang, A. Gurney, F. Martin, S. Fong, et al. “Immune Response in Silico (IRIS): Immune-Specific Genes Identified from a Compendium of Microarray Expression Data.” Genes and Immunity 6, no. 4 (June 2005): 319–31. https://doi.org/10.1038/sj.gene.6364173. 

- `TCGA_immune_classification.xlsx` - PanImmune Feature Matrix of Immune Characteristics. TCGA sample IDs, clinical characteristice, scores for key immune charactics, six immune signatures, individual cell types. [Table S1](https://www.cell.com/cms/10.1016/j.immuni.2018.03.023/attachment/49497a4e-cf22-46fe-b33e-27b91efae222/mmc2.xlsx)
    - Thorsson, Vésteinn, David L. Gibbs, Scott D. Brown, Denise Wolf, Dante S. Bortone, Tai-Hsien Ou Yang, Eduard Porta-Pardo, et al. “The Immune Landscape of Cancer.” Immunity 48, no. 4 (April 2018): 812-830.e14. https://doi.org/10.1016/j.immuni.2018.03.023.


## Misc Notes

- `Immuneatlas.org` - differences in immune cell composition and immune responses between humans, mice, and non-human primates. Example - mice have ~10-times fewer neutrophils than all primates, and ~10-times more B cells. Results can be explored at https://immuneatlas.org/#/, data at https://flowrepository.org/id/FR-FCM-Z2ZY
    - Bjornson-Hooper, Zachary B, Gabriela K Fragiadakis, Matthew H Spitzer, Deepthi Madhireddy, Dave McIlwain, and Garry P Nolan. “A Comprehensive Atlas of Immunological Differences between Humans, Mice and Non-Human Primates.” BioRxiv, March 11, 2019. https://doi.org/10.1101/574160.

- **T cell signature**: CD8A, CCL2, CCL3, CCL4, CXCL9, CXCL10, ICOS, GZMK, IRF1, HLA-DMA, HLA-DMB, HLA-DOA, and HLA-DOB
- **CTNNB1 score**: mean expression of TCF1, TCF12, MYC, EFNB3, VEGFA, and APC2, to be correlated with CD8b expression
    - Spranger, Stefani, Jason J. Luke, Riyue Bao, Yuanyuan Zha, Kyle M. Hernandez, Yan Li, Alexander P. Gajewski, Jorge Andrade, and Thomas F. Gajewski. “Density of Immunogenic Antigens Does Not Explain the Presence or Absence of the T-Cell-Inflamed Tumor Microenvironment in Melanoma.” Proceedings of the National Academy of Sciences of the United States of America 113, no. 48 (29 2016): E7759–68. https://doi.org/10.1073/pnas.1609376113.

- Markers used to type cells: NCAM1, NCR1, NKG2 (NK-cells), GNLY, PFN1, GZMA, GZMB, GMZM, GZMH (cytotoxic T, NK), FOXP3, CTLA4, TIGIT, TNFRSF4, LAG3, PDCD1 (Exhausted T cell, T-regulatory Cell), CD8, CD3, CD4 (T cells), IL7R (Naive T cells), CD19 (B cells), ENPP3, KIT (Mast cells), IL3RA, LILRA4 (plasmacytoid DC), HLA-DR, FCGR3A, CD68, ANPEP, ITGAX, CD14, ITGAM, CD33 (Monocytic Lineage).
    - Azizi, Elham, Ambrose J. Carr, George Plitas, Andrew E. Cornish, Catherine Konopacki, Sandhya Prabhakaran, Juozas Nainys, et al. “Single-Cell Map of Diverse Immune Phenotypes in the Breast Tumor Microenvironment.” Cell, June 2018. https://doi.org/10.1016/j.cell.2018.05.060.

