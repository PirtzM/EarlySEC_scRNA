# Data Processing and Analysis Pipeline

## Raw Data Alignment  
### CellRanger  (v7.1.0, 10X Genomics): [tutorial](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-ct)
Completed within a Linux terminal and with CellRanger downloaded with the proper path, run the `count` command to process raw scRNA-seq reads  
`cellranger count --id=mU29_align \
                 --transcriptome=/workdir/mgp73/MouseRef/RodortdT \
                 --fastqs=/workdir/mgp73/FASTQs/DiestrusMice_FASTQs/Sample_200521_AACLY5WM5_1_mU29 \
                 --localcores=8 \
                 --localmem=64`

Note: `cellranger count` is computationally intensive. For local PC users, it is recommended to start from the processed gene count matrices or prepared Seurat Objects

## Preparing RStudio
Install R (v4.1.1) and [Rstudio](https://rstudio-education.github.io/hopr/starting.html): approximately 15-30 min  

Install packages found in the 'sessionInfo.txt' to properly recreate the figures: approximately 60 min

## Preprocessing Seurat Objects
Refer to code in `EarlySEC_scRNA/DataPreprocessing/` to reproduce preprocessing results

### Seurat Objects for All Sample, Normal, Pre-dysplastic, and Dysplastic datasets
Use `SC_diestrus_metadirectory_final.csv` to identify sample file directories and additional sample level metadata.  
All samples were processed together according to `AllSamples_data_preprocess`. Here, you will find step by step details on how to remove ambient RNA signals with SoupX (v1.6.2), identifying doublets with DoubletFinder (v2.0.3), and preparing integrated/batch corrected objects with Seurat (v4.3.0) and harmony (v0.1.0).  

Each SEC specific dataset was subsetted from this All Sample dataset according to its stage. After subsetting, samples were re-scaled, re-integrated, and reclustered. Follow their individual preprocessing file for specific information.  

### Subsetting epithelial cells from each dataset
Follow files that include `EpiSubset` for specific workflow for subsetting epithelial cells. In general, full cell datasets of All Samples, Normal, Pre-dysplastic, and Dysplastic datasets were subsetted for their epithelial populations (identified by *Epcam* expression. Subsets were re-scaled, re-integrated, and reclustered. Follow their individual preprocessing file for specific information.  

### Figure Recreation
With available objects ([Dryad]()-active upon publication), you can recreate the figures from the manuscript using `EarlySEC_scRNA/Figure_Scripts/`.  
Each figure should only take about 5 minutes to recreate with proper objects and packages.
