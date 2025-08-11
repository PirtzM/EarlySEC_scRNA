# EarlySEC_scRNA


## Preprint Link
[Flesken-Nikitin, A.\*, Pirtz, M.G.\*, et al., *bioRxiv*, 2024](https://www.biorxiv.org/content/10.1101/2024.03.15.585274v1)


### Data Availability
`Seurat Objects`: data objects for quick analysis of cell types within the mouse uterus at different stages of SEC progression and CellChat analysis. Available through [Dryad](10.5061/dryad.59zw3r2h5) (DOI: 10.5061/dryad.59zw3r2h5). [active upon publication]  

`Raw Data`: raw and processed scRNA-seq and Visium data for the mouse uterus at normal, pre-dysplastic, and dysplastic stages of SEC. Available through [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269332) (GSE269332).

Note: All Visium analysis was completed using LoupeBrowser7. Any spots visually outside of the tissue boarders were ignored for analysis. Aligned files for analysis are available through [Dryad](10.5061/dryad.59zw3r2h5) [active upon publication]

### Subdirectories
`DataPreprocessing`: quality control and object preparation for all cells and epithelial subsets of combined datasets and normal, pre-dysplastic, and dysplastic datasets. Preprocessing is also including gene expression analysis and cluster renaming.

`Figure_Scripts`: code for recreating figures from the manuscript. 

`Data_Analysis`: code for recreating differential expression analysis used for identifying genes of interest in the manuscript.

`Reproducibility`: data alignment and R session information to recreate conclusions.  


## Contact
Matalin G. Pirtz ([mgp73@cornell.edu](mgp73@cornell.edu))  
Alexander Yu. Nikitin, PhD ([an58@cornell.edu](an58@cornell.edu))  
Benjamin D. Cosgrove, PhD ([bdc68@cornell.edu](bdc68@cornell.edu))  

