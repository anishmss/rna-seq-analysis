# Introduction
This repository contains scripts to reproduce the results presented in 
*Shrestha et al.: Comparative transcriptome profiling of heat stress response of the mangrove crab Scylla serrata
across sites of varying climate profiles*

# Data
The raw RNA-seq reads in fastq format can be downloaded from DDBJ Sequence Read Archive under the accession number DRA010977. Here is a direct [link](https://ddbj.nig.ac.jp/DRASearch/submission?acc=DRA010977). 
There are 29 RNA-seq samples: DRX242230 -- DRX242258. Each sample is from one of 3 sites: Bataan, Cagayan, Bicol; and one of 2 conditions: control or heat-stressed. Sample information can be obtained by clicking the link to each experiment.


# Assembly-mapping-quantification

## Software requirement
This pipeline requires requires Conda and [Snakemake](https://snakemake.readthedocs.io/en/stable/).
Other dependencies are handled by Snakemake automatically. 
Snakemake recommends installation via Conda:
```
$ conda install -c conda-forge mamba
$ mamba create -c conda-forge -c bioconda -n snakemake snakemake
```
This creates an isolated enviroment containing the latest Snakemake. To activate it:
```
$ conda activate snakemake
```
To test snakemake installation 
```
$ snakemake --help
```

## Running the pipeline
Download or clone this repository. In the *config.yaml* file, set (1) the absolute paths to the RNA-seq reads, 
(2) an absolute path to the directory where you want the outputs to be saved, and (3) the absolute paths to the databases for annotation.
Then at the top-level directory, run the following :
```
snakemake --configfile config.yaml --use-conda --cores all 
```

## Output
The main output is arranged into several folders, of which: *trinity_out_dir* which contains the assembly *Trinity.fasta*, 
*trinity_abundance* which contains the counts required for subsequent differential expression analysis.

The transciptome assembly we obtained can be downloaded [here](https://drive.google.com/file/d/1rw2rbbGzz2etlSLihpLCq1hU29PJBXGX/view?usp=sharing).
The quantification results we obtained can be downloaded [here](https://drive.google.com/file/d/1rw2rbbGzz2etlSLihpLCq1hU29PJBXGX/view?usp=sharing).


# Assembly annotation

## Data
For annotation, download the Swiss-Prot database and the fruit fly proteome UP000000803 .

## Software requirement
We used Dammit (v 1.2), which can be installed using bioconda:
```
$ mamba create -c conda-forge -c bioconda -n dammit dammit=1.2
```
Next activate the environment:
```
$conda activate dammit
```

## Running Dammit

With the Swiss-Prot dataset and fruit fly proteome, launch Dammit:
```
dammit databases --install --quick --busco-group metazoa 
dammit annotate trinity_out/Trinity.fasta --quick --user-databases Swiss-Prot.fa fruifly.fa -e 1e-10 --output-dir dammit_out
```

## Ouput
The output contains many files, which we have summarized as a data-frame required by downstream analysis [here](https://drive.google.com/file/d/1k5S_lzy4_NgnF7sXa5dPKTInMbE4fz9L/view?usp=sharing).


# Differential expression analysis
Differential expression analysis can be performed using the R Markdown notebooks \`1 - Explore and QC.Rmd\` and \`2 - DE Analysis.Rmd\`.

## Software requirement
The following steps were executed on RStudio (v.1.2.5042) and R (v. 4.0.3). 
RStudio should prompt you for all other library dependecies.

## Data
Place the following files in a `data` folder:

1. **Folder of RSEM Counts**
	
    Current folder: `counts_per_sample` [💾](https://drive.google.com/file/d/1rw2rbbGzz2etlSLihpLCq1hU29PJBXGX/view?usp=sharing) 

	Generated from RSEM after TRINITY assembly. Contains ".genes.results" file for each sample.  
	
	Filename format : "BAT_C_1.genes.results" - This refers to Bataan Control group sample 1


2. **Annotation Dataframe**
	
    Current file: `annotation.20200701.csv` [💾](https://drive.google.com/file/d/1k5S_lzy4_NgnF7sXa5dPKTInMbE4fz9L/view?usp=sharing) 
  
  	Dataframe should have columns "gene_id" and also for labeling and filtering.
    
	| gene_id | swiss_prot_acc | swiss_prot_stat_sig |
        |-|-|-|
        | TRINITY_DN10022_c0_g1 | Q9UQ80 | 1.2e-176 |
        | TRINITY_DN10031_c0_g1 | Q8VDW4 | 1.0e-99 |
        
3. **Descendant Genes Dataframe**

	Current file: `geneID_isDescendant.csv` [💾](https://docs.google.com/spreadsheets/d/1CEU09PBeoSlpwWcEvkLSivj9E2YCYtWVbQxASdg0QWM/edit?usp=sharing) 
    
    Generated from "supplementary notebooks/GO_stressResponse.Rmd". View [below](###GO_stressResponse.Rmd) for more details.
    

### DEG workflow

Run notebook \`**1 - Explore and QC.Rmd**\` for exploring data and perform quality checking

-   Chunk 1 : Prepare DESeq dataset using txImport

    -   Reads count data produced by Trinity as a DESeq2 dataset

-   Chunk 2 : Filter Genes and Samples

    -   Removes outlier samples

    -   Removes genes with low annotation significance and alignment

-   Chunk 3 : (optional) Save parameters of this experiment

-   Chunk 4 : Quality Check

    -   Graphs sample cluster with heatmap

    -   Graphs sample cluster dendrogram

    -   Graphs sample PCA plot

    -   Graphs gene heatmap

    -   Graphs Mean SD plot (for checking if variance is dependent on the mean)

-   Chunk 5 : Save filtered DESEq2 object to be loaded in \`2 - DE Analysis.Rmd\`

Run notebook \`**2 - DE Analysis.Rmd**\` for performing significance testing and graphing

-   Chunk 1 : Load the filtered DESeq2 dataset saved from \`1 - Explore and QC.Rmd\`

-   Chunk 2 : Perform significance testing

-   Chunk 3 : Export list of Differentially Expressed (DE) genes

-   Chunk 4 : Graph results from Main effects (the 3 sites: "CAG", "BAT", "BIC") testing

    -   Graphs expression of top most significant DE genes for each site

    -   Graphs venn diagram of the DE genes across the 3 site

    -   Graphs bar graph of DE gene counts for each site

    -   Graphs MA plot for each site

-   Chunk 5 : Graph results from Interaction effects testing

    -   Graphs venn diagram of intersecting DDE genes

    -   Graphs a count table of the interaction category for each site

    -   Graphs gene heatmaps of CAG vs BAT and CAG vs BIC

    -   Graphs interaction plots for each site

-   Chunk 6 : Combine plots for publication (Figure 7)

-   Chunk 7 : Combine plots for publication (Figure 8)

-   Chunk 8 : Plot figure of Top 20 read counts

# Supplementary Notebooks
### GO_stressResponse.Rmd
**Requirements :**

1. Annotation Dataframe with columns "gene_id" and its Uniprot Accession ID
    
	| gene_id | swiss_prot_acc |
        |-|-|
        | TRINITY_DN10022_c0_g1 | Q9UQ80 |
        | TRINITY_DN10031_c0_g1 | Q8VDW4 |
2. At least 1 hour stable internet

**Workflow :**

1. Get list of Accession IDs from specified column (currently "swiss_prot_acc")

2. Using Uniprot API, for each Accession ID, retrieve all GO terms related to Biological Process. 

3. Using QuickGO API, for each GO term of an Accession ID, retrieve its ancestors and look for Stress Response To Heat ("GO:0006950"). If so, isDescendant is True.

4. Map Gene IDs to isDecendant.
    
**Output :** This notebook will produce 3 CSV files :

1. **geneID_isDescendant.csv**
    
	| gene_id | isDescendant |
    |------------|--------------|
    | TRINITY_DN100016_c0_g1     | TRUE         |

2. uniprotID_isDescendant.csv
	* with columns : uniprot_id | isDescendant |

3. uniprotID_to_goID.csv
	* with columns : uniprot_id | GO_id      | isDescendant |
    
**Notes :**

1. Current DGE Analysis as of 08/24/2020 uses swiss_prot_acc to determine if gene is descendant. But the current filter criteria uses fly annotation. I tried using fly_acc for isDescendant gene and cross-referenced to the current isDescendant file, around 17% isDescendant did not match.

---
### climate_data_analysis.Rmd
**Select Input file :**

1. ERSSTv5 NetCDF files from [NOAA website](https://www1.ncdc.noaa.gov/pub/data/cmb/ersst/v5/netcdf/) (.zip here [💾](https://drive.google.com/file/d/1O9BcjtofW6ck03TohvgMoJp321fRa2vn/view?usp=sharing))
     * Use this if site coordinates are not yet final and you are still exploring.
	
2. `Climate Data.csv` [💾](https://docs.google.com/spreadsheets/d/11sglmWb_LP2KpaXw847h1F3LeH6vzHY7ha79iQW9EZ4/edit?usp=sharing)
    * Use this if site coordinates and date range is already final.
	* This file is generated in this notebook. Contains the SST and SSTA for all pre-specified sites within a date range. 
	     
		 | site | date      | SST | SSTA |
      |------------|------------|--------------|--------------|
      | CAG     | 185401 | 26.54        |1.19        |

**Workflow (Using NetCDF files) :**

1. For each site, define its 2x2 lon and lat grid.

2. Read all NetCDF files and extract SST and SSTAs from all site coordinates.

3. Generates the dataframe in Climate Data.csv, which can then be used for plotting.


**Output :** This notebook will produce 2 graphs :

1. Distribution of SSTs as a Violin plot

2. Time series of SSTA w/ yearly moving average


---
### pca_metadata.Rmd
**Requirements :**

1. Metadata CSV file

	* Current file : `NGS_morphodata_summary.csv`[💾](https://drive.google.com/file/d/1hKav0bmqGAm-sTRshKRQLrde60oxotpq/view?usp=sharing)

	* Columns should contain "samplenames" and other factors of that sample

	  | samplename | sex | weight_init | weight_final | carap_init | carap_final |
		|-|-|-|-|-|-|
		| BIC_E_1 | M | 100.8 | 116.4 | 96.58 | 96.66 |
		| BIC_E_2 | I | 88.8 | 97.4 | 91.92 | 92.21 |
    
    
**Workflow :**

1. PCA of samples can then be color coded by their characteristics.
