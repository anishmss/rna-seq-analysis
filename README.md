# Introduction
This repository contains scripts to reproduce the results presented in *Shrestha et al.: Comparative transcriptome profiling of heat stress response of the mangrove crab Scylla serrata across different sites.*

# Data
The raw RNA-seq reads in fastq format can be downloaded from DDBJ Sequence Read Archive under the accession number DRA010977. Here is a direct [link](https://ddbj.nig.ac.jp/DRASearch/submission?acc=DRA010977). 
There are 29 RNA-seq samples: DRX242230 -- DRX242258. Each sample is from one of 3 sites: Bataan, Cagayan, Bicol; and one of 2 conditions: control or heat-stressed. Sample information can be obtained by clicking the link to each experiment.

For annotation, download the Swiss-Prot database and the fruit fly proteome.

# Assembly-mapping-quantification and annotation pipeline

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
The main output is arranged into 3 folders: *trinity_out_dir* which contains the assembly, 
*trinity_abundance* which contains the counts required for subsequent differential expression analysis, 
and *dammit_out* which contains the output of annotation.

The transciptome assembly we obtained can be downloaded here.
The Bowtie2-RSEM quantification results we obtained can be downloaded here.


# Differential expression analysis
Differential expression analysis can be performed using the R Markdown notebook Main.Rmd.

## Software requirement
The following steps were executed on RStudio (v.1.2.5042) and R (v. 4.0.3). 
RStudio should prompt you for all other library dependecies.

## Data
Place the following files in a `data` folder:

1. **Folder of RSEM Counts**
	
    Current folder: `counts_per_sample` [💾](https://drive.google.com/file/d/1EuCITWPChEdc-SPXOluMmZNqYx81abG7/view?usp=sharing) 

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
    

### DEG workflow functions
1. Filtering
	* `filter.gene(dseq, df_filter)`
	* `filter.samples(dseq, removeSamples)`
	* `filter.printSummary(dseq)` 
2. QC
	* **Transform dseq to be homoskedastic** - `qcPrepare(dseq)`
	* **Saving QC plot outputs** - set `
qc.setSaveImg(TRUE)`, plots will be saved in _results/QC/QC (%timestamp%)_. To create directory in new timestamp, call `setTimestampNow()`.
	* **QC plots** 
		* `qc.plotGene_MeanSD()` - Diagnostic plot to check homoskedasticity of transformed data
		* `qc.plotGene_Heatmap()`
      * `qc.plotSample_ClusterHeatmap()`
      * `qc.plotSample_ClusterDendrograms()` - Same as the dendrograms of the clusters shown in the heatmap
      * `qc.plotSample_PCA()`
3. DEG
	* **Saving DEG plot outputs** - set `
deg.setSaveImg(TRUE)`, plots will be saved in _results/DE/DE (%timestamp%)_. To create directory in new timestamp, call `setTimestampNow()`.
	* **DEG plots** 
	* (Under construction)
	* (Where did testID come from?, what is its format?)
	* (How to perform signif test?)


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
