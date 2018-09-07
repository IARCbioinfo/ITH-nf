# ITH_nf #
## In development - do not use ##
![Image ITH](https://github.com/ImaneLboukili/ITH_nf/blob/master/ITH-nf.png)
## Description ##
Perform Intra-Tumor Heterogeneity (ITH) analysis.

## Dependencies ##
1. This pipeline is based on [nextflow](https://www.nextflow.io). As we have several nextflow pipelines, we have centralized the common information in the [IARC-nf](https://github.com/IARCbioinfo/IARC-nf) repository. Please read it carefully as it contains essential information for the installation, basic usage and configuration of nextflow and our pipelines.
2.Strelka2: see official installation [here](https://github.com/Illumina/strelka). )
3.Platypus: see official installation [here](http://www.well.ox.ac.uk/platypus).
4.Bcftools: see official installation [here](https://samtools.github.io/bcftools/bcftools.html).
5.Tabix: see official installation [here](http://www.htslib.org/doc/tabix.html).
6.Falcon: see official installation [here](https://omictools.com/falcon-3-tool).
7.Canopy: see official installation [here](https://github.com/yuchaojiang/Canopy).

## Input ##

**Name**         | **Description**
---------------  | -------------
--bam_folder     |  Folder containing BAM files
--output_folder  |  Path to output folder
--strelka        | Path to strelka installation dir
--correspondance | File containing the correspondance between the normal and two tumor samples and the sample
--ref            | Reference file
--regions        | Regions
--lib            | Path to libraries : falcon.output.R falcon.output.R falcon.getASCN.epsilon.R custom_canopy.plottree.R
--K              | Number of subclones to generate by Canopy
--tabix          |	Path to tabix installation dir
--platypus		   |	Path to platypus installation dir
--Rcodes 			   |	Path to folder containing R codes




### Flags ###

Flags are special parameters without value.

**Name**      | **Description**
------------- | -------------
--help        | Display help



## Usage ##

`nextflow run iarcbioinfo/pipeline.nf   --bam_folder path/to/bams/ --correspondance path/to/correpondance/csv/  --output_folder /path/to/output --strelka /path/to/trelka --bcftools  /path/to/bcftools --tabix /path/to/tabix --platypus /path/to/platypus --Rcodes --lib /path/to/R/libraries --K integer --ref /path/to/ref --regions path/to/regions`


## Output ##

**Name**                             | **Description**
-------------------------------------| -------------
 VCFs (and corresponding .tbi)       | Results of the germline and somatic callings, the somatic and germline coverage
 PDFs                                | Reports for Falcon and Canopy runs
 rda                                 | Falcon and Canopy coordinates
 txt                                 | Standard error for Falcon
 SVGs                                | Canopy trees
