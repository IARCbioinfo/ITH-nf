#!/usr/bin/env nextflow
// Copyright (C) 2017 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


params.input_BAM = null
params.input_ref = null
params.input_regions = null

log.info ""
log.info "----------------------------------------------------------------"
log.info " fastqc-0.11.3/  : Quality control with FastQC  and MultiQC     "
log.info "----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "-------------------QC-------------------------------"
    log.info "" 
    log.info "nextflow run iarcbioinfo/FastQC.nf   -j java --input_folder path/to/fasta/ --fastqc path/to/fastqc/ -o --output_folder /path/to/output"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--samtools              PATH               samtools installation dir"
    log.info "--sambamba              PATH               sambamba installation dir"
    log.info "--k8              PATH              k8 installation dir from bwakit-0.7.15" 
    log.info "--bwa-postalt              PATH              bwa-postalt.js installation dir from bwakit-0.7.15"
    log.info "--hs38DH.fa.alt              PATH            hs38DH.fa.alt installation dir from bwakit-0.7.15"
    log.info "--strelka              PATH                strelka installation dir"
    log.info "--platypus              PATH               platypus installation dir"
    log.info "--R              PATH               R installation dir"
    log.info "--vt              PATH               vt installation dir"
    log.info "--annovar              PATH               annovar installation dir"
    log.info "--tabix             PATH               tabix installation dir"
    log.info "--input_folder         FOLDER               Folder containing bam files and R scripts"
    log.info ""
    log.info "Optional arguments:"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=2)"
    log.info "--output_folder        PATH                 Output directory for html and zip files (default=fastqc_ouptut)"
    log.info "--config               FILE                 Use custom configuration file"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 1
} 


// BAM correspondance file
bams = file(params.input_BAM)

bams_repo = Channel.fromPath(bams)
  .splitCsv(header: ['patient_id', 'normal', 'tumor1', 'tumor2'], skip: 1 )
  .map {  row ->
  [row.patient_id ,
    [file(row.normal), file(row.normal.replace(".bam",".bam.bai"))],
    [file(row.tumor1), file(row.tumor1.replace(".bam",".bam.bai"))],
    [file(row.tumor2), file(row.tumor2.replace(".bam",".bam.bai"))]
  ] }

// Print repo
//repo = bams_repo.map {path -> path }
//repo.subscribe {println it}

/* Data access
patient_id : path[0]
normal array : path[1]
normal bam : path[1][0]
normal bai : path[1][1]
tumor1 array : path[2]
tumor2 array : path[3]
*/

bams_repo.into { bams_repo ; bam ; tumor1 ; tumor2 }
normal = bam.map { path -> path[1][0] }
tumor1 = tumor1.map { path -> path[2][0] }
tumor2 = tumor2.map { path -> path[3][0] }
all_bams = normal.mix(tumor1, tumor2)
patient_bams = bam.map { path -> [ path[1][0], path[2][0], path[3][0] ] }

/** PROCESSES **/


// Etape de post-alignement à finir
process post_alignment {
  input:
  file i from all_bams

  output:
  file("${i}.HEAD") into BAM_post_al

  shell:
  '''
samtools view -h ${BAM} | k8 bwa-postalt.js hs38DH.fa.alt | \
sambamba view -S -f bam -l 0 /dev/stdin | \
sambamba sort -t 8 -m 6G --tmpdir=!{BAM_post_al}/temp -o !{BAM_post_al}/!{i}_pa.bam /dev/stdin
  //head -n1 !{i} > !{i}.HEAD
  '''
}


process germline_calling {
  input:
  file i from normal
  file input_ref
  file input_regions

  output:
  file("${i}.germline.vcf") into VCF_germline

  shell:
  '''
  Platypus.py callVariants --bamFiles=!{i} --output=!{i}.germline.vcf --refFile=!{input_ref} --regions=!{input_regions} --badReadsThreshold=0 --qdThreshold=0 --rmsmqThreshold=20 --hapScoreThreshold=10 --scThreshold=0.99
  '''
}

process germline_calling_pass {
  input:
  file i from VCF_germline


  output:
  file("${i}.germline.vcf") into VCF_germline_pass

  shell:
  '''
cat !{i} | grep "PASS" > !{VCF_germline_pass}/!{i}.germline.pass.vcf
'''
}


process germline_AF {
  input :
  file i from VCF_germline_pass
  
  output :
  
  file ("${i}.svg") into germline_AF
  shell :
  '''
  Rscript --vanilla germline_fractions_alleliques.R VCF_germline_pass germline_AF 
  '''
}

process somatic_calling {
  input:
  file i from patient_bams
  file input_ref
  file input_regions

  output:
  file("${i}.somatic.vcf") into VCF_somatic

  shell:
  '''
  strelka/configureStrelkaSomaticWorkflow.py --tumorBam='$TUMORBAM' --normalBam=!{i}[1][0] --referenceFasta=!{input_ref} --callRegions=!{input_regions} --callMemMb=1024 --runDir=STRELKA/!{i}
  STRELKA/!{i}/runWorkflow.py -m local -j 28

  '''
}


process somatic_calling_pass {
  
  input:
  file i from VCF_somatic


  output:
  file("${i}.germline.vcf") into VCF_somatic_pass

  shell:
  '''
	cat !{i} | grep "PASS" > !{VCF_somatic_pass}/!{i}.somatic.pass.vcf
'''
}


process somatic_AF {
  input :
  file i from patient_bams

  output :
   file ("${i}.svg") into somatic_AF
  
  shell :
  '''
  Rscript --vanilla somatic_fractions_alleliques.R !{VCF_somatic_pass}/!{i}[2][0] !{VCF_somatic_pass}/!{i}[3][0] somatic_AF !{i}  !{i}[2][0] !{i}[3][0]
  '''
}

process somatic_germline_intersect {
  input :
  file i from patient_bams

  output :
   file ("${i}.svg") into somatic_germline_intersect
  
  shell :
  '''
  Rscript --vanilla somatic_germline_fractions_alleliques.R !{i}[1][0] somatic_germline_intersect !{i} !{i}[2][0] !{i}[3][0]
  '''
}


process venn_diagramm {
  input :
  file i from patient_bams

  output :
   file ("${i}.svg") into venn_diagramm

  shell :
  '''
  Rscript --vanilla somatic_overlap_venn.R !{VCF_somatic_pass}/!{i}[2][0] !{VCF_somatic_pass}/!{i}[2][0] venn_diagramm !{i} !{i}[2][0] !{i}[3][0]
  '''
}


process germline_tumor_coverage {
  input:
  file i from patient_bams
  file input_ref
  file input_regions
  file j from VCF_germline_pass
  
  output:
  file("${i}.coverage.germline.vcf") into coverage_germline
}
shell :
if(i==j){

Platypus.py callVariants --bamFiles=!{i}[1][0],!{i}[2][0],!{i}[3][0] --refFile=!{input_ref} --regions=!{input_regions} --nCPU=12 --output=!{coverage_germline}/!{i}.coverage.germline.vcf --source=!{VCF_germline_pass}/!{i}.germline.pass.vcf.gz --minPosterior=0 --getVariantsFromBAMs=0
}

process somatic_tumor_coverage {

   input:
  file i from patient_bams
  file input_ref
  file input_regions
  file j from VCF_somatic_pass
  
  output:
  file("${i}.coverage.germline.vcf") into coverage_somatic
  
  shell :
  '''
  if(i==j){
  Platypus.py callVariants --bamFiles=!{i}[1][0] --refFile=!{input_ref} --regions=!{input_regions} --nCPU=12 --output=!{coverage_somatic}/!{i}.coverage.somatic.vcf --source=!{VCF_somatic_pass}/!{i}.germline.pass.vcf --minPosterior=0 --getVariantsFromBAMs=0
  }
  '''
}




process split_VCF_somatic {

   input:
  file i from coverage_somatic
  
  output:
  file("${i}.vcf") into split_vcf_germline
  
  shell :
  

declare -a arr=("chr1" "chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY")

for( j in "${arr[@]}"){
do
   tabix -h $1 $i > $1.$j.vcf
done
}





process copy_numbers {
  // run falcon
}

process copy_numbers_epsilon {
  // run falcon_epsilon (to compute errors)
}

process compile_copy_numbers {
  // compile falcon results in a TSV file
}

process MCMC {
  // canopy pre-clustering & MCMC sampling
}

process tree {
  // canopy tree
}

process filter {
  // filter informative events
}