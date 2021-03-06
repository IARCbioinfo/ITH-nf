#!/usr/bin/env nextflow
// Copyright (C) 2018 IARC/WHO

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

params.help = null
params.fasta_ref = null
assert (params.fasta_ref != null) : "please specify --fasta_ref option to provide a reference genome"
params.regions = null
assert (params.regions != null) : "please specify --regions option to provide a BED file"
params.correspondance = null
assert (params.correspondance != null) : "please specify --correspondance option to provide a tumor/normal correspondance file"
params.bam_folder = null
assert (params.bam_folder != null) : "please specify --bam_folder option to provide path of BAM folder"
params.cpu  = 1
params.mem  = 4
params.strelka2  = null
assert (params.strelka2 != null) : "please specify --strelka2 option to provide the installation folder of strelka2"
params.bcftools   = "bcftools"
params.tabix  = "tabix"

log.info ""
log.info "----------------------------------------------------------------"
log.info "              Intra tumor heterogeneity pipeline                "
log.info "----------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "------------ IN DEVELOPMENT - DO NOT USE ---------------"
    log.info "--------------------------------------------------------"
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "----- Intra tumor heterogeneity nextflow pipeline ------"
    log.info ""
    log.info "nextflow run iarcbioinfo/pipeline.nf   --bam_folder path/to/bams/ --correspondance path/to/correpondance/csv/  --output_folder /path/to/output --strelka2 /path/to/trelka --bcftools  /path/to/bcftools --tabix /path/to/tabix --platypus /path/to/platypus --K integer --fasta_ref /path/to/ref --regions path/to/regions"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--strelka2            PATH          Path to strelka2 installation dir "
    log.info "--R                   PATH          R installation dir"
    log.info "--bam_folder          FOLDER        Folder containing bam files "
    log.info "--correspondance		  FILE				  File containing the correspondance between the normal and two tumor samples and the sample id"
    log.info "--fasta_ref           FILE          Reference file"
    log.info "--regions             FILE          Regions"
    log.info "--lib                 PATH  				Path to libraries : falcon.output.R falcon.output.R falcon.getASCN.epsilon.R custom_canopy.plottree.R"
    log.info "--K                   INTEGER				Number of subclones to generate by Canopy"
    log.info "--tabix 				      PATH  				Path to tabix installation dir"
    log.info "--platypus			      PATH		  		Path to platypus installation dir"
    log.info ""
    log.info "Optional arguments:"
    log.info "--cpu                 INTEGER       Number of cpu to use (default=28)"
    log.info "--output_folder       PATH          Output directory for html and zip files (default=fastqc_ouptut)"
    log.info "--config              FILE          Use custom configuration file"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 1
}


fasta_ref = file(params.fasta_ref)
fasta_ref_fai = file( params.fasta_ref+'.fai' )
regions = 	file(params.regions)
regions_tbi = file( params.regions+'.tbi' )
correspondance = file(params.correspondance)

bams_TTN= Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor1+'.bai'),file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.tumor2+'.bai'),file(params.bam_folder + "/" + row.normal),file(params.bam_folder + "/" + row.normal+'.bai') ]}

(bams_TT,bams_TT2) =  Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor1+'.bai'),file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.tumor2+'.bai') ]}.into(2)

bams_N = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.normal),file(params.bam_folder + "/" + row.normal+'.bai') ]}

bams_T1N = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor1+'.bai'),file(params.bam_folder + "/" + row.normal),file(params.bam_folder + "/" + row.normal+'.bai') ]}

bams_T2N = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.tumor2+'.bai'),file(params.bam_folder + "/" + row.normal),file(params.bam_folder + "/" + row.normal+'.bai') ]}

IDs_TTN = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,row.tumor1, row.tumor2, row.normal ]}

(IDs_TT,IDs_TT2) = Channel.fromPath(correspondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,row.tumor1, row.tumor2 ]}.into(2)

chromosomes = Channel.from(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,'X','Y')

strelka2_germline = params.strelka2 + '/bin/configureStrelkaGermlineWorkflow.py'


process germline_calling {
  tag { ID }

  publishDir params.output_folder, mode: 'copy'

  input:
  set val(ID), file(tumor1), file(tumor1_bai), file(tumor2), file(tumor2_bai), file(normal), file(normal_bai) from bams_TTN
  file fasta_ref
  file fasta_ref_fai
  file regions
  file regions_tbi

  output:
  set val("${ID}"),file("${ID}.germline.variants.vcf.gz") into VCF_germline
  set val("${ID}"), file("${ID}.germline.variants.vcf.gz.tbi") into TBI_Germline

  shell:
  '''
  runDir="results"
  !{strelka2_germline} --bam !{normal} --bam !{tumor1} --bam !{tumor2} --referenceFasta !{fasta_ref} --callRegions !{regions} --runDir .
  ./runWorkflow.py -m local -j !{params.cpu}
  sm=`samtools view -H !{normal} | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`
  idx=`bcftools query -l results/variants/variants.vcf.gz | nl -v 0 | grep $sm | cut -f1`
  echo "bcftools view -m2 -M2 -v snps -f PASS -i 'GT[$idx]=\"1/0\" | GT[$idx]=\"0/1\" | GT[$idx]=\"1/1\"' results/variants/variants.vcf.gz" | bash - | bgzip -c > !{ID}.germline.variants.vcf.gz
  tabix -p vcf !{ID}.germline.variants.vcf.gz
  '''
}

/*

strelka2_somatic= params.strelka2 + '/bin/configureStrelkaSomaticWorkflow.py'

process somatic_calling_T1 {

 publishDir params.output_folder, mode: 'copy'

  input:
set  val( ID) ,file (tumor1),file(tumor1_bai), file(normal), file(normal_bai) from bams_T1N
  file fasta_ref
  file regions

  output:
   set val (ID), file ('*.indels.vcf.gz') into VCF_somatic1_indels
   set val (ID), file ('*.snvs.vcf.gz') into (VCF_somatic1_snvs,VCF_somatic1_snvs2)
   set val (ID), file ('*.tbi') into TBI_somatic1
  shell:
  '''

 !{strelka2_somatic} --tumorBam=!{tumor1} --normalBam=!{normal} --referenceFasta=!{fasta_ref} --callRegions=!{params.regions} --callMemMb=1024   --runDir strelkaSomatic1/!{ID}
 cd strelkaSomatic1/!{ID}
     ./runWorkflow.py -m local -j 28

  cd results/variants

     !{params.bcftools} view -i'FILTER="PASS"' somatic.indels.vcf.gz > somatic.indels.vcf.gz
     mv somatic.indels.vcf.gz !{tumor1.baseName}.somatic.indels.vcf.gz

     !{params.bcftools} view -i'FILTER="PASS"' somatic.snvs.vcf.gz >  somatic.snvs.vcf.gz
     mv somatic.snvs.vcf.gz !{tumor1.baseName}.somatic.snvs.vcf.gz

     mv somatic.indels.vcf.gz.tbi !{tumor1.baseName}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{tumor1.baseName}.somatic.snvs.vcf.gz.tbi


  '''
}

process somatic_calling_T2 {

publishDir params.output_folder, mode: 'copy'

  input:
set  val( ID) ,file (tumor2), file(tumor2_bai), file(normal),file(normal_bai) from bams_T2N
  file fasta_ref
  file regions

  output:
   set val (ID), file ('*.indels.vcf.gz') into VCF_somatic2_indels
   set val (ID), file ('*.snvs.vcf.gz') into (VCF_somatic2_snvs,VCF_somatic2_snvs2)
   set val (ID), file ('*.tbi') into TBI_somatic2

  shell:
  '''

 !{strelka2_somatic} --tumorBam=!{tumor2} --normalBam=!{normal} --referenceFasta=!{fasta_ref} --callRegions=!{params.regions} --callMemMb=1024   --runDir strelkaSomatic2/!{ID}
 cd strelkaSomatic2/!{ID}
     ./runWorkflow.py -m local -j 28
cd results/variants

     !{params.bcftools} view -i'FILTER="PASS"' somatic.indels.vcf.gz > somatic.indels.vcf.gz
     mv somatic.indels.vcf.gz !{tumor2.baseName}.somatic.indels.vcf.gz

     !{params.bcftools} view -i'FILTER="PASS"' somatic.snvs.vcf.gz >  somatic.snvs.vcf.gz
     mv somatic.snvs.vcf.gz !{tumor2.baseName}.somatic.snvs.vcf.gz

     mv somatic.indels.vcf.gz.tbi !{tumor2.baseName}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{tumor2.baseName}.somatic.snvs.vcf.gz.tbi


  '''
}

VCF_somatic = VCF_somatic1_snvs.join(VCF_somatic2_snvs)
input_somaticCoverage = bams_TT.join(VCF_somatic)

VCF_somatic2 = VCF_somatic1_snvs2.join(VCF_somatic2_snvs2)


process somatic_tumor_coverage {

publishDir params.output_folder, mode: 'copy'

   input:
  set val(ID),file (bamtumor1),file(bamtumor2), file(somaticVCF1),file(somaticVCF2) from input_somaticCoverage
  file fasta_ref
  file regions


  output:
set val(ID),file("${ID}_covargeSomatic_T1.vcf.gz"),file("${ID}_covargeSomatic_T2.vcf.gz") into VCF_coverage_somatic
 set val(ID),file("${ID}_covargeSomatic_T1.vcf.gz.tbi"),file("${ID}_covargeSomatic_T2.vcf.gz.tbi") into TBI_coverage_somatic
 set val(ID),file("variants.vcf.gz"),file("variants.vcf.gz.tbi") into variants_coverage_somatic

  shell :
  '''
 !{strelka2_germline} --bam=!{bamtumor1} --bam !{bamtumor2} --forcedGT !{somaticVCF1} --forcedGT !{somaticVCF2}  --referenceFasta=!{fasta_ref}   --callRegions=!{params.regions} --runDir strelkaCoverageSomatic/!{ID}
 cd strelkaCoverageSomatic/!{ID}
     ./runWorkflow.py -m local -j 28

     cd results/variants

     mv genome.S1.vcf.gz !{ID}_covargeSomatic_T1.vcf.gz
     mv genome.S1.vcf.gz.tbi !{ID}_covargeSomatic_T1.vcf.gz.tbi
     mv genome.S2.vcf.gz !{ID}_covargeSomatic_T2.vcf.gz
     mv genome.S2.vcf.gz.tbi !{ID}_covargeSomatic_T2.vcf.gz.tbi

  '''
}


process split_into_chr {

	input :
	set val(ID),file(vcf) from coverage_germline
	set val(chr) from  chromosomes

	output :
	set val(ID), val(chr) file ("*.vcf.gz") into VCF_by_chr

	shell :

 """
   !{params.tabix} -p !{vcf}
   !{params.tabix} -h !{filevcf} chr${chr} > germline_chr${chr}.vcf
"""
}

input_Falcon = VCF_by_chr.join(IDs_TTN)

process Falcon {

publishDir params.output_folder, mode: 'copy'

 	input :
 	set val(ID), val(chr), file(vcf_splitted), val(T1_ID), val(T2_ID), val(N_ID) from input_Falcon

 	output :

 	set val(ID),  file('*.pdf') into Falcon_PDF_report
	set val(ID),  file('*.txt') into (Falcon_CNVs_txt,Falcon_CNVs_txt2)
	set val(ID), file('*.rda') into Falcon_CNVs_rda


 	shell :
 	 '''
 	 Rscript --vanilla !{baseDir}/Falcon.R !{vcf_splitted} !{ID} !{N_ID} !{T1_ID} !{T2_ID} !{chr} !{params.output_folder} !{params.baseDir}/libs}/falcon.output.R !{baseDir}/libs/falcon.qc.R
   '''
}


 input_Falcon_eps = Falcon_CNVs_txt.join(Falcon_CNVs_rda)

process Falcon_stderr {

publishDir params.output_folder, mode: 'copy'

	input :
	set val(ID), file(txt), file(coord1),file(coord2)  from input_Falcon_eps


	output :

	set val(ID), file('*.txt') into Falcon_stderr

	shell :
	 '''
Rscript --vanilla !{baseDir}/Falcon_epsilon.R   !{txt} !{coord1} !{coord2}   !{params.output_folder} !{baseDir}/libs/falcon.output.R !{baseDir}/libs/falcon.getASCN.epsilon.R
 '''
}


VCF_somatic_Canopy = VCF_somatic2.join(VCF_coverage_somatic)
input_Canopy = Falcon_CNVs_txt2.join(VCF_somatic_Canopy.join(IDs_TT2.join(Falcon_stderr)))


process Canopy {

publishDir params.output_folder, mode: 'copy'

	input :

	set val(ID), file(falcontxt), file(somatic1),file(somatic2), file(coveragesomatic1), file(coveragesomatic2) ,val(T1_ID), val(T2_ID),file(txt1),file(txt2) from input_Canopy


	output :

	set val(ID), file('*.bic') into Canopy
	set val(ID),file('*.svg'),file('*.pdf'),file('*.rda') into Canopy_reports

	shell :
  '''
  Rscript --vanilla !{baseDir}/Canopy.R !{falcontxt} !{ID} !{T1_ID}  !{T2_ID} !{somatic1} !{somatic2} !{coveragesomatic1} !{coveragesomatic2} !{params.output_folder} !{params.K} !{baseDir}/libs/custom_canopy.sample.cluster.R !{txt1} !{txt2}
  '''
}


process Canopy_tree {

publishDir params.output_folder, mode: 'copy'

	input :

	set val(ID), file(bic) from Canopy
	output :

	set val(ID), file('*.SVG'),file('*.pdf'),file('*.txt') into Canopy_trees


	shell :
	 '''
	Rscript --vanilla !{baseDir}/Canopy_tree.R !{ID}   !{bic}  !{baseDir}/libs/custom_canopy.plottree.R
 '''
}

*/
