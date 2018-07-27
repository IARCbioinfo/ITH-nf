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


params.help = null
params.ref = null
params.regions = null
params.cpu            		= "28"
params.mem           		 = "20"
params.strelka  = null
params.bcftools   = null
params.R   = null
params.tabix  = null




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
    log.info "                     USAGE                              "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "-------------------ITH-------------------------------"
    log.info "" 
    log.info "nextflow run iarcbioinfo/pipeline.nf   --bam_folder path/to/bams/ --correspondance path/to/correpondance/csv/  --output_folder /path/to/output --strelka --bcftools --R  --tabix --ref --regions path/to/"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--strelka             PATH                Path to strelka installation dir "
    log.info "--bcftools              PATH              bcftools installation dir"
    log.info "--R              PATH               R installation dir"
    log.info "--tabix             PATH               tabix installation dir"
    log.info "--falcon_qc             PATH            falcon.qc.r dir"
    log.info "--bam_folder         FOLDER               Folder containing bam files "
    log.info "--correspondance		FILE				File containing the correspondance between the normal and two tumor samples and the sample id"
    log.info "--ref                 FILE                Reference file"
    log.info "--regions             FILE                 Regions "
    log.info ""
    log.info "Optional arguments:"
    log.info "--cpu                  INTEGER              Number of cpu to use (default=28)"
    log.info "--output_folder        PATH                 Output directory for html and zip files (default=fastqc_ouptut)"
    log.info "--config               FILE                 Use custom configuration file"
    log.info ""
    log.info "Flags:"
    log.info "--help                                      Display this message"
    log.info ""
    exit 1
} 


ref = file(params.ref)
regions = 	file(params.regions)
correpondance = file(params.correpondance)



bams_TTN= Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.normal) ]}


(bams_TT,bams_TT2) =  Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor2) ]}.into(2)
		
bams_N = Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.normal) ]}
		
bams_T1N = Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.normal) ]}

bams_T2N = Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.normal) ]}
	




strelka_germline= params.strelka + '/bin/configureStrelkaGermlineWorkflow.py' 		
		
process germline_calling {
  input:
  set val(ID) ,file (normal) from bams_N
  file ref
  file regions

  output:
  set val(ID_tag),file("${normal.baseName}.vcf.gz") into VCF_germline
  set val(ID_tag),file("${normal.baseName}.variants.vcf.gz") into VCF_germlineVariants
  set val(ID_tag), file("${normal.baseName}.vcf.gz.tbi"), file("${normal.baseName}.variants.vcf.gz.tbi") into TBI_Germline

  shell:
  '''
  ID_tag=!{ID}  
  
  runDir="results/variants/"
  !{strelka_germline} --bam !{normal} --referenceFasta !{params.ref}   --callRegions !{params.regions} --runDir strelkaAnalysis
  cd strelkaAnalysis
  ./runWorkflow.py -m local -j !{params.cpu} 
  
  mv  genome.vcf.gz !{normal.baseName}.vcf.gz
  mv  variants.vcf.gz !{normal.baseName}.variants.vcf.gz
  mv  genome.vcf.gz.tbi !{normal.baseName}.vcf.gz.tbi
  mv  variants.vcf.gz.tbi !{normal.baseName}.variants.vcf.gz.tbi

{params.bcftools} norm  -f !{normal.baseName}.vcf.gz 
{params.bcftools} view -i'FILTER="PASS"' !{normal.baseName}.vcf.gz  > !{normal.baseName}.vcf.gz
!{params.tabix} -p vcf !{normal.baseName}.vcf.gz
  '''
}


input_germlineCoverage = bams_TTN.join(VCF_germline)
process germline_tumor_coverage {
  input:
  set val(ID),file (bamtumor1),file(bamtumor2),file(bamnormal),file(germlineVCF) from input_germlineCoverage
  file ref
  file regions

  
  output:
  set val(ID),file("${ID}.vcf") into coverage_germline

shell :
'''
 !{strelka_germline} --bam !{bamtumor1} --bam !{bamtumor2} --bam !{bamnormal} --forcedGT !{germlineVCF}  --referenceFasta=!{params.ref}   --callRegions=!{params.regions} --runDir strelkaAnalysis
 
'''
}



strelka_somatic= params.strelka + '/bin/configureStrelkaSomaticWorkflow.py' 

process somatic_calling_T1 {
  input:
  set val(ID) ,file (tumor1), file(normal) from bams_T1N
  file ref
  file regions

  output:
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.indels.vcf.gz' into VCF_somatic1_indels
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.snvs.vcf.gz' into VCF_somatic1_snvs
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.tbi' into TBI_somatic1
  shell:
  '''
  ID_tag=!{ID}
 !{strelka_somatic} --tumorBam=!{tumor1} --normalBam=!{normal} --referenceFasta=!{params.ref} --callRegions=!{params.regions} --callMemMb=1024   --runDir strelkaAnalysis
 cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd results/variants
     mv somatic.indels.vcf.gz !{tumor1.baseName}.somatic.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{tumor1.baseName}.somatic.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{tumor1.baseName}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{tumor1.baseName}.somatic.snvs.vcf.gz.tbi
     
     {params.bcftools} norm  -f !{tumor1.baseName}.somatic.indels.vcf.gz
    !{params.bcftools} view -i'FILTER="PASS"' !{tumor1.baseName}.somatic.indels.vcf.gz  > !{tumor1.baseName}.somatic.indels.vcf.gz
    !{params.tabix} -p vcf !{tumor1.baseName}.somatic.indels.vcf.gz
    
    {params.bcftools} norm  -f !{tumor1.baseName}.somatic.snvs.vcf.gz
    !{params.bcftools} view -i'FILTER="PASS"' !{tumor1.baseName}.somatic.snvs.vcf.gz  >  !{tumor1.baseName}.somatic.snvs.vcf.gz
    !{params.tabix} -p vcf !{tumor1.baseName}.somatic.snvs.vcf.gz
  '''
}


process somatic_calling_T2 {
  input:
  set val(ID) ,file(tumor2), file(normal) from bams_T2N
  file ref
  file regions

  output:
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.indels.vcf.gz' into VCF_somatic2_indels
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.snvs.vcf.gz' into VCF_somatic2_snvs
  set val(ID_tag), file 'strelkaAnalysis/results/variants/*.tbi' into TBI_somatic2
  shell:
  '''
  ID_tag=!{ID}
 !{strelka_somatic} --tumorBam=!{tumor2} --normalBam=!{normal} --referenceFasta=!{params.ref} --callRegions=!{params.regions} --callMemMb=1024  --runDir strelkaAnalysis
  cd strelkaAnalysis
     ./runWorkflow.py -m local -j !{params.cpu} -g !{params.mem}
     cd results/variants
     mv somatic.indels.vcf.gz !{tumor2.baseName}.somatic.indels.vcf.gz
     mv somatic.snvs.vcf.gz !{tumor2.baseName}.somatic.snvs.vcf.gz
     mv somatic.indels.vcf.gz.tbi !{tumor2.baseName}.somatic.indels.vcf.gz.tbi
     mv somatic.snvs.vcf.gz.tbi !{tumor2.baseName}.somatic.snvs.vcf.gz.tbi
   
   
    {params.bcftools} norm  -f !{tumor2.baseName}.somatic.indels.vcf.gz
    !{params.bcftools} view -i'FILTER="PASS"' !{tumor2.baseName}.somatic.indels.vcf.gz  > !{tumor2.baseName}.somatic.indels.vcf.gz
    !{params.tabix} -p vcf !{tumor2.baseName}.somatic.indels.vcf.gz
    
    {params.bcftools} norm  -f !{tumor2.baseName}.somatic.snvs.vcf.gz
    !{params.bcftools} view -i'FILTER="PASS"' !{tumor2.baseName}.somatic.snvs.vcf.gz  >  !{tumor2.baseName}.somatic.snvs.vcf.gz
    !{params.tabix} -p vcf !{tumor2.baseName}.somatic.snvs.vcf.gz
  '''
}

VCF_somatic = VCF_somatic1.join(VCF_somatic2)
input_somaticCoverage = bams_TT.join(VCF_somatic)


process somatic_tumor_coverage {
   input:
  set val(ID),file (bamtumor1),file(bamtumor2), file(somaticVCF1),file(somaticVCF2) from input_somaticCoverage
  file ref
  file regions
  
  
  output:
  set val(ID),file("${bamtumor1.baseName}.vcf"),file("${bamtumor2.baseName}.vcf") into coverage_somatic
  
  shell :
  '''
 !{strelka_germloine} --bam=!{bamtumor1} --bam !{bamtumor2} --forcedGT !{somaticVCF1} --forcedGT !{somaticVCF2}  --referenceFasta=!{params.ref}   --callRegions=!{params.regions} --runDir strelkaAnalysis

  '''
}
