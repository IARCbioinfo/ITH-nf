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
    log.info "--strelka              PATH                Path to strelka configureStrelkaSomaticWorkflow.py "
    log.info "--platypus              PATH               platypus installation dir"
    log.info "--R              PATH               R installation dir"
    log.info "--vt              PATH               vt installation dir"
    log.info "--annovar              PATH               annovar installation dir"
    log.info "--tabix             PATH               tabix installation dir"
    log.info "--falcon_qc             PATH            falcon.qc.r dir"
    log.info "--bam_folder         FOLDER               Folder containing bam files "
    log.info "--correspondance		FILE				File containing the correspondance between the normal and two tumor samples and the sample id"
    log.info "--ref                 FILE                Reference file"
    log.info "--regions             FILE                 Regions "
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


ref = file(params.ref)
regions = 	file(params.regions)
correpondance = file(params.correpondance)

(bams_TTN,bams_TTN2 )= Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor2),file(params.bam_folder + "/" + row.normal) ]}.into(2)

bams_TTN.subscribe { println "value: $it" }


(bams_TT,bams_TT2) =  Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.tumor1),file(params.bam_folder + "/" + row.tumor2) ]}.into(2)
		
bams_N = Channel.fromPath(correpondance).splitCsv(header: true, sep: '\t', strip: true)
		.map{row -> [ row.ID,file(params.bam_folder + "/" + row.normal) ]}
		
	
		
		
process germline_calling {
  input:
  set val(ID) ,file (normal) from bams_N
  file ref
  file regions

  output:
  set val(ID_tag),file("${normal.baseName}.vcf") into VCF_germline

  shell:
  '''
  ID_tag=!{ID}
  !{params.platypus} callVariants --bamFiles=!{normal} --output=!{normal.baseName}.vcf --refFile=!{params.ref} --regions=!{regions} --badReadsThreshold=0 --qdThreshold=0 --rmsmqThreshold=20 --hapScoreThreshold=10 --scThreshold=0.99
  '''
}

process germline_calling_pass {
  input:
  set val(ID) ,file(germlinevcf) from VCF_germline

  output:
  set val(ID_tag),file("${germlinevcf.baseName}.vcf") into (VCF_germline_passQC,VCF_germline_passCoverage)

  shell:
  '''
  ID_tag=!{ID}
cat !{germlinevcf} | grep "PASS" > !{germlinevcf.basename}.vcf
'''
}



process somatic_calling {
  input:
  set val(ID) ,file (tumor1), file(tumor2), file(normal) from bams_TTN
  file ref
  file regions

  output:
  set val(ID_tag),file("${tumor1.baseName}.vcf"),file("${tumor2.baseName}.vcf") into VCF_somatic

  shell:
  '''
  ID_tag=!{ID}
 !{params.strelka} --tumorBam=!{tumor1} --normalBam=!{normal} --referenceFasta=!{params.ref} --callRegions=!{params.regions} --callMemMb=1024  -m local -j 28
 !{params.strelka} --tumorBam=!{tumor2} --normalBam=!{normal} --referenceFasta=!{params.ref} --callRegions=!{params.regions} --callMemMb=1024  -m local -j 28
  '''
}

process somatic_calling_pass {
  
  input:
  set val(ID) ,file(somaticVCF1),file(somaticVCF2) from VCF_somatic

  output:
  set val(ID_tag),file("${somaticVCF1.baseName}.vcf"), file("${somaticVCF2.baseName}.vcf") into (VCF_somatic_passQC,VCF_somatic_passCoverage)

  shell:
  '''
  ID_tag=!{ID}
	cat !{somaticVCF1} | grep "PASS" > !{somaticVCF1}.vcf
'''
}

input_germlineCoverage = bams_TTN2.join(VCF_germline_passCoverage)
process germline_tumor_coverage {
  input:
  set val(ID),file (bamtumor1),file(bamtumor2),file(bamnormal),file(germlineVCF) from input_germlineCoverage
  file ref
  file regions

  
  output:
  set val(ID),file("${ID}.vcf") into coverage_germline
}
shell :
'''

Platypus.py callVariants --bamFiles=!{bamnormal},!{bamtumor1},!{bamtumor2} --refFile=!{input_ref} --regions=!{params.regions} --nCPU=12 --output=!{ID}.coverage.germline.vcf --source=!{germlineVCF} --minPosterior=0 --getVariantsFromBAMs=0
'''
}

input_somaticCoverage=bams_TT2.joint(VCF_somatic_passCoverage)
process somatic_tumor_coverage {
   input:
  set val(ID),file (bamtumor1),file(bamtumor2), file(somaticVCF1),file(somaticVCF2) from input_somaticCoverage
  file ref
  file regions
  
  
  output:
  set val(ID),file("${bamtumor1.baseName}.vcf"),file("${bamtumor2.baseName}.vcf") into coverage_somatic
  
  shell :
  '''
Platypus.py callVariants --bamFiles=!{bamtumor1} --refFile=!{params.ref} --regions=!{params.regions} --nCPU=12 --output=!{bamtumor1.baseName}.vcf --source=!{somaticVCF2} --minPosterior=0 --getVariantsFromBAMs=0
Platypus.py callVariants --bamFiles=!{bamtumor2} --refFile=!{params.ref} --regions=!{params.regions} --nCPU=12 --output=!{bamtumor2.baseName}.vcf --source=!{somaticVCF1} --minPosterior=0 --getVariantsFromBAMs=0

  '''
}
