#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl = 2

params.outD = "analysis_addLoci_2022_02_03/"
params.tartanIndex_projDir = "/clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/"
params.genome = 'hg19'

params.exBedFile = "bedFile_regionsToAdd.bed"
params.minibams_dir = "analysis_save/minibams/"


//parameters for wiggle room around locus
params.locus_band_indel = 10
params.locus_band_indel_sv = 350 //this example shows it must be 350
params.locus_band_del = 10 //wiggle room before and after a del 
//SJCBF100_E  21  36231771    -   60  8   93029591    -   44  CTX GENIC   RUNX1_RUNX1T1
//SJCBF100_E  21  36231771    -   4   8   93074937    -   9   CTX GENIC   RUNX1_RUNX1T1

params.locus_band_sv = 500

// random seed to generate a random subject id
params.random_seed = 901


//workflows

// generate bed files
include { generateBedFilesFromExBed } from './modules/utilities.nf' params(genome:params.genome, locus_band_indel:params.locus_band_indel,  locus_band_indel_sv:params.locus_band_indel_sv, locus_band_sv:params.locus_band_sv, locus_band_del:params.locus_band_del, tartanIndex_projDir:params.tartanIndex_projDir, outD:params.outD) 

// generate subbams
include { generateSubBams as generateSubBams_RNA_tumor } from './modules/get_minibams' params( outD:params.outD)

include { generateSubBams as generateSubBams_WGS_tumor } from './modules/get_minibams' params( outD:params.outD)
include { generateSubBams as generateSubBams_WGS_gl } from './modules/get_minibams' params( outD:params.outD)

include { generateSubBams as generateSubBams_WES_tumor } from './modules/get_minibams' params( outD:params.outD)
include { generateSubBams as generateSubBams_WES_gl } from './modules/get_minibams' params( outD:params.outD)


// merge subbams to minibams
include { mergeToMinibams as mergeToMinibams_RNA_tumor } from './modules/get_minibams' params(genome:params.genome, random_seed:params.random_seed, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WGS_tumor } from './modules/get_minibams' params(genome:params.genome, random_seed:params.random_seed, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WGS_gl } from './modules/get_minibams' params(genome:params.genome, random_seed:params.random_seed, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WES_tumor } from './modules/get_minibams' params(genome:params.genome, random_seed:params.random_seed, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WES_gl } from './modules/get_minibams' params(genome:params.genome, random_seed:params.random_seed, outD:params.outD)

workflow {
    main:

    // read external bed file
    extBedFile_ch = channel.fromPath(params.exBedFile)
    
    // read already generated minibams using main 
    minibam_RNA_tumor_ch = channel.fromPath(params.minibams_dir + '/TRANSCRIPTOME_' + params.genome + '/*_D1.bam')

    minibam_WGS_tumor_ch = channel.fromPath(params.minibams_dir + '/WHOLE_GENOME_' + params.genome + '/*_D1.bam')
    minibam_WGS_gl_ch = channel.fromPath(params.minibams_dir + '/WHOLE_GENOME_' + params.genome + '/*_G1.bam')

    minibam_WES_tumor_ch = channel.fromPath(params.minibams_dir + '/EXOME_' + params.genome + '/*_D1.bam')
    minibam_WES_gl_ch = channel.fromPath(params.minibams_dir + '/EXOME_' + params.genome + '/*_G1.bam')

  
    // generate bed files for regions to be extracted in a next step 
    generateBedFilesFromExBed(extBedFile_ch)

    //generateBedFilesFromExBed.out.map{it->it[1]}.view()
    // generate subBams for merging one subBam for each bed file passed
    generateSubBams_RNA_tumor(generateBedFilesFromExBed.out.map{it -> it[0]}.flatten(), channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    generateSubBams_WGS_tumor(generateBedFilesFromExBed.out.map{it -> it[1]}.flatten(), channel.value('WHOLE_GENOME'), channel.value('tumor'))
    generateSubBams_WGS_gl(generateBedFilesFromExBed.out.map{it -> it[2]}.flatten(), channel.value('WHOLE_GENOME'), channel.value('germline'))

    generateSubBams_WES_tumor(generateBedFilesFromExBed.out.map{it -> it[3]}.flatten(), channel.value('EXOME'), channel.value('tumor'))
    generateSubBams_WES_gl(generateBedFilesFromExBed.out.map{it -> it[4]}.flatten(), channel.value('EXOME'), channel.value('germline'))


    // collect all subBams mix in previously generated minibams and merge to new minibams 
    mergeToMinibams_RNA_tumor(generateSubBams_RNA_tumor.out.concat(minibam_RNA_tumor_ch).collect(), channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    mergeToMinibams_WGS_tumor(generateSubBams_WGS_tumor.out.concat(minibam_WGS_tumor_ch).collect(), channel.value('WHOLE_GENOME'), channel.value('tumor'))
    mergeToMinibams_WGS_gl(generateSubBams_WGS_gl.out.concat(minibam_WGS_gl_ch).collect(), channel.value('WHOLE_GENOME'), channel.value('germline'))

    mergeToMinibams_WES_tumor(generateSubBams_WES_tumor.out.concat(minibam_WES_tumor_ch).collect(), channel.value('EXOME'), channel.value('tumor'))
    mergeToMinibams_WES_gl(generateSubBams_WES_gl.out.concat(minibam_WES_gl_ch).collect(), channel.value('EXOME'), channel.value('germline'))
}
