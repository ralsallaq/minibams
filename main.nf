#!/usr/bin/env nextflow

// enable DSL2
nextflow.enable.dsl = 2

params.outD = "analysis/"
params.bams_base = "/research_jude/rgs01_jude/resgen/dev/tartan/index/data/"
params.abnormal_eventsFile = "${workflow.projectDir}/locii_input_file.tsv"
params.lookupDircs = ['/clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/','/clinical/cgs01/clingen/prod/tartan/index/data/Clinical/2021']
params.genome = 'hg19'
params.get_clean_manifest = false
params.clean_manifest_file = "analysis/manifest_all_wRegions.tsv"

if (! params.get_clean_manifest) {
    if (params.clean_manifest_file == null) {
        log.error("if params.get_clean_manifest is false, params.clean_manifest_file is expected to be a clean manifest path")
        exit 1
       }
 } else if (params.abnormal_eventsFile == null) {
     log.error("if params.get_clean_manifest is true, params.abnormal_eventsFile is expected to be a path to file with events to be collected from tartan")
     exit 1
} else {
     log.error("unrecognized value for params.get_clean_manifest")
     exit 1
}


//parameters for wiggle room around locus
params.locus_band_indel = 10
params.locus_band_indel_sv = 350 //this example shows it must be 350
params.locus_band_del = 10 //wiggle room before and after a del 
//SJCBF100_E  21  36231771    -   60  8   93029591    -   44  CTX GENIC   RUNX1_RUNX1T1
//SJCBF100_E  21  36231771    -   4   8   93074937    -   9   CTX GENIC   RUNX1_RUNX1T1

params.locus_band_sv = 500

#### random seed to generate a random subject id
params.random_seed = 900


//workflows
//include { get_homog_manifests } from './modules/get_homog_manifests' params(genome:params.genome, lookupDircs:params.lookupDircs, outD:params.outD, locus_band_indel:params.locus_band_indel,  locus_band_indel_sv:params.locus_band_indel_sv, locus_band_sv:params.locus_band_sv) 

// generate ground-truth manifests

include { get_cytoband_loci as get_cytoband_loci1 } from './modules/get_cytoband_loci' params(genome:params.genome, lookupDircs:params.lookupDircs[0])
include { get_cytoband_loci as get_cytoband_loci2 } from './modules/get_cytoband_loci' params(genome:params.genome, lookupDircs:params.lookupDircs[1])

include { get_clean_manifest as get_clean_manifest } from './modules/get_minibams' params(genome:params.genome, locus_band_indel:params.locus_band_indel,  locus_band_indel_sv:params.locus_band_indel_sv, locus_band_sv:params.locus_band_sv, locus_band_del:params.locus_band_del, outD:params.outD)

include { validate_manifest as validate_manifest } from './modules/get_minibams' params(genome:params.genome, outD:params.outD)

include { get_GT_manifest as get_GT_manifest_RNA_tumor } from './modules/get_minibams' params( outD:params.outD) 

include { get_GT_manifest as get_GT_manifest_WGS_tumor } from './modules/get_minibams' params( outD:params.outD) 
include { get_GT_manifest as get_GT_manifest_WGS_gl } from './modules/get_minibams' params( outD:params.outD) 

include { get_GT_manifest as get_GT_manifest_WES_tumor } from './modules/get_minibams' params( outD:params.outD) 
include { get_GT_manifest as get_GT_manifest_WES_gl } from './modules/get_minibams' params( outD:params.outD) 

// generate bed files

include { generateBEDFiles as generateBEDFiles_RNA_tumor } from './modules/get_minibams' params( outD:params.outD)

include { generateBEDFiles as generateBEDFiles_WGS_tumor } from './modules/get_minibams' params( outD:params.outD)
include { generateBEDFiles as generateBEDFiles_WGS_gl } from './modules/get_minibams' params( outD:params.outD)

include { generateBEDFiles as generateBEDFiles_WES_tumor } from './modules/get_minibams' params( outD:params.outD)
include { generateBEDFiles as generateBEDFiles_WES_gl } from './modules/get_minibams' params( outD:params.outD)


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



    if (params.get_clean_manifest) {

        inputFile_ch = channel.fromPath(params.abnormal_eventsFile).splitCsv(header: true, sep: "\t", strip: true)
        inputFile_ch.map{ event -> [event.sample, event.gene, event.cytoLocus, event.type] }.set {abnormal_events_ch}

        //abnormal_events_ch.view()

        get_cytoband_loci1(abnormal_events_ch)
        get_cytoband_loci2(abnormal_events_ch)
        events_capture_ch = get_cytoband_loci1.out.map{it->it[4]}.concat(get_cytoband_loci2.out.map{it->it[4]}).collectFile(name:'sample_gene_locus_raw.txt',storeDir:params.outD, newLine:true)

        // get clean manifest for downstream analyses
        get_clean_manifest(events_capture_ch)

        // validate the manifest
        validate_manifest(get_clean_manifest.out)

        // get ground truth manifests by analysis-type and sample type
        get_GT_manifest_RNA_tumor(get_clean_manifest.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
        get_GT_manifest_WGS_tumor(get_clean_manifest.out, channel.value('WHOLE_GENOME'), channel.value('tumor'))
        get_GT_manifest_WGS_gl(get_clean_manifest.out, channel.value('WHOLE_GENOME'), channel.value('germline'))
    
        get_GT_manifest_WES_tumor(get_clean_manifest.out, channel.value('EXOME'), channel.value('tumor'))
        get_GT_manifest_WES_gl(get_clean_manifest.out, channel.value('EXOME'), channel.value('germline'))

    } else {

        get_clean_manifest_ch =  channel.fromPath(params.clean_manifest_file) //.splitCsv(header: true, sep: "\t", strip: true)

        // validate the manifest
        //validate_manifest(get_clean_manifest_ch)

        // get ground truth manifests by analysis-type and sample type
    
        get_GT_manifest_RNA_tumor(get_clean_manifest_ch, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
        get_GT_manifest_WGS_tumor(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('tumor'))
        get_GT_manifest_WGS_gl(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('germline'))
    
        get_GT_manifest_WES_tumor(get_clean_manifest_ch, channel.value('EXOME'), channel.value('tumor'))
        get_GT_manifest_WES_gl(get_clean_manifest_ch, channel.value('EXOME'), channel.value('germline'))

    }



    // generate bed files for regions to be extracted in a next step 
    generateBEDFiles_RNA_tumor(get_GT_manifest_RNA_tumor.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    generateBEDFiles_WGS_tumor(get_GT_manifest_WGS_tumor.out, channel.value('WHOLE_GENOME'), channel.value('tumor'))
    generateBEDFiles_WGS_gl(get_GT_manifest_WGS_gl.out, channel.value('WHOLE_GENOME'), channel.value('germline'))

    generateBEDFiles_WES_tumor(get_GT_manifest_WES_tumor.out, channel.value('EXOME'), channel.value('tumor'))
    generateBEDFiles_WES_gl(get_GT_manifest_WES_gl.out, channel.value('EXOME'), channel.value('germline'))

   
    // generate subBams for merging one subBam for each bed file passed
    generateSubBams_RNA_tumor(generateBEDFiles_RNA_tumor.out.flatten(), channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    generateSubBams_WGS_tumor(generateBEDFiles_WGS_tumor.out.flatten(), channel.value('WHOLE_GENOME'), channel.value('tumor'))
    generateSubBams_WGS_gl(generateBEDFiles_WGS_gl.out.flatten(), channel.value('WHOLE_GENOME'), channel.value('germline'))

    generateSubBams_WES_tumor(generateBEDFiles_WES_tumor.out.flatten(), channel.value('EXOME'), channel.value('tumor'))
    generateSubBams_WES_gl(generateBEDFiles_WES_gl.out.flatten(), channel.value('EXOME'), channel.value('germline'))


    // collect all subBams and merge to minibams 
    mergeToMinibams_RNA_tumor(generateSubBams_RNA_tumor.out.collect(), channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    mergeToMinibams_WGS_tumor(generateSubBams_WGS_tumor.out.collect(), channel.value('WHOLE_GENOME'), channel.value('tumor'))
    mergeToMinibams_WGS_gl(generateSubBams_WGS_gl.out.collect(), channel.value('WHOLE_GENOME'), channel.value('germline'))

    mergeToMinibams_WES_tumor(generateSubBams_WES_tumor.out.collect(), channel.value('EXOME'), channel.value('tumor'))
    mergeToMinibams_WES_gl(generateSubBams_WES_gl.out.collect(), channel.value('EXOME'), channel.value('germline'))

}
