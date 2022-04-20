#!/usr/bin/env nextflow

// Enable DSL2
nextflow.enable.dsl = 2

params.outD = "analysis/"
params.abnormal_eventsFile = "${workflow.projectDir}/locii_input_file.tsv"
params.lookupDircs = '/clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/,/clinical/cgs01/clingen/prod/tartan/index/data/Clinical/2021'
lookupDircs_list = params.lookupDircs?.tokenize(',') - null
params.genome = 'hg19'
// Implies get_clean_manifest when null (default is null)
params.use_manifest_file = null 
// By default the recover mode is false use --recover to set true
params.recover = false 
// By default the validation on the manifest is off use --validate_manifest turn it on
params.validate_manifest = false

//Parameters for wiggle room around locus
params.locus_band_indel = 10
params.locus_band_indel_sv = 350 //this example shows it must be 350
params.locus_band_del = 10 //wiggle room before and after a del 
//SJCBF100_E  21  36231771    -   60  8   93029591    -   44  CTX GENIC   RUNX1_RUNX1T1
//SJCBF100_E  21  36231771    -   4   8   93074937    -   9   CTX GENIC   RUNX1_RUNX1T1

params.locus_band_sv = 500

// random seed to generate a random subject id
params.random_seed = 900

// Parameters for addLociToMinibam_wf
// By default addLociToExisitngMB is false; use --addLociToExisitngMB to set it true
params.addLociToExisitngMB = false
params.tartanIndex_projDir = "/clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/"
params.exBedFile = "bedFile_regionsToAdd.bed"
params.minibams_dir = "analysis_save/minibams/"

if (params.use_manifest_file == null & params.abnormal_eventsFile == null){
        log.error("if --use_manifest_file is null --abnormal_eventsFile should point to a file")
        exit 1
}

if (params.recover & params.abnormal_eventsFile == null & params.validate_manifest) {
    log.error("for --recover option --abnormal_eventsFile is required and --validate_manifest is not allowed")
    exit 1
}

if (params.addLociToExisitngMB) {
    if (params.tartanIndex_projDir == null | params.exBedFile == null | params.minibams_dir == null) {
        log.error("with --addLociToExisitngMB, --exBedFile and --minibams_dir must not be null")
        exit 1
        }
    }



//import addLociToMinibam workflow
include { addLociToMinibam_wf } from './modules/addLociToMinibam_wf' params(outD:'analysis/add_loci/', \
                 tartanIndex_projDir:'/clinical/cgs01/clingen/prod/tartan/index/data/ClinicalPilot/ClinicalPilot/', \
                 exBedFile:'bedFile_regionsToAdd.bed', minibams_dir:'analysis/minibams/', \
                 genome:params.genome, locus_band_indel:params.locus_band_indel, locus_band_indel_sv:params.locus_band_indel_sv, \
                 locus_band_sv:params.locus_band_sv, random_seed:params.random_seed)

// generate ground-truth manifests

include { get_cytoband_loci as get_cytoband_loci1 } from './modules/get_cytoband_loci' params(genome:params.genome, lookupDircs:lookupDircs_list[0])
include { get_cytoband_loci as get_cytoband_loci2 } from './modules/get_cytoband_loci' params(genome:params.genome, lookupDircs:lookupDircs_list[1])

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

    //If this run is for adding loci to already existing minibams

    if (params.addLociToExisitngMB) {

        addLociToMinibam_wf()
     
    } else {


        if (params.use_manifest_file == null) {
  
            // Get a manifest 
            inputFile_ch = channel.fromPath(params.abnormal_eventsFile).splitCsv(header: true, sep: "\t", strip: true)
            inputFile_ch.map{ event -> [event.sample, event.gene, event.cytoLocus, event.type] }.set {abnormal_events_ch}
    
            //abnormal_events_ch.view()
    
            get_cytoband_loci1(abnormal_events_ch)
            get_cytoband_loci2(abnormal_events_ch)
            events_capture_ch = get_cytoband_loci1.out.map{it->it[4]}.concat(get_cytoband_loci2.out.map{it->it[4]}).collectFile(name:'sample_gene_locus_raw.txt',storeDir:params.outD, newLine:true)
    
            // Get clean manifest for downstream analyses
            get_clean_manifest(events_capture_ch)
    
            // Validate the manifest if asked to
            if (params.validate_manifest) {
                validate_manifest(get_clean_manifest.out[0])
            }
    
            // Get ground truth manifests by analysis-type and sample type
            get_GT_manifest_RNA_tumor(get_clean_manifest.out[0], channel.value('TRANSCRIPTOME'), channel.value('tumor'))
        
            get_GT_manifest_WGS_tumor(get_clean_manifest.out[0], channel.value('WHOLE_GENOME'), channel.value('tumor'))
            get_GT_manifest_WGS_gl(get_clean_manifest.out[0], channel.value('WHOLE_GENOME'), channel.value('germline'))
        
            get_GT_manifest_WES_tumor(get_clean_manifest.out[0], channel.value('EXOME'), channel.value('tumor'))
            get_GT_manifest_WES_gl(get_clean_manifest.out[0], channel.value('EXOME'), channel.value('germline'))
    
        } else {
       
            // Use existing manifest
    
            get_clean_manifest_ch =  channel.fromPath(params.use_manifest_file) //.splitCsv(header: true, sep: "\t", strip: true)
    
            // Validate the manifest if asked to
            if (params.validate_manifest) {
                validate_manifest(get_clean_manifest.out)
            }
    
            // Get ground truth manifests by analysis-type and sample type
        
            get_GT_manifest_RNA_tumor(get_clean_manifest_ch, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
        
            get_GT_manifest_WGS_tumor(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('tumor'))
            get_GT_manifest_WGS_gl(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('germline'))
        
            get_GT_manifest_WES_tumor(get_clean_manifest_ch, channel.value('EXOME'), channel.value('tumor'))
            get_GT_manifest_WES_gl(get_clean_manifest_ch, channel.value('EXOME'), channel.value('germline'))
    
        }
    
    
        // When not in recover mode process the manifest to minibams
        if (! params.recover) {
    
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

      }
}
