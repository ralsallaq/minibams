/*=================== workflow to get a manifest for each target ======================*/
//defaults
params.get_clean_manifest = false
params.genome = null
params.lookupDircs = null 

params.locus_band_indel = null
params.locus_band_snv = null
params.locus_band_indel_sv = null
params.locus_band_sv = null
params.locus_band_del = 10 //wiggle room before and after a del 

// generate ground-truth manifests

include { get_cytoband_loci as get_cytoband_loci1 } from './get_cytoband_loci' params(genome:params.genome, lookupDircs:params.lookupDircs[0])
include { get_cytoband_loci as get_cytoband_loci2 } from './get_cytoband_loci' params(genome:params.genome, lookupDircs:params.lookupDircs[1])

include { get_clean_manifest as get_clean_manifest } from './get_minibams' params(genome:params.genome, locus_band_indel:params.locus_band_indel,  locus_band_indel_sv:params.locus_band_indel_sv, locus_band_sv:params.locus_band_sv, locus_band_del:params.locus_band_del, outD:params.outD)

include { get_GT_manifest as get_GT_manifest_RNA_tumor } from './get_minibams' params( outD:params.outD) 

include { get_GT_manifest as get_GT_manifest_WGS_tumor } from './get_minibams' params( outD:params.outD) 
include { get_GT_manifest as get_GT_manifest_WGS_gl } from './get_minibams' params( outD:params.outD) 

include { get_GT_manifest as get_GT_manifest_WES_tumor } from './get_minibams' params( outD:params.outD) 
include { get_GT_manifest as get_GT_manifest_WES_gl } from './get_minibams' params( outD:params.outD) 

// generate subbams
include { generateSubBams as generateSubBams_RNA_tumor } from './get_minibams' params( outD:params.outD)

include { generateSubBams as generateSubBams_WGS_tumor } from './get_minibams' params( outD:params.outD)
include { generateSubBams as generateSubBams_WGS_gl } from './get_minibams' params( outD:params.outD)

include { generateSubBams as generateSubBams_WES_tumor } from './get_minibams' params( outD:params.outD)
include { generateSubBams as generateSubBams_WES_gl } from './get_minibams' params( outD:params.outD)


// merge subbams to minibams
include { mergeToMinibams as mergeToMinibams_RNA_tumor } from './get_minibams' params(genome:params.genome, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WGS_tumor } from './get_minibams' params(genome:params.genome, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WGS_gl } from './get_minibams' params(genome:params.genome, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WES_tumor } from './get_minibams' params(genome:params.genome, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WES_gl } from './get_minibams' params(genome:params.genome, outD:params.outD)

workflow get_homog_manifests {
    take:abnormal_events_ch

    main:
    if (params.get_clean_manifest) {

        get_cytoband_loci1(abnormal_events_ch)
        get_cytoband_loci2(abnormal_events_ch)
        events_capture_ch = get_cytoband_loci1.out.map{it->it[4]}.concat(get_cytoband_loci2.out.map{it->it[4]}).collectFile(name:'sample_gene_locus_raw.txt',storeDir:params.outD, newLine:true)
        // get clean manifest for downstream analyses
        get_clean_manifest(events_capture_ch)

        // get ground truth manifests by analysis-type and sample type
    
        get_GT_manifest_RNA_tumor(get_clean_manifest.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
        get_GT_manifest_WGS_tumor(get_clean_manifest.out, channel.value('WHOLE_GENOME'), channel.value('tumor'))
        get_GT_manifest_WGS_gl(get_clean_manifest.out, channel.value('WHOLE_GENOME'), channel.value('germline'))
    
        get_GT_manifest_WES_tumor(get_clean_manifest.out, channel.value('EXOME'), channel.value('tumor'))
        get_GT_manifest_WES_gl(get_clean_manifest.out, channel.value('EXOME'), channel.value('germline'))

    } else {
        get_clean_manifest_ch =  channel.fromPath(params.clean_manifest_file).splitCsv(header: true, sep: "\t", strip: true)

        // get ground truth manifests by analysis-type and sample type
    
        get_GT_manifest_RNA_tumor(get_clean_manifest_ch, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
        get_GT_manifest_WGS_tumor(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('tumor'))
        get_GT_manifest_WGS_gl(get_clean_manifest_ch, channel.value('WHOLE_GENOME'), channel.value('germline'))
    
        get_GT_manifest_WES_tumor(get_clean_manifest_ch, channel.value('EXOME'), channel.value('tumor'))
        get_GT_manifest_WES_gl(get_clean_manifest_ch, channel.value('EXOME'), channel.value('germline'))

    }



   
    // generate subBams for merging 
    generateSubBams_RNA_tumor(get_GT_manifest_RNA_tumor.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    generateSubBams_WGS_tumor(get_GT_manifest_WGS_tumor.out, channel.value('WHOLE_GENOME'), channel.value('tumor'))
    generateSubBams_WGS_gl(get_GT_manifest_WGS_gl.out, channel.value('WHOLE_GENOME'), channel.value('germline'))

    generateSubBams_WES_tumor(get_GT_manifest_WES_tumor.out, channel.value('EXOME'), channel.value('tumor'))
    generateSubBams_WES_gl(get_GT_manifest_WES_gl.out, channel.value('EXOME'), channel.value('germline'))


    // merge suBams to minibams 
    mergeToMinibams_RNA_tumor(generateSubBams_RNA_tumor.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'))
    
    mergeToMinibams_WGS_tumor(generateSubBams_WGS_tumor.out, channel.value('WHOLE_GENOME'), channel.value('tumor'))
    mergeToMinibams_WGS_gl(generateSubBams_WGS_gl.out, channel.value('WHOLE_GENOME'), channel.value('germline'))

    mergeToMinibams_WES_tumor(generateSubBams_WES_tumor.out, channel.value('EXOME'), channel.value('tumor'))
    mergeToMinibams_WES_gl(generateSubBams_WES_gl.out, channel.value('EXOME'), channel.value('germline'))

    //emit:
    // should emit one homogeneous manifest file at a time
    //get_minibams.out
} //end of workflow
