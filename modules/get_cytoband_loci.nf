/*=================== workflow to get cytobands of localized abnormal events ======================*/

params.genome = null
params.lookupDircs = null 

process get_cytoband_loci {
    label 'io_limited'
    //publishDir "${params.outD}", mode: 'copy'

    input:
    tuple val(sample), val(gene), val(cytoLocus), val(type)

    output:
    tuple  val(sample), val(gene), val(cytoLocus), val(type), file ("events_locii_fixed.tsv"), emit:cytoband_loci_ch
    

    """

    #!/usr/bin/bash
    module load parallel/20201222

    
    if [[ ${params.genome} == "hg38"  || ${params.genome} == "GRCh38" ]]; then
       export official_genome="GRCh38"
    elif [[ ${params.genome} == "hg19" ||  ${params.genome} == "GRCh37-lite" ]]; then
       export official_genome="GRCh37-lite"
    else
       echo "the genome ${params.genome} is not recognized!"
       exit 1
    fi

    echo -e "official genome is \$official_genome"
    echo -e "lookupDir is ${params.lookupDircs}"

    function getRecordsForSamples() {
        #ingore errors of not finding something
        ### looks under two directories in clingen for the sample data
        ### and collects and saves a list of records for the sample
        #### the types for the samples

        mkdir -p sample_records

        echo getting records for ${sample}

        recs=\$(find "${params.lookupDircs}"/"${sample}"*) || recs="" 

        echo \$recs|tr " " "\\n" >> sample_records/"${sample}".records

        unset recordsFound
        recordsFound=\$(cat sample_records/"${sample}".records|wc -l) 

   }

    function getEvents() {
        ### goes through cytoLocus events for samples and saves a file encompassing the event

        rm -f events_capture.info
    
        recFile="sample_records/${sample}.records"
        echo ${type}

        outputFile=\$(echo "${sample}_${gene}_${cytoLocus}_${type}.records"|tr "[,\\ \\/\\>]" "_")
        rm -f sample_records/"\$outputFile"
        touch sample_records/"\$outputFile"

        if [ "${type}" == "Somatic SV / CNA" ]; then
            
            genes=\$(echo ${gene}|tr "\\/\\ " "|")

            ### this covers deletions
            cat \$recFile | grep -w 'conserting_crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "DEL" >> sample_records/"'\$outputFile'" '  
            ### this seems to cover amplifications
            cat \$recFile | grep -w 'conserting_crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "ITX" >> sample_records/"'\$outputFile'" '  
            cat \$recFile | grep -w 'conserting_crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "CTX" >> sample_records/"'\$outputFile'" '  
            ### this seems to cover duplications
            cat \$recFile | grep -w 'conserting_crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "INS" >> sample_records/"'\$outputFile'" '  

            if [[ \$(cat sample_records/"\$outputFile"|wc -l) -eq 0 ]]; then
                
                echo "could not find any in conserting_crest-post will look into crest-post"
                cat \$recFile | grep -w 'crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "DEL" >> sample_records/"'\$outputFile'" '  
                ### this seems to cover amplifications
                cat \$recFile | grep -w 'crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "ITX" >> sample_records/"'\$outputFile'" '  
                cat \$recFile | grep -w 'crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "CTX" >> sample_records/"'\$outputFile'" '  
                ### this seems to cover duplications
                cat \$recFile | grep -w 'crest-post' |parallel ' grep -iwE "'\$genes'" {}/*-event_fusion.txt | grep -iE "INS" >> sample_records/"'\$outputFile'" '  
            fi

        elif [ "${type}" == "Somatic INDEL/SV" ]; then

            grep -iE "Somatic" \$recFile | grep 'sv-somatic-classified' | parallel ' grep -w "'${gene}'" {}/*.txt| grep -iE "'${cytoLocus}'" >> sample_records/"'\$outputFile'" ' 

        elif [ "${type}" == "Somatic CNA" ] && [ "${cytoLocus}" == "Deletion or Disruption" ]; then

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"|sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*.txt | grep -iE "Del" >> sample_records/"'\$outputFile'" '


        elif [ "${type}" == "Germline CNA" ] && [ "${cytoLocus}" == "Deletion, Intragenic" ]; then
            ## example SJMB030020_G1_gl_conserting_unpaired_smcls.txt

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"|sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*.txt | grep -iE "Del" >> sample_records/"'\$outputFile'" '


        elif [ "${type}" == "Somatic SV" ] && ([ "${cytoLocus}" == "Fusion" ] || [ "${cytoLocus}" == "Missed Fusion" ] || [ "${cytoLocus}" == "Position Effect" ] || [ "${cytoLocus}" == "Amplification / Disruption" ]); then

            fusion=\$(echo ${gene}|tr "-" "_")
            echo "fusion is \$fusion"

            #grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}" | sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'\$fusion'" {}/*_predicted_sv.txt >> sample_records/"'\$outputFile'" '
            #### this will look for the fusion in many places: cicero-itd-post, crest-post,
            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}" | sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'\$fusion'" {}/*.txt >> sample_records/"'\$outputFile'" '


        elif [ "${type}" == "Somatic SV" ] && [ "${cytoLocus}" == "Deletion or Disruption" ]; then

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"|sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*.txt | grep -iE "Del" >> sample_records/"'\$outputFile'" '


        elif [ "${type}" == "Somatic SNV" ]; then

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"| sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*tier1* | grep -iE "'${cytoLocus}'" >> sample_records/"'\$outputFile'" '


        elif [ "${type}" == "Somatic INDEL" ]; then

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"| sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*_putative_exon_somatic_indel_mutation-multi_target.txt | grep -iE "'${cytoLocus}'" >> sample_records/"'\$outputFile'" '


        else

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"|sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*.txt | grep -iE "'${cytoLocus}'" >> sample_records/"'\$outputFile'" '


        fi

        unset eventsFound
        eventsFound=\$(cat sample_records/"\$outputFile"|wc -l) 

        wait
        echo "getting bam paths regardless to if bams exist or not under these paths"
        cat \$recFile | grep "\\/bam" | grep -E "WHOLE_GENOME|EXOME|TRANSCRIPTOME" | parallel 'echo -e \$(basename {//}) "'\t'" \$(basename \$(dirname {//})) "'\t'" {}' >   list_of_bamPaths.info_unflitered
        echo "collecting info on \$outputFile"
        echo -e "${params.lookupDircs}\t"  "${sample}\t" "${gene}\t" "${cytoLocus}\t" "${type}\t" "\$outputFile\t" "\$(cat sample_records/\$outputFile|wc -l)" > events_capture.info


        }


    function getAvailableBams() {
            nbams=\$(cat list_of_bamPaths.info_unflitered|wc -l)
            rm -f list_of_bamPaths.info
            rm -f list_of_missing_bamPaths.info
            if [[ \$nbams -eq 0 ]]; then
                ### create empty file
                touch list_of_bamPaths.info
            else
                 #### subset for uniq and available bams
                 cat list_of_bamPaths.info_unflitered | sort | uniq > list_of_bamPaths.info_unflitered_uniq
                 unset nlines
                 nlines=\$(cat list_of_bamPaths.info_unflitered_uniq|wc -l)
                 ### check bam one by one and collect those exit as well as those do not exist
                 for ((i=1;i<=\$nlines;i++)); do 
                     line=\$(cat list_of_bamPaths.info_unflitered_uniq|head -\$i|tail -1)
                     bdir=\$(cat list_of_bamPaths.info_unflitered_uniq|head -\$i|tail -1|awk '{print \$3}')
                     if [ -f \$bdir/*.bam ]; then 
                         echo \$line >> list_of_bamPaths.info
                     else
                         echo \$line >> list_of_missing_bamPaths.info
                     fi
                 done
            fi
    }
 
    function prepManifest() {
        rm -f events_locii.tsv
        event=\$(cat events_capture.info)
        eventFile=\$(cat events_capture.info|cut -f 6|cut -f 2 -d " ")
        nDuplicates=\$(cat events_capture.info|cut -f 7)
        echo -e \$eventFile \t \$nDuplicates
        if [[ "\$nDuplicates" -gt 0 ]]; then
            echo "there records for this event"
            echo -e "\$eventFile\t\$nDuplicates"
            nbams=\$(cat list_of_bamPaths.info|wc -l)
            ### make each collected event available for each target
            for ((i=1;i<="\$nDuplicates";i++));do
                echo -e "\$event\t\$(cat "sample_records/\$eventFile"|sed 's/.*txt://g' |cut -f1-10|sort|uniq|head -\$i|tail -1)"  > events_locii.tsv_temp
                if [ \$nbams -lt 10 ] || [ \$nbams != "" ]; then
                    paste -d"\\t" <(cat \$(yes events_locii.tsv_temp|head -\$nbams)) <(cat list_of_bamPaths.info) >> events_locii.tsv
                else
                    echo "something is wrong with number of bams \$nbams"
                fi
            done
        ### drop duplicated events
        cat events_locii.tsv|sort|uniq > temp_file
        mv temp_file events_locii.tsv
        else
            echo "no records for the event found ... creating empty file"
            touch events_locii.tsv
        fi
    

        }

    function fix_multiplicity() {

        module load python/3.7.0
        python ${workflow.projectDir}/bin/fix_multiplicity.py -i events_locii.tsv -o events_locii_fixed.tsv

        }


    getRecordsForSamples
    echo "done getting \$recordsFound records"
    getEvents
    echo "done getting \$eventsFound events"
    getAvailableBams
    echo "done getting available bams"
    prepManifest
    echo "done preparing raw manifest encompassing available bams"
    echo "assuring that the number of lines in the manifest is either 0 or multiples of 5"

    if [[ "\$(cat events_locii.tsv|wc -l)" -gt 0 ]]; then

        echo 'fixing multiplicity'
        fix_multiplicity

    else
        echo 'no need to fix anything empty events'
        touch events_locii_fixed.tsv
    fi



    echo "Done"

    
    """
}
