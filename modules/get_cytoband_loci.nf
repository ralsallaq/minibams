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

            #### this will look for the fusion in many places: cicero-itd-post, crest-post, is Somatic SV is in the name of the file
            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}" | sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'\$fusion'" {}/*.txt >> sample_records/"'\$outputFile'" '

            if [[ \$(cat sample_records/"\$outputFile"|wc -l) -eq 0 ]]; then
            
                echo "could not find records with somatic SV in their name; will look into crest-post"
                cat \$recFile | grep -w 'crest-post'| parallel ' grep -iwE "'\$fusion'" {}/*-event_fusion.txt  >> sample_records/"'\$outputFile'" '  
            fi


        elif [ "${type}" == "Somatic SV" ] && [ "${cytoLocus}" == "Deletion or Disruption" ]; then

            grep -iE \$(echo "${type}"|cut -f 1 -d " ") \$recFile | grep -iE \$(echo "${type}"|sed "s/CNA/CNV/g"|cut -f 2 -d " ") | parallel ' grep -w "'${gene}'" {}/*.txt | grep -iE "DEL" >> sample_records/"'\$outputFile'" '

            if [[ \$(cat sample_records/"\$outputFile"|wc -l) -eq 0 ]]; then

                echo "could not find records with somatic SV in their name; will look into all text files"
                cat \$recFile | grep -w 'crest-post'| parallel ' grep -iwE "'${gene}'" {}/*-event_fusion.txt | grep -iE "DEL"  >> sample_records/"'\$outputFile'" '  
            fi


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
                 cat list_of_bamPaths.info_unflitered | sort | uniq |while read target sample bdir
                 do 
                     if [ -f "\$bdir/\$sample.bam" ]; then 
                         echo -e "\$target \\t \$sample \\t \$bdir" >> list_of_bamPaths.info
                     else
                         echo \$bdir >> list_of_missing_bamPaths.info
                     fi
                 done
            fi
    }
 

    function prepManifest() {

        module load python/3.7.0
        eventsFileFound=false
        while [[ "\$eventsFileFound" == "false" ]]
            do
                
                if [ -f events_capture.info ]; then
                    eventsFileFound=true
                    echo "python ${workflow.projectDir}/bin/prepManifest.py -i events_capture.info -b list_of_bamPaths.info -o events_locii_fixed.tsv"
                    python ${workflow.projectDir}/bin/prepManifest.py -i events_capture.info -b list_of_bamPaths.info -o events_locii_fixed.tsv
                else
                    >&2 echo "events_capture.info is not done yet or does not exist"
                    sleep 1
                fi
            done

        }


    getRecordsForSamples
    echo "done getting \$recordsFound records"
    getEvents
    echo "done getting \$eventsFound events"
    getAvailableBams
    echo "done getting available bams"

    if [[ "\$(cat events_capture.info|cut -f 7| sed 's/^[[:space:]]*//g')" == "0" ]]; then

        echo 'no need to fix anything empty events'
        touch events_locii_fixed.tsv

    else
        prepManifest

    fi

    echo "done preparing raw manifest encompassing available bams"
    echo "assuring that the number of lines in the manifest is either 0 or multiples of 5"


    echo "Done"

    
    """
}
