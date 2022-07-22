
process splitMinibamByRGs {
    tag "split minibams by RGs"
    cpus 1
    memory 8.GB

    input:
    path(minibam)
    val(analysisType)
    val(sampleType)

    output:
    path("bam_splits/*.bam"), emit: bam_splits_ch
    
    """
    #!/usr/bin/bash 

    which samtools
    if [ \$? != 0 ]
    then
       echo "samtools is not found in path" >&2
       exit 1
    fi

    mkdir bam_splits
    samtools split -@ "$task.cpus" -u bam_splits/reads_with_no_RG -f bam_splits/'%*_%!.%.' ${minibam} 

    if [[ \$(ls -p bam_splits/*.bam|wc -l) == 0 ]]
    then
      >&2 echo "no split bams found"
      exit 1
    fi


    # Make sure there are no reads without RGs; exit if those are found
    nreads_noRGs=\$(samtools view bam_splits/reads_with_no_RG | wc -l)
    if [[ "\$nreads_noRGs" -gt "0" ]]
    then
       echo "there are reads that does not belong to any RGs!! ..exiting" >&2
       exit 1
    fi

    """

    
} 


process bam2fastq {
    tag "convert coord sorted bam to fastqs"
    label 'io_mem'

    input:
    path(bamFile)
    val(analysisType)
    val(sampleType)

    output:
    tuple env(rg), path('RG.txt'), path("*.fastq.gz"), emit: fastqs_from_bam_ch
    
    """
    #!/usr/bin/bash 
    module load biobambam/2.0.87
    which bamtofastq
    if [ \$? != 0 ]
    then
       echo "bamtofastq is not found in path" >&2
       exit 1
    fi

    echo "making sure that the bam file is sorted by coord" >&2
    grep -w 'SO:coordinate' <(samtools view -H ${bamFile} |head -1)
    if [ \$? != 0 ]
    then
       echo "the bam is not sorted by coord ... exiting" >&2
       exit 1
    fi

    # Keep going but save exit code
    EXIT_CODE=0
    samtools quickcheck ${bamFile} || EXIT_CODE=\$?

    echo "quickcheck returned exit code \$EXIT_CODE on ${bamFile}" >&2

    # The assumed for of split is like SJOther992967_D1_RGVALUE.bam 
    outputPrefix=\$(basename ${bamFile} .bam)
    samFile=\${outputPrefix}.sam
    baseName=\$(echo \${outputPrefix} | cut -f 1-2 -d "_")
    #rg=\$(echo \${outputPrefix} | sed "s/'\$baseName'\\_//g")

    
    if [[ \$EXIT_CODE -eq 0 ]]
    then
       echo "saving RG to a text file" >&2
       samtools view -H  ${bamFile} | grep '^@RG' > RG.txt
       rg=\$(cat RG.txt|sed -r 's/.*(ID:)(.*\\s{1})PL.*/\\2/g'|xargs)

       echo "processing ${bamFile} to fastqs \${outputPrefix}.R1.\${rg}.fastq.gz, \${outputPrefix}.R2.\${rg}.fastq.gz and \${outputPrefix}.O1.\${rg}.fastq.gz, \${outputPrefix}.O2.\${rg}.fastq.gz" >&2
       bamtofastq inputformat=bam disablevalidation=1 gz=1 filename="${bamFile}" F="\${outputPrefix}.R1.\${rg}.fastq.gz" F2="\${outputPrefix}.R2.\${rg}.fastq.gz" O="\${outputPrefix}.O1.\${rg}.fastq.gz" O2="\${outputPrefix}.O2.\${rg}.fastq.gz"

       if [[ \$? != 0 ]]; then echo "error in converting splitBam ${bamFile} to fastq" >&2 ; exit 1; fi
       
    else
       echo "quickcheck failed; processing \$samFile to fastqs" >&2
       exit 1
   fi
       
    """
}

process bwaMEM {
    tag "map fastqs to human genome"
    label 'io_mem'

    input:
    tuple val(rg), path(RGtxt), path(fastqs)
    val(analysisType)
    val(sampleType)
    val(genome)

    output:
    tuple val("${rg}"), path("${rg}.bam"), emit: subbam_ch
    
    """
    #!/usr/bin/bash 
    module load biobambam/2.0.87

    . import_config.sh genome ${genome}
    echo \${WG_BWA_DB} 

    if [ -z \${WG_BWA_DB} ]
    then
       >&2 echo "The bwa index directory is not found ..exiting"
       exit 1
    fi

    R1=`ls -p *.fastq.gz|grep 'R1'`
    R2=`ls -p *.fastq.gz|grep 'R2'`
    O1=`ls -p *.fastq.gz|grep 'O1'`
    O2=`ls -p *.fastq.gz|grep 'O2'`

    # Align paired fastqs
    bwa mem -t "$task.cpus" -M \${WG_BWA_DB} \$R1 \$R2 > paired.${rg}.sam 
    if [[ \$? != 0 ]]; then echo "error in mapping paired reads to bam" >&2 ; exit 1; fi

    # Align unpaired fastqs
    zcat \$O1 \$O2 > unpaired.${rg}.fq
    bwa mem -t "$task.cpus" -aM \${WG_BWA_DB} unpaired.${rg}.fq > unpaired.${rg}.sam 
    if [[ \$? != 0 ]]; then echo "error in mapping unpaired reads to bam" >&2 ; exit 1; fi

    # Sort by coord

    samtools sort -o paired.${rg}.bam -@ "$task.cpus" paired.${rg}.sam
    if [[ \$? != 0 ]]; then echo "error in sorting paired sam to bam" >&2 ; exit 1; fi
    samtools sort -o unpaired.${rg}.bam -@ "$task.cpus" unpaired.${rg}.sam
    if [[ \$? != 0 ]]; then echo "error in sorting unpaired sam to bam" >&2 ; exit 1; fi

    nbams_RG=`ls -p *.bam|wc -l`
    if [[ "\$nbams_RG" -gt "2" ]]
    then
       >&2 echo "Expected 2 bams per RG ${rg} got \$nbams_RG ..exiting"
       exit 1
    fi

    # Merge sorted files to give sorted output for the RG

    echo "output bam is ${rg}.bam do not attach RG tag when merging" >&2
    samtools merge -@  "$task.cpus"  -r -c -p -f --output-fmt BAM  ${rg}_noRG.bam paired.${rg}.bam unpaired.${rg}.bam

    # Prep the header
    echo adding RG line to the bam >&2
    samtools view -H ${rg}_noRG.bam > ${rg}.fixedH
    cat "${RGtxt}" >> ${rg}.fixedH 
 
    # Fix RG in the alignment part
    samtools view ${rg}_noRG.bam | sed 's/\\(paired\\|unpaired\\)\\.//g' > ${rg}_wRG.sam

    samtools view -S -b <(cat ${rg}.fixedH ${rg}_wRG.sam) > ${rg}.bam
    

    """
}
