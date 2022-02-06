/*=================== workflow to get minibams of localized abnormal events for specific target and sample type======================*/
params.genome = null
params.locus_band_indel = null
params.locus_band_indel_sv = null
params.locus_band_sv = null
params.locus_band_del = null
params.random_seed = 900


process get_clean_manifest {
    tag "generate a manifest for ground truth events"
    label 'io_limited'
    publishDir "${params.outD}/", mode: 'copy'

    input:
    path(manifestFile)

    output:
    path("manifest_all_wRegions.tsv"), emit:manifest_all_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0

echo -e "number of lines in the raw file = \\t" \$(cat "${manifestFile}"|wc -l)
echo "delete empty lines"
cat "${manifestFile}" | sed '/^\$/d' > manifest.tsv
echo -e "number of lines in the raw file after removing empty lines = \\t" \$(cat  manifest.tsv|wc -l)

python -c "import pandas as pd
import numpy as np
import os
import sys


locus_band_indel = int("${params.locus_band_indel}")
locus_band_indel_sv = int("${params.locus_band_indel_sv}")
locus_band_sv = int("${params.locus_band_sv}")
locus_band_del = int("${params.locus_band_del}")


def setChromosomeName(chr, genome='${params.genome}'):
    if genome in ['hg38','GRCh38']:
        chr_ = chr if len(chr.split('chr'))>1 else 'chr'+chr
        if not chr_ in ['chrX','chrY']:
            try:
                chr_ = int(chr_.split('chr')[1])
            except:
                raise  ValueError

    elif genome in ['hg19','GRCh37-lite']:
        chr_ = chr if len(chr.split('chr'))==1 else chr.split('chr')[1]
        if not chr_ in ['X','Y']:
            try:
                chr_ = int(float(chr_))
            except:
                raise  ValueError
    return str(chr_)


def getRegion(row):
    if row[4]=='Somatic SNV' or row[4]=='Germline SNV' or row[4]=='Somatic INDEL' or row[4]=='Germline INDEL':
        try:
            chrA = setChromosomeName(row[9])
            posA = int(row[10]) 
        except ValueError:
            chrA = setChromosomeName(row[10])
            posA = int(row[11]) 
        return [[chrA, posA-locus_band_indel, posA+locus_band_indel ]]

    elif row[4]=='Somatic INDEL/SV':
        chrA = setChromosomeName(row[8])
        posA = int(row[9])
        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        return [[chrA, posA-locus_band_indel_sv, posA+locus_band_indel_sv ] , \
                [chrB, posB-locus_band_indel_sv, posB+locus_band_indel_sv ] ]

    elif row[4]=='Somatic SV / CNA' and row[16] == 'DEL':
        chrA = setChromosomeName(row[8])
        posA = int(row[9])
        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        assert(chrB==chrA),'deletions cannon occur on different chromosomes'
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        return [[chrA, posA-locus_band_del, posB+locus_band_del]]

    elif row[4]=='Somatic SV / CNA' and row[3]=='Deletion or Disruption' and \
           setChromosomeName(row[8]) == setChromosomeName(row[12]): #deletion as well
        chrA = setChromosomeName(row[8])
        posA = int(row[9])
        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        assert(chrB==chrA),'deletions cannon occur on different chromosomes'
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        return [[chrA, posA-locus_band_del, posB+locus_band_del]]

    #elif row[4]=='Somatic SV / CNA' and row[16] == 'INS':
    #    return str(row[8])+':'+str(int(row[9])+'-'+str(int(row[13])

    elif row[4] == 'Somatic SV / CNA': #more like a disruption
        chrA = setChromosomeName(row[8])
        posA = int(row[9])
        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        return [[chrA, posA-locus_band_sv, posA+locus_band_sv ] , \
                [chrB, posB-locus_band_sv, posB+locus_band_sv ] ]

    elif row[4]=='Somatic SV' and row[16] == 'DEL': 
        try:
            chrA = setChromosomeName(row[8])
            posA = int(row[9]) 
        except ValueError:
            chrA = setChromosomeName(row[9])
            posA = int(row[10]) 

        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        assert(chrB==chrA),'deletions cannon occur on different chromosomes'
        return [[chrA, posA-locus_band_del, posB+locus_band_del]]

    elif row[4]=='Somatic SV': 
        try:
            chrA = setChromosomeName(row[8])
            posA = int(row[9]) 
        except ValueError:
            chrA = setChromosomeName(row[9])
            posA = int(row[10]) 

        chrB = setChromosomeName(row[12])
        posB = int(float(row[13]))
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)
        return [[chrA, posA-locus_band_sv, posA+locus_band_sv ] , \
                [chrB, posB-locus_band_sv, posB+locus_band_sv ] ]

    elif (row[4]=='Germline CNA' and row[3]=='Deletion, Intragenic') or (row[4]=='Somatic CNA' and row[3]=='Deletion or Disruption'):
        try:
            #chrom   loc.start       loc.end
            chrA = setChromosomeName(row[7])
            posA = int(row[8]) 
            chrB = chrA
            posB = int(float(row[9]))
        except ValueError:
            #crest like output
            chrA = setChromosomeName(row[8])
            posA = int(row[9]) 
            chrB = setChromosomeName(row[12])
            posB = int(float(row[13]))

        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)

        assert(chrB==chrA),'deletions cannon occur on different chromosomes'
        return [[chrA, posA-locus_band_del, posB+locus_band_del]]

##### fix the raw manifest such that we get the event and the coordinates in the same place regardless to event
dfn = pd.read_csv('manifest.tsv', sep='\t', header=None)
print('here the assumption is that the manifest has 21 columns')  
assert(dfn.columns.shape[0] == 21),'unexpected number of columns in manifest'

assert(dfn.loc[:,19].str.contains('.bam').sum() == dfn.shape[0]),'some entries in the manifest do not have bams in column 19'

print('after this step any null bams would make problems')

assert(dfn[17].isin(['WHOLE_GENOME', 'TRANSCRIPTOME', 'EXOME']).sum() == dfn.shape[0]),'unknown analysis types or manifest is not in expected format'

##### now we fix the manifest to a standard shape: baseSampleName, gene, cytolocus, type, chrA, posA1,posA2, chrB, posB1,posB2, target, sampleName, bam, sampleType
##### note that we add wiggle room to posA and posB when it applies

dfn = dfn[[1,2,3,4,7,8,9,10,11,12,13,16,17,18,19,20]].drop_duplicates()


dfn.index = range(dfn.shape[0])

n_events = dfn.shape[0]

print('number of events after naive deduplication = ',n_events)

##### now we fix the manifest to a standard shape: baseSampleName, gene, cytolocus, type, chrA, posA, chrB, posB, target, sampleName, bam, sampleType
manifest_wregions = pd.DataFrame()

for i,row in dfn.iterrows():
    regions = getRegion(row)
    short_row = row.loc[[1,2,3,4,17,18,19,20]]
    assert(sum([len(kk) for kk in regions])>0),'missing regions for event {}'.format(row)

    print(short_row,short_row.values,row[[4,7,8,9,10,11,12,13,16]])
    for cc in range(short_row.shape[0]):
        manifest_wregions.loc[i,cc] = short_row.values[cc] 

    if len(regions)==1: #one interval

        manifest_wregions.loc[i,'chrA'] = regions[0][0] 
        manifest_wregions.loc[i,'posA1'] = regions[0][1] 
        manifest_wregions.loc[i,'posA2'] = regions[0][2] 
        manifest_wregions.loc[i,'chrB'] = None 
        manifest_wregions.loc[i,'posB1'] = None 
        manifest_wregions.loc[i,'posB2'] = None 

    elif len(regions)==2: #two intervals

        manifest_wregions.loc[i,'chrA'] = regions[0][0] 
        manifest_wregions.loc[i,'posA1'] = regions[0][1] 
        manifest_wregions.loc[i,'posA2'] = regions[0][2] 
        manifest_wregions.loc[i,'chrB'] = regions[1][0] 
        manifest_wregions.loc[i,'posB1'] = regions[1][1] 
        manifest_wregions.loc[i,'posB2'] = regions[1][2] 

    else:
       
        print('unknown output for regions')
        sys.exit(1)


assert(manifest_wregions.shape[0]==n_events),'not all events have regions, some events are missing!!!'

manifest_wregions.columns = ['basename','gene','cytolocus','type','target','sample','bamPath','sampleType','chrA','posA1','posA2','chrB','posB1','posB2']
manifest_wregions.loc[:,'genome'] = '${params.genome}'

print(manifest_wregions.shape)
print(manifest_wregions.head())
print(manifest_wregions.columns)

#### Now after getting all in one row, we do the real deduplication
manifest_wregions = manifest_wregions.drop_duplicates()
print('size of manifest after real deduplication = ',manifest_wregions.shape[0])


manifest_wregions = manifest_wregions.convert_dtypes()

manifest_wregions.to_csv('manifest_all_wRegions.tsv', sep='\\t',index=False)


"
echo "Done"
    """     
}



process validate_manifest {

    tag "validates the manifest"
    label 'io_limited'
    publishDir "${params.outD}/", mode: 'copy'

    input:
    path(manifestWRegionsFile)

    output:
    path("manifest_validation.txt"), emit:manifest_all_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0

python -c "import pandas as pd
import numpy as np
import os
import sys

manifest_wregions = pd.read_csv('${manifestWRegionsFile}', sep='\\t')

#### check that we got the same events across target/sampleTypes

gpbyAnlyTypeSampleType = manifest_wregions.groupby(['target','sampleType'], sort=True)

number_of_uniq_events = []

for g, dfg in gpbyAnlyTypeSampleType:

   assert(dfg['sample'].unique().shape[0]==dfg['bamPath'].unique().shape[0]),'the number of bams are different than the number of samples for group {}'.format(g) 
   number_of_uniq_events.append(dfg.drop_duplicates().shape[0])
   analysisType = g[0]
   sampleType = g[1]
   print('saving a file ','manifest_'+analysisType+'_'+sampleType+'.tsv',' this is really for checking')
   dfg = dfg.drop_duplicates().convert_dtypes()
   dfg.to_csv('manifest_'+analysisType+'_'+sampleType+'.tsv',sep='\\t', index=False)


#assert(number_of_uniq_events.count(number_of_uniq_events[0]) == len(number_of_uniq_events)),'the events are not the same across target/sampleType'

df_valid = pd.DataFrame()

    
if number_of_uniq_events.count(number_of_uniq_events[0]) == len(number_of_uniq_events):

    df_valid.loc['VALID1','STATUS'] = 'PASS'
    df_valid.loc['VALID1','NOTES'] = 'The manifest passed the first validation: the events are the same across target/sample-type combinations'

else:

    df_valid.loc['VALID1','STATUS'] = 'FAIL'
    df_valid.loc['VALID1','NOTES'] = 'The manifest failed the first validation: the events are not the same across target/sample-type combinations'


#Now we do the coverage validation:

hg19_chr=list(map(str,list(range(1,23))))
hg38_chr = ['chr'+i for i in hg19_chr]

missing_chr = []

if '${params.genome}' == 'hg19':

    for chrom in hg19_chr:
        if manifest_wregions['chrA'].isin([chrom]).sum()==0 and manifest_wregions['chrB'].isin([chrom]).sum()==0:
            missing_chr.append(chrom)

elif '${params.genome}' == 'hg38':

    for chrom in hg38_chr:
        if manifest_wregions['chrA'].isin([chrom]).sum()==0 and manifest_wregions['chrB'].isin([chrom]).sum()==0:
            missing_chr.append(chrom)

else:
    print('unknown genome')
    sys.exit(1)
    

if len(missing_chr) == 0:

    df_valid.loc['VALID2','STATUS'] = 'PASS'
    df_valid.loc['VALID2','NOTES'] = 'The manifest passed the chromosome voverage test'

else:

    df_valid.loc['VALID2','STATUS'] = 'FAIL'
    df_valid.loc['VALID2','NOTES'] = 'The manifest failed the second validation: there are chromosomes with 0 coverage this will break the pipeline, please fill the following chromosomes with some events: {}'.format(missing_chr)

df_valid.to_csv('manifest_validation.txt', sep='\\t',index=False)

if df_valid['STATUS'].isin(['FAIL']).sum() > 0 :
    print('Manifest failed validation, see below for more information')
    print(df_valid)
    sys.exit(1)

print('Done')



"

"""
}




process get_GT_manifest {
    tag "generate a manifest for ground truth events for each anlysis/sample types"
    label 'io_limited'
    publishDir "${params.outD}/", mode: 'copy'

    input:
    path(manifestCleanFile)
    val(analysisType)
    val(sampleType)

    output:
    path ("manifest_${analysisType}_${sampleType}.tsv"), emit:minibam_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0

echo -e "number of lines in the raw file = \\t" \$(cat "${manifestCleanFile}"|wc -l)

python -c "import pandas as pd
import numpy as np
import sys

##### fix the raw manifest
dfn = pd.read_csv('${manifestCleanFile}', sep='\\t')

cols = ['basename','gene','cytolocus','type','target','sample','bamPath','sampleType','chrA','posA1','posA2','chrB','posB1','posB2','genome']

assert(dfn.columns.isin(cols).sum() == len(cols)),'error reading the file ${manifestCleanFile}'

#### eventually all samples all genes all loci under the same analysis type (WGS/WES/RNA) and the same sample type (tumor/germline) will be merged

gpbyAnlyTypeSampleType = dfn.groupby(['target','sampleType'], sort=True)
print(gpbyAnlyTypeSampleType.groups.keys())

#### this is what we want

dfx = gpbyAnlyTypeSampleType.get_group(('${analysisType}', '${sampleType}')).drop_duplicates()
assert(dfx['sample'].unique().shape[0]==dfx['bamPath'].unique().shape[0]),'the number of bams are different than the number of samples'

print('saving a file ','manifest_'+'${analysisType}'+'_'+'${sampleType}'+'.tsv')

dfx = dfx.convert_dtypes()
dfx.to_csv('manifest_'+'${analysisType}'+'_'+'${sampleType}'+'.tsv',sep='\\t', index=False)

print(dfx)


"
echo "Done"
    """     
}


process generateBEDFiles {

    tag "generate bed files"
    label 'io_limited'
    
    input:
    path(manifestFile)
    val(analysisType)
    val(sampleType)

    output:
    path ("bam_beds/*.bed"), emit:bed_files_ch

    """
    #!/usr/bin/bash
    module load python/3.7.0
    mkdir bam_beds

    python -c "import pandas as pd
import os
import sys
import numpy as np


dfx = pd.read_csv('${manifestFile}', sep = '\\t')
print(dfx.head(3))
#### group by bams
print(dfx.columns)
#### columns ['basename','gene','cytolocus','type','target','sample','bamPath','sampleType','chrA','posA1','posA2','chrB','posB1','posB2','genome']
assert(dfx.columns.shape[0] == 15),'unexpected number of columns in manifest'
gpbyBams = dfx.groupby('bamPath')


for bamP, dfg in gpbyBams:
    assert(dfg['sample'].drop_duplicates().shape[0]==1),'sample name is not 1-to-1 correspondence with bam path'
    sampleName = dfg['sample'].drop_duplicates().iloc[0]
    bamFilePath = bamP+'/'+sampleName+'.bam'
    print('sample:',sampleName,'\\t','bamP:',bamFilePath)

    ### keeping what matters at this point --> different regions for the bam/sample
    temp = dfg.loc[:,['sample','bamPath','chrA','posA1','posA2','chrB','posB1','posB2']].drop_duplicates()
    
    ### at least each event must have chrA, posA1, posA2
    assert(temp[['chrA','posA1','posA2']].isnull().sum().sum() == 0),'missing regions for bam path {}'.format(bamP)

    #### generate bed file for the current bam
    colsA = {'chrA':'chr','posA1':'pos1','posA2':'pos2'}
    colsB= {'chrB':'chr','posB1':'pos1','posB2':'pos2'}
    bam_regions = [temp[['chrA','posA1','posA2']].convert_dtypes().rename(columns=colsA), temp[['chrB','posB1','posB2']].dropna().convert_dtypes().rename(columns=colsB)]
    print(bam_regions)
    bam_bed = pd.concat(bam_regions, axis=0)
    print(bam_bed.head())
    bam_bed.loc[:,'sample'] = sampleName
    bam_bed.loc[:,'bamPath'] = bamFilePath
    bam_bed = bam_bed.drop_duplicates().dropna()
    bam_bed.columns = ['chr','pos1','pos2','sample','bamPath']
    bam_bed = bam_bed.convert_dtypes()
    bam_bed.to_csv('bam_beds/'+sampleName+'_wheader.bed', sep='\\t',index=False) 
    

print('Done')

"

    """
}



process generateSubBams {
    tag "small bams from $analysisType $sampleType"
    //container "stjudecloud/samtools:branch-chipseq-1.0.2"
    label 'multithread'
    //publishDir "${params.outD}/", mode: 'copy'
    
    // do not run if results exist
    // storeDir "${params.outD}/runDir"



    input:
    path(bedFile)
    val(analysisType)
    val(sampleType)

    output:
    //path ("bamsToCombine/*.bam"), emit:minibamToCombine_ch
    path ("*.bam"), emit:minibamToCombine_ch
    
    """

    module load samtools/1.12
    
    function getSmallBam() {

        outputBam=\$(echo \$(basename \$1 _wheader.bed).bam)
        inputBam=\$(cat \$1|sed 1d|cut -f 5|sort|uniq)
        inputSam=\$(echo "\$(basename \$inputBam .bam)".sam)

        #keep going but save exit code
        EXIT_CODE=0

        samtools quickcheck \$inputBam || EXIT_CODE=\$?


        if [[ \$EXIT_CODE -eq 0 ]]; then
           echo "processing \$inputBam to \$outputBam"
           samtools view -@ "$task.cpus" -L <(cat \$1 |sed 1d|cut -f 1-3)  -o \$outputBam \$inputBam
        else
           echo "processing \$inputSam to \$outputBam"
           samtools view -h \$inputBam > \$inputSam || true 
           samtools view -@ "$task.cpus"  -b -L <(cat \$1 |sed 1d|cut -f 1-3) -o \$outputBam  \$inputSam
        fi
    }

    export -f getSmallBam
        
        
    getSmallBam "${bedFile}"


    """
}



process mergeToMinibams {
    tag "merge to minibams"
    //container "stjudecloud/samtools:branch-chipseq-1.0.2"
    label 'multithread'
    publishDir "${params.outD}/minibams", mode: 'copy'

    input:
    path ("*.bam")
    val(analysisType)
    val(sampleType)

    output:
    //path ("${analysisType}_${sampleType}_${params.genome}.*"), emit:minibam_ch
    path ("${analysisType}_${params.genome}/*"), emit:minibam_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0
    module load samtools/1.12

    mkdir -p ${analysisType}_${params.genome}

    find *.bam > bamlistfile 

    echo "merge small bams into a minibam"
    samtools merge -@ "$task.cpus" -r -c -p -f --output-fmt BAM -b bamlistfile ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam


    echo "changing SM in the header to be unique single for a single minibam file"

    #### this produces the same number 92966; change the seed to get a different number
    random6Digit=\$(echo 9\$(python -c "import numpy as np; np.random.seed (int('${params.random_seed}')); print(np.random.randint(90000,99999,1)[0])"))
   
    if [ "${sampleType}" == "tumor" ]; then

        echo "sample type is tumor"
       
        SM="SJOther"\$random6Digit"_D1"

        echo "\$SM"

    elif [ "${sampleType}" == "germline" ]; then
  
        echo "sample type is germline"
        
        SM="SJOther"\$random6Digit"_G1"
  
        echo "\$SM"

    else

        echo "unknown sample type"

    fi

    
    #### this could produce duplicate lines, no need to sort to get unique lines
    samtools view -H ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam | sed -r '/^@RG/s/SM:(\\S*)(\\s*)/SM:'\$SM'\\2/p' | uniq > ${analysisType}_${sampleType}_${params.genome}.notsorted.fixedH 

    samtools reheader ${analysisType}_${sampleType}_${params.genome}.notsorted.fixedH ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam > ${analysisType}_${sampleType}_${params.genome}.notsorted.bam 
    
    echo "sort by coordinates"
    samtools sort -@ "$task.cpus" ${analysisType}_${sampleType}_${params.genome}.notsorted.bam -o ${analysisType}_${params.genome}/"\$SM".bam 

    cd ${analysisType}_${params.genome}/

    echo "index to ready for testing"
    samtools index -@ "$task.cpus" "\$SM".bam


    echo "done"

    
    """
}
