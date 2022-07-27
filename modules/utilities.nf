/*=================== utility scripts ======================*/
params.genome = null
params.locus_band_indel = null
params.locus_band_indel_sv = null
params.locus_band_sv = null
params.locus_band_del = null
params.tartanIndex_projDir = null

process multiLoci2Bed {
    tag "from multiloci (e.g. crest) to bed"
    label 'io_limited'
    publishDir "${params.outD}/tmp/", mode: 'copy'

    input:
    path(simpleCrestFile)

    output:
    path("simpleCrest.bed"), multiLoci2Bed_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0
    python -c "import pandas as pd
import sys

#### a crest file with header that has the following columns
####  disease_sample  chr_a      pos_a  chr_b      pos_b type
#### where type = 'DEL', 'ITX', ..etc

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



def get_interval(row):
    if row['type']=='DEL':
        row['pos_a'] = int(row['pos_a'])-locus_band_del
        row['pos_b'] = int(row['pos_b'])+locus_band_del
        return row.loc[['disease_sample','chr_a','pos_a','pos_b']]
    else:
        s1=row.copy()
        s2=row.copy()
        posA = int(row['pos_a'])
        posB = int(row['pos_b'])
        s1['pos_a'] = posA-locus_band_sv
        s1['pos_b'] = posA+locus_band_sv
        s2['pos_a'] = posB-locus_band_sv
        s2['pos_b'] = posB+locus_band_sv
        s2['chr_a'] = s2['chr_b']
        s1 = s1.loc[['disease_sample','chr_a','pos_a','pos_b']]
        s2 = s2.loc[['disease_sample','chr_a','pos_a','pos_b']]
        return pd.DataFrame([s1,s2])

df_crest = pd.read_csv('${simpleCrestFile}',sep='\\t')

df_crest.loc[:,'chr_a'] = df_crest['chr_a'].apply(lambda r: setChromosomeName(r))

df_crest.loc[:,'chr_b'] = df_crest['chr_b'].apply(lambda r: setChromosomeName(r))

print(df_crest.head())
                                         
temp_df1 = df_crest.apply(lambda r: get_interval(r), axis=1)

temp_df2 = pd.DataFrame()

for i, row in temp_df1.iteritems():
    if pd.api.types.is_string_dtype(row):
        temp_df2 = temp_df2.append(pd.DataFrame(row).T)
    else:
        temp_df2 = temp_df2.append(row)

df_crest = temp_df2

df_crest.columns = ['disease_sample', 'chr', 'posA', 'posB']

df_crest.to_csv('${simpleCrest.bed}', sep='\\t', index=False)

print('Done')

"
        """
}


process generateBedFilesFromExBed {
    tag "generate bed files w/ headers"
    label 'io_limited'
    publishDir "${params.outD}/tmp/", mode: 'copy'

    input:
    path(bedFileWSampleName)

    output:
    tuple path("TRANSCRIPTOME_tumor/*.bed"), path("WHOLE_GENOME_tumor/*.bed"), path("WHOLE_GENOME_germline/*.bed"), path("EXOME_tumor/*.bed"), path("EXOME_germline/*.bed"), emit: ext_bed_files_ch
    
    """
    #!/usr/bin/bash
    module load python/3.7.0
python -c "import numpy as np
import pandas as pd
import glob
import sys
import os
from pathlib import Path


from collections import defaultdict

bedFile = pd.read_csv('${bedFileWSampleName}', sep='\\t')

assert(bedFile.columns.isin(['Sample', 'chr', 'posA', 'posB']).sum() == 4), 'unexpected columns in {}'.format('${bedFileWSampleName}')

# This takes care if the samples provided are SJMB030021_D1, SJMB030021_G1 or just SJMB030021
samples_ = bedFile['Sample'].apply(lambda r:r.split('_')[0]+'_')

bams = defaultdict(set)

for s in samples_:
    bams[s.split('_')[0]] = bams[s.split('_')[0]].union(glob.glob('${params.tartanIndex_projDir}'+'/'+s+'*/*/bam/*.bam'))

print(bams.keys())
final_df = pd.DataFrame()

idx = 0

for i, row in bedFile.iterrows():
    bams_i = list(bams[row['Sample'].split('_')[0]]) #list of 5bams
    print(row['Sample'], bams_i)
    assert(len(bams_i)==5),'not 5 bams!!'
    temp = pd.DataFrame({'bamPath':bams_i})
    temp.loc[:,'chr'] = row['chr']
    temp.loc[:,'posA'] = row['posA']
    temp.loc[:,'posB'] = row['posB']
    final_df = final_df.append(temp)
    idx += 1


assert(bedFile.shape[0]*5 == final_df.shape[0]),'bed events failed to duplicate by 5'

final_df.loc[:,'sample'] = final_df['bamPath'].apply(lambda r:os.path.basename(r).split('.bam')[0])

final_df.loc[:,'target'] = final_df['bamPath'].apply(lambda r:os.path.basename(Path(r).parent.parent).strip())

final_df.loc[:,'sampleType'] = np.where((final_df['sample'].str.contains('_D') | final_df['sample'].str.contains('_A') | final_df['sample'].str.contains('_E') | final_df['sample'].str.contains('_M')),'tumor','germline')

gpby = final_df.groupby(['target','sampleType'], as_index=False)

for i,row in final_df[['target','sampleType']].drop_duplicates().iterrows():
    os.makedirs(row['target']+'_'+row['sampleType'], exist_ok=True)
    dfx = gpby.get_group((row['target'],row['sampleType']))
    #### save a tsv file for each target/sampleType
    dfx.to_csv(row['target']+'_'+row['sampleType']+'/bedFile_wheader.tsv',sep='\\t', index=False)

    
    #### save beds for each sample in target/sampleType group
    gpbyBams = dfx.groupby('bamPath')

    for bamP, dfg in gpbyBams:
        assert(dfg['sample'].drop_duplicates().shape[0]==1),'sample name is not 1-to-1 correspondence with bam path'
        sampleName = dfg['sample'].drop_duplicates().iloc[0]
    
        ### keeping what matters at this point --> different regions for the bam/sample
        bam_bed_s = dfg.loc[:,['sample','bamPath','chr','posA','posB']].drop_duplicates()
        
        ### at least each event must have chr, posA, posB
        assert(bam_bed_s[['chr','posA','posB']].isnull().sum().sum() == 0),'missing regions for bam path {}'.format(bamP)
    
        print(bam_bed_s.head())
        bam_bed_s = bam_bed_s[['chr','posA','posB','sample','bamPath']]
        bam_bed_s = bam_bed_s.convert_dtypes()
        bam_bed_s.to_csv(row['target']+'_'+row['sampleType']+'/'+sampleName+'_wheader.bed', sep='\\t',index=False) 
    

print('Done')

"
        """
}
