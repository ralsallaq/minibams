/*=================== workflow to get minibams of localized abnormal events for specific target and sample type======================*/
/* 
   Inputs:
          genome: from which to determine chr and tag results
          locus_band_indel: wiggle room for indels loci
          locus_band_indel_sv: wiggle room for SV indels loci
          locus_band_sv: wiggle room for SV loci 
          locus_band_del: wiggle room before and after deletions
   Outputs:
          varies by process
*/

params.outD = null
params.genome = null
params.locus_band_indel = null
params.locus_band_indel_sv = null
params.locus_band_sv = null
params.locus_band_del = null
params.random_seed = null

process get_clean_manifest {
    tag "generate a manifest for ground truth events"
    cpus 1
    memory 5.GB
    publishDir "${params.outD}/", mode: 'copy'

    input:
    path(manifestFile)

    output:
    tuple path("manifest_all_wRegions.tsv"), path("bed_like_events.tsv"), path("bedpe_like_events.tsv"),
          path("bed_events_wheader.bed"), path("bedpe_events_wheader.bedpe"), emit:manifest_all_ch
    
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
locus_band_sv = int("${params.locus_band_sv}")
locus_band_del = int("${params.locus_band_del}")


def setChromosomeName(chr, genome='${params.genome}'):
#       Takes chr (e.g. chrX or X) and genome
#       name (hg38 or GRCh38, hg19 or GRCh37-lite)
#       and returns the correct format 
#       for the chromosome (X for hg19 and chrX for hg38)
#       Args:
#            chr either chrX or X
#            genome hg19 (or GRCh37-lite) or hg38 (or GRCh38)
#      Output:
#            string (e.g. chrX for hg38)
 
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


def reformatRaw(rawD):
#        Take a data frame for the 
#        raw collected events
#        reformats the raw dataframe into two
#        one is bed-like for events with one
#        loci and the other is bedpe-like
#        for events with two loci
#        Args:
#            rawD: a data-frame of events as collected 
#        Output: 
#            bedpe_like data-frame defining two-loci events
#            bed_like data-frame defining one-loci events
#            bedpe a BEDPE file
#            bed a BED file

    clskeep = rawD.columns[[1, 2, 3, 4]+list(range(7,rawD.shape[1]))]
    # We first split to SVs and non-SVs
    # Germline CNA and Somatic CNA should be SVs
    nonSVs = rawD[(~rawD[4].str.contains('SV')) & (~rawD[4].str.contains('CNA'))].iloc[:,clskeep]
    SVs = rawD[rawD[4].str.contains('SV') | rawD[4].str.contains('CNA')].iloc[:,clskeep]

    #format SVs to bedpe format and nonSVs to bed format
    str_chr = [str(k) for k in range(1,23)] + ['X', 'Y']
    def classify_svs(row):
        if row[4] == 'Somatic SV / CNA':
            #return [8, 9, 10, 12, 13, 14, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]
            return 't1' 
        if row[4] == 'Somatic SV':
            # This type has NaN pair and score
            if row[7].split('chr')[-1] in str_chr: 
                #return [7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 1, 2, 3, 4]
                return 't0'
            elif row[8].split('chr')[-1] in str_chr: 
                #return [8, 9, 10, 12, 13, 14, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]
                return 't1' 
            else:
                #return [9, 10, 11, 12, 13, 14, 8, 15, 7, 17, 18, 19, 20, 1, 2, 3, 4]
                return 't2' 
        if row[4] == 'Somatic INDEL/SV':
            #return [8, 9, 10, 12, 13, 14, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4] 
            return 't1' 
        if row[4] in ['Somatic CNA', 'Germline CNA']: 
            if row[7].split('chr')[-1] in str_chr:
                #return [7, 8, 11, 7, 9, 11, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]
                return 't3' 
            else:
                #return [8, 9, 10, 12, 13, 14, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]
                return 't1' 

    if SVs.shape[0] > 0:
        SVs.loc[:,'type'] = SVs.apply(classify_svs, axis=1)
        # Three types of rows in the SVs that are easy to convert to bedpe:
        # We need to have something with at least the columns:
        cols=['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','score']
    
        bedpe_like0 = SVs.loc[SVs['type']=='t0', [7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 19, 20, 1, 2, 3, 4]]
        bedpe_like0.columns = ['ChrA_', 'PosA', 'OrtA', 'ChrB_', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                          'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event']
        bedpe_like1 = SVs.loc[SVs['type']=='t1', [8, 9, 10, 12, 13, 14, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]]
        bedpe_like1.columns = ['ChrA_', 'PosA', 'OrtA', 'ChrB_', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                          'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event']
    
        bedpe_like2 = SVs.loc[SVs['type']=='t2', [9, 10, 11, 12, 13, 14, 8, 15, 7, 17, 18, 19, 20, 1, 2, 3, 4]]
        bedpe_like2.columns = ['ChrA_', 'PosA', 'OrtA', 'ChrB_', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                          'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event']
    
        bedpe_like3 = SVs.loc[SVs['type']=='t3', [7, 8, 11, 7, 9, 11, 16, 11, 7, 17, 18, 19, 20, 1, 2, 3, 4]]
        bedpe_like3.columns = ['ChrA_', 'PosA', 'OrtA', 'ChrB_', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                          'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event']
    
    
        bedpe_like = bedpe_like0.append(bedpe_like1)
        bedpe_like = bedpe_like.append(bedpe_like1)
        bedpe_like = bedpe_like.append(bedpe_like2)
        bedpe_like = bedpe_like.append(bedpe_like3)
    
        # Drop those with Null type as they are usually duplicates
        # This also means do not supply events with null Type
        test_numeric = pd.to_numeric(bedpe_like['Type'], errors='coerce')
        bedpe_like = bedpe_like.loc[test_numeric.isnull()]
        bedpe_like.loc[:,'PosA'] = bedpe_like['PosA'].astype(float).astype(int)
        bedpe_like.loc[:,'PosB'] = bedpe_like['PosB'].astype(float).astype(int)
    else:
        bedpe_like = pd.DataFrame(columns=['ChrA', 'PosA', 'OrtA', 'ChrB', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                               'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event', 'ChrA_', 'ChrB_'])

    # Next are the nonSVs
    def classify_nonSvs(row):
        if row[10].startswith('chr'):
            return 't1'
        else:
            return 't2'

    if nonSVs.shape[0] > 0:
        nonSVs.loc[:,'type'] = nonSVs.apply(classify_nonSvs, axis=1)
        bed_like1 = nonSVs.loc[nonSVs['type']=='t1', [10, 11, 7, 12, 9,  17, 18, 19, 20, 1, 2, 3, 4]]
        bed_like1.columns = ['ChrA_', 'PosA', 'gene', 'type', 'pair', 'target', 'sample', 'bam', 'sampleType', \
                        'sname', 'gene', 'cytolocus', 'event']
    
        bed_like2 = nonSVs.loc[nonSVs['type']=='t2', [9, 10, 7, 11, 8, 17, 18, 19, 20, 1, 2, 3, 4]]
        bed_like2.columns = ['ChrA_', 'PosA', 'gene', 'type', 'pair', 'target', 'sample', 'bam', 'sampleType', \
                        'sname', 'gene', 'cytolocus', 'event']
    
        bed_like = bed_like1.append(bed_like2)
        bed_like.loc[:,'PosA'] = bed_like['PosA'].astype(float).astype(int)
    else:
        bed_like = pd.DataFrame(columns=['ChrA', 'PosA', 'type', 'pair', 'target', 'sample', 'bam', 'sampleType', \
                                 'sname', 'gene', 'cytolocus', 'event', 'ChrA_'])



    bedpe_like.loc[:,'ChrA'] = bedpe_like['ChrA_'].apply(lambda r: setChromosomeName(r))
    bedpe_like.loc[:,'ChrB'] = bedpe_like['ChrB_'].apply(lambda r: setChromosomeName(r))
    bed_like.loc[:,'ChrA'] = bed_like['ChrA_'].apply(lambda r: setChromosomeName(r))

    bedpe_like.drop(['ChrA_', 'ChrB_'], axis=1, inplace=True)
    bed_like.drop('ChrA_', axis=1, inplace=True)

    bedpe_like = bedpe_like[['ChrA', 'PosA', 'OrtA', 'ChrB', 'PosB', 'OrtB', 'Type', 'score', 'pair', 'target', \
                   'sample', 'bam', 'sampleType', 'sname', 'gene', 'cytolocus', 'event']].drop_duplicates()
    bed_like = bed_like[['ChrA', 'PosA', 'type', 'pair', 'target', 'sample', 'bam', 'sampleType', \
                   'sname', 'gene', 'cytolocus', 'event']].drop_duplicates()

    bedpe_like = bedpe_like.reset_index().drop('index', axis=1)
    bed_like = bed_like.reset_index().drop('index', axis=1)

    ## get bedpe and bed format to use with bedtools
    bedpe = bedpe_like.copy()
    bedpe.loc[:,'PosA1'] = bedpe['PosA'].apply(lambda r: int(r)-1 if pd.notna(r) else r).values
    bedpe.loc[:,'PosB1'] = bedpe['PosB'].apply(lambda r: int(r)-1 if pd.notna(r) else r).values

    bedpe = bedpe[['ChrA', 'PosA1', 'PosA', 'ChrB', 'PosB1', 'PosB', 'Type', 'score', 'OrtA', 'OrtB', \
                   'gene', 'cytolocus', 'event']].drop_duplicates()
  
    bed = bed_like.copy()
    bed.loc[:, 'PosA1'] = bed['PosA'].apply(lambda r: int(r)-1 if pd.notna(r) else r).values

    bed = bed[['ChrA', 'PosA1', 'PosA', 'type', \
               'gene', 'cytolocus', 'event']].drop_duplicates()

    return bedpe_like, bed_like, bedpe, bed



def getRegion(row, dftype='bedpe'):
#        Given a row from a data frame which
#        is either a bed-like or a bedpe-like
#        data frame, this function defines the
#        regions for the row using the characteristics
#        of the event
#        the output is a list of lists encompassing:
#        [[chrA, posA1, posA2], [chrB, posB1, posB2]] or
#        [[chrA, posA1, posA2]] 
# 

    if dftype == 'bed':
        if row['event']=='Somatic SNV' or row['event']=='Germline SNV' or row['event']=='Somatic INDEL' or row['event']=='Germline INDEL':
            return [[ row['ChrA'], int(row['PosA'])-locus_band_indel, int(row['PosA'])+locus_band_indel ]]
        
    elif dftype == 'bedpe':

        chrA = row['ChrA']; posA = int(row['PosA']); chrB = row['ChrB']; posB = int(row['PosB'])
        posA,posB = (posA,posB) if posA<posB else (posB,posA)
        chrA,chrB = (chrA,chrB) if posA<posB else (chrB,chrA)

        if (row['event'] == 'Germline CNA' and row['cytolocus'] == 'Deletion, Intragenic') or (row['event'] == 'Somatic CNA' and row['cytolocus'] == 'Deletion or Disruption'):
             assert(chrB==chrA),'deletions cannon occur on different chromosomes'
             return [[chrA, posA-locus_band_del, posB+locus_band_del]]

        elif row['event'] == 'Somatic INDEL/SV':
             # This is either a deletion or insertion
             return [[chrA, posA-locus_band_del, posB+locus_band_del]]
            
        elif row['event'] == 'Somatic SV / CNA' and row['Type'] == 'DEL': 
            assert(chrB==chrA),'deletions cannon occur on different chromosomes'
            return [[chrA, posA-locus_band_del, posB+locus_band_del]]
       
        elif row['event'] == 'Somatic SV / CNA' and row['cytolocus'] == 'Deletion or Disruption' and \
                      chrB==chrA: #deletion
            return [[chrA, posA-locus_band_del, posB+locus_band_del]] 
       
        elif row['event'] == 'Somatic SV / CNA': #more like a disruption
            return [[chrA, posA-locus_band_sv, posA+locus_band_sv ] , \
                   [chrB, posB-locus_band_sv, posB+locus_band_sv ]]

        elif row['event'] == 'Somatic SV' and not pd.isna(row['Type']) and row['Type'] == 'DEL':
            assert(chrB==chrA),'deletions cannot occur on different chromosomes'
            return [[chrA, posA-locus_band_del, posB+locus_band_del]]

        elif row['event'] == 'Somatic SV' and not pd.isna(row['Type']) and row['Type'] == 'INS':
            assert(chrB==chrA),'insertions might not occur on different chromosomes'
            return [[chrA, posA-locus_band_del, posB+locus_band_del]]

        elif  row['event'] == 'Somatic SV':
            # SV event that manifest on two loci
            return [[chrA, posA-locus_band_sv, posA+locus_band_sv ] , \
                   [chrB, posB-locus_band_sv, posB+locus_band_sv ] ]

           

# Fix the raw manifest such that we get the event and the coordinates in the same place regardless to event
dfn = pd.read_csv('manifest.tsv', sep='\\t', header=None)
bedpe_like, bed_like, bedpe, bed = reformatRaw(dfn)

# Save to be used with liftover when needed
bedpe_like.to_csv('bedpe_like_events.tsv', sep='\\t', index=False)
bed_like.to_csv('bed_like_events.tsv', sep='\\t', index=False)

# Save to use with bedtools when needed
bedpe.to_csv('bedpe_events_wheader.bedpe', sep='\\t', index=False)
bed.to_csv('bed_events_wheader.bed', sep='\\t', index=False)

assert(bed_like.loc[:,'bam'].str.contains('.bam').sum() == bed_like.shape[0]),'some entries in the bed_like file do not have bams'
assert(bedpe_like.loc[:,'bam'].str.contains('.bam').sum() == bedpe_like.shape[0]),'some entries in the bedpe_like file do not have bams'

print('after this step any null bams would make problems', file=sys.stderr)

# Now we fix the manifest to a standard shape: baseSampleName, gene, cytolocus, type, chrA, posA1,posA2, chrB, posB1,posB2, target, sampleName, bam, sampleType
# Note that we add wiggle room to posA and posB when it applies

n_events = bedpe_like.shape[0] + bed_like.shape[0]

print('number of events after deduplication = ', n_events, file=sys.stderr)

# Now we fix the manifest to a standard shape: baseSampleName, gene, cytolocus, type, chrA, posA, chrB, posB, target, sampleName, bam, sampleType
manifest_wregions = pd.DataFrame()

for i,row in bedpe_like.iterrows():
    # If no enough information skip
    if (~row.loc[['ChrA', 'PosA', 'ChrB', 'PosB', 'Type', 'cytolocus']].isnull()).sum() < 6:
        print('skipping event {}'.format(row), file = sys.stderr)
        continue

    regions = getRegion(row, dftype='bedpe')
    assert(sum([len(kk) for kk in regions])>0),'missing regions for event {}'.format(row)
    short_row = row.loc[['sname', 'gene', 'cytolocus', 'event', 'target', 'sample', 'bam', 'sampleType']]

    print(short_row, short_row.values, file=sys.stderr)
    for cc in range(short_row.shape[0]):
        manifest_wregions.loc[i,short_row.index[cc]] = short_row.values[cc] 

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
       
        print('unknown output for regions', file=sys.stderr)
        sys.exit(1)

for j, row in bed_like.iterrows():
    # If no enough information skip
    if (~row.loc[['ChrA', 'PosA', 'type', 'cytolocus']].isnull()).sum() < 4:
        print('skipping event {}'.format(row), file = sys.stderr)
        continue

    regions = getRegion(row, dftype='bed')
    assert(sum([len(kk) for kk in regions])>0),'missing regions for event {}'.format(row)
    short_row = row.loc[['sname', 'gene', 'cytolocus', 'event', 'target', 'sample', 'bam', 'sampleType']]

    print(short_row, short_row.values, file=sys.stderr)
    for cc in range(short_row.shape[0]):
        manifest_wregions.loc[i+j+1,short_row.index[cc]] = short_row.values[cc] 

    if len(regions)==1: #one interval

        manifest_wregions.loc[i+j+1,'chrA'] = regions[0][0] 
        manifest_wregions.loc[i+j+1,'posA1'] = regions[0][1] 
        manifest_wregions.loc[i+j+1,'posA2'] = regions[0][2] 
        manifest_wregions.loc[i+j+1,'chrB'] = None 
        manifest_wregions.loc[i+j+1,'posB1'] = None 
        manifest_wregions.loc[i+j+1,'posB2'] = None 

    elif len(regions)==2: #two intervals

        manifest_wregions.loc[i+j+1,'chrA'] = regions[0][0] 
        manifest_wregions.loc[i+j+1,'posA1'] = regions[0][1] 
        manifest_wregions.loc[i+j+1,'posA2'] = regions[0][2] 
        manifest_wregions.loc[i+j+1,'chrB'] = regions[1][0] 
        manifest_wregions.loc[i+j+1,'posB1'] = regions[1][1] 
        manifest_wregions.loc[i+j+1,'posB2'] = regions[1][2] 

    else:
       
        print('unknown output for regions', file=sys.stderr)
        sys.exit(1)


manifest_wregions.columns = ['basename','gene','cytolocus','type','target','sample','bamPath','sampleType','chrA','posA1','posA2','chrB','posB1','posB2']
manifest_wregions.loc[:,'genome'] = '${params.genome}'

print(manifest_wregions.shape, file=sys.stderr)
print(manifest_wregions.head(), file=sys.stderr)
print(manifest_wregions.columns, file=sys.stderr)

#### Now after getting all in one row, we do the real deduplication
manifest_wregions = manifest_wregions.drop_duplicates()
print('size of manifest after real deduplication = ',manifest_wregions.shape[0], file=sys.stderr)


manifest_wregions = manifest_wregions.convert_dtypes()

manifest_wregions.to_csv('manifest_all_wRegions.tsv', sep='\\t',index=False)


"
echo "Done"
    """     
}



process validate_manifest {

    tag "validates the manifest"
    cpus 1
    memory 5.GB
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
   print('saving a file ','manifest_'+analysisType+'_'+sampleType+'.tsv',' this is really for checking', file=sys.stderr)
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


# Now we do the coverage validation:

hg19_chr=list(map(str,list(range(1,23)))) + ['MT','Y','X']
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
    print('unknown genome', file=sys.stderr)
    sys.exit(1)
    

if len(missing_chr) == 0:

    df_valid.loc['VALID2','STATUS'] = 'PASS'
    df_valid.loc['VALID2','NOTES'] = 'The manifest passed the chromosome voverage test'

else:

    df_valid.loc['VALID2','STATUS'] = 'FAIL'
    df_valid.loc['VALID2','NOTES'] = 'The manifest failed the second validation: there are chromosomes with 0 coverage this will break the pipeline, please fill the following chromosomes with some events: {}'.format(missing_chr)

df_valid.to_csv('manifest_validation.txt', sep='\\t',index=False)

if df_valid['STATUS'].isin(['FAIL']).sum() > 0 :
    print('Manifest failed validation, see below for more information', file=sys.stderr)
    print(df_valid, file=sys.stderr)
    with open('exit_status.txt', 'w') as f:
        f.write('1')
else:
    with open('exit_status.txt', 'w') as f:
        f.write('0')

"
exit_status=\$(echo exit_status.txt)
echo \$exit_status >&2

"""
}




process get_GT_manifest {
    tag "generate a manifest for ground truth events for each anlysis/sample types"
    cpus 1
    memory 5.GB
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

echo -e "number of lines in the raw file with regions = \\t" \$(cat manifest_all_wRegions.tsv|wc -l)

python -c "import pandas as pd
import numpy as np
import sys

##### fix the raw manifest
dfn = pd.read_csv('manifest_all_wRegions.tsv', sep='\\t')

cols = ['basename','gene','cytolocus','type','target','sample','bamPath','sampleType','chrA','posA1','posA2','chrB','posB1','posB2','genome']

assert(dfn.columns.isin(cols).sum() == len(cols)),'error reading the file manifest_all_wRegions.tsv'

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
    cpus 1
    memory 5.GB
    
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
    # Convert it to 0-based as it should be for a bed file expected by samtools
    bam_bed = bam_bed.convert_dtypes()
    bam_bed.loc[:,'pos1'] = bam_bed['pos1'].values - 1
    bam_bed.to_csv('bam_beds/'+sampleName+'_wheader.bed', sep='\\t',index=False) 
    

print('Done')

"

    """
}



process generateSubBams {
    tag "small bams from $analysisType $sampleType"
    //container "stjudecloud/samtools:branch-chipseq-1.0.2"
    cpus 8
    memory 80.GB

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

    module load samtools/1.15.1
    
    function getSmallBam() {

        outputBam="\$(basename \$1 _wheader.bed)".bam
        inputBed="\$(basename \$outputBam .bam)".bed
        # Sort as it seems to help speed up samtools view
        cat \$1 | sed 1d| sort -k1,1V -k2,2n > bed_noheader_sorted.txt
        inputBam=\$(cat \$1|sed 1d|cut -f 5|sort|uniq)
        inputSam="\$(basename \$inputBam .bam)".sam
        inputBam2="\$(basename \$inputBam .bam)"_orig.bam
        cat bed_noheader_sorted.txt| cut -f 1-3  > \$inputBed

        #keep going but save exit code
        EXIT_CODE=0

        samtools quickcheck \$inputBam || EXIT_CODE=\$?


        if [[ \$EXIT_CODE -eq 0 ]]; then
           echo "processing \$inputBam to \$outputBam"
           # This would not extract mates if they are mapped outside the region
           # and will end up with singleton reads (mapped reads with no mapped mates)
           # samtools view -@ "$task.cpus" -L \$inputBed  -o \$outputBam \$inputBam
           #java.sh -XX:ParallelGCThreads=$task.cpus -jar \$(which samviewwithmate.jar) --bed \$inputBed --samoutputformat BAM -o \$outputBam \$inputBam
           # Reads in the region and their mates wherever they are (note the -P option available since 1.15)
           samtools view -@ "$task.cpus" -P -L \$inputBed  -o \$outputBam \$inputBam 
 
        else
           echo "processing \$inputSam to \$outputBam"
           samtools view -h \$inputBam > \$inputSam || true 
           # This would not extract mates if they are mapped outside the region
           # and will end up with singleton reads (mapped reads with no mapped mates)
           # samtools view -@ "$task.cpus"  -b -L \$inputBed  -o \$outputBam  \$inputSam
           # This will output the reads within the regions and their mates wherever they are
           # java.sh -XX:ParallelGCThreads=$task.cpus -jar \$(which samviewwithmate.jar) --bed \$inputBed  --samoutputformat BAM -o \$outputBam \$inputSam
           # sort and index the sam file
           samtools sort -@ "$task.cpus" \$inputSam -o \$inputBam2
           samtools index \$inputBam2
           # Reads in the region and their mates wherever they are (note the -P option available since 1.15)
           samtools view -@ "$task.cpus" -P -L \$inputBed  -o \$outputBam \$inputBam2
        fi
    }

    export -f getSmallBam
        
        
    getSmallBam "${bedFile}"


    """
}



process mergeToMinibams {
    tag "merge to minibams"
    //container "stjudecloud/samtools:branch-chipseq-1.0.2"
    cpus 8
    memory 80.GB
    publishDir "${params.outD}", mode: 'copy'

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
    module load samtools/1.15.1

    mkdir -p ${analysisType}_${params.genome}
   
    echo "firt split each bam by read group" >&2
    mkdir -p bams_split_by_rgs

    find *.bam > bams_split_by_rgs/bamlistfile_raw 

    cd bams_split_by_rgs

    for bam in \$(cat bamlistfile_raw | xargs)
    do
      samtools split -@ "$task.cpus" --no-PG ../\${bam}
    done

    echo "merging \$(ls *.bam|wc -l) bams" >&2
  
    find *.bam > bamlistfile 
    # Sort the splits by coordinates for merging
    for bam in \$(cat bamlistfile|xargs)
    do
      samtools sort -@ "$task.cpus" \$bam -o \${bam}.sorted
      if [ "\$?" == "0" ]; then
          mv \${bam}.sorted \$bam
      fi   
    done

    echo "merge read groups sorted bams into a minibam" >&2
    samtools merge -@ "$task.cpus" -r -c -p -f --no-PG --output-fmt BAM -b bamlistfile  ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam

    echo "changing SM in the header to be unique single for a single minibam file" >&2

    # This produces the same number 92966; change the seed to get a different number
    random6Digit=\$(echo 9\$(python -c "import numpy as np; np.random.seed (int('${params.random_seed}')); print(np.random.randint(90000,99999,1)[0])"))
   
    if [ "${sampleType}" == "tumor" ]; then

        >&2 echo "sample type is tumor"
       
        SM="SJOther"\$random6Digit"_D1"

        >&2 echo "\$SM"

    elif [ "${sampleType}" == "germline" ]; then
  
        >&2 echo "sample type is germline"
        
        SM="SJOther"\$random6Digit"_G1"
  
        >&2 echo "\$SM"

    else

        >&2 echo "unknown sample type"

    fi

    
    #### this could produce duplicate lines, no need to sort to get unique lines
    samtools view -H ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam | sed -r '/^@RG/s/SM:(\\S*)(\\s*)/SM:'\$SM'\\2/p' | uniq > ${analysisType}_${sampleType}_${params.genome}.notsorted.fixedH 

    samtools reheader ${analysisType}_${sampleType}_${params.genome}.notsorted.fixedH ${analysisType}_${sampleType}_${params.genome}.notsorted.rawH.bam > ${analysisType}_${sampleType}_${params.genome}.notsorted.bam 
    
    echo "sort by coordinates" >&2
    samtools sort -@ "$task.cpus" ${analysisType}_${sampleType}_${params.genome}.notsorted.bam -o ../${analysisType}_${params.genome}/"\$SM".bam 

    cd ../${analysisType}_${params.genome}/

    echo "index to ready for testing" >&2
    samtools index -@ "$task.cpus" "\$SM".bam


    echo "done" >&2

    
    """
}
