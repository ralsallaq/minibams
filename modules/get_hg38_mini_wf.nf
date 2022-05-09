/* ========================= workflow to get hg38 minibams  =================================*/
params.outD = null
params.events_bed_file = null
params.events_bedpe_file = null
params.minibams_dir = null 
params.target_genome = null
params.genome = null
params.random_seed = null

// Workflows


// Split minibams by RGs
include { splitMinibamByRGs as splitMinibamByRGs_RNA_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { splitMinibamByRGs as splitMinibamByRGs_WGS_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { splitMinibamByRGs as splitMinibamByRGs_WGS_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { splitMinibamByRGs as splitMinibamByRGs_WES_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { splitMinibamByRGs as splitMinibamByRGs_WES_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)


// Convert bams to fastqs
include { bam2fastq as bam2fastq_RNA_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { bam2fastq as bam2fastq_WGS_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { bam2fastq as bam2fastq_WGS_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { bam2fastq as bam2fastq_WES_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { bam2fastq as bam2fastq_WES_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

// Map fastqs to target genome 
include { bwaMEM as bwaMEM_RNA_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { bwaMEM as bwaMEM_WGS_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { bwaMEM as bwaMEM_WGS_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

include { bwaMEM as bwaMEM_WES_tumor } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)
include { bwaMEM as bwaMEM_WES_gl } from './split_map_minibam' params(genome:params.target_genome, outD:params.outD)

// Merge subbams to minibams
include { mergeToMinibams as mergeToMinibams_RNA_tumor } from './get_minibams' params(genome:params.target_genome, random_seed:params.random_seed, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WGS_tumor } from './get_minibams' params(genome:params.target_genome, random_seed:params.random_seed, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WGS_gl } from './get_minibams' params(genome:params.target_genome, random_seed:params.random_seed, outD:params.outD)

include { mergeToMinibams as mergeToMinibams_WES_tumor } from './get_minibams' params(genome:params.target_genome, random_seed:params.random_seed, outD:params.outD)
include { mergeToMinibams as mergeToMinibams_WES_gl } from './get_minibams' params(genome:params.target_genome, random_seed:params.random_seed, outD:params.outD)

workflow get_hg38_mini_wf { 
    main:

        liftoverSVs(channel.value('GRCh37-lite'), channel.value('GRCh38'), channel.fromPath(params.events_bedpe_file)) 
        liftoverBED(channel.value('GRCh37-lite'), channel.value('GRCh38'), channel.fromPath(params.events_bed_file))

        // read already generated minibams using main 
        minibam_RNA_tumor_ch = channel.fromPath(params.minibams_dir + '/TRANSCRIPTOME_' + params.genome + '/*_D1.bam')
    
        minibam_WGS_tumor_ch = channel.fromPath(params.minibams_dir + '/WHOLE_GENOME_' + params.genome + '/*_D1.bam')
        minibam_WGS_gl_ch = channel.fromPath(params.minibams_dir + '/WHOLE_GENOME_' + params.genome + '/*_G1.bam')
    
        minibam_WES_tumor_ch = channel.fromPath(params.minibams_dir + '/EXOME_' + params.genome + '/*_D1.bam')
        minibam_WES_gl_ch = channel.fromPath(params.minibams_dir + '/EXOME_' + params.genome + '/*_G1.bam')

        // Split minibams by RGs
        splitMinibamByRGs_RNA_tumor(minibam_RNA_tumor_ch, channel.value('TRANSCRIPTOME'), channel.value('tumor')) 
    
        splitMinibamByRGs_WGS_tumor(minibam_WGS_tumor_ch, channel.value('WHOLE_GENOME'), channel.value('tumor')) 
        splitMinibamByRGs_WGS_gl(minibam_WGS_gl_ch, channel.value('WHOLE_GENOME'), channel.value('germline')) 
    
        splitMinibamByRGs_WES_tumor(minibam_WES_tumor_ch, channel.value('EXOME'), channel.value('tumor'))
        splitMinibamByRGs_WES_gl(minibam_WES_gl_ch, channel.value('EXOME'), channel.value('germline')) 

        // Get fastqs from each bam split
        bam2fastq_RNA_tumor(splitMinibamByRGs_RNA_tumor.out.flatten(), channel.value('TRANSCRIPTOME'), channel.value('tumor')) 

        bam2fastq_WGS_tumor(splitMinibamByRGs_WGS_tumor.out.flatten(), channel.value('WHOLE_GENOME'), channel.value('tumor')) 
        bam2fastq_WGS_gl(splitMinibamByRGs_WGS_gl.out.flatten(), channel.value('WHOLE_GENOME'), channel.value('germline')) 
    
        bam2fastq_WES_tumor(splitMinibamByRGs_WES_tumor.out.flatten(), channel.value('EXOME'), channel.value('tumor'))
        bam2fastq_WES_gl(splitMinibamByRGs_WES_gl.out.flatten(), channel.value('EXOME'), channel.value('germline')) 

        // Map fastqs to human genome
        bwaMEM_RNA_tumor(bam2fastq_RNA_tumor.out, channel.value('TRANSCRIPTOME'), channel.value('tumor'), channel.value('GRCh38'))
    
        bwaMEM_WGS_tumor(bam2fastq_WGS_tumor.out, channel.value('WHOLE_GENOME'), channel.value('tumor'), channel.value('GRCh38'))
        bwaMEM_WGS_gl(bam2fastq_WGS_gl.out, channel.value('WHOLE_GENOME'), channel.value('germline'), channel.value('GRCh38'))
    
        bwaMEM_WES_tumor(bam2fastq_WES_tumor.out, channel.value('EXOME'), channel.value('tumor'), channel.value('GRCh38'))
        bwaMEM_WES_gl(bam2fastq_WES_gl.out, channel.value('EXOME'), channel.value('germline'), channel.value('GRCh38'))


        mergeToMinibams_RNA_tumor(bwaMEM_RNA_tumor.out.map{it->it[1]}.collect(), channel.value('TRANSCRIPTOME'), channel.value('tumor'))
        mergeToMinibams_WGS_tumor(bwaMEM_WGS_tumor.out.map{it->it[1]}.collect(), channel.value('WHOLE_GENOME'), channel.value('tumor'))
        mergeToMinibams_WGS_gl(bwaMEM_WGS_gl.out.map{it->it[1]}.collect(), channel.value('WHOLE_GENOME'), channel.value('germline'))
    
        mergeToMinibams_WES_tumor(bwaMEM_WES_tumor.out.map{it->it[1]}.collect(), channel.value('EXOME'), channel.value('tumor'))
        mergeToMinibams_WES_gl(bwaMEM_WES_gl.out.map{it->it[1]}.collect(), channel.value('EXOME'), channel.value('germline'))
        

}

process liftoverSVs {
    tag 'liftover segments'
    label 'io_mem'
    publishDir "${params.outD}/", mode: 'copy'
 
    input:
    val(genomeFrom)
    val(genomeTo)
    path(bedpe_events)

    output:
    path("bedpe_events.bedpe"), emit: liftoverBEDPE_ch

    """
    #!/usr/bin/bash
    module load python/3.7.0
    echo "you need to set up the environment for liftover_flatfile.pl and other programs by: \
    cbload configs \
    cbload common-scripts-internal \
    module load perl/5.10.1 #sjcb perl 5.30.0 did not work
    module load sjcb/ucsc \
    module load sjcb/vcftools \
    module load sjcb/blat/35 \
      \
    cbload seq_anls_ops \
    cbload snv-annovar \
    "
    if [ ! \$(echo ${genomeFrom} | grep 'GRCh') ] || [ ! \$(echo ${genomeTo} | grep 'GRCh') ]
    then
       echo "official names of genome from and to are needed for liftover_sv.sh" >&2
       exit 1
    fi

    echo "processing liftover from ${genomeFrom} to ${genomeTo}"

    echo "prep the raw SV calls file for liftover_sv.sh script"
    python -c "import pandas as pd
import sys
import os
import pandas as pd
outputFile = 'reformat_bedpe.txt' 
df = pd.read_csv('${bedpe_events}', sep='\t', header=None) 
cols_idx=[0,2,3,5,6,7,8,9] + list(range(10,df.shape[1])) 
cols = ['ChrA','PosA','ChrB','PosB','Type','score','OrtA','OrtB', 'gene', 'cytolocus', 'eventType']
# No change the shape of df
df = df.loc[:,cols_idx]
df.columns = cols

# Drop null lines
dx = (~df['ChrA'].isnull()) & (~df['ChrB'].isnull()) & (~df['PosA'].isnull()) & (~df['PosB'].isnull())

# make sure no non-numeric values
ix1 = pd.to_numeric(df['PosA'], errors='coerce').isnull()
df = df.loc[~ix1]
assert(df.loc[ix1].shape[0]==0),'non-numeric values detected in the swab file in PosA'

ix2 = pd.to_numeric(df['PosB'], errors='coerce').isnull()
df = df.loc[~ix2]
assert(df.loc[ix2].shape[0]==0),'non-numeric values detected in the swab file in PosB'

## change the name of chromosomes from chr1 to 1 
def addChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==2 else 'chr'+row.strip()
    else:
        new_r = str(int(row)).strip()
        return 'chr'+new_r

def removeChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==1 else row.split('chr')[1].strip()
    else:
        new_r = str(int(row)).strip()
        return new_r
def formatChrom(row, fromChromosomeFormat):
    if len(fromChromosomeFormat.split('chr'))>1: #from format is chrX, ..etc
        return addChr(row)
    else:
        return removeChr(row)

    
if ('${genomeFrom}' == 'GRCh37-lite'):
    fromChromosomeFormat = 'X'
elif ('${genomeFrom}' == 'GRCh38'):
    fromChromosomeFormat = 'chrX'

df.loc[:,'ChrA'] = df['ChrA'].apply(lambda r: formatChrom(r,fromChromosomeFormat)).astype(pd.StringDtype())

#### If from format is X (i.e. we are lifting over from hg19, chromosome names such as Un_g1000224, ..etc are not acceptable and will cause liftover_sv.sh to fail

if not fromChromosomeFormat.startswith('chr'): #lifting from hg19
    ## detecting when chromosome name is not X or Y and is not numeric
    idx1 = (~df['ChrA'].isin(['X','Y']) ) & (~df['ChrA'].str.isnumeric()) 
    df = df.loc[~idx1]
    idx2 = (~df['ChrB'].isin(['X','Y']) ) & (~df['ChrB'].str.isnumeric()) 
    df = df.loc[~idx2]

#### If from format is chrX (i.e. we are lifting over from hg38, chromosome names with alt in them (e.g. chr7_KI270803v1_alt will cause liftover_sv.sh to fail b/c these are not in the reference /research/rgs01/resgen/ref/tartan/runs/ad_hoc/NAME-4tmOgtgc/output/sequence/fasta/GRCh38_no_alt.fa
if fromChromosomeFormat.startswith('chr'): #lifting from hg38
    ## detecting when chromosoe name has _alt in it
    idx1 = df['ChrA'].apply(lambda r: True if len(r.split('_alt'))>1 else False)
    df = df.loc[~idx1]
    idx2 = df['ChrB'].apply(lambda r: True if len(r.split('_alt'))>1 else False)
    df = df.loc[~idx2]

### Make sure that positions are intgers
df.loc[:,'PosA'] = df['PosA'].astype(int)
df.loc[:,'PosB'] = df['PosB'].astype(int)

df.head()
df.to_csv(outputFile, sep='\t', index=False)
" 
# End of python script

    # apply leftover 
    echo "applying liftover_sv.sh which produces the same *.predSV file with extra columns"
    if [ ! grep -w 'GRCh' <(echo ${genomeFrom}) ]; then echo "official names of genomes must be supplied" >&2; exit 1; fi
    if [ ! grep -w 'GRCh' <(echo ${genomeTo}) ]; then echo "official names of genomes must be supplied" >&2; exit 1; fi
    liftover_sv.sh -w tempDir  ${genomeFrom} ${genomeTo}  reformat_bedpe.txt  bedpe_events_raw
    wait
    echo "save a bedpe file"

python -c "import numpy as np
import pandas as pd
import os
import sys
import time
inputFWpath = 'bedpe_events_raw'
outputFWpath = os.path.basename(inputFWpath).split('_raw')[0]+'.bedpe'
dflo = pd.read_csv(inputFWpath, sep='\t')
print(dflo.head())
assert(dflo.columns.isin('ChrA\tPosA\tOrtA\tChrB\tPosB\tOrtB\tType\tgene\tcytolocus\teventType\tliftover_ok_PosA\tliftover_chr_PosA\tliftover_base_PosA\tliftover_ort_PosA\tliftover_ok_PosB\tliftover_chr_PosB\tliftover_base_PosB\tliftover_ort_PosB\tliftover_min_match\tliftover_interval_length_delta\tliftover_interval_length_delta_fraction'.split('\t')).sum()==21),'some columns are missing! check the dataframe columns{}'.format(dflo.columns.values)
# Take only the liftover columns and assign them simple names
dflo = dflo[['liftover_chr_PosA','liftover_base_PosA','liftover_ort_PosA','liftover_chr_PosB','liftover_base_PosB','liftover_ort_PosB','Type','score', 'gene', 'cytolocus', 'eventType']]
dflo.columns = ['ChrA','PosA','OrtA','ChrB','PosB','OrtB','Type','score', 'gene', 'cytolocus', 'eventType']

dflo.loc[:, 'PosA1'] = dflo['PosA'].astype(int) - 1
dflo.loc[:, 'PosB1'] = dflo['PosB'].astype(int) - 1

dflo = dflo[['ChrA', 'PosA1', 'PosA', 'OrtA', 'ChrB', 'PosB1', 'PosB', \
           'OrtB', 'Type', 'score', 'gene', 'cytolocus', 'eventType']]

print('saving a liftover file')
dflo.to_csv(outputFWpath, sep='\t', index=False)
"
# End of python script
echo "Done"

        """
}

process liftoverBED {
    tag 'liftover loci'
    label 'io_mem'
    publishDir "${params.outD}/", mode: 'copy'
 
    input:
    val(genomeFrom)
    val(genomeTo)
    path(bed_events)

    output:
    path("bed_events.bed"), emit: liftoverBED_ch

    """
    #!/usr/bin/bash
    module load python/3.7.0
    echo "you need to set up the environment for liftover_flatfile.pl and other programs by: \
    cbload configs \
    cbload common-scripts-internal \
    module load perl/5.10.1 #sjcb perl 5.30.0 did not work
    module load sjcb/ucsc \
    module load sjcb/vcftools \
    module load sjcb/blat/35 \
      \
    cbload seq_anls_ops \
    cbload snv-annovar \
    "

    if [ \$(echo ${genomeFrom} | grep 'GRCh38') ] 
    then
        genomeFrom="hg38"
    elif [ \$(echo ${genomeFrom} | grep 'GRCh37') ]
    then
        genomeFrom="hg19"
    else
       echo "The official name ${genomeFrom} is not recognized" >&2
       exit 1
    fi

    if [ \$(echo ${genomeTo} | grep 'GRCh38') ] 
    then
        genomeTo="hg38"
    elif [ \$(echo ${genomeTo} | grep 'GRCh37') ]
    then
        genomeTo="hg19"
    else
       echo "The official name ${genomeTo} is not recognized" >&2
       exit 1
    fi

    echo "processing liftover from ${genomeFrom} to ${genomeTo}"

    echo "prep the bed file for liftover_flatfile.pl script"
    python -c "import pandas as pd
import sys
import os
import pandas as pd
outputFile = 'reformat_bed.txt' 
df = pd.read_csv('${bed_events}', sep='\t', header=None) 
cols_idx=[0,2,3] + list(range(5,df.shape[1])) 
cols = ['ChrA','PosA','Type', 'gene', 'cytolocus', 'eventType']
# Now change the shape of df
df = df.loc[:,cols_idx]
df.columns = cols

# Drop null lines
dx = (~df['ChrA'].isnull()) & (~df['PosA'].isnull()) 

# make sure no non-numeric values
ix1 = pd.to_numeric(df['PosA'], errors='coerce').isnull()
df = df.loc[~ix1]
assert(df.loc[ix1].shape[0]==0),'non-numeric values detected in the swab file in PosA'

## change the name of chromosomes from chr1 to 1 
def addChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==2 else 'chr'+row.strip()
    else:
        new_r = str(int(row)).strip()
        return 'chr'+new_r

def removeChr(row):
    if type(row) is str:
        return row.strip() if len(row.split('chr'))==1 else row.split('chr')[1].strip()
    else:
        new_r = str(int(row)).strip()
        return new_r
def formatChrom(row, fromChromosomeFormat):
    if len(fromChromosomeFormat.split('chr'))>1: #from format is chrX, ..etc
        return addChr(row)
    else:
        return removeChr(row)

    
if ('${genomeFrom}' == 'GRCh37-lite'):
    fromChromosomeFormat = 'X'
elif ('${genomeFrom}' == 'GRCh38'):
    fromChromosomeFormat = 'chrX'

df.loc[:,'ChrA'] = df['ChrA'].apply(lambda r: formatChrom(r,fromChromosomeFormat)).astype(pd.StringDtype())


if not fromChromosomeFormat.startswith('chr'): #lifting from hg19
    ## detecting when chromosome name is not X or Y and is not numeric
    idx1 = (~df['ChrA'].isin(['X','Y']) ) & (~df['ChrA'].str.isnumeric()) 
    df = df.loc[~idx1]

if fromChromosomeFormat.startswith('chr'): #lifting from hg38
    ## detecting when chromosoe name has _alt in it
    idx1 = df['ChrA'].apply(lambda r: True if len(r.split('_alt'))>1 else False)
    df = df.loc[~idx1]

### Make sure that positions are intgers
df.loc[:,'PosA'] = df['PosA'].astype(int)

df.head()
df.to_csv(outputFile, sep='\t', index=False)
" 
# End of python script

    # apply leftover 
    echo "applying liftover_flatfile.pl"
    if [ ! grep -w 'GRCh' <(echo ${genomeFrom}) ]; then echo "official names of genomes must be supplied" >&2; exit 1; fi
    if [ ! grep -w 'GRCh' <(echo ${genomeTo}) ]; then echo "official names of genomes must be supplied" >&2; exit 1; fi
    liftover_flatfile.pl -genome-from  \$genomeFrom -genome-to \$genomeTo  -file reformat_bed.txt -field-chr ChrA -field-pos PosA -out bed_events_raw 
    wait
    echo "save a bed file"

python -c "import numpy as np
import pandas as pd
import os
import sys

inputFWpath = 'bed_events_raw'
outputFWpath = os.path.basename(inputFWpath).split('_raw')[0]+'.bed'
dflo = pd.read_csv(inputFWpath, sep='\t')
print(dflo.head())
assert(dflo.columns.isin('ChrA\tPosA\tType\tgene\tcytolocus\teventType\tliftover_ok\tliftover_chr\tliftover_base\tliftover_min_match\tliftover_interval_length_delta\tliftover_interval_length_delta_fraction'.split('\t')).sum()==12),'some columns are missing! check the dataframe columns{}'.format(dflo.columns.values)
# Take only the liftover columns and assign them simple names
dflo = dflo[['liftover_chr', 'liftover_base', 'Type', 'gene', 'cytolocus', 'eventType']]
dflo.columns = ['ChrA', 'PosA', 'Type', 'gene', 'cytolocus', 'eventType']

dflo.loc[:,'PosA1'] = dflo['PosA'].astype(int) - 1

dflo = dflo[['ChrA', 'PosA1', 'PosA', 'Type', 'gene', 'cytolocus', 'eventType']]


print('saving a liftover file')
dflo.to_csv(outputFWpath, sep='\t', index=False)
"
# End of python script
echo "Done"

        """
}
