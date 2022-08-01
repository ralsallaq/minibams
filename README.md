# FrankenBam minibams
These are bams built from regions of interest that have been collected from the bams of several real samples. The goal is to be able to use these in pipelines testing and validation.

## Getting Started
Get nextflow and test it
```
module load sjcb/openjdk/17.0.0_35
curl -s https://get.nextflow.io | bash

./nextflow run nextflow-io/hello

```

You would need a nextflow config file in your run directory. Here a suggested one to put in a file nextflow.config

```

profiles{

    lsfCluster {
        process {
            executor = 'lsf'
            queue = 'compbio_auto'

            withLabel: 'io_limited' {
                cpus = 1
                memory = 4.GB
            }
            withLabel: 'io_mem' {
                cpus = 1
                memory = 8.GB
            }
            withLabel: 'multithread' {
                cpus = 8
                memory = 40.GB
            }
        }

    }

}

```

### Get frankenbams from events file
If you have an event file that looks like this
```
cat events_input_file.tsv
sample          gene           abnormality       type
SJALCL030011    NPM1-ALK       Fusion       Somatic SV
SJAML001417     GATA2          R337fs       Somatic INDEL
SJAML001417     SETD2          S1572fs      Somatic INDEL
SJAML001417     WT1            R380fs       Somatic INDEL
SJAML001417     FLT3           ITD          Somatic INDEL/SV
SJAML001417     BCL11A-GRIP2   Fusion       Somatic SV
SJAML001417     NPM1-CCDC88A   Fusion       Somatic SV

```
Then you can build five frankenbams using the following command which executes
on lsf. Change the random seed to change the integer part in the subject id (e.g. SJ992966).
The subject id is generated randomly and always will be in the range [90000, 99999].
The lookupDircs are the path to tartan index project/subproject directory where minibams workflow
would be looking for the samples and events in the events_input_file.tsv.

```
./nextflow run ralsallaq/minibams main.nf \
      --random_seed 902 \
      --outD /path/to/output/dir \
      --abnormal_eventsFile events_input_file.tsv \
      --lookupDircs /path/to/tartan/index/project/subproject/dir
      --genome hg19
      -w /path/to/work/dir \
      -c nextflow.config \
      -profile lsfCluster \
      -resume"

```

### Add regions from a bed file
Chances are you have already built five minibams and now you have a project under tartan
and you want to add regions from samples under that project given a bed file,
which a 0-based tab-delimited text file that have the following required format
```
Sample  chr     posA    posB
SJHYPER011    MT      12223   13223
SJHYPER011    Y       3130987 3131987
SJMLL002     14      22564394        22565394
```

Then you can run the minibam workflow like so

```
./nextflow run ralsallaq/minibams main.nf \
               --addLociToExisitngMB \
               --tartanIndex_projDir /path/to/project/subproject/dir/under/tartan \
               --exBedFile /path/to/the/external/bed/file \
               --genome hg19 \
               --outD path/to/output/dir \
               --random_seed 902 \
               --minibams_dir /path/to/exisiting/minibams/dir \
               -c nextflow.config \
               -profile lsfCluster \
               -w /path/to/work/dir \
               -resume
```
Note here that I used the random seed 902 as I do not want to change the subject id after
adding regions to the exisiting bams. In case you need to change the subject id after adding
regions, then change the random seed to a different integer than the one used to generate the
exisiting minibams.


### Recover events from a post-analysis run index directory
Now that you have run the frankenbams through pipeline(s) and you want to recover the events, how do you do that? 
Here, we are going to use the recover mode. First edit the events_input_file.tsv to replace sample names by the frankenbam name:
```
cat events_input_file.tsv | awk '{NR !=1 ?$1="SJOther9XXXXX":$1=$1; print;}' > events_input_file_recover.tsv

# Then run:

./nextflow run ralsallaq/minibams main.nf \
      --recover \
      --random_seed 902 \
      --outD /path/to/recover/output/dir \
      --abnormal_eventsFile events_input_file_recover.tsv \
      --lookupDircs /path/to/tartan/index/project/subproject/dir/for/frankenbam/
      --genome hg19
      -w /path/to/work/dir \
      -c nextflow.config \
      -profile lsfCluster \
      -resume"
```

The output in the recovery mode will encompass bed-like and bedpe-like files that can be used to compare with 
the corresponding files produced when generating the frankenbams
 
