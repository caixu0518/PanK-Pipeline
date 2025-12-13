
# Here, we provide a test dataset in the file testdata.zip (1.90 GB), which includes all input files, intermediate outputs from each step, and the command lines used in the workflow.


## Test dataset download
The testdata.zip file can be downloaded directly from the shared Google Drive link or the shared Baidu Cloud Netdisk.

From Google Drive (https://drive.google.com/file/d/1vocLkzop2grbC2qVCUrMEN0jjSQlF3D-/view?usp=sharing)

From Baidu Netdisk (https://pan.baidu.com/s/191eTuPq2gpQRilfgmkIXjg, extraction code: 1111)

## How to Use the Test Dataset: A Step-by-Step Workflow

'runme.sh' is a shell script containing the commands to run the pipeline. Please update the absolute path to the script according to the actual installation directory of PanK-Pipeline. The scripts Generate_PolymorphicKmers.pl, Generate_RepresentativeKmers.pl, and PopKmerGenotypesToVCF.pl can be copied directly to the current working directory and run.

```
unzip testdata.zip
cd testdata

##- Step1. Pipeline for identifying Pan-genome polymorphic k-mers
perl  Generate_PolymorphyicKmers.pl    -species  rapa  -ksize 17   -pangenome   Pangenome.txt  -PipelinePath    /mydata/caix/PanK-Pipeline

##- Step2: Pipeline for identifying Pan-genome representative k-mers
perl  Generate_RepresentativeKmers.pl  -species  Brapa  -ksize  17  -pangenome  Pangenome.txt    -PipelinePath   /mydata/caix/PanK-Pipeline   -PolymorphicKmer  rapa.merged.kmer.k17.list.Polymorphic_kmers.List

##- Application of Pan-genome representative k-mers for population structure analysis
cd genotyping
perl PopKmerGenotypesToVCF.pl    -Sam   samid.txt   -KmerGTDir  /mydata/caix/PanK-Pipeline/testdata/genotyping/countresults   -KmerList rapa.merged.kmer.k17.list.Polymorphic_kmers.List.representative.list   -threads  20

```
