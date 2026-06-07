#!/bin/bash 
set -u

##- Author:  Xu Cai
##- Email:   caixu@caas.cn 
##- Date:    06/08/2026
##----------------------------------------------------------------------------------------------------------------
##- The current script was developed to query the frequencies of given k-mers in resequencing data. 
##----------------------------------------------------------------------------------------------------------------

ulimit -n 20480
##-- inputs --------------------------------------------------
samid=$1                   ##- samid
representativekmerList=$2  ##- merged.kmers.list.Polymorphic_kmers.List.representative.list
ksize=31                   ##- ksize

##----- check the TEMP DIR ----------------------------------
if [ -d "logs" ];then 
     sleep 0.0001
else 
     mkdir logs
fi

if [ -f "./logs/${samid}.ok" ];then
     exit 0
else
     sleep  0.0001
fi

if [ -d "kmc_tmp" ];then
     sleep 0.0001
else  
     mkdir kmc_tmp
fi

if [ -d "tmp_sort" ]; then
     sleep 0.0001
else
     mkdir tmp_sort
fi

##-- peocessing --------------------------------------------------
sortedkmerList=${representativekmerList}.sorted
kmerdb=${representativekmerList}.kmcdb

if [ -n "${representativekmerList}" ];then
     if [ -f "${sortedkmerList}" ];then
	  sleep 0.0001
     else
          cut -f 1 ${representativekmerList} | LC_ALL=C sort --parallel=50 --buffer-size=50% --temporary-directory=./tmp_sort -k1,1   > ${representativekmerList}.sorted
     fi	  

     if [ -f "${kmerdb}.kmc_suf" ];then
          sleep 0.0001
     else
          awk '{print ">"NR"\n"$0}' ${sortedkmerList} > ${representativekmerList}.sorted.fasta
          kinfo="-k"${ksize}
	  kmc ${kinfo}  -t50  -m80 -ci1  -fa  ${representativekmerList}.sorted.fasta   ${representativekmerList}.kmcdb  ./kmc_tmp
	  rm ${representativekmerList}.sorted.fasta
     fi 
else
    echo -e "Error: Please check the inputs ... "	
fi

##--- start processing ---------------------------------------------------------------------------------------
querykmcdb=${kmerdb}                     ##- merged.Bra.kmc_suf  merged.Bra.kmc_pre 
querydumpkmersorted=${sortedkmerList}    ##- Bra_42genomes.31mers.list.sorted

currentTime=$(date)
echo -e "${currentTime} ---- Start: ${samid} ----"

##- file extension
tmptypeA=${samid}"_1.fq.ft.gz"
tmptypeB=${samid}"_1.fastq.gz"

if [ -f "${tmptypeA}" ];then
     leftread=${samid}"_1.fq.ft.gz"
     rightread=${samid}"_2.fq.ft.gz"
fi

if [ -f "${tmptypeB}" ];then
     leftread=${samid}"_1.fastq.gz"
     rightread=${samid}"_2.fastq.gz"
fi

querydbfile=${querykmcdb}".kmc_suf"
querydbindex=${querykmcdb}".kmc_pre"

if [ -f ${querydbfile} ] && [ -f ${querydbindex} ];then
     sleep 0.0001

else
     echo "Querykmcdb(i.e. merged.Bra.kmc_suf  merged.Bra.kmc_pre) not exists !!! "
     exit 0
fi


###---------- start ------------------------------------------

if [ -f "${leftread}" ] && [ -f "${rightread}" ]; then

     filename=${samid}.filename.txt
     ls ${leftread} ${rightread} > ${filename}     
     #tmpFastq=${samid}".fastq"
     #pigz -dc ${leftread}  -p 10  >  ${tmpFastq}
     #pigz -dc ${rightread} -p 10  >> ${tmpFastq}

     ##- generate kmcdb
     currentTime=$(date)
     echo -e "${currentTime} ---- generate ${samid} kmcdb ----"
     kmcdb=${samid}".kmcdb.kmc_suf"
     kmcindex=${samid}".kmcdb.kmc_pre"
     kinfo="-k"${ksize}
     kmc  ${kinfo} -t50 -m80  -ci1 -cx10000 -fq  @${filename}  ${samid}".kmcdb"   kmc_tmp       ##- it depends_   m: RAM: 80G; t: threads: 50  ci min/ cx max

     ##- overlap with query
     currentTime=$(date)
     echo -e "${currentTime} ---- get overlaped kmcdb ----"
     overlapfile=${samid}"_overlap"
     kmc_tools  simple  ${samid}".kmcdb"  ${querykmcdb}  intersect  ${overlapfile}   -ci1   -ocleft 
   
     ##- kmc_dump
     currentTime=$(date)
     echo -e "${currentTime} ---- dump overlaped kmcdb ----"
     ovelapkmerdump=${samid}"_overlap.dump.list"
     kmc_dump  ${overlapfile} ${ovelapkmerdump}

     ##- sort kmc_dump list
     currentTime=$(date)
     echo -e "${currentTime} ---- sort overlap kmers ----"
     ovelapkmerdumpsorted=${samid}"_overlap.dump.sorted.list"
     LC_ALL=C sort --parallel=50 --buffer-size=50% --temporary-directory=./kmc_tmp  -k1,1  ${ovelapkmerdump} >  ${ovelapkmerdumpsorted}

     ##- generate kmer counts
     currentTime=$(date)
     echo -e "${currentTime} ---- run join to generate ${samid}  kmer counts ----" 
     kcountFile=${samid}".k${ksize}.counts"
     LC_ALL=C  join -a 1 -e 0 -o  2.2  ${querydumpkmersorted}   ${ovelapkmerdumpsorted}  >  ${kcountFile} 
   
     ##- pigz 
     pigz ${kcountFile} 
      
     touch ${samid}".ok"
     mv ${samid}".ok"  logs

     ##- clean
     rm -rf ${filename}  ${kmcdb} ${kmcindex}  ${overlapfile} ${ovelapkmerdump} ${ovelapkmerdumpsorted}   ${overlapfile}".kmc_suf"   ${overlapfile}".kmc_pre"  ${leftread}  ${rightread}
     currentTime=$(date)
     echo -e "${currentTime} ---- Finished: ${samid} ----"    
      
else
    sleep 0.0001
fi
