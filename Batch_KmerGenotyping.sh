#!/bin/bash  


##- A shell script to call CountKmersInReads.pl for batch genotyping of k-mers in each resequenced sample.
   
    ##- Bug reports: Xu Cai
    ##- caixu@caas.cn       
    ##- 18,8,2025


kmerList=$1      ##- *.list.Polymorphic_kmers.List.representative.list
samidfile=$2     ##-   Brapasamids.txt
ksize=$3         ##-   k-mer size  17
jfFileDir=$4     ##-   the dir contains jf files   "/mydata/caix/CC_k17_analysis/jffiles"
scriptdir=$5     ##-   /mydata/caix/PanK-Pipeline/scripts



if [ -d "countresults" ];then
     sleep 0.001
else
     mkdir countresults
fi


for eachsam in  `cat  ${samidfile}`

do
	countresult=${eachsam}.gz
        jfFile=${eachsam}.jf

	if [ -f "./countresults/${countresult}" ]; then
		sleep 0.0001
        else
               
                if [ -f "${jfFileDir}/${jfFile}" ]; then   
		    ln -s ${jfFileDir}/${jfFile}   .

		    echo "[Current system time: $(date)] Start: ${eachsam} ..."
		    perl  ${scriptdir}/CountKmersInReads.pl   ${kmerList}    ${eachsam}  ${ksize}  ${eachsam}   ${scriptdir}
                    echo "[Current system time: $(date)] Finish: ${eachsam} ..."   

		    mv  ${countresult}   countresults
                    rm -rf ${jfFile} 		 
                else
                     echo "${jfFileDir}/${jfFile} not exists"
		fi
	fi

done
