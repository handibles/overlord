#!/bin/bash

#SBATCH --job-name=fqc20
#SBATCH --output=fqc20.%j.slrm
#SBATCH --ntasks=20

INTXT=$1
OUT=$2


mkdir -p $OUT/fwd $OUT/rev ${OUT}_multi

if [[ $? = 0 ]] && [[ -s $INTXT ]] ;
then
  echo " + + +   running fastqc   --- - - -- - -   "
  module load fastqc multiqc
  fastqc -t 20 $(cat $INTXT) -o $OUT   #  just about works!!
else
  echo " < ! >   something awry with the fastq.gz given to FQC/MQC, or the output <dirs> :
     < ! >   check <dir> structure and file manifest : $INTXT"
  # exit 1
fi

mv $OUT/*R1* $OUT/fwd/
mv $OUT/*R2* $OUT/rev/
echo " + + +   running multiqc   -- -- --- -- - - --- "
multiqc $OUT/fwd -o ${OUT}_multi/multi_fwd
multiqc $OUT/rev -o ${OUT}_multi/multi_rev
