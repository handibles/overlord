#!/bin/bash

#SBATCH --job-name=ommfire_samp
#SBATCH --output=ommfire_samp.%j.slrm
#SBATCH --ntasks=15
#SBATCH --time=90:00


##   i n f r a s t r u c t u r e   -------------------------


# just give sample and WRKing dir
SAMPLE=$1
PPATH=$2
RAWD=$3

          # # ##test stuff
          # SAMPLE=fhi__tester__FHI112_S78
          # PPATH=$ODIR
          # RAWD=$RAWDIR


# for all steps, check matches SLURM
OMMTHREADS=15

PROJECT=${PPATH##*/}
PRSA=${PROJECT}__${SAMPLE}
QC=${PPATH}/1__qc
FILT=${PPATH}/2__filt
DECO=${PPATH}/3__deco
KRAK=${PPATH}/4__krak2
MAT=${PPATH}/Materials

if [[ ! -d $PPATH ]] || [[ ! -d $MAT ]]; then
  echo " < ! >   main dir for sample-fire $SAMP (dir: $PPATH) not sane - check"
  # exit 1
fi

for i in $QC $FILT $DECO $KRAK \
  $QC/${PROJECT}_trim $QC/${PROJECT}_trim_multi \
  $QC/${PROJECT}_deco $QC/${PROJECT}_deco_multi ;
do 
  mkdir -p $i ;
done


##   f a s t q c   t r i m m    p o i n t s   -------------------------

# $MSDAT/$GENU/1__qc/fhi__tester_raw_multi/multi_fwd/multiqc_report
## unify this as will all be the same
OMM_FW5=$( head -1 $QC/${PROJECT}_raw_multi/multi_fwd/multiqc_report_MQC-trim-points.txt )
OMM_FW3=$( tail -n 1 $QC/${PROJECT}_raw_multi/multi_fwd/multiqc_report_MQC-trim-points.txt )
OMM_RE5=$( head -1 $QC/${PROJECT}_raw_multi/multi_rev/multiqc_report_MQC-trim-points.txt )
OMM_RE3=$( tail -n 1 $QC/${PROJECT}_raw_multi/multi_rev/multiqc_report_MQC-trim-points.txt )


##   t r i m m o m a t i c   -------------------------

module load trimmomatic
echo " + + +   trimmo reads together  -   read order unfixed   <!>  ---- "
echo " + + +   trimmo :: ILLUMINACLIP:fa:2:30:10:5, CROP: ${OMM_FW3}  --- "
echo " + + +   trimmo :: HEADCROP: ${OMM_FW5}, SLIDINGWINDOW:6:15, MINLEN:100  --- "
trimmomatic PE \
  $RAWD/${SAMPLE}_R1.fastq.gz \
  $RAWD/${SAMPLE}_R2.fastq.gz \
  $FILT/${PRSA}_trimm_R1.fastq.gz \
  $FILT/${PRSA}_trimtemp1_R1_unpaired.fastq.gz \
  $FILT/${PRSA}_trimm_R2.fastq.gz \
  $FILT/${PRSA}_trimtemp1_R2_unpaired.fastq.gz \
  ILLUMINACLIP:$MAT/fqc_trimmo_ill_ref.fa:2:30:10:5 \
  CROP:$OMM_FW3 \
  HEADCROP:$OMM_FW5 \
  SLIDINGWINDOW:6:15 \
  MINLEN:100 \
  -threads $OMMTHREADS


echo " + + +   H.sap BT2 deco on ${PRSA}:   -------- "

# load necess
module unload samtools
module load bowtie2 samtools/1.10
BT2_DB=$( grep BT2_DB $MAT/ommfire__0__vars.tsv | cut -f 2 -d ' ' )

if [[ ! -z $BT2_DB ]];
then 
  ## stamp for the output
  DAT_TIM=$(date +"%T %F")
  echo " + + +   --------------------------------------------
  + + +   $DAT_TIM
  + + +   sampl: $PRSA
  + + +   input: $FILT
  + + +   outpt: $DECO
  + + +   to DB: $BT2_DB/hum_deco
  + + +   --------------------------------------------"
else
  echo " < ! >   couldn't find BT2 database ( BT2_DB = )" ; exit 1
fi

echo " + + +   align   -------------------------------------"
time bowtie2 -p $OMMTHREADS \
  -x $BT2_DB/hum_deco \
  -1 $FILT/${PRSA}_trimm_R1.fastq.gz \
  -2 $FILT/${PRSA}_trimm_R2.fastq.gz \
  -S $DECO/${PRSA}_bt2_refmapped.sam   #> $MAT/ommfire__3__deco__${PRSA}.log


if [[ $? = 0 ]] && [[ -s $DECO/${PRSA}_bt2_refmapped.sam ]]
then
  echo " + + +   --   change format..  ------------------------------------"
  samtools view -bS $DECO/${PRSA}_bt2_refmapped.sam > $DECO/${PRSA}_bt2_refmapped.bam
  
  echo " + + +   ----------------   ..remove matches..  --------------------"
  samtools view -b -f 12 -F 256 $DECO/${PRSA}_bt2_refmapped.bam > $DECO/${PRSA}_bt2_decontamd.bam
  
  ## obviated - not enough memory, and not essential for kraken & co.
  # time samtools sort -n -m 5G -@ $OMMTHREADS $DECO/${PRSA}_bt2_decontamd.bam -o $DECO/${PRSA}_bt2_decontamd_sort.bam
  echo " + + +     ( WARNING: DIDNT re-organise reads, for K2/mem sake )  --
  + + +   --------------------------  ...make fastq...   ------------"
  samtools fastq -@ $OMMTHREADS $DECO/${PRSA}_bt2_decontamd.bam -1 $DECO/${PRSA}_bt2decon_R1.fastq.gz -2 $DECO/${PRSA}_bt2decon_R2.fastq.gz -0 /dev/null -s /dev/null -n
else
  echo " < ! >   ---   BT2 align failed   ::   $SAMPLE   ---   < ! >"
  # exit 1
fi


## shrink intermediates
if [[ $? -eq 0 ]] && [[ -s $DECO/${PRSA}_bt2decon_R1.fastq.gz ]] && [[ -s $DECO/${PRSA}_bt2decon_R2.fastq.gz ]] ;
then 
  echo " + + +   -----------------------------  ..shrink/rm *AMs.. ---"
  rm $DECO/${PRSA}_*sam &
  gzip $DECO/${PRSA}_*bam &
  echo " + + +   --------------------------------...  ${PRSA} trimmed.
  +-+-+
  +-+-+"
fi


echo " + + +   Kraken on 2 ${PRSA}   ------------------------------ "

K2_REF=$( grep REF $MAT/ommfire__0__vars.tsv | cut -f 2 -d ' ')

if [[ -s $KRAK/${i}*_kraken_report ]] 
then
  echo " + + +   krak2 output for $i already in $KRAK - stopping   ---
  +-+-+
  +-+-+   "
elif [[ -s  $DECO/${PRSA}_bt2decon_R1.fastq.gz ]] && [[  -s $DECO/${PRSA}_bt2decon_R2.fastq.gz ]] && [[ -s $K2_REF ]]
then
  echo "  + + +   krak2: firing.. .   .    . 
  + + +   to DB: ${K2_REF} .. .  .  . 
  + + +   outpt: ${KRAK}.. .
  + + +   unass: to /dev/null..."
  module load kraken2/2.1.1
  time kraken2 --db $K2_REF \
    $DECO/${PRSA}_bt2decon_R1.fastq.gz \
    $DECO/${PRSA}_bt2decon_R2.fastq.gz \
    --paired \
    --confidence 0.1 \
    --minimum-hit-groups 5 \
    --minimum-base-quality 20 \
    --threads $OMMTHREADS \
    --gzip-compressed \
    --report $KRAK/${PRSA}_kraken_report \
    --output $KRAK/${PRSA}_kraken_output \
    --unclassified-out /dev/null/${PRSA}_kraken_unclass#
  echo " +-+-+   
  + + +   "
else
  echo " < ! >   missing pieces for kraken2 - check deco step or database   --- "
fi


echo " + + +   Beat ${PRSA} back with Bracken   ------------------------------ "
echo " < ! >   Warning: Bracken is using fixed kmer length of 120   ---------- "
if [[ ! -s $KRAK/${PRSA}_kraken_report ]] ;
then 
  echo " < ! >   Kraken2 report ( $KRAK/${PRSA}_kraken_report ) not found "
  # exit 1
else
  BR_r=120  # 3trim-5trim=120 
  BR_l=S
  BR_t=10   # counts! not threads
  module unload kraken2
  module load braken
  bracken -d $K2_REF/ -i $KRAK/${PRSA}_kraken_report -o $KRAK/${PRSA}.bracken -r $BR_r -l $BR_l -t $BR_t
fi

echo " + + +   Bracken output: should be $KRAK/${PRSA}.bracken --------------
+-+-+
+-+-+
+-+-+"


##   w r a p   u p   t h e   m a t t e r 
if [[ $? -eq 0 ]] && \
  [[ -s $FILT/${PRSA}_trimm_R1.fastq.gz ]] && \
  [[ -s $FILT/${PRSA}_trimm_R2.fastq.gz ]] && \
  [[ -s $DECO/${PRSA}_bt2decon_R1.fastq.gz ]] && \
  [[ -s $DECO/${PRSA}_bt2decon_R2.fastq.gz ]] && \
  [[ -s $KRAK/${PRSA}_kraken_report ]] && \
  [[ -s $KRAK/${PRSA}.bracken ]] ; then

  echo " +-+-+   ${PRSA} completed successfully.   =================
 +-+-+   ========================================================================
 
 "
  # these get large
  gzip $KRAK/${PRSA}_kraken_output &

  echo "$SAMPLE \
    $FILT/${PRSA}_trimm_R1.fastq.gz \
    $FILT/${PRSA}_trimm_R2.fastq.gz \
    $DECO/${PRSA}_bt2decon_R1.fastq.gz \
    $DECO/${PRSA}_bt2decon_R2.fastq.gz \
    $KRAK/${PRSA}_kraken_report \
    $KRAK/${PRSA}.bracken" >> $MAT/ommfire__0__samp_complete.tsv
  
  
else
  echo " < ! >   something went wrong: ======================="
  for i in $FILT/${PRSA}_trimm_R1.fastq.gz \
    $FILT/${PRSA}_trimm_R2.fastq.gz \
    $DECO/${PRSA}_bt2decon_R1.fastq.gz \
    $DECO/${PRSA}_bt2decon_R2.fastq.gz \
    $KRAK/${PRSA}_kraken_report \
    $KRAK/${PRSA}.bracken ;
  do
    if [[ ! -s $i ]]; then echo " < ! >   ${i##*/} is missing   --------------" ; fi
  done
  echo " \ ^ /   ______________________________________________________________"
  echo " < ! >   ERROR :: do not pass go, do not collect 200 slurm"
  echo " / v \   --------------------------------------------------------------"
  # exit 1
fi  


# :) 
exit 0
