#!/bin/bash
    # set -e
    # set -u
    # set -o pipefail
    
    
  ## 0.4 - feb 2024  -  need checkpoints!! see wbartonisms
  ## 0.5 - dec 2024  -  remove hallucinations & dross for AMUKHERJEE. Same as mmoverlord__chunk__over0.5.sh
  

  ## consider
  # terminate if 1; no run on unset VAR; allow finishing on nonzero(?)
  # https://data-skills.github.io/unix-and-bash/03-bash-scripts/index.html
  # use =() to set RAWF to an array, then count ; echo
  # see https://mywiki.wooledge.org/BashFAQ/004


## ** NB ** :: assume that wherever you're at thats where you're from. Change if you see fit.
OBIN=$( dirname $0 )


##  ____________________________________________________________________________
##   S E T U P   ===============================================================

# fhi__omm_starter.sh $RAWFQGLOB OUTDIR NSAMPS
RAWF=$1   # raw fastq: a glob for matching some/not others, or a raw dir
ODIR=$2   # desired project output dir
RUNNSAMP=3 # $3  n samples in parallel  -  potentially disastrous - using "sbatch -W", and multiply by N tasks


echo " + + +   RAWF = $RAWF
 + + +   ODIR = $ODIR"  


## check a) that the list of globbed fastq > 0,  OR that it at least ends with "fastq.gz"
if [[ ! $( shopt -s nullglob ; FQGZ=($RAWF) ; echo ${#FQGZ[@]} ) -gt 0 ]] || [[ ! $( echo $RAWF | sed -r "s/.*(fastq.gz)/\1/" ) = fastq.gz ]] ;
then 
  echo " < ! >    inputs failed first check - arg1 doesn't glob to any *fastq.gz"
  # exit 1
fi


# fire out
PROJ=${ODIR##*/}
ROOT=${ODIR%/*}
RAWDIR=${RAWF%/*}
mkdir -p $ODIR

# references
REF=$ROOT/ref       # !!! echo this, as hardcoded
BT2_DB=$REF/bt2

# materiel
MAT=${ODIR}/Materials
mkdir -p $MAT

# file listing samples to run on
find $RAWF | grep 'fastq.gz' > $MAT/ommfire__fastq_files.tsv
find $RAWF | grep 'fastq.gz' | sed -r "s/.*\/(.*)_R[12].*\.fastq\.gz/${PROJ}__\1/g" | sort | uniq > $MAT/ommfire__samples.tsv


echo "PROJ $PROJ
ODIR $ODIR
RAWDIR $RAWDIR
ROOT $ROOT
MAT $MAT
REF $REF
BT2_DB $BT2_DB" > $MAT/ommfire__0__vars.tsv


# reference seqs for trimmomatic
echo '>Ampli_Tru_Seq_adapter__MOST_IMPORTANT_FEEEL_ME
CTGTCTCTTATACACATCT
>Transposase_Adap__for_tagmentation_1
TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_2
GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
>PCR_primer_index_1
CAAGCAGAAGACGGCATACGAGATNNNNNNNGTCTCGTGGGCTCGG
>PCR_primer_index_2
AATGATACGGCGACCACCGAGATCTACACNNNNNTCGTCGGCAGCGTC
>TruSeq_single_index_LT_CD_HT__1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq_single_index_LT_CD_HT__2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
>PCR-Free_Prep__Tagm__additional_seq
ATGTGTATAAGAGACA
>Ampli_Tru_Seq_adapter__MOST_IMPORTANT_FEEEL_ME__mate_RC
AGATGTGTATAAGAGACAG
>Transposase_Adap__for_tagmentation_1_RC:
CTGTCTCTTATACACATCTGACGCTGCCGACGA
>Transposase_Adap__for_tagmentation_2_RC:
CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
>PCR_primer_index_1_RC:
CCGAGCCCACGAGACNNNNNNNATCTCGTATGCCGTCTTCTGCTTG
>PCR_primer_index_2_RC:
GACGCTGCCGACGANNNNNGTGTAGATCTCGGTGGTCGCCGTATCATT
>TruSeq_single_index_LT_CD_HT__1_RC:
TGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
>TruSeq_single_index_LT_CD_HT__2_RC:
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>PCR-Free_Prep__Tagm__additional_seq_RC:
TGTCTCTTATACACAT
' > $MAT/fqc_trimmo_ill_ref.fa


##  ____________________________________________________________________________
##   c h e c k s   =============================================================

# flush_mat() {
#     echo " \ ^ /   ______________________________________________________________"
#     echo " < ! >   interrupted: try restarting from the caught samples"
#     echo " / v \   --------------------------------------------------------------"
# }
# trap "flush_mat;# exit 1" SIGINT SIGTERM # SIGERR   not sure if thats helping or no


if [[ -z $PROJ ]] || [[ ! -d $ODIR ]] || [[ ! -d $ROOT ]] || [[ ! -d $RAWDIR ]] ; then
    echo " < ! >   input-args (project, work dir, raw fastq dir) not sane - check"
    exit 1
fi


if [[ ! -d $REF ]] || [[ ! -d $BT2_DB ]] ; then
    echo " < ! >   ref-args (\$HGR_BT, \$BT2_DB, \$REF) not sane - check"
    exit 1
fi


if [[ ! -s $MAT/ommfire__fastq_files.tsv ]] ; then
    echo " < ! >   \$MAT/ommfire__0__vars.tsv not found - check"
    exit 1
fi


if [[ ! -s $MAT/ommfire__samples.tsv ]] || [[ -z $MAT/ommfire__samples.tsv ]]; then
    echo " < ! >   \$MAT/ommfire__samples.tsv empty/not found - check"
    exit 1
fi


if  [ ! $( echo $( wc -l $MAT/ommfire__fastq_files.tsv | cut -f 1 -d ' ' )/2 | bc ) -eq $( wc -l $MAT/ommfire__samples.tsv | cut -f 1 -d ' ' ) ] ; then
    echo " < ! >   n_samp != n_fastq/2  -  check \$MAT manifests"
    exit 1
fi  

for i in $( tail -6 $MAT/ommfire__0__vars.tsv | cut -f 2 -d ' ' ) ;
do
  if [ -z $i ] || [ ! -d $i ] ; then
    echo " < ! >   var \$$i not found - check (terrible check - this will be empty if VAR is empty)"
    exit 1
  fi
done


##  ____________________________________________________________________________
##   l o g   s e t u p   =======================================================

echo "
 +-+-+ 
 + + +   reference data expected in:   -----------------------------
 + + +   $REF
 + + +   setting parameters as below: view them in: ----------------
 + + +   $MAT 
 +++++   -----------------------------------------------------------
 +-+-+
$( sed -E 's/^/ + + + > /' $MAT/ommfire__0__vars.tsv)
 + + + 
 +-+-+ 
 +++++   -----------------------------------------------------------
 +-+-+ "



##  ____________________________________________________________________________
##   N O W   R U N  !  =========================================================
  
    # 3 scripts per sample, outputs in dirs specified above

      
##  ____________________________________________________________________________
##   f q c   /   m q c    ======================================================

    # use the sample list through $( cat fq.txt ) subshell
      sbatch -W $OBIN/omm_fqc20.sh \
        $MAT/ommfire__fastq_files.tsv \
        $ODIR/1__qc/${PROJ}_raw # INDIR OUTDIR

      wait   # suspect this is redundant, considering "sbatch -W" above


##  ____________________________________________________________________________
##  s u r v e y o . R   ::  get overall trimming spots    ======================

      # R 4.0 env - see surveyo.R for libs - just stringr and reshape2
      # source /data/Food/analysis/R0602_microsupport/jamie.fitzgerald/bin/miniforge3/bin/activate jfg_r4
      Rscript $OBIN/omm_surveyo.R $ODIR/1__qc/${PROJ}_raw_multi/multi_fwd/multiqc_report.html &&
      Rscript $OBIN/omm_surveyo.R $ODIR/1__qc/${PROJ}_raw_multi/multi_rev/multiqc_report.html
      conda deactivate


##  ____________________________________________________________________________
##   q u i c k w a s h   ::  per sample      ===================================

    ### 0.4+  -  use parallel -W to manage +1 samples at a time  
    module load parallel
    
    sed -r "s/.*${PROJ}__//" $ODIR/Materials/ommfire__samples.tsv | \
      parallel -j $RUNNSAMP sbatch -W $OBIN/omm_quickwa.sh {} $ODIR $RAWDIR    #   NSAMP   SCRIPT   SAMPLE   BASEDIR   RAWDIR


      echo " +++++   completed $PROJ run for  $RUNSAMP :  ===============

      ## consider adding run stats when time lotto


