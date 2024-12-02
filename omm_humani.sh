
## same as fhi__mmo__humann_0.2, to AMUKHERJEE


## setup  ===================

RAW=    # no  -  raw  -  till dir ln
WRK=/data/Food/analysis/R0602_microsupport/jamie.fitzgerald/fhi__omm_SC40
REF=$MSDAT/ref
H3REF=$REF/h3_ref
MAT=$WRK/Materials


# outputs etc
QC=$WRK/1__qc        # QC data
FILT=$WRK/2__filt    # fitlered & trimmed (trimmo)
DECO=$WRK/3__deco    # host-removed (BT2)
H3=$WRK/7__h3        # aspirant output

#!# set a sample-name to test things on
TEST=fhi__omm_SC40__FHI101_S67

mkdir $H3
mkdir -p $H3REF



## excised  ==========================================================

    # from v 0.1 - messing, checking, wrangling
    

## humann setup  ==========================================================
      # 
      # screen -R hum
      # module load humann/3.8 
      # 
      # # check - ok
      # humann_test
      # 
      # # demo run - fail
      # find /install/software/restart/py3/humann_3.5/* | grep demo.fastq | less
      # # > [...]/pkgs/humann-3.6-pyh7cba7a3_0/site-packages/humann/tests/data/demo.fastq
      # H3DEMO=/install/software/restart/py3/humann_3.5/pkgs/humann-3.6-pyh7cba7a3_0/site-packages/humann/tests/data
      # humann -i $H3DEMO/demo.fastq -o $H3/sample_results
      # 
      # ## CRITICAL ERROR: The directory provided for ChocoPhlAn contains files
      #   # - ( g__Phaeobacter.s__Phaeobacter_sp_CECT_7735.centroids.v296_v201901b.ffn.gz ) that are not
      #   # - of the expected version. Please install the latest version of the database: v201901_v31
      # 
      # ## install and update database locations
      # humann_databases --download chocophlan full $H3REF
      # humann_databases --download uniref uniref90_diamond $H3REF
      # 
      # ## don't have the permissions to update this (where is that cached?!) so specify at runtime
      # # humann_config --update database_folders nucleotide $H3REF
      # # humann_config --update database_folders protein /data/databases/Humann3/uniref


## H3 slurm  =============================================

echo '#!/bin/bash

#SBATCH --job-name=humanatee3
#SBATCH --output=humanatee3.%j.slrm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=400:00

SAMPLE=$1
IN=$2
OUT=$3
DB1=$4
DB2=$5

NTHREADS=10 # as above so below
module load humann/3.8

# echo "## Sys.time"

# -----------------------------------------------------

if [ ! -f $IN/${SAMPLE}*_bt2decon_R1.fastq.gz ] | [ ! -f $IN/${SAMPLE}*_bt2decon_R2.fastq.gz ]
then
  echo "## UH-OH :: inputs not found  -  stopping slurm_hum.sh" ; exit 1
fi


# shout  ----------------------------------------------

echo "## HUMAnN3 on ${SAMPLE}:
echo "##    - -  :: F.f.gz - $( ls $IN/${SAMPLE}*R1.fastq.gz )"
echo "##    - -  :: R.f.gz - $( ls $IN/${SAMPLE}*R2.fastq.gz )"
echo "##    - -  :: choco  - ${DB1}"
echo "##    - -  :: uniref - ${DB2}"
echo "##    - -  :: output - $OUT"


# operate samplewise, referencing list   -------------

if [ -s $IN/${SAMPLE}_bt2decon_R1R2.fastq.gz ]
then
  echo "##  - -  :: fastq input for ${SAMPLE} R1+R2 already catd in ${OUT} - skipping cat step"
else
  echo "##  - -  :: cat ${SAMPLE} R1+R2 -> ${SAMPLE}_bt2decon_R1R2.fastq"
  zcat $IN/${SAMPLE}*_bt2decon_R1.fastq.gz $IN/${SAMPLE}*_bt2decon_R2.fastq.gz > $IN/${SAMPLE}_bt2decon_R1R2.fastq
fi


if [ -s $OUT/${SAMPLE}_hum3/${SAMPLE}*_genefamilies.tsv ]
then
  echo "##    - -  :: Humann3 gene_fam_output.tsv for ${SAMPLE} already in ${OUT} - stopping"
else
  echo "##    - -  :: firing Humann3 on ${SAMPLE}"


  time humann \
  -i $IN/${SAMPLE}_bt2decon_R1R2.fastq \
  -o $OUT/${SAMPLE}_hum3 \
  --threads $NTHREADS \
  --evalue 1.0 \
  --output-basename ${SAMPLE}_hum3 \
  --prescreen-threshold 0.01 \
  --nucleotide-database $DB1 \
  --protein-database $DB2

  # something wrong in call to tar
  if [ $? == 0 ] && [ -s $OUT/${SAMPLE}_hum3/${SAMPLE}*_genefamilies.tsv ]
  then
    rm $IN/${SAMPLE}_bt2decon_R1R2.fastq
    tar -czf $OUT/${SAMPLE}_hum3/${SAMPLE}_hum3_temp.tar.gz --remove-files $OUT/${SAMPLE}_hum3/${SAMPLE}_hum3_humann_temp &
    echo "## "
    echo "## normalising to counts-per-million..."
    humann_renorm_table --input $OUT/${SAMPLE}_hum3/${SAMPLE}_hum3_genefamilies.tsv \
    --output $OUT/${SAMPLE}_hum3/${SAMPLE}_hum3_genefamilies-cpm.tsv \
    --units cpm --update-snames
    echo "## "
    echo "## "
    echo "## HUMAnN3 on ${SAMPLE}   +++ complete +++ "
  else
    echo "<!>   HUMAnN3 exit code != 0  /  no genefam table - ERROR! <!> "
    # exit 1
  fi


fi' > $MAT/slurm_hum3.sh  # SAMPLE IN OUT CHOCDB DIAMDB


## total set =====================================================
  
$TEST
$DECO
$H3

## CHOC and DIAM set above - reenter if necessary
CHOCDB=$H3REF/chocophlan
DIAMDB=/data/databases/humann3/manualdownload/uniref90 
file $CHOCDB
file $DIAMDB


# check - 11:30 - 14:05 - 2 hours
# sbatch $MAT/slurm_hum3.sh  $TEST $DECO $H3 $CHOCDB $DIAMDB

# where $MAT/*subdat.tsv is a list of file (sample) names
for S in $( cut  $MAT/*subdat.tsv -f 2 | sed -r 's/\"//g' | sed -r 's/^e[0-9]*_//g' ) ;
do
  sbatch -w $MAT/slurm_hum3.sh  fhi__omm_SC40_${S} $DECO $H3 $CHOCDB $DIAMDB
done > less -S 


## get outputs   ---------------------------------------------------------------

  cat $MAT/samples_doneH3
  module load humann3/3.6 parallel


  ## EC PWY are the def\ult output for H3, though shoulfd renorm (CPM default)
  
  ## get entire mapping set! 
  # humann_databases --download utility_mapping full $REF/
  #   # map_ec_name.txt.gz          map_go_name.txt.gz      map_ko_uniref50.txt.gz        map_pfam_name.txt.gz       map_uniref50_uniref90.txt.gz
  #   # map_eggnog_name.txt.gz      map_go_uniref50.txt.gz  map_ko_uniref90.txt.gz        map_pfam_uniref50.txt.gz   map_uniref90_name.txt.bz2
  #   # map_eggnog_uniref50.txt.gz  map_go_uniref90.txt.gz  map_level4ec_uniref50.txt.gz  map_pfam_uniref90.txt.gz   uniref50-tol-lca.dat.bz2
  #   # map_eggnog_uniref90.txt.gz  map_ko_name.txt.gz      map_level4ec_uniref90.txt.gz  map_uniref50_name.txt.bz2  uniref90-tol-lca.dat.bz2


  ## normalise outputs 
  
  # - GENEFAM cpm
    cat $MAT/samples_doneH3 | parallel -j 25 humann_renorm_table \
    --input $H3/{}_hum3/{}_hum3_genefamilies.tsv \
    --output $H3/{}_hum3/{}_hum3_genefamilies-cpm.tsv \
    --units cpm --update-snames > $MAT/genefam_cpm.log

  # - PATHWAY cpm
    cat $MAT/samples_doneH3 | parallel -j 25 humann_renorm_table \
    --input $H3/{}_hum3/{}_hum3_pathabundance.tsv \
    --output $H3/{}_hum3/{}_hum3_pathabundance-cpm.tsv \
    --units cpm --update-snames > $MAT/genefam_pathw_cpm.log


  # regroup gene abundances to metabolic RXNs - MetaCyc

    cat $MAT/samples_doneH3 | parallel -j 25 humann_regroup_table \
    --input $H3/{}_hum3/{}_hum3_genefamilies-cpm.tsv \
    --output $H3/{}_hum3/{}_hum3_genefamilies-MetaCyc.tsv --groups uniref90_rxn > $MAT/genefam_RXNs.log

  # regroup gene abundances to ECs- MetaCyc - for Deyan_55, 11-12%

    cat $MAT/samples_doneH3 | parallel -j 25 humann_regroup_table \
    --input $H3/{}_hum3/{}_hum3_genefamilies-cpm.tsv \
    --output $H3/{}_hum3/{}_hum3_genefamilies-MetaCyc-EC.tsv \
    -c $REF/utility_mapping/map_level4ec_uniref90.txt.gz > $MAT/genefam_ECs.log


  # regroup GENE ABUNDANCES to Protein FAmilies - Pfam! - for Dyen_55, 50-65%

    cat $MAT/samples_doneH3 | parallel -j 25 humann_regroup_table \
    --input $H3/{}_hum3/{}_hum3_genefamilies-cpm.tsv \
    --output $H3/{}_hum3/{}_hum3_genefamilies-pfam.tsv \
    -c $REF/utility_mapping/map_pfam_uniref90.txt.gz  > $MAT/genefam_PFAMs.log


  ## join tables

  # group PATHWAY abundances
  
    mkdir $MAT/fhi__mmoverlord__H3_pathabundance-cpm ; cp $H3/*/*_hum3_pathabundance-cpm.tsv $MAT/fhi__mmoverlord__H3_pathabundance-cpm
    humann_join_tables -i $MAT/fhi__mmoverlord__H3_pathabundance-cpm -o $MAT/fhi__mmoverlord_pathabundance-cpm.tsv --file_name pathabundance-cpm


  # group GENE abundances
  
    mkdir $MAT/fhi__mmoverlord__H3_genefam-cpm ; cp $H3/*/*_hum3_genefamilies-cpm.tsv $MAT/fhi__mmoverlord__H3_genefam-cpm
    humann_join_tables -i $MAT/fhi__mmoverlord__H3_genefam-cpm -o $MAT/fhi__mmoverlord_genefam-cpm.tsv --file_name genefamilies-cpm

    mkdir $MAT/fhi__mmoverlord__H3_MetaCyc ; cp $H3/*/*_hum3_genefamilies-MetaCyc.tsv $MAT/fhi__mmoverlord__H3_MetaCyc
    humann_join_tables -i $MAT/fhi__mmoverlord__H3_MetaCyc -o $MAT/fhi__mmoverlord_MetaCyc.tsv --file_name MetaCyc

    mkdir $MAT/fhi__mmoverlord__H3_MetaCyc-EC ; cp $H3/*/*_hum3_genefamilies-MetaCyc-EC.tsv $MAT/fhi__mmoverlord__H3_MetaCyc-EC
    humann_join_tables -i $MAT/fhi__mmoverlord__H3_MetaCyc-EC -o $MAT/fhi__mmoverlord_MetaCyc-EC.tsv --file_name MetaCyc-EC

    mkdir $MAT/fhi__mmoverlord__H3_pfam ; cp $H3/*/*_hum3_genefamilies-pfam.tsv $MAT/fhi__mmoverlord__H3_pfam
    humann_join_tables -i $MAT/fhi__mmoverlord__H3_pfam -o $MAT/fhi__mmoverlord_pfam.tsv --file_name pfam


##  j u n k   =======================================================================================================================

  ### count things
  # if [ -s $OUT/${SAMPLE}*_kaiju ]
  # then
  #   KAI_C=$(grep -cE '^C' $OUT/${SAMPLE}_kaiju ) &
  #   KAI_U=$(grep -cE '^U' $OUT/${SAMPLE}_kaiju ) 
  #   KAI_TOT=$(echo "$KAI_C + $KAI_U" | bc)
  #   KAI_PC=$(echo "scale=2; ($KAI_C / $KAI_TOT)*100" | bc)
  #   echo "## kaiju sample processed: ${SAMPLE} : total classified: ${KAI_PC}% (total: $KAI_TOT read-pairs)  ------"
  # else
  #   echo "## no output - kaiju for sample ${SAMPLE} failed"
  # fi