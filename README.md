
# about

`mmo / micromatrix overlord` is a custom pipe to consistently process large volumes of microbiome `fastq` data from different projects, and congeal them to some sort of unified output. It won't rationalise bizarre and discordant metadata, but it will try and stabilise data quality to give comparable microbiome communities across runs/projects. 

Note that as the pipeline is designed to standardise the quality parameters across runs, different runs should be called to `overlord` separately (i.e. **spawn more overlords**). 


## structure

Given:

 - `/a/glob/for/all/*fastq.gz`	# a path to your `fastq.gz` ; should work with a dir too, but
 - `./omm_overlord.sh`		# setup & run 
 - `./omm_fqc20.sh`		# `FastQC / MultiQC `script
 - `./omm_surveyo.R`		# `R` script for consistent `trimmomatic` decisions
 - `./omm_quickwa.sh`		# runs `trimmomatic, bowtie2, kraken2`
 - `output/name`			# path and basename for output 
 
<p>

 - `./omm_humani.sh`  # separate, accessory code to concatenate R1+R2 and run `HUMANN3` 


## operation

When you run `omm_overlord.sh`, it:  

 - sets up folders, parameters, QC file
 - calls fqc20 to get quality data
 - calls surveyor to pick trimming lengths
 - calls quickwash to filter, trimm, remove host, and run kraken2
 - should put out a load of annoying comments, so consider sinking to a logfile

There is also a secondary workflow for running `HUMAnN3`, which tries to account for having pre-processed your `fastq` (for consistency). It's not a script (yet), but is consistent with `overlord`. `HUMAnN3` can take **>10 hours per sample and uses a lot of resources** - thus `overlord` does not run it automatically.


### pooling runs for analysis

The above steps allow multiple separate sequencing runs to be combined, even if the quality & QC requirements differ. All the `kraken2/Bracken` outputs generated are pooled, and combined with an `XLSX` file of sample sheets prepared manually. 
Two post-processing scripts then anneal the run metadata with the microbiome output, and analyse the microbial communities hidden within the `FASTQ`:

#### `omm_post_spreadsheeting`:

 - load & format `kraken2/Bracken` abundance estimates/taxonomy into feature/hierarchy tables
 - load metadata (prepared manually in `XLSX`) for all shotgun samples
 - create alpha and beta diversity metrics across dataset
 - perform centre-log ratio (CLR) standardisation of microbiome counts, with count-zero multiplicative (CZM) treatment of zero-values
 - pick some basic-bacillus colours, theme options for `R`
 
#### `omm_post_strainanalysis.R`:

 - source spreadsheeting above, thereby loading required data (samples, metadata)
 - using ``dbug` variable, subset data to only samples (test strain) of interest, plus all FFM / RSM-ve samples from those runs
 - plot relative abundance, alpha & beta diversity
 - run differential testing (linear mixed-effect model, with _post hoc_ estimated marginal means (LMM-EMM)), plot outcomes
 - save all output plots, tables, and workspace as RDS, RDS, and R Workspaces respectively.
 
These saved plots/tables/workspaces allow the user to store/access specific parts of the workflow later. 


## by example

Note that the `globbable` path (e.g. `/home/user/downloads/project/project*fastq.gz`) is **quoted**.

```
# presuming everything in your path and executable... 

./omm_overlord.sh \
	"/a/globbable/path/to/all/relevant/*fastq.gz" \
	output/name \
	> ./otuput.log \
	2> ./output.err &
	
```

Formatting of metadata and unifying of all `Bracken` outputs is up to the user. Only the post processing scripts are run, something like this would check them and save all the figures (or tables, etc.):

```
# in R 

# the data will have already been subset using debug in overlord_post_strainanalysis.R
debug 

# this RDS file is created by overlord_post_strainanalysis.R
if( file.exists( paste0("output/mm_", debug, "/fhi__mmoverlord__", debug, "__plots.RDS") )    ){
  print( debug)
  omm_plots <- readRDS( paste0("output/mm_", debug, "/fhi__mmoverlord__", debug, "__plots.RDS") )
  # ggsave( omm_plots[[1]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__relab_full.png"), width = 12, height = 10)
  # ggsave( omm_plots[[2]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__relab_cool.png"), width = 12, height = 6)
  ggsave( omm_plots[[3]] & 
            theme(
              axis.text.x = element_text( size = 14),
              legend.text = element_text( face = "italic", size = 14),
              NULL
            ) & 
            guides(fill = guide_legend(ncol =2)), 
          file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__relab_pool.png"), width = 12, height = 10)
  # ggsave( omm_plots[[4]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__alpha.png"), width = 10, height = 6)
  # ggsave( omm_plots[[5]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__beta-full-CCA.png"), width = 16, height = 12)
  # ggsave( omm_plots[[6]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__beta-48-CCA.png"), width = 16, height = 12)
  # ggsave( omm_plots[[7]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__beta-48-CCA.png"), width = 16, height = 12)
  # ggsave( omm_plots[[8]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__LMM-EMM_", debug, "-RSMve.png"), width = 8, height = 12.5)
  # ggsave( omm_plots[[9]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__LMM-EMM_", debug, "-FFM.png"), width = 8, height = 20)
  # ggsave( omm_plots[[10]], file = paste0("vis/mm_", debug, "/fhi__mmoverlord__", debug, "__LMM-EMM_FFM-RSMve.png"), width = 8, height = 25)
  
}else{ 
  paste0( "no plots for you - ", debug)
}
})

```


<br/>


### apologies 

Sorry for the generally amateur code, and the bugs that are no doubt hiding herein. Sorry also about the (lack of) naming convention in the scripts etc..


### dependencies / catastrophic assumptions

 - assumes that you have many, many things **installed**:
 - - standard :: `mamba create -n mmo -c bioconda multiqc fastqc trimmomatic bowtie2 kraken2 -y`
 - - HUMAniNg :: `mamba create -n hum -c biobakery humann -y`
 - - R        :: `install.packages(c("ggplot2", "patchwork", "vegan", "ggvegan", "lme4", "emmeans", "parallel"))`
 - assumes that you've built the **databases** for `kraken2`, `bracken`, and `bowtie2`! (see params echo'd to `$MAT/ommfire__0__vars.tsv` by `overlord.sh`, line 62)
 - assumes that all four scripts are in, and run from, the same folder - to change this behaviour, modify the `OBIN` parameter when set at the start of `omm_overlord.sh` (line 19 of `overlord.sh`).
 - makes _very_ specific assumptions about the naming convention of your `fastq` - specifically, that they match `some/path/file_R*.fastq.gz` (line 59 of `overlord.sh`). 

My databases as used here are located at:

 - `bt2`  : `/../analysis/R***/user/ref/bt2/hum_deco`
 - `kr2`  : `/../analysis/R***/user/ref/kraken..`, probably
 - `brkn` : `/../analysis/R***/user/ref/kraken...`, probably also...

My conda envs are at `/../analysis/R***/user/bin/miniforge3/envs/`, likely 'twas the `base` env used for the above (though thats inadvisable).
 
Last tested on Teagasc HPC-1.



##### jfg dec 2024
