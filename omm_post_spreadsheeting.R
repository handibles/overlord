

## fhi__mmoverlord__spreadsheeting :: copies of copies of copies of copies of OVERLORDS OF

      # rm(list=ls())    # View( mgdat )
      set.seed(2205)


## libs  =========

      library(vegan)
      library(ggvegan)
      library(ggplot2)
      library(patchwork)


#   only run once   ====================================================================================================================
# 
#         # Kraken2 - this is the pooled output across all available runs
#             krak2_brack_file="input/fhi__omm__krakenStnd_abundances__Feb2024.tsv"
#             krak2_mpa_file="input/fhi__omm__krakenStnd_taxonomy__Feb2024.tsv"
# 
# 
#           ## Re-arrange data   -----------------------------------
# 
#             k2dat <- read.table(file = krak2_brack_file, header=TRUE, quote = "", sep = '\t')
#             k2dat[1:10, 1:10]
# 
#             prenames_k2dat <- paste0("k2_", stringr::str_pad(k2dat[,2], width=7, pad=0))
#             ## non-unique names -  assume only one dupe of each
#             prenames_k2dat[ duplicated(k2dat[,2])] <- gsub('k2_', 'k2-A_', prenames_k2dat[ duplicated(k2dat[,2])])
#             rownames(k2dat) <- prenames_k2dat
# 
#             # counts
#             k2_counts <- k2dat[ , grep('_num', colnames(k2dat))]
#             colnames(k2_counts) <- gsub( "_S\\d*.bracken_num", "", colnames(k2_counts), perl = TRUE )
# 
#             # relab%
#             k2_relab <- k2dat[ , grep('_frac', colnames(k2dat))]
#             colnames(k2_relab) <- gsub( "_S\\d*.bracken_frac", "", colnames(k2_relab), perl = TRUE )
# 
# 
#           ## start with data from the abundance table names  --------------------------------
# 
#             tax_a <- data.frame(taxon=k2dat[ ,1], k2_id=k2dat[ ,2], row.names = rownames(k2dat) , stringsAsFactors = FALSE)
#             tax_a[,1] <- gsub('\\/','', tax_a[,1])
#             tax_a[,1] <- gsub('ALOs_','ALOs-', tax_a[,1])
#             head(tax_a)
# 
#             ## NOTE use quote="" to avoid nightmare of special characters like '"`/,| etc in names
#             tax_b <- read.table( krak2_mpa_file, sep='\t', header=FALSE, fill = TRUE, stringsAsFactors = FALSE, quote="")
#             head(tax_b)
# 
# 
#           ## regularise taxonomic ranks (not all have K:P:C:O:F:G:S etc.)   -----------------
# 
#             ranks <- c("k__", "p__", "c__", "o__" , "f__", "g__", "s__")
#             tax_b <- t(apply(tax_b, 1, function(aa){  # aa <- tax_b[1,]
#               sapply(ranks, function(aaa){ # aaa <- "c__"
#                 bbb <- grep(aaa, unlist(aa), value=TRUE)
#                 ifelse( length(bbb) == 0, "unkn.", bbb)
#               })
#             }))
#             dim(tax_b) ; str(tax_b)
# 
# 
#           ## more tax issues due to names: doubled ranks:
# 
#             prior_index <- (unlist( lapply(1:nrow(tax_b), function(aa){ if( sum( grepl("unkn.", tax_b[aa,])) == 6){aa} }) ) - 1)
#             tax_b[ prior_index , "s__"] <- tax_b[ (prior_index+1) , "s__"]
#             # row is only useful if has a species: cannot have a species, and hit 6, And be worth keeping
#             tax_b <- tax_b[ !apply(tax_b, 1, function(aa) sum(grepl("unkn.", unlist(aa)) ) == 6) , ]
# 
#             colnames(tax_b) <- c("d", "p", "c", "o", "f", "g", "s")
# 
#             tax_a[ , "taxon"] <- gsub(" ", "_", paste0("s__", tax_a[ , "taxon"]))
#             head(tax_a) ; head(tax_a[ , "taxon"]) ;  dim(tax_a)
#             head(tax_b) ;   head(tax_b[ , "s"]) ;  dim(tax_b)
# 
# 
#           ## stitch  &  name   ----------------------------------
# 
#             tax_c <- merge(tax_b, tax_a, by.x="s", by.y="taxon")[, c(2:7,1,8)]
#             head(tax_c) ; dim(tax_c)
#             prenames_tax <- paste0("k2_", stringr::str_pad(tax_c[,"k2_id"], width=7, pad=0))
#             rownames(tax_c) <- prenames_tax
# 
# 
#           ## check all on the same page   -----------------------
# 
#             all(rownames(k2dat) == rownames(k2_counts) & rownames(k2dat) == rownames(k2_relab))
#             all(rownames(k2dat) %in% rownames(tax_c))
# 
#             kraken_ids <- sort(rownames(k2dat))
#             k2dat <- k2dat [kraken_ids , ]
#             k2_counts <- k2_counts [kraken_ids , ]
#             k2_relab <- k2_relab [kraken_ids , ]
#             tax_c <- tax_c[kraken_ids , ]
# 
# 
#           ## remove if no values (kraken2 set to give all taxa)
# 
#             pos_abund <- rownames(k2_relab)[rowSums(k2_relab) > 0]
#             pos_abund_c <- colnames(k2_relab)[colSums(k2_relab) > 0]
#             k2_relab <- k2_relab[ pos_abund , pos_abund_c ]
#             k2_counts <- k2_counts[ pos_abund , pos_abund_c ]
#             k2_tax <- tax_c[ pos_abund , ]
# 
# 
#       ## add metadata from spreadsheet   --------------------------------
# 
#         ##   from xls-process-2 (v.3)  ----------------------------------------------------------
#             ( fhidat <- readRDS("output/fhi__mmoverlord__FHIdat_1611-8.RDS") )[1:10,]
# 
#           # match?
#             dim(fhidat)
#             dim(k2_counts)
#             all( rownames(fhidat) %in% colnames(k2_counts))
#             all( colnames(k2_counts) %in% rownames(fhidat))
# 
# 
#           ## fixes    -------------------
# 
#             fhidat$desc <- gsub("^T\\d*", "", fhidat$desc, perl = TRUE)
#             fhidat$desc <- gsub("Ha119", "HA119", fhidat$desc, perl = TRUE)
#             fhidat$desc <- gsub("Sc6A", "SC6A", fhidat$desc, perl = TRUE)
# 
#             substrates <- unique( c("HA136", "LH88","HA119", "R0033Anaerobic", "R0033Aerobic",
#                             "SC40", "HA179","HA110","R1058", "HA118",
#                             "ST765","ST906","SC6","R0389",
#                             "SC17","DPC1817","DHA118",
#                             "P200622","R0070Anaerobic","R0070Aerobic","R0419","FFM",
#                             "HA110","RSM-ve"),
#                             # not included elsewhere
#                             c("HA108", "SC6A","R0215","R0418","R0422","SC15","SC18","SC76","SC71","SL119")
#             )
# 
#           ## retroactive - from email, and PROTEAN-HEAD  : the top 19 things
#             # length(strains_19 <- c("HA136", "DHA118", "DPC1817", "HA110", "HA119", "HA179",
#             #                        "LH88", "P200622", "R0033Aerobic", "R0033Anaerobic", "R0389",
#             #                        "R0419", "R1058", "SC17", "SC40", "SC6",
#             #                        "ST765", "ST906",
#             #                        "R0070Aerobic", "R0070Anaerobic",
#             #                        "HA118"))
#             # strains_19[ !(strains_19 %in% substrates) ]
#             # substrates[ !(substrates %in% strains_19) ]
# 
# 
#             fhidat$strain <- sapply( fhidat$desc, function(aa){   # aa <- "SC40T1SU"
#                 bb_v <- substrates[ sapply( substrates, function(aaa){ grepl( aaa, aa) }) ] # aaa <- "SC40
#                 if( length(bb_v) > 1){
#                   paste0( sort(bb_v), collapse = "/")
#                 }else if( length(bb_v) == 0 ){
#                   "none"
#                 }else{
#                   bb_v
#                 }
# 
#               })
#             fhidat[ grep("R0033Genome", fhidat$desc), "strain"] <- "R0033"
# 
#             fhidat$fsi <- sapply( fhidat$exp, function(aa){
#               bb_log <- as.numeric( gsub("e", "", aa))
#               if( is.na(bb_log)){
#                   "none"
#                 }else if( bb_log < 18){
#                   "FSI1"
#                 }else if( bb_log  %in% 18:23 ){
#                     "FSI2"
#                 }else if( bb_log %in% 24:26 ){
#                     "FSI3"
#                 }else if( bb_log > 26 ){
#                     "FSI4"
#                 } })
# 
# 
#           ## alpha   -------------------
# 
#             fhidat$invs <- vegan::diversity( t(k2_counts), index = "invs")
#             fhidat$shan <- vegan::diversity( t(k2_counts), index = "shan")
#             fhidat$rich <- apply( t(k2_counts), 1, function(aa){ sum(aa > 0) })
#             fhidat$chao1 <- apply( t(k2_counts), 1, function(aa){ vegan::estimateR(aa) })[2,]
#             fhidat$SeqDepth <- rowSums( t(k2_counts) )
# 
# 
#           ## dissims - here for rownames and dissim   ------------------
# 
#             k2_ra_bcdist <- as.dist(vegdist( t(k2_relab[ , rownames(fhidat)]), method = "bray"))
#             k2_ra_jacdist <- as.dist(vegdist( t(k2_relab[ , rownames(fhidat)]), method = "jacc"))
# 
#           ## revised route from paedad, after inexplciable difficulties getting order from old method...
#             fhidat$bc_wardD2 <- factor( rownames(fhidat), levels = colnames(k2_relab)[hclust( k2_ra_bcdist, method = "ward.D2")$order ]  )
#             fhidat$jacc_wardD2 <- factor( rownames(fhidat), levels = colnames(k2_relab)[hclust( k2_ra_jacdist, method = "ward.D2")$order ]  )
#             fhidat$bc_seriate <- factor( rownames(fhidat), levels = colnames(k2_relab)[
#                   seriation::get_order( seriation::seriate( k2_ra_bcdist, method = "TSP") )
#             ]  )
# 
#           ## check - do all match?   -----------------------
# 
#             dim(k2_counts)    # 730 1611  ~~711 1267
#             dim(k2_relab)     # 730 1611  ~~711 1267
#             dim(k2_tax)       # 730 8     ~~711  8
#             # dim(mgdat)        # 1611 7    ~~1267 6
#             dim(fhidat)        # 1611 11    ~~1267 6
# 
#             # look
#             k2_counts[1:10, 1:10]
#             k2_relab[1:10, 1:10]
# 
#             # View(k2_counts)
#             # hist( log10(k2_counts ))
# 
# 
#       ## transformations etc   ---------------------------------------------------
# 
#           # CLR    -------------
# 
#             ##  CZM/CLR Sanity check - omm checked   --------------------------
#                 ## from ?cmultRepl::
#                   # > "If method="user", user ... *** for each count vector (row) in X ***.
#                 ## zCompositions assumes samples are ~~COLUMNS~~ ROWS!!!
#                 # zPatterns(k2_counts, 0)
#                 # heatmap( apply(k2_counts[ , mgdat$bc_wardD2 ] > 0, 1, as.numeric), Rowv = NULL )
# 
# 
#               dim(k2_counts) ; dim(fhidat)
#                 # > [1]  730 1611
#                 # > [1] 1611   18
#             ## zCompositions assumes samples are ROWS!!!
#               k2_cmult <- zCompositions::cmultRepl( t(k2_counts), method = "CZM")
#               dim(k2_cmult)
#                 # > [1] 1611  730
#               k2_clr <- t( apply( k2_cmult, 1, function(aa){ log(aa) - mean(log(aa)) }))   # formula from vegan
#               dim(k2_clr)
#                 # > [1] 1611  730
# 
# 
# #       # ## write out   ---------------------------------------------------
# # 
# #           ## general metadata
# #             # dim(mgdat)
# #             dim(fhidat)        # 1611 11    ~~1267 6
# #             # saveRDS( mgdat, "output/fhi__mmoverlord__dat-1611-7.RDS")
# #             saveRDS( fhidat, "output/fhi__mmoverlord__FHIdat-1611-16.RDS")
# # 
# #           ## tsv output
# #             dir.create("output/tsv")
# #             # write.table(mgdat, "output/tsv/fhi__mmoverlord__dat-1611-7.tsv", sep = "\t")
# #             write.table(fhidat, "output/tsv/fhi__mmoverlord__dat-1611-16.tsv", sep = "\t")
# #             write.table(k2_counts, "output/tsv/fhi__mmoverlord__kraken.10.5.20.50-1611__num.tsv", sep = "\t")
# #             write.table(k2_relab, "output/tsv/fhi__mmoverlord__kraken.10.5.20.50-1611__frac.tsv", sep = "\t")
# #             write.table(k2_tax, "output/tsv/fhi__mmoverlord__kraken.10.5.20.50-1611__tax.tsv", sep = "\t")
# # 
# #           ## RDS for R output - note samples are ROWS
# #             saveRDS( t(k2_counts), "output/fhi__mmoverlord__kraken.10.5.20.50-1611__num.RDS")
# #             saveRDS( t(k2_relab), "output/fhi__mmoverlord__kraken.10.5.20.50-1611__frac.RDS")
# #             saveRDS( t(k2_ra_bcdist), "output/fhi__mmoverlord__kraken.10.5.20.50-1611__BCdist.RDS")
# #             saveRDS( t(k2_ra_jacdist), "output/fhi__mmoverlord__kraken.10.5.20.50-1611__Jaccdist.RDS")
# #             saveRDS( (k2_clr), "output/fhi__mmoverlord__kraken.10.5.20.50-1611__clr.RDS")
# #             saveRDS( k2_tax, "output/fhi__mmoverlord__kraken.10.5.20.50-1611__tax.RDS")
#        
      
## read in   ===================================================================
      
      # mgdat <- readRDS("output/fhi__mmoverlord__dat-1611-7.RDS")   # basicb
      # mgdat <- readRDS( "output/fhi__mmoverlord__unif_expdat-1611-10.RDS")
      print("  + + +   < ! >   l o a d i n g   t h e   F H I - d a t . . . .       + + +")
      mgdat <- readRDS( "output/fhi__mmoverlord__FHIdat-1611-16.RDS")
      
      mgfeat <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__num.RDS")
      mgfeat_ra <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__frac.RDS")
      mgtax <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__tax.RDS")
      
      mgfeat_clr <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__clr.RDS")
      mgfeat_dist_bray <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__BCdist.RDS")
      mgfeat_dist_jac <- readRDS("output/fhi__mmoverlord__kraken.10.5.20.50-1611__Jaccdist.RDS")
      
    # missed this - HA118, HA110 not properly represented
      mgdat$strain <- ifelse( mgdat$desc == "HA118", "HA118", mgdat$strain)
      mgdat$strain <- ifelse( mgdat$desc == "HA110", "HA110", mgdat$strain)
      
      
    ## plotting info  =========
      
      # colours, factors, whathavveyou  
      length(nr_col <- c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#FF1493",
                         "#b15928", "#737f33", "#8B008B", "#32fbd8", "#fdbf6f",
                         RColorBrewer::brewer.pal(9, "Spectral")[],
                         "#b2df8a", "#fb9a99", "#d9e627", "#EE82EE", "#DEB887",
                         "#a6cee3" ))

    ## you want colours?
      # library("vapoRwave")
      # 
      # length(unique(va_col <- sort(c(
      #   vapoRwave_pal()(n=12),
      #   seaPunk_pal()(n=12),
      #   sunSet_pal()(n=12),
      #   newRetro_pal()(n=12),
      #   mallSoft_pal()(n=12),
      #   macPlus_pal()(n=12),
      #   jwz_pal()(n=12)
      # ))))
      length(va_cols <- c("#20DE8B"
                          ,"#CCDE8B"
                          ,"#FFDE8B"
                          ,"#FFA88B"
                          ,"#FF6A8B"
                          ,"#FF6AD5"
                          ,"#C874AA"
                          ,"#C774E7"
                          ,"#AD8CFF"
                          ,"#966BFF"
                          ,"#90CFFF"
                          ,"#296656"
                          ,"#569874"
                          ,"#7EC488"
                          ,"#A997AB"
                          ,"#532E57"
                          ,"#F9897B"
                          ,"#D7509F"
                          ,"#F9247E"
                          ,"#AE1357"
                          ,"#661246"
                          ,"#9239F6"
                          ,"#903495"
                          ,"#6F3460"
                          ,"#4A354F"
                          ,"#D20076"
                          ,"#FF0076"
                          ,"#FF4373"
                          ,"#FF6B58"
                          ,"#F8B660"
                          ,"#7FD4C1"
                          ,"#30BFDD"
                          ,"#8690FF"
                          ,"#ACD0F4"
                          ,"#F7C0BB"
                          ,"#FBCFF3"
                          ,"#65323E"
                          ,"#FE7F9D"
                          ,"#FFC0CB"
                          ,"#75D8D5"
                          ,"#09979B"
                          ,"#063B41"
                          ,"#7B556C"
                          ,"#86486F"
                          ,"#F1956E"
                          ,"#EB9B60"
                          ,"#864A42"
                          ,"#FAD36A"
                          ,"#F5FFBF"
                          ,"#EFDB9C" ))
      
      
      ## theme variables    -------------------------------------------------------------------
      
      print("  + + +   applying theme_update in spreadsheeting_shortcuts   + + +    ")
      theme_set( theme_minimal())
      theme_update(
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.line = element_line(colour = "grey80", linewidth = 0.2),
        #
        legend.box = NULL,  # controls the direction of multiple legend boxes
        legend.position = "right",
        legend.text = element_text(size = 9),
        legend.key.size = unit(0.75, "cm"),
        #
        panel.border = element_rect(colour = "grey20", fill=NA, linewidth=0.3),
        panel.spacing.y = unit(2, "lines"),
        
        # -----
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        # -----
        
        # panel.grid.major.x = element_line(colour = "lightskyblue3", linewidth = 0.01),
        # panel.grid.major.y = element_line(colour = "lightskyblue3", linewidth = 0.01),
        # panel.grid.minor.x = element_line(colour = "lightskyblue3", linewidth = 0.01),
        # panel.grid.minor.y = element_line(colour = "lightskyblue3", linewidth = 0.01),
        panel.spacing.x = unit(1, "lines"),
        #
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),  # if plotting class using super
        plot.tag = element_text(size = 25, face = "bold", margin = margin()),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.title = element_text(size = 18, face = "bold", hjust = 0),  #  margin = margin(l = 10)),  # if plotting class using super
        plot.subtitle = element_text(size=14, hjust = 0.5),
        plot.background = element_rect(fill = "white", colour = "white"),
        #
        strip.background = element_blank(),
        strip.text.x = element_text(angle = 0, size = 16),   # Top
        strip.text.y = element_text(angle = 0, size =14),
        
        text = element_text("sans", colour = "grey20"),
        NULL
      )
      
      
##  F U N C T I O N S   ===============================
      
      ## k_A...
      k_A <- function(
        counts,     # count table, taxa are ROWS
        k=0,        # raw counts
        A=0         # prevalence
          ){
            # simplify
            if( any( is.na(counts)) ){   
              stop( "filthy NAs in your count table. Fix and try try again")
            }
            print("k_A: assuming samples are COLUMNS, and na.rm = TRUE")
            # parachute
            data <- counts     
            
            # to % as necessary
            if(k < 1 & any(colSums(counts)>1)){ counts <- apply(counts, 2, function(a) a/sum(a))}   
            
            if(A < 1){                                                                    # x <- counts[20,]
              filt <- apply(counts, 1, function(x) sum(x>k)/ncol(counts) > A )          # above 0.5% in 0.1 of samples
            }else if(A > 1){
              filt <- apply(counts, 1, function(x) sum(x > k) > A )                     # above 0.5% in 10 samples
            }else{
              filt <- apply(counts, 1, function(x) (sum(x)/sum(counts)) > k )           # above 0.5% across all samples
            }
            return(data[ filt, ])
            
        }
      
      
      ## colour the taxa in the manner of our forebearers
      shade_ranks2 <- function(...){
        print( paste0("   < ! >   don't use shade_ranks2\n   < ! >   use the reduced and transparent paragraph\n   < ! >   written for paedad:
      }
      
      
      uniq_g <- unique( gsub(\" .*\", \"\", fhi_ag$plot_var))
      cc_cols <- unlist(sapply( 1:length(uniq_g), function(aaa){  # aaa <- 9
        
        bb_s_in_g <- unique(
          mgtax_filt[ , \"s\" ][ grep( uniq_g[aaa], mgtax_filt[ , \"g\" ]) ]
        )
        
        cc_colours <- colorRampPalette( c(va_cols[aaa], \"white\") )( length(bb_s_in_g) + 3 )[1:length(bb_s_in_g)]
        names(cc_colours) <- gsub(\"_\", \" \", bb_s_in_g)
        cc_colours
        
      } ) )
      names(cc_cols)
      # scales::show_col(cc_cols)
      
      "))
      }

      