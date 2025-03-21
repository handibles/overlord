
## per - strain evaluation

    ## protean starts with an empty env
      # rm(list = ls())
      
      source("omm_post_spreadsheeting.R")
  
      # dbug set in the env above
      # if( !all(grepl( "dbug", ls())) ){
      #   stop( paste(" + + +   dbug:", dbug, " not set to anything"))
      # }
      
    ## defines the bug you look at - for example:
      dbug <- "R0033Aer"

    ## need to replace the / in strain names 
      mgdat$desc <- gsub("/", "-", mgdat$desc)
      mgdat$desc <- gsub("-ve", "", mgdat$desc)
      mgdat$strain <- gsub("/", "-", mgdat$strain)
      
      mgdat$desc <- gsub("^T\\d*", "", mgdat$desc, perl = TRUE)
      # View(mgdat)
  
      system( paste0("mkdir output/mm_", dbug, " vis/mm_", dbug ))

        
# subset
    
    ## ha! likely, we wish to compare like with like where possible. 
    ## only take FFM and RSM from runs where they were alongside SC40
    
      target_runs <- unique( dplyr::filter( mgdat, grepl( paste0("^", dbug, "$"), desc, perl = TRUE))$exp ) 
      if( length(target_runs) == 0){
        stop( " + + +   dbug:", dbug, " not recruiting any experiments")
      }
      dim( subdat <- dplyr::filter( mgdat, grepl( paste0("^", dbug, "$|^FFM$|^RSM$"), mgdat$desc, perl = TRUE), exp %in% target_runs ))
      
      dim(subfeat <- mgfeat[ rownames(subdat) , colSums(mgfeat[ rownames(subdat) , ]) > 0 ])
      subfeat_ra <- mgfeat_ra[ rownames(subdat) ,  colSums(mgfeat[ rownames(subdat) , ]) > 0]
      subfeat_clr <- (mgfeat_clr)[ rownames(subdat) , colnames(subfeat_ra) ]
      
      
      tp_order <- c("T-16", "T0", "T24", "T48",  "T72", "T144", "EB2",  "none", "EB1", NA)
      exp_order <- c("e0", "e1", "e9", "e11", "e12", "e13", "e14", "e15", "e16", "e17", "e18", "e19", "e2122",
                     "e21", "e22", "e23", "e24", "e25", "e27", "e28", "e29", "e30",
                     "EB2")

    ## unify taxa subset
      dim(subfeat_ra)[[1]]*0.05
      dim( suka_lesser <- t(k_A( t(subfeat_ra), k = 0.01, A = 0.05)))
      handy_tax_g <- unique( mgtax[ colnames(suka_lesser) , "g" ])
      handy_tax_c <- unique( mgtax[ colnames(suka_lesser) , "c" ])

    ## great cols
      uniq_krak <- colnames(suka_lesser)
      # rank to group to
      uniq_rank <- sort(unique( mgtax[ uniq_krak , "g"]))
      # for each entry in (.e.g.) genus, get all relevant (e.g.) species 
      tax_cols <- unlist( sapply( 1:length(uniq_rank), function(aa){    # aa <- 12 "c__Bacteroidia"
        bb_taxa_in_rank <- mgtax[ uniq_krak , "s"][ grep( uniq_rank[aa], mgtax[ uniq_krak , "g"])]
        cc_colours <- colorRampPalette( c( va_cols[aa], "white"))( length(bb_taxa_in_rank)+3)[1: length(bb_taxa_in_rank)]
        names(cc_colours) <- sort(bb_taxa_in_rank)
        cc_colours
      }))
      # scales::show_col(tax_cols)
      names(tax_cols) <- gsub("s__", "", names(tax_cols))
      tax_cols <- c( tax_cols, " other" = "grey95")
    
    ## for genus - could be more generic I supppose...
      tax_cols_higher <- sapply( unique( gsub("_.*", "", names(tax_cols))), function(aa){ tax_cols[ grep( aa, names(tax_cols))[1] ] })
      names(tax_cols_higher) <- unique( gsub("_.*", "", names(tax_cols)))
      
      
      # scales::show_col( 
        # desc_cols <- c("#00BFFF", "#009ca0","#EEE7D3", "#E4B777")
      desc_cols <- c("#00BFFF", "#009ca0","#D7E5EA", "#E4B777")
        # )   # "#BF5B17","#984EA3", 
      names( desc_cols ) <- c(NA, dbug, "FFM", "RSM")

      ## does not work with alpha,not diff enough
      # scales::show_col( 
        # time_cols <- colorRampPalette(colors = c("#CCE6E6", "#00647A"))(3)
        time_cols <- c("#A4243B", "#D8973C", "#00647A")
      # )
      names(time_cols) <- c("T0", "T24", "T48")

      cca_shapes <- c("FFM" = 23, "RSM" = 24, NULL = 22, "species" = 21)
      names(cca_shapes)[3] <- dbug

            
      # write.table( subdat, 
      #              paste0("output/mm_", dbug, "/fhi__mmoverlord__", dbug, "__subdat.tsv"),
      #              sep = "\t")
      # write.table( matrix( nrow = 5, c("exp", target_runs)), 
      #              paste0("output/mm_", dbug, "/fhi__mmoverlord__", dbug, "__subdat_target_runs.tsv"),
      #              sep = "\t")
      
      
      
##   relab   =====================================================================
      
    ## group samples by exp-timepint-desc   -------------------------------------
      
      # suka_lesser defined at opening section - unified subset
      suka <- cbind( suka_lesser, " other" = c(1 - rowSums(suka_lesser)) )
      head( sukam <- as.data.frame(reshape2::melt( suka ), stringsAsFactors = FALSE))
      sukam$Var1 <- as.character( sukam$Var1 )
      sukam$Var2 <- as.character( sukam$Var2 )
      
      sukam$etd <- paste0( mgdat[ sukam$Var1 , "exp"], "_", mgdat[ sukam$Var1 , "timepoint"], "_", mgdat[ sukam$Var1 , "desc"])
      head(suag <- aggregate( value ~ Var2 + etd, FUN = mean, sukam))
      head(suag_ra <- do.call("rbind", lapply( unique(suag$etd), function(aa){
        dd_df <- dplyr::filter( suag, etd == aa)
        dd_df$value <- (dd_df$value / sum(dd_df$value))*100
        dd_df
      })))
      
      suag_ra$exp <- gsub("^(.*)_.*_.*", "\\1", suag_ra$etd)
      suag_ra$timepoint <- gsub(".*_(.*)_.*", "\\1", suag_ra$etd)
      suag_ra$desc <- gsub(".*_.*_(.*)", "\\1", suag_ra$etd)
      suag_ra$plot_var <- gsub(".__", "", mgtax[ suag_ra$Var2 , "s" ])
      suag_ra$plot_var <- ifelse( is.na( suag_ra$plot_var), " other", suag_ra$plot_var)
      suag_ra$facet_var <- gsub(".__", "", mgtax[ suag_ra$Var2 , "p" ])
      suag_ra$fsi <- subdat[ match( suag_ra$exp, subdat$exp), "fsi" ]
      head( suag_ra) # View( suag_ra)
      
      
      ## great colours created at top when subsetting  
      ## great colours
      # taxa to plot (as KrakIDs)
      
      ## sort legend, indirectly
      # NOTE this doesn;t seem to be working - its still alphabetical
      suag_ra$plot_var <- factor( gsub("s__", "", suag_ra$plot_var), levels = names(tax_cols))
      
      tax_cols[ grep("coli", names(tax_cols)) ] <- "grey"
      
      suag_ra_plots <- lapply( c("T0", "T24", "T48"), function(aa){
        bb_df <- dplyr::filter( suag_ra, timepoint==aa)        
        cc_plot <- ggplot( bb_df, aes( fill = gsub("s__", "", plot_var), x = gsub("_.*", "", etd), y = value))+
          facet_grid( . ~ desc, space = "free_x", scales = "free") +
          geom_col(position = "stack", colour = "grey30", linewidth = 0.2) + 
          scale_fill_manual( values = tax_cols, "species", drop = FALSE, na.value = "red") +
          theme(
            legend.position = "bottom",
            legend.text = element_text(face = "italic", size = 11),
            legend.title = element_text(size = 14),
            axis.text.x = element_text(angle = 90, size = 12),
            axis.text.y = element_blank(),
            panel.spacing.y = unit(0.4, "lines"),
            panel.spacing.x = unit(0.4, "lines"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.y = element_blank(),
            NULL) +
          labs( x = "", y = "", title = aa, subtitle = "_______________________________" )
        
        if( aa == "T48"){
          cc_plot + theme(
            strip.text.y = element_text(),
            NULL)
        }else if( aa == "T0"){
          cc_plot + theme(
            axis.text.y = element_text(),
            NULL)
        }else{
          cc_plot
        }
        
      })
      
      plots_suag_ra <- suag_ra_plots[[1]] + suag_ra_plots[[2]] + suag_ra_plots[[3]] + plot_layout(guides = "collect")  &
        theme(legend.position = "bottom") &
        plot_annotation(
          caption = paste0("features filtered to retain if abundant above >1% in >5% of samples. Replicates combined within experiments. Only FFM, RSM, and ", dbug, " shown"),
          title = paste0( dbug, " abundance versus FFM and RSM ")
        )
      

    ## aggregated to one column per condition, ie. 9 bars total   -----------------------------------------
      
      sukam$td <- paste0( mgdat[ sukam$Var1 , "timepoint"], "_", mgdat[ sukam$Var1 , "desc"])
      head(reag <- aggregate( value ~ Var2 + td, FUN = mean, sukam))
      head(reag_ra <- do.call("rbind", lapply( unique(reag$td), function(aa){
        dd_df <- dplyr::filter( reag, td == aa)
        dd_df$value <- (dd_df$value / sum(dd_df$value))*100
        dd_df
      })))
      reag_ra$timepoint <- gsub("(.*)_.*", "\\1", reag_ra$td)
      reag_ra$desc <- gsub(".*_(.*)", "\\1", reag_ra$td)
      reag_ra$plot_var <- gsub(".__", "", mgtax[ reag_ra$Var2 , "s" ])
      reag_ra$plot_var <- ifelse( is.na( reag_ra$plot_var), " other", reag_ra$plot_var)
      reag_ra$facet_var <- gsub(".__", "", mgtax[ reag_ra$Var2 , "p" ])
      head( reag_ra) # View( reag_ra)

      reag_ra_plots <- lapply( c("T0", "T24", "T48"), function(aa){
        bb_df <- dplyr::filter( reag_ra, timepoint==aa)        
        cc_plot <- ggplot( bb_df, aes( fill = gsub("s__", "", plot_var), x = desc, y = value))+
          geom_col(position = "stack", colour = "grey30", linewidth = 0.2) + 
          scale_fill_manual( values = tax_cols, "species", drop = FALSE, na.value = "red") +
          theme(
            axis.ticks.length=unit(.15, "cm"),
            axis.ticks.x.bottom = element_line(linewidth = 0.6, colour = "lightskyblue3"),
            legend.position = "left",
            legend.text = element_text(face = "italic", size = 11),
            legend.title = element_text(size = 14),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
            axis.text.y = element_blank(),
            panel.spacing.y = unit(0.4, "lines"),
            panel.spacing.x = unit(0.4, "lines"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.y = element_blank(),
            strip.text.x = element_text(angle = 90, hjust = 0),
            NULL) +
          labs( x = "", y = "", title = aa) #, subtitle = "______________________" )
        
        if( aa == "T48"){
          cc_plot + theme(
            strip.text.y = element_text(),
            NULL)
        }else if( aa == "T0"){
          cc_plot + theme(
            axis.ticks.y.left = element_line(linewidth = 0.6, colour = "lightskyblue3"),
            axis.text.y = element_text(),
            NULL)
        }else{
          cc_plot
        }
        
      })
      
     ( plots_reag_ra <- reag_ra_plots[[1]] + reag_ra_plots[[2]] + reag_ra_plots[[3]] + plot_layout(guides = "collect")  &
        theme(legend.position = "left") &
        plot_annotation(
          caption = paste0("features filtered to retain if abundant above >1% in >5% of samples. Replicates and experiments combined. Only FFM, RSM, and ", dbug, " shown"),
          title = paste0( dbug, " abundance versus FFM and RSM ")
        ) & 
        guides( fill = guide_legend(ncol = 2)))
      
      
    ## mean levels of bif / lac / akk   -----------------------------------------
      
      # of_interest <- sapply( c("Bifidobacterium", "Lactobacillus", "Akkermansia", "Faecalibacterium"), function(aa){ rownames(mgtax)[ grep( aa, mgtax[ , "s"]) ] })
      of_interest <- sapply( c("Bifidobacterium", "Lactobacillus", "Akkermansia", "Faecalibacterium"), function(aa){ unique(grep( aa, suag_ra$plot_var, value = TRUE)) })
      
      indiv_plots <- lapply( c("T0", "T24", "T48"), function(aa){
        bb_df <- dplyr::filter( suag_ra, timepoint==aa)        
        cc_plot <- ggplot( 
          dplyr::filter( bb_df, plot_var %in% unlist( of_interest)),
          aes( fill = gsub("s__", "", plot_var), x = etd, y = value))+
          facet_grid( facet_var ~ desc, space = "free_x", scales = "free") +
          geom_col(position = "stack", colour = "grey30", linewidth = 0.2) + 
          scale_fill_manual( values = tax_cols, "species of\ninterest", drop = FALSE, na.value = "red") +
          theme(
            legend.position = "bottom",
            legend.text = element_text(face = "italic", size = 11),
            legend.title = element_text(size = 14),
            axis.text.x = element_text(angle = 90, size = 12),
            axis.text.y = element_blank(),
            panel.spacing.y = unit(0.4, "lines"),
            panel.spacing.x = unit(0.4, "lines"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.x = element_text(face = "italic"),
            strip.text.y = element_blank(),
            NULL) +
          labs( x = "", y = "", title = aa, subtitle = "_______________________________" )
        
        if( aa == "T48"){
          cc_plot + theme(
            strip.text.y = element_text(),
            NULL)
        }else if( aa == "T0"){
          cc_plot + theme(
            axis.text.y = element_text(),
            NULL)
        }else{
          cc_plot
        }
        
      })
      
      plots_indiv <- indiv_plots[[1]] + indiv_plots[[2]] + indiv_plots[[3]] + plot_layout(guides = "collect")  & 
        theme(legend.position = "bottom") &
        plot_annotation(
          caption = paste0("features filtered to retain if abundant above >1% in >5% of samples. Replicates combined within experiments. Only FFM, RSM, and ", dbug,  " shown"),
          title = paste0( dbug, " abundance versus FFM and RSM for select taxa")
        )
      
      
    ## relatively specific :: nine columns only for the cool stuff.    -------------------
      
      of_interest_2 <- sapply( c("Bifidobacterium", "Lactobacillus", "Faecalibacterium"), function(aa){ unique(grep( aa, suag_ra$plot_var, value = TRUE)) })
      scales::show_col(ran_cols <- sample(va_cols, 10, FALSE))
      # [1] "#FFDE8B" "#903495" "#8690FF" "#4A354F" "#FF4373" "#F7C0BB" "#FF6AD5" "#AE1357" "#CCDE8B"
      # [10] "#532E57"
      
      specifc_plots <- lapply( c("T0", "T24", "T48"), function(aa){
        bb_df <- dplyr::filter( reag_ra, timepoint==aa)        
        cc_plot <- ggplot( 
          dplyr::filter( bb_df, plot_var %in% unlist( of_interest_2)),
          aes( fill = gsub("s__", "", plot_var), x = td, y = value))+
          facet_grid( facet_var ~ desc, space = "free_x", scales = "free") +
          geom_col(position = "stack", colour = "grey30", linewidth = 0.2) + 
          # scale_fill_manual( values = tax_cols, "species of interest", drop = FALSE, na.value = "red") +
          scale_fill_manual( values = ran_cols, "species of interest", drop = FALSE, na.value = "red") +
          theme(
            legend.position = "right",
            legend.text = element_text(face = "italic", size = 11),
            legend.title = element_text(size = 14),
            axis.text.x = element_blank(), # text(angle = 90, size = 12),
            axis.text.y = element_blank(),
            panel.spacing.y = unit(0.4, "lines"),
            panel.spacing.x = unit(0.4, "lines"),
            plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            strip.text.y = element_blank(),
            strip.text.x = element_text(face = "italic"),
            NULL) +
          labs( x = "", y = "", title = aa, subtitle = "___________________" )
        
        if( aa == "T48"){
          cc_plot + theme(
            strip.text.y = element_text(),
            NULL)
        }else if( aa == "T0"){
          cc_plot + theme(
            axis.text.y = element_text(),
            NULL)
        }else{
          cc_plot
        }
        
      })
      
      plots_specif <- specifc_plots[[1]] + specifc_plots[[2]] + specifc_plots[[3]] + plot_layout(guides = "collect")  & 
        theme(legend.position = "right") &
        plot_annotation(
          caption = paste0("features filtered to retain if abundant above >1% in >5% of samples. Replicates combined within experiments. Only FFM, RSM, and ", dbug,  " shown"),
          title = paste0( dbug, " abundance versus FFM and RSM for select taxa")
        )
      
      
## alpha   ==========================================================================
    
      divs <- c("invs", "shan", "chao1")
  
          
  ## alpha beta testing

    ## strain more tham FFM
      
      wtest1 <- smooth_out(do.call("rbind", lapply( divs, function(aa){   # aa <- "chao1"
        wilcox.test( subdat[ dplyr::filter(subdat, desc == "FFM", timepoint == "T48")$fhi_samp, aa] , subdat[ dplyr::filter(subdat, desc == dbug, timepoint == "T48")$fhi_samp, aa], data = subdat)
      }))[ , c(1,3,5:6)] )
      rownames(wtest1) <- paste0( divs, "_FFM_", dbug)
      
      wtest2 <- smooth_out(do.call("rbind", lapply( divs, function(aa){   # aa <- "chao1"
        wilcox.test( subdat[ dplyr::filter(subdat, desc == "RSM", timepoint == "T48")$fhi_samp, aa] , subdat[ dplyr::filter(subdat, desc == dbug, timepoint == "T48")$fhi_samp, aa], data = subdat)
      }))[ , c(1,3,5:6)] )
      rownames(wtest2) <- paste0( divs, "_RSM_", dbug)

      wts <- rbind( wtest1, wtest2)
      wts$FDR <- p.adjust( wts[ , "p.value"], method = "fdr")
      table_alpha <- wts[ , c(1,2,5)]
      

    ## plots ---
      al_dat <- dplyr::filter(subdat, timepoint != "none")
      al_plots <- lapply( divs, function(aa){    # aa <- "shan"
        bb_plot <- ggplot( al_dat, aes_string(y = aa)) +
          # facet_grid( . ~ fsi) + 
          scale_fill_manual( values = desc_cols, "condition") +
          geom_boxplot(aes(x = timepoint, fill = desc), outlier.shape = 21, alpha = 0.6) +
          # geom_point(aes(x = timepoint, fill = desc), shape = 21, alpha = 0.6, size = 2.5, position = position_jitterdodge(jitter.width = 0.05)) + #position_dodge(width = 0.75)) +
          # scale_shape_manual( values = c(21,24,23,25))
          theme(
            plot.title = element_text(size = 12, face = "plain"),
            NULL) +
          labs(
            x = ""
          )
            
        if( aa == "chao1"){
          bb_plot + labs(title = "n species (Chao1)", y = "n species")
        }else if(aa == "shan"){
          bb_plot + 
            labs(title = "Shannon's H index", y = "H")
        }else if(aa == "invs"){
          bb_plot + labs(title = "Inverse  Simpson's index", y = "1/Simp")
        }else{ bb_plot}
        
        })
      
      plots_alpha <- (al_plots[[3]] + al_plots[[1]] + al_plots[[2]] + plot_layout(guides = "collect")) & plot_annotation(title = "alpha diversity measures")
    
      
# beta 
      
    ## beta testing - heterosced??
    
      # desc  
        bc_beta <- betadisper( vegdist( subfeat_ra, method = "bray"), group = subdat[rownames(subfeat_ra) , "desc"] )
        table_het1_desc <- anova( bc_beta)
          # >  Analysis of Variance Table
          # >  
          # >  Response: Distances
          # >  Df  Sum Sq   Mean Sq F value Pr(>F)
          # >  Groups      2 0.01498 0.0074879  0.3723   0.69
          # >  Residuals 111 2.23240 0.0201118           
      # fsi  
      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
        bc_beta_fsi <- betadisper( vegdist( subfeat_ra, method = "bray"), group = subdat[rownames(subfeat_ra) , "fsi"] )
        table_het2_fsi <- anova( bc_beta_fsi)
      }else{
        table_het2_fsi <- NULL
      }
          # >  Analysis of Variance Table
          # >  
          # >  Response: Distances
          # >  Df  Sum Sq   Mean Sq F value Pr(>F)
          # >  Groups      2 0.02966 0.014830  1.7898 0.1718
          # >  Residuals 111 0.91975 0.008286                  

      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
      table_adonis_full <- adonis2( 
        vegdist( subfeat_ra, method = "bray") ~
          desc + timepoint + exp + fsi, data = subdat, by = "margin")
      }else{
      table_adonis_full <- adonis2( 
        vegdist( subfeat_ra, method = "bray") ~
          desc + timepoint, data = subdat, by = "margin")
        
      }
          # >  Permutation test for adonis under reduced model
          # >  Marginal effects of terms
          # >  Permutation: free
          # >  Number of permutations: 999
          # >  
          # >  adonis2(formula = vegdist(subfeat_ra, method = "bray") ~ desc + timepoint + exp + fsi, data = subdat, by = "margin")
          # >  Df SumOfSqs      R2       F Pr(>F)    
          # >  desc        2   2.4663 0.08528 13.4492  0.001 ***
          # >  timepoint   2  13.5839 0.46970 74.0748  0.001 ***
          # >  exp         1   0.4537 0.01569  4.9484  0.002 ** 
          # >  fsi         0   0.0000 0.00000    -Inf           
          # >  Residual  106   9.7192 0.33607                   
          # >  Total     113  28.9205 1.00000                   
          # >  ---
          # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
      
    ## beta div pairwise
      compars <- combn( c("RSM", "FFM", dbug), 2)
      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
      table_adonis_comb <- do.call( "rbind", apply( compars, 2, function(aa){   # aa <- compars [, 3]
                bb_compar_samps <- rownames( dplyr::filter( subdat, desc == aa[1] | desc == aa[2] ))
                cc_df <- as.data.frame( adonis2( vegdist( subfeat_ra[ bb_compar_samps, ], method = "bray") ~ desc + timepoint + (exp) + (fsi), data = subdat[ bb_compar_samps, ]))
                data.frame( "combo" = rep(paste0(aa, collapse = "-"), 3),
                            cc_df[ c(1:3) , ])
      }))
      }else{
      table_adonis_comb <- do.call( "rbind", apply( compars, 2, function(aa){   # aa <- compars [, 3]
                bb_compar_samps <- rownames( dplyr::filter( subdat, desc == aa[1] | desc == aa[2] ))
                cc_df <- as.data.frame( adonis2( vegdist( subfeat_ra[ bb_compar_samps, ], method = "bray") ~ desc + timepoint, data = subdat[ bb_compar_samps, ]))
                data.frame( "combo" = rep(paste0(aa, collapse = "-"), 3),
                            cc_df[ c(1:3) , ])
      }))
      }
      
      
    # different FSI induce significantly different structures
      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
        table_adonis_fsi <- adonis2( vegdist( subfeat_ra, method = "bray") ~ fsi, data = subdat)
      }else{
        table_adonis_fsi <- NULL
      }
      
          # >  Permutation test for adonis under reduced model
          # >  Terms added sequentially (first to last)
          # >  Permutation: free
          # >  Number of permutations: 999
          # >  
          # >  adonis2(formula = vegdist(subfeat_ra, method = "bray") ~ fsi, data = subdat)
          # >  Df SumOfSqs      R2      F Pr(>F)    
          # >  fsi        2   2.6974 0.09327 5.7088  0.001 ***
          # >  Residual 111  26.2231 0.90673                  
          # >  Total    113  28.9205 1.00000                  
          # >  ---
          # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
          # >
      
      
##   diversity orientation   ====================================================
      
      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
        sub_cca <- cca( decostand( subfeat_ra, method = "hel") ~ Condition(fsi) + desc + timepoint + Condition(exp), data = subdat)
      }else{
        sub_cca <- cca( decostand( subfeat_ra, method = "hel") ~ desc + timepoint, data = subdat)
      }
      fort <- fortify(sub_cca)
      fort$Label <- fort$label
      fort$Score <- fort$score
      fort$exp <- subdat[ fort$Label, "exp"]
      fort$desc <- subdat[ fort$Label, "desc"]
      fort$timepoint <- subdat[ fort$Label, "timepoint"]
      fort$fsi <- subdat[ fort$Label, "fsi"]
      fort$plot_var <- sapply( mgtax[ fort$Label , "g"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_g ){ aa }else{"other"} ) })
      fort$facet_var <- sapply( mgtax[ fort$Label , "c"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_c ){ aa }else{"other"} ) })
      
    # abundances - prob a waste of your time
      fort$abundo <- ifelse( fort$Label %in% colnames(subfeat_ra), colMeans(subfeat_ra)[fort$Label], NA)
      cca_spec_centr <- data.frame( 
        aggregate( abundo ~ plot_var, FUN = sum, data = fort),
        "CCA1" = aggregate( CCA1 ~ plot_var, FUN = mean, data = fort)[,2],
        "CCA2" = aggregate( CCA2 ~ plot_var, FUN = mean, data = fort)[,2],
        stringsAsFactors = FALSE
      )      
      
      head(fort)
      
      plots_beta_full <- ggplot( dplyr::filter(fort, Score == "sites")) +
        coord_fixed(ratio = 1) +
        geom_hline(yintercept = 0, linewidth = 0.1) + geom_vline(xintercept = 0, linewidth = 0.1) +
        
        geom_point( size = 5, alpha = 0.2, aes( fill = desc, shape = timepoint, x = CCA1, y = CCA2), colour = "black") +
        geom_label(
          data = dplyr::filter(fort, Score == "centroids"),
          aes( label = gsub("timepoint|desc", "", Label), x = CCA1, y = CCA2), colour = "grey50") +
        geom_point(
          data = cca_spec_centr, shape = 21,
          aes( size = (abundo), x = CCA1, y = CCA2), colour = "black") +
        geom_point(
          data = cca_spec_centr, shape = 20,
          aes( colour = plot_var, size = (abundo), x = CCA1, y = CCA2)) +
        
        scale_colour_manual(values = tax_cols_higher) +
        scale_fill_manual(values = desc_cols, "microMatrix timepoint:") +
        scale_shape_manual(values = unname(cca_shapes), drop = FALSE, "condition:" ) +
        theme( legend.spacing = unit(units = "cm", 0.01)) +
        guides(
          size = "none",
          colour = guide_legend(override.aes = list(shape = 20, size = 8), 
                                ncol = 1, "average position of Genus:", direction = "horizontal", title.position = "top", order = 3),
          fill = guide_legend(override.aes = list(shape = 22, size = 8, alpha =1), direction = "horizontal", title.position = "top", order= 1),
          shape = guide_legend(override.aes = list(colour = "black", alpha = 1), direction = "horizontal", title.position = "top", order = 2)
        )
                                                                            
      
    ##   diversity & T48? -------------------------------------------------------
      
      t48_samps <- dplyr::filter( subdat, timepoint == "T48" )$fhi_samp
      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
        sub_cca_t48 <- cca( decostand( subfeat_ra[t48_samps,], method = "hel") ~ Condition(fsi) + desc + shan + Condition(exp), data = subdat[t48_samps,])
      }else{
        sub_cca_t48 <- cca( decostand( subfeat_ra[t48_samps,], method = "hel") ~ desc + shan, data = subdat[t48_samps,])
      }
        
      forte8 <- fortify(sub_cca_t48)
      forte8$Score <- forte8$score
      forte8$Label <- forte8$label
      forte8$exp <- subdat[ forte8$Label, "exp"]
      forte8$desc <- subdat[ forte8$Label, "desc"]
      forte8$timepoint <- subdat[ forte8$Label, "timepoint"]
      forte8$fsi <- subdat[ forte8$Label, "fsi"]
      forte8$plot_var <- sapply( mgtax[ forte8$Label , "g"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_g ){ aa }else{"other"} ) })
      forte8$facet_var <- sapply( mgtax[ forte8$Label , "c"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_c ){ aa }else{"other"} ) })
    # abundances - prob a waste of your time
      forte8$abundo <- ifelse( forte8$Label %in% colnames(subfeat_ra), colMeans(subfeat_ra)[forte8$Label], NA)
      # abundances - prob a waste of your time
      fort$abundo <- ifelse( fort$Label %in% colnames(subfeat_ra), colMeans(subfeat_ra)[fort$Label], NA)
      cca_spec_centr_48 <- data.frame( 
        aggregate( abundo ~ plot_var, FUN = sum, data = forte8),
        "CCA1" = aggregate( CCA1 ~ plot_var, FUN = mean, data = forte8)[,2],
        "CCA2" = aggregate( CCA2 ~ plot_var, FUN = mean, data = forte8)[,2],
        stringsAsFactors = FALSE
      )      
      
      plots_beta_48 <- ggplot(forte8) +
        coord_fixed(ratio = 1) +
        geom_hline(yintercept = 0, linewidth = 0.1) + geom_vline(xintercept = 0, linewidth = 0.1) +
        geom_point( data = subset(forte8, Score == "sites"), size = 5, alpha = 0.4, aes( fill = desc, x = CCA1, y = CCA2, shape = Scores), colour = "grey40", shape = 23) +
        geom_point(
          data = cca_spec_centr_48, shape = 21,
          aes( size = (abundo), x = CCA1, y = CCA2), colour = "black") +
        geom_point(
          data = cca_spec_centr_48, shape = 20,
          aes( colour = plot_var, size = (abundo), x = CCA1, y = CCA2)) +
        geom_segment(
          data = dplyr::filter( forte8, Label == "shan"),
          aes(x = 0, y = 0, xend = CCA1, yend = CCA2),
          colour = "grey20", arrow = arrow(angle = 15, length = unit(0.2, "cm"),
                                           ends = "last", type = "closed"),  )  + # xend = CCA1, yend = CCA2) +
        ggrepel::geom_text_repel(
          data = dplyr::filter( forte8, Label == "shan"),
          aes( label = "Shannon\ndiversity", x = CCA1, y = CCA2), force = 0.2, size = 5, ylim = -1.6, segment.colour = NA ) +

        
        
        scale_colour_manual(values = tax_cols_higher) +
        scale_fill_manual(values = desc_cols, "fermentate:") +
        scale_shape_manual(values = c("species" = 21, "sites" = 23), drop = FALSE) +
        guides(
          size = "none",
          colour = guide_legend(override.aes = list(shape = 20, size = 8), ncol = 1, "average position of Genus:"),
          fill = guide_legend(override.aes = list(shape = 22, size = 8, alpha =1)),
          shape = guide_legend(override.aes = list(colour = "black", alpha = 1))
        )

    table_cca <- anova(sub_cca, by = "margin")
      # >  Permutation test for cca under reduced model
      # >  Marginal effects of terms
      # >  Permutation: free
      # >  Number of permutations: 999
      # >  
      # >  Model: cca(formula = decostand(subfeat_ra, method = "hel") ~ Condition(fsi) + desc + timepoint + Condition(exp), data = subdat)
      # >  Df ChiSquare       F Pr(>F)    
      # >  desc        2   0.13416  5.7617  0.001 ***
      # >  timepoint   2   0.36516 15.6819  0.001 ***
      # >  Residual  106   1.23412                   
      # >  ---
      # >  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

    
##   L M M   n u m   n u m   n u m s   ===========================================

    ## differential abundance: For T48 - T24 of less interest, T0 of accademic interest so far.

    # ## check clr
      # - checked spreadsheet_shortcuts:
      ## >    dim(k2_counts) ; dim(fhidat)
      ## >>     [1]  730 1611
      ## >>     [1] 1611   18
      ## >  ## zCompositions assumes samples are ROWS!!!
      ## >    k2_cmult <- zCompositions::cmultRepl( t(k2_counts), method = "CZM")
      ## >    dim(k2_cmult)
      ## >>     [1] 1611   730
      ## >    k2_clr <- t( apply( k2_cmult, 1, function(aa){ log(aa) - log(mean(aa)) }))   # formula from vegan
      ## >    dim(k2_clr)
      ## >>     [1] 1611   730


    ## subset correctly
      # - not k_A'd - simply subset to taxa seen in this subset etc.
      # now subset to FFM, RSM, DBUG samples
      dim(subfeat_clr)
      rownames(subfeat_clr)
      # hist( subfeat_clr, breaks = 50)


    ## from fhi__deyan__chunk__LMER-test.R, now in fhi__mmoverlord__chunk__LMER-EMM-outliers.R
      
      # - see original for thinking on outliers, LASSO, log-ratio LASSO etc.
      library(emmeans)
      library(lme4)

      
    ## outlier testing outcome:   ------------------------------------------
      
      # large number of outliers, but there is nothing to motivate their REMOVAL in this case
      # we must try stay cognisant of the data we have
      
      
  ##   feat test ---------------------------------------------------------------------
      
    ## "for the final works" a continuous ranef for age of fsi would be great (or a cont. semi-dummy ranef for age of each)
    ## still stuck in univariate methods. need to go multivariate. 
      
      prot_cluster <- parallel::makeCluster(spec = 7, type = "PSOCK")
      parallel::clusterExport( cl = prot_cluster, c("subfeat_clr", "subdat", "mgtax"), envir = .GlobalEnv )
      
      feat_test <- do.call ("rbind",
                            # parallel::mclapply(
                            parallel::parLapply( prot_cluster,
                            # lapply(
                              1:ncol(subfeat_clr),
                              function(aa){  # aa <- 12    # aa <- 173 reuteri     aa_kid <- "k2_0001598"
                                
                                library(lme4)
                                library(emmeans)
                                
                                aa_kid <- colnames(subfeat_clr)[aa]
                                gsub("s__", "", mgtax[ aa_kid, "s"])
                                head(bb_df <- data.frame(
                                  "var" = subfeat_clr[ , aa_kid ],
                                  "exp" = subdat[ rownames(subfeat_clr) , "exp"],
                                  "fsi" = subdat[ rownames(subfeat_clr) , "fsi"],
                                  "desc" = subdat[ rownames(subfeat_clr) , "desc"],
                                  "timepoint" = subdat[ rownames(subfeat_clr) , "timepoint"],
                                  stringsAsFactors = FALSE
                                ))
                                
                                if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
                                  cc_lmer <- lmer( var ~ timepoint * desc + (1|exp), data = bb_df)
                                }else{
                                  cc_lmer <- lm( var ~ timepoint * desc, data = bb_df)
                                }
                                
                                # ## resid diagnostic :: qqplot
                                #     # cc_res <- data.frame(
                                #     #   "sample" = names(resid( cc_lmer )),
                                #     #   "resid" = resid( cc_lmer ),
                                #     #   "exp" = subdat[ names(resid( cc_lmer )) , "exp"],
                                #     #   "fsi" = subdat[ names(resid( cc_lmer )) , "fsi"], stringsAsFactors = FALSE )
                                #     # ggplot( cc_res) +    # , aes( fill = desc, shape = timepoint)
                                #     #   facet_wrap( exp ~ .) +
                                #     #   geom_qq(aes(sample = resid, colour = exp), shape = 21) +
                                #     #   geom_qq_line(aes(sample = resid, colour = exp)) +
                                #     #   labs( title = paste0( "qq plot for ", gsub("s__", "", mgtax[aa_kid, "s"])), subtitle =paste0("nrows = ", nrow(cc_res)),
                                #     #         x = "theoretical normal quantiles", y = "observed resid quantiles") +
                                #     #   theme(legend.direction = "horizontal", legend.position = "bottom", strip.text.x = element_text(size = 8)) +
                                #     #   guides( fill = guide_legend(override.aes = list(shape = 22, size = 10)) )
                                
                                
                                ## create all the estimates first, then do the pairwise comparisons
                                ## note:  | sorts out a grouping factor - compare WITHIN timepoints, and not across (* behaviour)
                                (dd_d_tp <- emmeans(cc_lmer, specs = ~ desc | timepoint))
                                ## doublecheck that the values are unique, then subset
                                feat_test <- data.frame( 
                                  "kid" = aa_kid, 
                                  "tax" = mgtax[ aa_kid, "s"],
                                  as.data.frame(pairs(dd_d_tp, adjust = "none" ))[ 7:9, ],
                                  stringsAsFactors = FALSE
                                )
                                
                              # }, mc.cores = 1, mc.cleanup = TRUE, mc.preschedule = TRUE) )
                              }) )
      
      
      feat_test$FDR <- p.adjust(feat_test$p.value, method = "fdr")
      # View(dplyr::filter(feat_test, FDR < 0.05) )
      
      length(unique( dplyr::filter(feat_test, FDR < 0.05)$kid))   # 227 for SC40, similar for LH88
      
      
  ## visualise in some manner ------------------------------------------------------------------------
      
      dim(subfeat_clr)      
      # not sparing much here
      dim(dplyr::filter(feat_test, FDR < 0.05, !grepl( dbug, contrast) ) )    # SC40-175
      # at least these are the important ones
      dim(dplyr::filter(feat_test, FDR < 0.05, !grepl( "RSM", contrast) ))        # SC40-154
      # and this one in particular
      # dim(dplyr::filter(feat_test, FDR < 0.05, contrast == paste0( "(RSM) - ", dbug) ))   # SC40-0!!!!!!!!!!!!
      dim(dplyr::filter(feat_test, FDR < 0.05, !grepl( "FFM", contrast)))   # SC40-0!!!!!!!!!!!!
      
      
  ## plots   -----------------------------------------------------------------------------------------
      
      theme_update(
        text = element_text(family = "sans"), 
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.ticks.length = unit(0.4, "cm"),
        axis.text.y = element_text(family = "sans", face = "italic"),
        plot.title = element_text( hjust = 0)
      )
      
      
      ## 1 - SC40-RSM
      
      head( sr_plot <- as.data.frame(t(sapply( unique(dplyr::filter(feat_test, !grepl( "FFM", contrast), FDR < 0.05)$kid),
                                               function(aa){   # aa <- "k2_0001354"
                                                 c("kid" = aa,
                                                   "lfc" = log2( mean( subfeat[ dplyr::filter(subdat, desc == dbug)$fhi_samp , aa ]) + 1) - log2( mean( subfeat[ dplyr::filter(subdat, desc == "RSM")$fhi_samp , aa ]) + 1),
                                                   "s" = gsub("s__", "", mgtax[ aa , 7]),
                                                   "g" = gsub("g__", "", mgtax[ aa , 6]),
                                                   "f" = gsub("f__", "", mgtax[ aa , 5]))
                                                 
                                               })), stringsAsFactors = FALSE) )
      sr_plot[ , "lfc"] <- as.numeric( sr_plot[ , "lfc"])
      sr_plot$outcome <- ifelse( sr_plot$lfc > 0, dbug, "RSM")
      
      ## needs to be fixed so that legend shows up - currently spinning two different dataframes - melt together instead, then plot
      if( nrow(sr_plot) != 0){
        pl1 <- ggplot( sr_plot, aes(x = reorder(gsub("_", " ", sr_plot$s), lfc, sum, decreasing = TRUE))) +
            coord_flip() + 
            geom_hline(yintercept = 0, colour = "black", size = 0.9) +
            geom_col( aes(y = lfc, fill = outcome), position = position_dodge(), colour = "black", size = 0.8, alpha = 0.55, width = 0.65) +
            labs(
              x = "" , y = "log2-fold change (LFC)", title = paste0("taxa - Effect of ", dbug, " w.r.t RSM")
            ) +
            scale_fill_manual( values = desc_cols, "higher in:") +
            theme( legend.position = c(0.75, 0.75),
                   legend.text = element_text(size = 12),
                   legend.title = element_text(size = 14),
                   # plot.title = element_text( hjust = 1.5)
                   NULL
            )+
          guides(
            fill = guide_legend(keywidth = unit(units = "cm", 1), keyheight = unit(units = "cm", 1))
          )
      }else{ 
        pl1 <- NULL
        }
      
      
      
      head( sf_plot <- as.data.frame(t(sapply( unique(dplyr::filter(feat_test, !grepl( "RSM", contrast), FDR < 0.05)$kid),
                                               function(aa){   # aa <- "k2_0001354"
                                                 c("kid" = aa,
                                                   "lfc" = log2( mean( subfeat[ dplyr::filter(subdat, desc == dbug)$fhi_samp , aa ]) + 1) - log2( mean( subfeat[ dplyr::filter(subdat, desc == "FFM")$fhi_samp , aa ]) + 1),
                                                   "s" = gsub("s__", "", mgtax[ aa , 7]),
                                                   "g" = gsub("g__", "", mgtax[ aa , 6]),
                                                   "f" = gsub("f__", "", mgtax[ aa , 5]))
                                                 
                                               })), stringsAsFactors = FALSE) )
      sf_plot[ , "lfc"] <- as.numeric( sf_plot[ , "lfc"])
      sf_plot$outcome <- ifelse( sf_plot$lfc > 0, dbug, "FFM")
      
      ## needs to be fixed so that legend shows up - currently spinning two different dataframes - melt together instead, then plot
      pl2 <- ggplot( sf_plot, aes(x = reorder(gsub("_", " ", sf_plot$s), lfc, sum, decreasing = TRUE))) +
          coord_flip() + 
          geom_hline(yintercept = 0, colour = "black", size = 0.9) +
          geom_col( aes(y = lfc, fill = outcome), position = position_dodge(), colour = "black", size = 0.8, alpha = 0.55, width = 0.65) +
          labs(
            x = "" , y = "log2-fold change (LFC)", title = paste0("taxa - Effect of ", dbug, " w.r.t. FFM")
          ) +
          scale_fill_manual( values = desc_cols, "higher in:") +
        theme( legend.position = c(0.75, 0.75),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 14),
               # plot.title = element_text( hjust = 1.5)
               NULL
        )+
        guides(
          fill = guide_legend(keywidth = unit(units = "cm", 1), keyheight = unit(units = "cm", 1))
        )
      
      
      head( rf_plot <- as.data.frame(t(sapply( unique(dplyr::filter(feat_test, !grepl( dbug, contrast), FDR < 0.05)$kid),
                                               function(aa){   # aa <- "k2_0001354"
                                                 c("kid" = aa,
                                                   "lfc" = log2( mean( subfeat[ dplyr::filter(subdat, desc == "RSM")$fhi_samp , aa ]) + 1) - log2( mean( subfeat[ dplyr::filter(subdat, desc == "FFM")$fhi_samp , aa ]) + 1),
                                                   "s" = gsub("s__", "", mgtax[ aa , 7]),
                                                   "g" = gsub("g__", "", mgtax[ aa , 6]),
                                                   "f" = gsub("f__", "", mgtax[ aa , 5]))
                                                 
                                               })), stringsAsFactors = FALSE) )
      rf_plot[ , "lfc"] <- as.numeric( rf_plot[ , "lfc"])
      rf_plot$outcome <- ifelse( rf_plot$lfc > 0, "RSM", "FFM")
      
      ## needs to be fixed so that legend shows up - currently spinning two different dataframes - melt together instead, then plot
      pl3 <- ggplot( rf_plot, aes(x = reorder(gsub("_", " ", rf_plot$s), lfc, sum, decreasing = TRUE))) +
          coord_flip() + 
          geom_hline(yintercept = 0, colour = "black", size = 0.9) +
          geom_col( aes(y = lfc, fill = outcome), position = position_dodge(), colour = "black", size = 0.8, alpha = 0.55, width = 0.65) +
          labs(
            x = "" , y = "log2-fold change (LFC)", title = "taxa - difference in controls (RSM:FFM)"
          ) +
          scale_fill_manual( values = desc_cols, "higher in:") +
        theme( legend.position = c(0.75, 0.75),
               legend.text = element_text(size = 12),
               legend.title = element_text(size = 14),
               # plot.title = element_text( hjust = 1.5)
               NULL
        )+
        guides(
          fill = guide_legend(keywidth = unit(units = "cm", 1), keyheight = unit(units = "cm", 1))
        )
      
      ## check a specific bacterium/phage   -----------------------------
        # "k2_0001584" %in% suag_ra$Var2
        # ggplot( dplyr::filter(suag_ra, Var2 == "k2_0001584")) +
        #   geom_boxplot(aes(x = desc, y = value)) +
        #   geom_jitter(aes(x = desc, y = value), width = 0.2, height = 0.2 ) +
        #   labs(
        #     title = gsub("_", "\n", gsub("s__", "", mgtax[ "k2_0001584" , "s"]))
        #   )

  
  ## plot just the important stuff in an ilustrative CCA

      if( length(unique(subdat$fsi)) >1 & length(unique(subdat$exp)) >1 ){
        lmm_cca <- cca( decostand( subfeat_ra[ , unique(dplyr::filter(feat_test, FDR < 0.01)$kid) ], method = "hel") ~ desc + timepoint + Condition(fsi) + Condition(exp), data = subdat)
      }else{
        lmm_cca <- cca( decostand( subfeat_ra[ , unique(dplyr::filter(feat_test, FDR < 0.01)$kid) ], method = "hel") ~ desc + timepoint, data = subdat)
      }
      fort_cca <- fortify(lmm_cca)
      fort_cca$Label <- fort_cca$label
      fort_cca$Score <- fort_cca$score
      fort_cca$exp <- subdat[ fort_cca$Label, "exp"]
      fort_cca$desc <- subdat[ fort_cca$Label, "desc"]
      fort_cca$timepoint <- subdat[ fort_cca$Label, "timepoint"]
      fort_cca$fsi <- subdat[ fort_cca$Label, "fsi"]
      fort_cca$plot_var <- sapply( mgtax[ fort_cca$Label , "g"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_g ){ aa }else{"other"} ) })
      fort_cca$facet_var <- sapply( mgtax[ fort_cca$Label , "c"], function(aa){
        gsub(".__", "", 
             if( is.na(aa)){ NA
             }else if( aa %in% handy_tax_c ){ aa }else{"other"} ) })
      
      # abundances - prob a waste of your time
      fort_cca$abundo <- ifelse( fort_cca$Label %in% colnames(subfeat_ra), colMeans(subfeat_ra)[fort_cca$Label], NA)
      cca_spec_centr <- data.frame( 
        aggregate( abundo ~ plot_var, FUN = sum, data = fort_cca),
        "CCA1" = aggregate( CCA1 ~ plot_var, FUN = mean, data = fort_cca)[,2],
        "CCA2" = aggregate( CCA2 ~ plot_var, FUN = mean, data = fort_cca)[,2],
        stringsAsFactors = FALSE
      )      
      
      head(fort_cca)
      
      plots_beta_lmm <- ggplot( dplyr::filter(fort_cca, Score == "sites")) +
        coord_fixed(ratio = 1) +
        geom_hline(yintercept = 0, linewidth = 0.1) + geom_vline(xintercept = 0, linewidth = 0.1) +
        
        geom_point( size = 5, alpha = 0.2, aes( fill = desc, shape = timepoint, x = CCA1, y = CCA2), colour = "black") +
        geom_label(
          data = dplyr::filter(fort_cca, Score == "centroids"),
          aes( label = gsub("timepoint|desc", "", Label), x = CCA1, y = CCA2), colour = "grey50") +
        geom_point(
          data = cca_spec_centr, shape = 21,
          aes( size = (abundo), x = CCA1, y = CCA2), colour = "black") +
        geom_point(
          data = cca_spec_centr, shape = 20,
          aes( colour = plot_var, size = (abundo), x = CCA1, y = CCA2)) +
        
        scale_colour_manual(values = tax_cols_higher) +
        scale_fill_manual(values = desc_cols, "microMatrix timepoint:") +
        scale_shape_manual(values = unname(cca_shapes), drop = FALSE, "condition:" ) +
        theme( legend.spacing = unit(units = "cm", 0.01)) +
        guides(
          size = "none",
          colour = guide_legend(override.aes = list(shape = 20, size = 8), 
                                ncol = 1, "average position of Genus:", direction = "horizontal", title.position = "top", order = 3),
          fill = guide_legend(override.aes = list(shape = 22, size = 8, alpha =1), direction = "horizontal", title.position = "top", order= 1),
          shape = guide_legend(override.aes = list(colour = "black", alpha = 1), direction = "horizontal", title.position = "top", order = 2)
        )
  
            
##   h e l l o   w o r l d   ===========================================================
      
      # plots
        plots_suag_ra + plots_indiv + plots_reag_ra
        plots_alpha + plots_beta_full + plots_beta_48 + plots_beta_lmm
        pl1 + pl2 + pl3

        
      ## lines are all different colours
        
        
      ## lfc  > 1, FDR < 0.05        
        pl1_bugset_gt1 <- ggplot( dplyr::filter(sr_plot, lfc > 1),
                aes(x = reorder( gsub("_", " ", dplyr::filter(sr_plot, lfc > 1)$s), lfc, sum, decreasing = TRUE))) +
          # coord_flip() + 
          geom_hline(yintercept = 0, colour = "black", size = 0.9) +
          geom_col( aes(y = lfc, fill = outcome), position = position_dodge(), colour = "black", size = 0.8, alpha = 0.55, width = 0.65) +    # fill = outcome
          labs(
            x = "" , y = "LFC", title = paste0("taxa - Effect of ", dbug, " w.r.t RSM")
          ) +
          scale_fill_manual( values = desc_cols, "higher in:") +
          theme( axis.text.x = element_text(angle = 90, face = "italic", hjust = 1, vjust = 0.5),
                 legend.position = "right", # c(0.75, 0.5),
                 legend.text = element_text(size = 12),
                 legend.title = element_text(size = 14, hjust = 0.5),
                 plot.title = element_text( size = 14, face = "bold"),
                 NULL
          ) +
          labs(
            title = "log2-fold change (LFC)\nw.r.t. RSM"
          ) +
          guides( fill = guide_legend(keywidth = unit(units = "cm", 1), keyheight = unit(units = "cm", 1)) ) 

        
      ##  ------------------------------------        
        
        ( 
          (reag_ra_plots[[1]] + reag_ra_plots[[2]] + reag_ra_plots[[3]] +
            plot_layout(guides = "collect", ncol = 3) &
            theme(
              legend.position = "left"
              ) &
            guides(
              fill = guide_legend(ncol = 2)) ) | plot_spacer() |
        
        (
          pl1_bugset_gt1 /
           (al_plots[[1]] + 
              labs(y = "", title = "Inverse Simpson's Index") +
              theme(
                 legend.title = element_text( hjust = 0.5, size = 14),
                 legend.text = element_text( hjust = 0.5, size = 12),
                 plot.title = element_text( hjust = 0.5, size = 14, face = "bold")
              ) )
          ) + plot_layout(heights = c(0.25,0.75))) +
              plot_layout(widths = c(0.68, 0.04, 0.28))

        
      ##  ------------------------------------        
        
      # tables
        knitr::kable(table_adonis_comb)
        knitr::kable(table_het1_desc)
        knitr::kable(table_het2_fsi)
        knitr::kable(table_adonis_fsi)
        knitr::kable(table_adonis_full)
        knitr::kable(anova(sub_cca, by = "margin"))
        knitr::kable(table_alpha)
        knitr::kable( dplyr::filter(feat_test, FDR < 0.05))

        
      saveRDS( list( plots_suag_ra, plots_indiv, plots_reag_ra, plots_alpha, plots_beta_full, plots_beta_48, plots_beta_lmm, pl1, pl2, pl3),
               paste0("output/mm_", dbug, "/fhi__mmoverlord__", dbug, "__plots.RDS") )

      saveRDS( list( feat_test, dplyr::filter(feat_test, FDR < 0.05), table_adonis_comb, table_het1_desc, table_het2_fsi, table_adonis_fsi,
                     table_adonis_full, anova(sub_cca, by = "margin"), table_alpha),
                     paste0("output/mm_", dbug, "/fhi__mmoverlord__", dbug, "__tables.RDS") )

      save.image( paste0("output/mm_", dbug, "/fhi__mmoverlord__", dbug, "__workspace.RData") )
      
