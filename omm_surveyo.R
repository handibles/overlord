#!/usr/bin/env Rscript

## ==   surveyor  - take a 3-point slope from variation in aggregated MultiQC reports
## + + + + + + + + would need to be very different to use rLib FastQC, but likely doable.
## + + + + + + + + relies on the "fastqc_seq_content_data" field in MQC output. Found circa line ~5720.
## +
## + jfg_init:0ct_23   ::   any chance of an aul docopt?
## + jfg_cont:Dec_24   ::   to AMUKHERJEE ; same as surveyo_0.2. 
## +


## note: assumes libs stringr & reshape2 available

## params   ------------------------------------------------------------

bash_args <- commandArgs(trailingOnly = TRUE)
html_path <- bash_args[1]
## test
# html_path <- "C:/Users/Jamie.FitzGerald.TEAGASCAD/Desktop/fhi__mm2630_fwd_multiqc_raw.html"

print( paste0( "+++++++  working on it :: ", html_path))
## always needs to be some value yr staying below. This one chosen based on appearances
varDiff <- 1.5  
varRat <- 1.5
## lowest possible value for 5prime trim
fiveprime_limit <- 15
## lowest possible value for 3prime trim
threeprime_limit <- 140
## generic term to separate text for each sample
split_term <- '}}, \\"'
## catch here if get nothing, or too much!! handy for testing setup if nothing else
if( !file.exists(html_path) ){ stop(c(" < ! >   the MultiQC Report.html file is a SHAM, and a FAKE, and.. \n\n.. and a SHAM. Check input and try again."))}else{print(c(" + + +   wait two secs while Beast eat Chris?"))}


## wrangle   -----------------------------------------------------------

# I dont know HTML well enough to consider premature EOL significant. plx mail if a bad
xh <-  readLines( html_path, warn = FALSE)[ which( grepl( split_term, readLines( html_path, warn = FALSE))) ]

## belete info from negative controls
xhs_pre <- unlist(strsplit( xh, split_term, perl = TRUE) )
xhs <- xhs_pre[ !grepl( "Neg", xhs_pre, ignore.case = TRUE)]

## break samples down into rows per base-position, based on field match (& drop first mangled instance)
xhs_l1 <- lapply( xhs, function(aa){  # aa <- xhs[3]
  strsplit(aa, '\\"\\d*\\"', perl = TRUE )[[1]][-1]
})

## break rows to df. lazily not keeping the individual sample names.
xhs_l2 <- as.data.frame( do.call( "rbind", lapply( xhs_l1,
                                                   function(aaa){  # aaa <- xhs_l1[65]
                                                     bbb <- unlist( aaa, recursive = TRUE )
                                                     do.call("rbind",
                                                             lapply( bbb, function(cccc){   # cccc <- bbb[65]
                                                               c( 
                                                                 "base" = stringr::str_pad( width = 3, pad = 0,
                                                                                            gsub("(\\d)\\.0", "\\1", gsub('"\\d*-(\\d*)"', '\\1', 
                                                                                                                          gsub('.*"base": (.*?),.*', '\\1', cccc, perl = TRUE),
                                                                                                                          perl = TRUE), perl = TRUE)),
                                                                 "A" = gsub('.*"a": (\\d*\\.\\d*),.*', '\\1', cccc, perl = TRUE),
                                                                 "T" = gsub('.*"t": (\\d*\\.\\d*),.*', '\\1', cccc, perl = TRUE),
                                                                 "G" = gsub('.*"g": (\\d*\\.\\d*).*', '\\1', cccc, perl = TRUE),
                                                                 "C" = gsub('.*"c": (\\d*\\.\\d*).*', '\\1', cccc, perl = TRUE) )
                                                             }))
                                                   })))
## for processing
xhs_l2[,1] <- paste0("base_", stringr::str_pad( width = 3, pad = 0,
                                                gsub( "(\\d)\\.0", "\\1", gsub( '"\\d*-(\\d*)"', '\\1', gsub('.*"base": (.*), "t.*', '\\1',
                                                                                                             xhs_l2[,1],
                                                                                                             perl = TRUE), perl = TRUE), perl = TRUE)))
colnames(xhs_l2) <- c("base", "A", "C", "T", "G")


## get clear drop in variation
agg_val <- aggregate( value ~ base, reshape2::melt(xhs_l2, id.vars = c("base")), FUN = var)


## process   -----------------------------------------------------------

## JE: get the RATIO between points -  might still be the way to go :: abs(agg_val[aa-1, "value"] / agg_val[aa, "value"]) < 1
je <- sapply( 1:nrow(agg_val), function(aa){
  if( aa == nrow(agg_val)){ Inf }else{
    agg_val[aa, "value"] / agg_val[aa+1, "value"]
  }
})


## instead get differences within a threesome: aa-1, aa, aa+1              
variance_ratio_bin <- sapply( 1:nrow(agg_val), function(aa){
  abs(agg_val[aa-1, "value"] / agg_val[aa, "value"]) < varRat && 
    abs(agg_val[aa, "value"] / agg_val[aa+1, "value"]) < varRat &&
    abs(agg_val[aa-1, "value"] / agg_val[aa+1, "value"]) < varRat
})
## ratios perform well - and as expected on trimmed data. Stick?
variance_diff_bin <- sapply( 1:nrow(agg_val), function(aa){
  abs(agg_val[aa-1, "value"] - agg_val[aa, "value"]) < varDiff && 
    abs(agg_val[aa, "value"] - agg_val[aa+1, "value"]) < varDiff &&
    abs(agg_val[aa-1, "value"] - agg_val[aa+1, "value"]) < varDiff
})

## need to get the longest sequential stretch - not beautiful but appears robust
seq_lengs <- do.call("rbind", sapply( 1:length(variance_ratio_bin), function(aa){ # aa <- 10
  do.call( "rbind",
           lapply( aa:length(variance_ratio_bin), function(aaa){   # aaa <- 38
             if( any( is.na(c(variance_ratio_bin[ aa], variance_ratio_bin[aaa]))) ){
               c("start" = aa, "end" = aaa, "length" = 0)
             }else if( all( variance_ratio_bin[ aa : aaa]) ){
               c("start" = aa, "end" = aaa, "length" = c(aaa-aa))
             } })
  )
})
)

# remember, these are the row numbers of agg_val - need base positions:
longest_seq <- seq_lengs[ which( seq_lengs[ , "length"] == max( seq_lengs[,"length"])) ,  ]
(slopes_out <- as.numeric( gsub("base_0*", "", c(agg_val[ longest_seq[1], 1], agg_val[ longest_seq[2], 1]) )) + c(0,-1))    # a slight caution


## diagnose / outputs   -----------------------------------------------

## stout fellow

print( paste0("-+-+-+-"))
if( any( slopes_out > 152)){ 
  warning(" < ! >   reads >150bp detected - we ASSUME this is 150bp illumina, this script should NOT be used for amplicon data")
}else if( any(is.na(slopes_out)) ){ 
  error( paste0( " < ! >   some trimming points were NA, possibly because of bad data munging. Suggested trims: ", paste0(slopes_out, collapse = ", ")))
}

if( slopes_out[1] < fiveprime_limit ){
  print( paste0(" + + +   suggested 5prime trim is too short: ", slopes_out[1], "bp; increasing to ", fiveprime_limit, "bp."))
  slopes_out[1] <- fiveprime_limit
}else if( slopes_out[2] > threeprime_limit ){
  print( paste0(" + + +   suggested 3prime trim is too short: ", slopes_out[2], "bp; increasing to ", threeprime_limit, "bp."))
  slopes_out[2] <- threeprime_limit
}

if( all(sapply( slopes_out,function(aa){ is.numeric(aa) && !is.null(aa) && !is.na(aa)})) && slopes_out[1] < slopes_out[2] ){
  print( paste0(" + + +   5prime trim is ", slopes_out[1], "; 3prime trim is ", slopes_out[2], "."))
  print( paste0(" + + +   (feel free to check the PNG output and decide otherwise)"))
}else{
  stop(" < ! >   5prime/3prime trim points were in a strange (string/NA/Null/error) format: check wrangling of DF and subsequent munging")  # ~~modelling~~
}
print( paste0("-+-+-+-"))


## write to small file : apologies if too basic for yr
try( write.table( x = slopes_out, file = gsub( ".html", "_MQC-trim-points.txt", html_path), sep = "\t", row.names = FALSE, col.names = FALSE), silent = FALSE)

## retrognostic diagram
try_plot <- function(){
  png(gsub(".html", "_MQC-diagnostic.png", html_path), width = 900, height = 1200)
  par(mfrow= c(3,2))
  ## plot some outputs.... 
  plot(je, type = "b", main = "JE ratio of variance in seq comp")
  plot( type = "b", agg_val[,2], main = "aggregated variance in seq comp", xlab = NA, xaxt ="n" )
  axis(1,at=1:nrow(agg_val), labels=agg_val[,1], las=2)
  
  plot(unlist(sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[(aa-1), "value"] - agg_val[aa, "value"]) })), main = "difference between a/a+1/a-1")
  points(unlist( sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[aa, "value"] - agg_val[(aa+1), "value"]) })), col = "red")
  points(unlist( sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[(aa-1), "value"] - agg_val[(aa+1), "value"]) })), col = "blue")
  
  plot(unlist(sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[(aa-1), "value"] / agg_val[aa, "value"]) })), main = "ratio between a/a+1/a-1")
  points(unlist( sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[aa, "value"] / agg_val[(aa+1), "value"]) })), col = "red")
  points(unlist( sapply( 1:nrow(agg_val), function(aa){ abs(agg_val[(aa-1), "value"] / agg_val[(aa+1), "value"]) })), col = "blue")
  # primary
  plot( unlist( sapply( 1:nrow(agg_val), function(aa){ agg_val[ aa, 2] - agg_val[ aa+1, 2] }) ), type = "b", xaxt = "n", xlab = "base position", ylab = "variance in base composition",
        main = paste0("trim along dotted line, at 5prime: ", slopes_out[1], " and 3prime: ", slopes_out[2],
                      ",\nbased on \"flat\" stretch of var RATIOs.\n(n.b.: script allowed a min 5prime trim of ", fiveprime_limit,  "bp)"))  # aa <- 10
  axis(1,at=1:nrow(agg_val), labels= gsub("base_", "", agg_val[,1]), las=2)
  abline(lm( 0~1), lty= 2)
  segments( x0 = longest_seq[1], x1 = longest_seq[2],
            y0 = 0, y1 = 0, lty= 2, col = "green", lwd = 4)
  dev.off()
}
try( try_plot(), silent = FALSE )


slopes_out


