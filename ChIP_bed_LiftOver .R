if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("liftOver")

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = nrow(TFdata), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar

#parameter setup

all_path = paste0(getwd(), "/")
TFdata <- read.csv(paste0(all_path,"82_TF_MOTIF.csv"), header = TRUE)
#
for( i in 1:nrow(TFdata)){
  tryCatch({
  TF_No <- i #enter TF number
  tfbs.path = paste0(TFdata[TF_No, 4])
  BED <- fread(tfbs.path,showProgress=FALSE)
  if(sum(grep("hg19",BED$V4)) == 0){
    gr <- makeGRangesFromDataFrame(BED, TRUE, seqnames.field="V1" ,start.field="V2", end.field="V3", strand.field = "V6")
    
    seqlevelsStyle(gr) = "UCSC"
    gr19 = liftOver(gr, ch)
    class(gr19)
    
    gr19 = unlist(gr19)
    genome(gr19) = "hg19"
    export.bed(gr,con=paste0(TFdata[TF_No, 3],'.bed'))
  }else{
    export.bed(BED,con=paste0(TFdata[TF_No, 3],'.bed'))
  }
  }, 
  error=function(e){})
  setTxtProgressBar(pb, i)
}
close(pb)


