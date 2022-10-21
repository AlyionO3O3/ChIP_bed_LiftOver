if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("liftOver")

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch

pb <- txtProgressBar(min = 0,
                     max = nrow(TFdata),
                     style = 3,
                     width = 50,
                     char = "=") 

#parameter setup

all_path = paste0(getwd(), "/")
TFdata <- read.csv(paste0(all_path,"82_TF_MOTIF.csv"), header = TRUE)

#liftover

trans_list <- list()
origi_list <- list()
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
    trans_list <- append(trans_list,TFdata[TF_No, 3])
    export.bed(gr,con=paste0(TFdata[TF_No, 3],'.bed'))
  }else{
    origi_list <- append(origi_list,TFdata[TF_No, 3])
    write.table(BED, paste0(TFdata[TF_No, 3],'.bed'), sep="\t", col.names=FALSE, row.names = FALSE, append = TRUE, quote = FALSE) 
  }
  }, 
  error=function(e){})
  setTxtProgressBar(pb, i)
}
close(pb)

trans_list
origi_list
