###################################################
## Project: genomics                             ##
## Script Purpose: collect the function I wrote  ##
## Data: 2022.05.22                              ##
## Author: Yiming Sun                            ##
###################################################

#' find one-by-one correspondence for peaksets between species.
#' @param ori_GRanges the original peakset to be lifted over, must be a GRanges object.
#' @param UCSC_liftOver_path the path of UCSC liftOver tools.
#' @param chain_file the path of UCSC genome chain file.
#' @param liftOver_mismatch maximum ratio of bases that fail to remap by USCS liftOver.
#' @param length_filter whether do the length filter, default is TRUE.
#' @param length_mismatch maximum ratio of the length difference between the new peak and the original peak.
#' @param chr_filter whether do the chr filter, default is TRUE.
#' @param mapped_chr among which chromosomes should the original peakset be lifted over to?
#' @param overlap_filter whether do the overlap filter, default is TRUE.
#' @param tmp_path the temp file path.
my_unique_peakset_liftover <- function(ori_GRanges,
                                       UCSC_liftOver_path,
                                       chain_file,
                                       liftOver_mismatch = 0.1,
                                       length_filter = TRUE,
                                       length_mismatch = 0.1,
                                       chr_filter = TRUE,
                                       mapped_chr,
                                       overlap_filter = TRUE,
                                       tmp_path){
  
  #load library
  require(rtracklayer)
  require(GenomicRanges)
  require(S4Vectors)
  require(ggplot2)
  require(cowplot)
  
  #check input
  if(class(ori_GRanges) != 'GRanges'){
    stop('ori_GRanges has to be GRanges object!')
  }
  if(!(file.exists(UCSC_liftOver_path))){
    stop('UCSC_liftOver_path do not exist!')
  }
  if(!(file.exists(chain_file))){
    stop('chain_file do not exist!')
  }
  if(!(file.exists(tmp_path))){
    stop('tmp_path do not exists')
  }
  if(sum(GenomicRanges::countOverlaps(query = ori_GRanges,subject = ori_GRanges) > 1)){
    stop('find overlap within original peakset!')
  }
  
  #rename ori_GRanges
  ori_GRanges <- GRanges(seqnames = ori_GRanges@seqnames,ranges = ori_GRanges@ranges,strand = ori_GRanges@strand)
  temp <- paste(ori_GRanges@seqnames,as.character(ori_GRanges@ranges),sep = '-')
  names(ori_GRanges) <- temp
  char <- paste(tmp_path,'ori_peakset.bed',sep = '/')
  rtracklayer::export.bed(object = ori_GRanges,con = char)
  
  #liftOver
  mapped_path <- paste(tmp_path,'lifted_peakset.bed',sep = '/')
  unmapped_path <- paste(tmp_path,'unlifted_peakset.bed',sep = '/')
  cmd <- paste(UCSC_liftOver_path,char,chain_file,mapped_path,unmapped_path,paste0('-minMatch=',(1-liftOver_mismatch)),sep = ' ')
  system(cmd)
  mapped_GRanges <- rtracklayer::import.bed(con = mapped_path)
  if(sum(duplicated(mapped_GRanges$name))){
    stop('ununique liftover occured!')
  }
  print(paste0('ori peakset length: ',length(ori_GRanges),', mapped peakset length: ',length(mapped_GRanges),', ratio: ',round(length(mapped_GRanges)/length(ori_GRanges)*100,digits = 2),'%'))
  
  #chr filter
  if(chr_filter){
    print('')
    print('chr filter...')
    mapped_GRanges <- mapped_GRanges[as.character(mapped_GRanges@seqnames) %in% mapped_chr]
    print(paste0('ori peakset length: ',length(ori_GRanges),', mapped peakset length: ',length(mapped_GRanges),', ratio: ',round(length(mapped_GRanges)/length(ori_GRanges)*100,digits = 2),'%'))
  }
  
  #length filter
  if(length_filter){
    print('')
    print('length mismatch filter...')
    temp <- ori_GRanges[mapped_GRanges$name]
    delta_length <- mapped_GRanges@ranges@width - temp@ranges@width
    delta_length <- delta_length/(temp@ranges@width)
    
    temp <- data.frame(ratio=delta_length)
    temp <- ggplot(temp,aes(x=ratio)) + 
      geom_density() + 
      xlim(c(-1*length_mismatch,1*length_mismatch)) + 
      theme_cowplot() + 
      theme(aspect.ratio = 0.75)
    print(temp)
    delta_length <- c(abs(delta_length) <= length_mismatch)
    mapped_GRanges <- mapped_GRanges[delta_length]
    print(paste0('ori peakset length: ',length(ori_GRanges),', mapped peakset length: ',length(mapped_GRanges),', ratio: ',round(length(mapped_GRanges)/length(ori_GRanges)*100,digits = 2),'%'))
  }
  
  #overlap filter
  if(overlap_filter){
    print('')
    print('overlapped peak filter...')
    mapped_GRanges <- mapped_GRanges[GenomicRanges::countOverlaps(query = mapped_GRanges,subject = mapped_GRanges) == 1]
    print(paste0('ori peakset length: ',length(ori_GRanges),', mapped peakset length: ',length(mapped_GRanges),', ratio: ',round(length(mapped_GRanges)/length(ori_GRanges)*100,digits = 2),'%'))
  }
  
  #return data
  temp <- S4Vectors::SimpleList(
    ori = ori_GRanges[mapped_GRanges$name],
    mapped = mapped_GRanges
  )
  return(temp)
}

#' merge peaksets using bedtools.
#' @param peakset_x the first peakset to be merged, must be a GRanges object.
#' @param peakset_y the second peakset to be merged, must be a GRanges object.
#' @param bedtools_path the bedtools software path.
#' @param d maximum distance between features allowed for features to be merged.
#' @param bedtools_param more parameter for bedtools, must be a single string.
#' @param tmp_path path for temp files while operating.
my_bedtools_merge <- function(peakset_x,
                              peakset_y,
                              bedtools_path,
                              d = 0,
                              bedtools_param = NULL,
                              tmp_path){
  
  #load library
  require(rtracklayer)
  require(GenomicRanges)
  require(ggplot2)
  require(cowplot)
  
  #check input
  if(class(peakset_x) != 'GRanges'){
    stop('peakset_x is not a GRanges object!')
  }
  if(class(peakset_y) != 'GRanges'){
    stop('peakset_y is not a GRanges object!')
  }
  if(!file.exists(bedtools_path)){
    stop('can not find bedtools!')
  }
  if(!file.exists(tmp_path)){
    stop('tmp_path do not exist!')
  }
  
  #append peaksets
  peakset_x <- GenomicRanges::GRanges(seqnames = peakset_x@seqnames,ranges = peakset_x@ranges,strand = peakset_x@strand)
  peakset_y <- GenomicRanges::GRanges(seqnames = peakset_y@seqnames,ranges = peakset_y@ranges,strand = peakset_y@strand)
  peakset <- base::append(peakset_x,peakset_y)
  char <- paste(tmp_path,'append_peakset.bed',sep = '/')
  rtracklayer::export.bed(object = peakset,con = char)
  
  #sort
  temp <- paste(tmp_path,'sorted_append_peakset.bed',sep = '/')
  char <- paste('sort -k1,1 -k2,2n',char,'>',temp,sep = ' ')
  system(char)
  
  #bedtools merge
  if(!(is.null(bedtools_param))){
    char <- paste0(bedtools_path,' merge -d ',d,' ',as.character(bedtools_param))
  }else{
    char <- paste0(bedtools_path,' merge -d ',d)
  }
  char <- paste(char,'-i',temp,'>',paste(tmp_path,'merged_peakset.bed',sep = '/'),sep = ' ')
  system(char)
  
  #return data
  peakset <- rtracklayer::import.bed(con = paste(tmp_path,'merged_peakset.bed',sep = '/'))
  names(peakset) <- paste(peakset@seqnames,as.character(peakset@ranges),sep = '-')
  print(summary(peakset@ranges@width))
  return(peakset)
}

#' a reimplementation of rtracklayer liftOver.
#' @param ori_GRanges the original peakset to be lifted over, must be a GRanges object.
#' @param chain_file the path of UCSC genome chain file.
#' @param merge whether merge the lifted peak pieces into a single one peak?
#' @param workers How many cores to be used for computation? Default use all cores.
#' @param future.globals.maxSize The max size of objects used in paralleled computation, default is 2GB.
my_rtracklayer_liftOver <- function(ori_GRanges,
                                    chain_file,
                                    merge = TRUE,
                                    workers = NULL,
                                    future.globals.maxSize = 2*(1024^3)){
  
  #load library
  require(rtracklayer)
  require(GenomicRanges)
  require(future.apply)
  require(future)
  
  #check input
  if(class(ori_GRanges) != 'GRanges'){
    stop('ori_GRanges has to be GRanges object!')
  }
  if(!(file.exists(chain_file))){
    stop('chain_file do not exists')
  }
  
  options(future.globals.maxSize = future.globals.maxSize)
  if(is.null(workers)){
    workers <- future::availableCores()
  }
  
  #rename ori_GRanges
  ori_GRanges <- GenomicRanges::GRanges(seqnames = ori_GRanges@seqnames,ranges = ori_GRanges@ranges,strand = ori_GRanges@strand)
  names(ori_GRanges) <- paste(ori_GRanges@seqnames,as.character(ori_GRanges@ranges),sep = '-')
  
  #load chain_file
  chain_file <- rtracklayer::import.chain(con = chain_file)
  
  #lift over
  if(merge){
    mapped_GRanges <- rtracklayer::liftOver(x = ori_GRanges,chain = chain_file)
    ori_peak <- names(mapped_GRanges)
    plan(multisession,workers = workers)
    mapped_GRanges <- future.apply::future_lapply(X = ori_peak,FUN = function(x){
      temp <- mapped_GRanges[[x]]
      
      #whether fail to remap
      if(length(temp) == 0){
        return(NULL)
      }
      
      temp <- rtracklayer::as.data.frame(temp)
      
      #whether remap on different chromosomes
      if(length(unique(as.character(temp$seqnames))) != 1){
        return(NULL)
      }
      
      start_idx <- min(c(temp$start,temp$end))
      end_idx <- max(c(temp$start,temp$end))
      temp <- paste0(unique(as.character(temp$seqnames)),':',start_idx,'-',end_idx)
      return(temp)
    })
    plan(sequential)
    names(mapped_GRanges) <- ori_peak
    mapped_GRanges <- base::unlist(x = mapped_GRanges,use.names = TRUE)
    mapped_GRanges <- as(mapped_GRanges,'GRanges')
    mapped_GRanges$ori_peak <- names(mapped_GRanges)
    names(mapped_GRanges) <- NULL
    gc()
    return(mapped_GRanges)
  }else{
    mapped_GRanges <- rtracklayer::liftOver(x = ori_GRanges,chain = chain_file)
    mapped_GRanges <- base::unlist(x = mapped_GRanges,use.names = TRUE)
    mapped_GRanges$ori_peak <- names(mapped_GRanges)
    names(mapped_GRanges) <- NULL
    return(mapped_GRanges)
  }
}
