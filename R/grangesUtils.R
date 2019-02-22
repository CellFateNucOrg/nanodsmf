# general useful functions for granges

#' Convert DNAStringSet to GenomicRanges
#' @param DNAss A DNAStringSet
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' DNAstringsetToGR(Biostrings::DNAStringSet(x=c(chr1="agctgtagct",chr2="agagagagttt"))
#' @export
DNAStringSetToGR<-function(DNAss,defaultStrand="*") {
  #extract only the first word in fasta '>' row to be the seq name
  chrNames<-sapply(strsplit(names(DNAss)," "),'[[',1)
  gr<-GenomicRanges::GRanges(seqnames=chrNames,
                             ranges=IRanges::IRanges(start=1,width=Biostrings::width(DNAss)),
                             strand=defaultStrand)
  return(gr)
}



#' Convert BSgenome to DNAStringSet
#' @param BSgnm A BSgenome object
#' @return A DNAStringSet object
#' @examples
#' BSgenomeToDNAStringSet(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToDNAStringSet<-function(BSgnm){
  listOfChr<-sapply(seqnames(BSgnm),function(x){BSgnm[[x]]})
  DNAss<-as(listOfChr,"DNAStringSet")
  return(DNAss)
}



#' Convert BSgenome to GenomicRanges
#' @param BSgnm A BSgenome object
#' @param defaultStrand One of "+", "-" or "*" for the strand from which the sequences come from (default="*")
#' @return A GenomicRanges object
#' @examples
#' BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)
#' @export
BSgenomeToGR<-function(BSgnm,defaultStrand="*") {
  gr<-GenomicRanges::GRanges(seqnames=BSgenome::seqnames(BSgnm),
                             ranges=IRanges::IRanges(start=1,end=GenomeInfoDb::seqlengths(BSgnm)),
                             strand="*")
  return(gr)
}



#' Convert UCSC to wormbase chromosome names
#' @param ucscGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans))
#' @export
ucscToWbGR<-function(ucscGR) {
  wbGR<-ucscGR
  GenomeInfoDb::seqlevels(wbGR)<-gsub("chr","",GenomeInfoDb::seqlevels(wbGR))
  GenomeInfoDb::seqlevels(wbGR)<-gsub("M","MtDNA",GenomeInfoDb::seqlevels(wbGR))
  return(wbGR)
}



#' Convert wormbase to ucsc chromosome names
#' @param wbGR A GRanges object with ucsc chomosome names ("chrI"..."chrM")
#' @return A GenomicRanges object
#' @examples
#' wbToUcscGR(ucscToWbGR(BSgenomeToGR(BSgenome.Celegans.UCSC.ce11::Celegans)))
#' @export
wbToUcscGR<-function(wbGR) {
  ucscGR<-wbGR
  GenomeInfoDb::seqlevels(ucscGR)<-gsub("MtDNA","M",GenomeInfoDb::seqlevels(ucscGR))
  GenomeInfoDb::seqlevels(ucscGR)<-past0("chr",GenomeInfoDb::seqlevels(ucscGR))
  return(ucscGR)
}



#' Convert MIndex object to Genomic Ranges object
#' @param mIdx A MIndex object
#' @return A GenomicRanges object
#' @examples
#' mIdxToGR(Biostrings::vmatchPattern("ATTTAGGGTTTTAGAATACTGCCATTAATTAAAAAT",genomeFile))
#' @export
mIdxToGR<-function(mIdx) {
  #function to convert Mindex object (obtained matching patterns on DNAstringset) to genomic ranges
  allGR<-GRanges()
  seqlevels(allGR)<-names(mIdx)
  for (n in names(mIdx)) {
    if (length(mIdx[[n]])>0) {
      gr<-GRanges(seqnames=Rle(c(n),length(mIdx[[n]])),
                     ranges=mIdx[[n]], strand="*")
      allGR<-append(allGR,gr)
    }
  }
  return(allGR)
}
