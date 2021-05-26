#' Get internal exons from a ballgown object
#'
#' Extract the coordinates of internal exons and the flanking
#' introns from a ballgown object
#'
#' @param bg a ballgown object
#' @param genome a DNAStringSet object of genomic sequences
#' @param min.junction set the minimual junction reads for the flanking introns 
#' @param samples a subset of samples
#' @return a list of exons and introns.

getex<-function(bg, genome, min.junction=5,samples=sampleNames(bg)){
    library(ballgown)
    library(GenomicRanges)
    library(Biostrings)
    library(BSgenome)
    library(matrixStats)
    exon<-structure(bg)$exon
    intron<-structure(bg)$intron
    ie<-subset(iexpr(bg,'rcount'),select=sampleNames(bg) %in% samples)
    int<-intron[rowMaxs(ie)>=min.junction]
    int_seq<-genome[int]
    int<-int[(subseq(int_seq,1,2)=='GT' | subseq(int_seq,1,2)=='GC') &
                 subseq(int_seq,width(int_seq)-1,width(int_seq))=='AG']
    names(int)<-int$id
    
    ex<-exon
    ie<-rowMeans(ie[match(int$id,rownames(ie)),])
    l<-findOverlaps(int,GenomicRanges::shift(ex,-width(ex)),type = 'end')
    r<-findOverlaps(int,GenomicRanges::shift(ex,width(ex)), type = 'start')
    l<-sapply(split(l@from,l@to), function(x) x[which.max(ie[x])])
    r<-sapply(split(r@from,r@to), function(x) x[which.max(ie[x])])
    
    ex$int1<-l[match(seq(ex),names(l))]
    ex$int2<-r[match(seq(ex),names(r))]
    ex<-ex[!is.na(ex$int1) & !is.na(ex$int2)]
    ex$int1<-int[ex$int1]$id
    ex$int2<-int[ex$int2]$id
    names(ex)<-ex$id
    ex$gene<-geneIDs(bg)[sapply(ex$transcripts,function(x) {
        eval(parse(text=x))[1]
    })]
    int$gene<-geneIDs(bg)[sapply(int$transcripts,function(x) {
        eval(parse(text=x))[1]
    })]
    list(int=int,ex=ex)
}