######### find microexons in a plant genome
MEPer<-function(genome, min.score='80%',include.intronLoss=TRUE,
                span=20000, min.intron=20, max.intron=10000, cores=1) {
    source('mapPWM.R')
    load('mex_pattern_blast.rda')
    suppressMessages(library(Biostrings))
    suppressMessages(library(GenomicRanges))
    suppressMessages(library(parallel))
    cat(paste(t0<-Sys.time(),'..... started run\n'))
    cat(paste(Sys.time(),'..... loading genome\n'))
    if (class(genome)!='DNAStringSet') {
        if (!is.character(genome)) {
            stop("'genome' must be the path(s) to the fasta file(s) or
                 a 'DNAStringSet' class" )
        }
        genome<-readDNAStringSet(genome)
        names(genome)<-sub("^(\\S+)\\s+.*","\\1",names(genome))
    }
    genome<-genome[width(genome)>10000]
    prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = T)
    cat(paste(Sys.time(),'..... finished loading genome\n'))
    res<-GRanges()
    for(i in seq(nrow(mex_pattern))) {
        cat(sprintf('%s ..... Cluster %d (size: %d; phase: %d; motif: %s)\n',
            Sys.time(),i,mex_pattern[i,1],mex_pattern[i,2],mex_pattern[i,3]))
        cons<- consensusMatrix(unique(c(mex_pattern$nt[[i]],
                                        mex_pattern$blast[[i]])))[DNA_BASES,]
        exons<-mex_pattern$blocks[[i]]
        res_i<-mapPWM(cons,exons,genome,focus=mex_pattern$me_order[i],
                      min.score,include.intronLoss,span, min.intron, max.intron,
                      prior.params=prior.params,cores = cores,check.params = FALSE)
        if (length(res_i)==0) next
        res_i$block.starts<-do.call(c,apply(res_i$block.starts, 1, IntegerList))
        res_i$block.sizes<-do.call(c,apply(res_i$block.sizes, 1, IntegerList))
        res_i$cluster<-i
        suppressWarnings(res<-c(res,res_i))
    }
    cat(paste(t1<-Sys.time(),'..... finished successfully\n'))
    cat(sprintf('### total regions found:\t%d\n', length(res)))
    cat(sprintf('### total time used:\t%s mins\n', 
                round(difftime(t1,t0,units = 'mins'),2)))
    res
}

