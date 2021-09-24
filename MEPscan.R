######### scan microexons in a plant genome
MEPscan<-function(genome, min.score='80%',include.intronLoss=TRUE,
                span=20000, min.intron=20, max.intron=10000, cores=1) {
    source('mapPWM.R')
    load('MEPdata.RData')
    suppressPackageStartupMessages({
        library(Biostrings)
        library(GenomicRanges)
        library(parallel)
    })
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
    for(i in seq(nrow(MEPdata))) {
        cat(sprintf('%s ..... Cluster %d (size: %d; phase: %d; motif: %s)\n',
                    Sys.time(),MEPdata$cluster[i],MEPdata$size[i],MEPdata$phase[i],
                    MEPdata$motif[i]))
        cons<- consensusMatrix(unique(c(MEPdata$nt[[i]],
                                        MEPdata$blast[[i]])))[DNA_BASES,]
        exons<-MEPdata$blocks[[i]]
        res_i<-mapPWM(cons,exons,genome,focus=MEPdata$me_order[i],
                      min.score,include.intronLoss,span, min.intron, max.intron,
                      prior.params=prior.params,cores = cores,check.params = FALSE)
        if (length(res_i)==0) next
        res_i$block.starts<-do.call(c,apply(res_i$block.starts, 1, IntegerList))
        res_i$block.sizes<-do.call(c,apply(res_i$block.sizes, 1, IntegerList))
        res_i$cluster<-MEPdata$cluster[i]
        suppressWarnings(res<-c(res,res_i))
    }
    cat(paste(t1<-Sys.time(),'..... finished successfully\n'))
    cat(sprintf('### total regions found:\t%d\n', length(res)))
    cat(sprintf('### total time used:\t%s mins\n', 
                round(difftime(t1,t0,units = 'mins'),2)))
    res
}

