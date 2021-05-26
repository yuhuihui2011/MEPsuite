### translate transcript sequences to proteins
trans2prot<-function(bg, genome, mex, minsize=50, mincov=1, MSTRG='At'){
    require(BSgenome)
    require(GenomicRanges)
    require(Biostrings)
    require(ORFik)
    txNames<-function(x) sub("_\\d+$","",names(x)) # rewrite txNames
    
    trans<-structure(bg)$trans
    names(trans)<-transcriptNames(bg)
    cov<-texpr(bg,'cov')
    trans<-trans[rowMaxs(cov)>=mincov,]
    trans<-GenomicRanges::sort(trans)
    minus<-which(unlist(runValue(strand(trans)=='-')))
    trans[minus] <- revElements(trans[minus])
    names(trans)<-sub('MSTRG',MSTRG,names(trans))
    
    orfs<-findMapORFs(trans,genome[trans],startCodon = 'ATG',
                      minimumLength = minsize, groupByTx = FALSE)
    size<-sum(width(orfs))
    idx<-tapply(size, txNames(orfs), function(x) names(x[which.max(x)]))
    orfs<-orfs[names(orfs)%in%idx]
    names(orfs)<-txNames(orfs)
    
    idx<-match(names(orfs@unlistData),sub('^MSTRG',MSTRG,transcriptNames(bg)))
    orfs@unlistData$gene<-geneIDs(bg)[idx]
    orfs@unlistData$gene<-sub('^MSTRG',MSTRG,orfs@unlistData$gene)
    orfs@unlistData$mex<-overlapsAny(orfs@unlistData,mex,type='equal')
    
    mex_orfs<-orfs[names(orfs) %in% names(orfs@unlistData[orfs@unlistData$mex]) ]
    prot<-Biostrings::translate(genome[orfs],if.fuzzy.codon='solve')
    list(prot=prot,trans=trans,orfs=orfs,mex_orfs=mex_orfs)
}

