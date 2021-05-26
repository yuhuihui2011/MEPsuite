####### get the position and the phase of miroexons in CDS and proteins
tag_size<-108
options(stringsAsFactors = F)
sci<-c('C.reinhardtii','P.patens','S.moellendorffii','P.somniferum',
       'A.thaliana','G.max','V.vinifera','H.annuus','O.sativa','Z.mays')
sn<-sub("^(\\w).(\\w)\\w+","\\1\\2",sci)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(parallel)
mex_CDS<-mclapply(seq_along(sci), function(i) {
    genome<-readDNAStringSet(paste0('~/igv/genomes/',sci[i],'.fa'))
    names(genome)<-gsub("^(\\S+)\\s+.*","\\1",names(genome))
    load(paste0('../proteins/',sn[i],'_prot.rda'))
    orfs<-ORFik::sortPerGroup(get(paste0(sn[i],'_mex_orfs')))
    mex<-lapply(orfs,function(x) {
        end<-cumsum(width(x))
        start<-c(0,head(end,-1))+1
        x$phase<- (start -1) %% 3
        x$cds_start<-start
        x$aa_start<- (start -1) %/% 3 +1
        x[x$mex]
    })
    mex<-unique(GRangesList(mex)@unlistData)
    mex$species<-sci[i]
    mex$cds<-DNAStringSet(as.character(genome[orfs[names(mex)]]))
    mex$prot<-Biostrings::translate(mex$cds)
    mex$mex_aa<-subseq(mex$prot,mex$aa_start,
                       mex$aa_start+ceiling((width(mex)+mex$phase)/3)-1)
    ir<-resize(IRanges(mex$cds_start,width=width(mex)),tag_size,fix = 'center')
    shift<-(start(ir)-1)%%3
    #shift<-(width(mex)+2*mex$phase)%/%2%%3
    shift[shift==1]<- -1; shift[shift==2]<-1
    ir<-shift(ir,shift)
    ir<-IRanges(pmax(1,start(ir)),pmin(width(mex$cds),end(ir)))
    mex$nt<-DNAStringSet(substr(mex$cds,start(ir),end(ir)))
    names(mex$nt)<-paste0(names(mex$nt),':',start(ir),'-',end(ir))
    mex$aa<-Biostrings::translate(mex$nt)
    
    mex$up<- mex$cds_start-start(ir)
    mex$down<- end(ir)-(mex$cds_start+width(mex)-1)
    mex$blocks<-IntegerList(lapply(seq_along(mex),function(j) {
        x<-mex[j]
        y<-orfs[[names(x)]]
        z<-IRanges(end=cumsum(width(y)),width = width(y))
        mcols(z)$id<-seq_along(z)
        z<-restrict(z,start(ir)[j],end(ir)[j])
        width(z)[width(z)>0]
    }))
    cat(sci[i],'\n')
    mex
},mc.cores = 3)
names(mex_CDS)<-sci
mex_CDS<-suppressWarnings(GRangesList(mex_CDS))
mex_prot<-mex_CDS@unlistData$prot
writeXStringSet(mex_prot[!duplicated(names(mex_prot))],'all_mex_prot.fa')
save(mex_CDS,file='mex_CDS.rda')

###### search protein motifs (https://www.genome.jp/tools/motif/)
motif<-read.table('all_mex_prot.motif',as.is = T)
colnames(motif)<-c('name','database','motif','pos','query','subject','score','evalue')
idx<-grep(',',motif$evalue)
motif1<-motif[-idx,]
motif2<-do.call(rbind,lapply(idx, function(i) {
    x<-motif[i,]
    data.frame(name=x$name,database=x$database,motif=x$motif,
               pos=unlist(strsplit(x$pos,',')),
               query=unlist(strsplit(x$query,',')),
               subject=unlist(strsplit(x$subject,',')),
               score=unlist(strsplit(x$score,',')),
               evalue=unlist(strsplit(x$evalue,',')))
}))
motif<-rbind(motif1,motif2)
motif$start<-as.numeric(sapply(strsplit(motif$pos,'\\.\\.'),'[',1))
motif$end<-as.numeric(sapply(strsplit(motif$pos,'\\.\\.'),'[',2))
motif$score<-as.numeric(motif$score)
motif$evalue<-as.numeric(motif$evalue)
motif<-split(motif,motif$name)
save(motif,file='all_mex_prot_motif.rda')

##### extract motifs encoded by microexons
load('mex_CDS.rda')
load('all_mex_prot_motif.rda')
library(parallel)
mex_motif<-motif[match(names(mex_CDS@unlistData),names(motif))]
mex_motif<-do.call(rbind,mclapply(seq_along(mex_motif), function(i) {
    x<-mex_motif[[i]]
    if (is.null(x)) {
        return(data.frame(motif='Unknown',evalue=NA,start=NA, end=NA))
    }
    y<-mex_CDS@unlistData[i]
    x<-x[x$start-5<=y$aa_start & x$end+5>=y$aa_start+nchar(y$mex_aa),,drop=FALSE]
    if (nrow(x)==0) {
        return(data.frame(motif='Unknown',evalue=NA,start=NA, end=NA))
    }
    idx<-which.min(x$evalue)
    x<-x[idx,c('motif','evalue','start','end'),drop=FALSE]
    return(x)
},mc.cores=3))
mex_CDS@unlistData$motif<-mex_motif$motif
mex_CDS@unlistData$evalue<-mex_motif$evalue
mex_CDS@unlistData$motif.start<-mex_motif$start
mex_CDS@unlistData$motif.end<-mex_motif$end
save(mex_CDS,file='mex_CDS_motif.rda')

##### classify microexons according to the size and phase, and the miroexon-tags
load('mex_CDS_motif.rda')
library(Biostrings)
library(GenomicRanges)

mex<-mex_CDS[-1]@unlistData #remove C.reinhardtii
mex$size<-width(mex)
mex$id<-seq(mex)
mex$evalue[is.na(mex$evalue)]<-Inf
mex_class<-data.frame(size=mex$size,phase=mex$phase)
mex_class<-unique(mex_class)
mex_class<-mex_class[order(mex_class$size,mex_class$phase),]
rownames(mex_class)<-seq_len(nrow(mex_class))

mex_pattern<-NULL;
min.species<-3; 
cut<-50
for(i in seq(nrow(mex_class))){
    x<-mex[mex$size==mex_class$size[i] & mex$phase==mex_class$phase[i]]
    if (length(unique(x$species))<min.species) next
    s<-outer(translate(x$nt),translate(x$nt),function(x1,x2) {
        pairwiseAlignment(x1,x2,substitutionMatrix='BLOSUM62',
                          scoreOnly=T,gapOpening=100)
    })
    # plot(hclust(as.dist(-s)))
    cluster<-cutree(hclust(as.dist(-s)),h=-cut)
    for(j in unique(cluster)){
        y<-x[cluster==j]
        if (length(unique(width(y$aa)))>1) {
            size.mode<-names(sort(-table(width(y$aa))))[1]
            y<-y[width(y$aa)==size.mode]
        }
        if (length(unique(y$species))<min.species) next
        cons<-DNAString(consensusString(y$nt,threshold=1e-10))
        pos<-sort(unique(unlist(cumsum(y$blocks))))
        pattern<-cons
        for(n in seq_along(head(pos,-1))) {
            pattern<-append(pattern,".",pos[n]+n-1)
        }
        blocks<-diff(c(0,pos))
        res<-data.frame(size=mex_class[i,1],phase=mex_class[i,2],
                        motif=y$motif[which.min(y$evalue)],
                        species=length(unique(y$species)),cases=length(y),
                        exons=length(pos),blocks=NA,
                        me_order=which(pos==y$up[1])+1,
                        pattern=DNAStringSet(pattern))
        res$blocks<-IntegerList(blocks)
        res$aa<-list(y$aa)
        res$nt<-list(y$nt)
        res$rep<-paste0(DNA_BASES[max.col(t(
            consensusMatrix(y$nt)[DNA_BASES,]))],collapse = '')
        mex_pattern<-rbind(mex_pattern,res)
    }
    cat('\r',i)
}
sum(mex_pattern$cases)

mex_pattern<-mex_pattern[-46, ]  # remove low complexity pattern
rownames(mex_pattern)<-seq(nrow(mex_pattern))

seq<-DNAStringSet(mex_pattern$rep)
names(seq)<-seq_along(seq)
writeXStringSet(seq,'mex_pattern_rep.fa')
save(mex_pattern,file='mex_pattern.rda')


