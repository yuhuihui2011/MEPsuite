### map gapped PWM to a genome
mapPWM.0<-function(cons,exons,genome,min.score="80%",include.intronLoss=TRUE,
                   span=20000, min.intron=20, max.intron=10000, 
                   prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = T),
                   cores=1){
    require(Biostrings)
    require(GenomicRanges)
    pwm<-PWM(cons,prior.params=prior.params)
    pwm1<-PWM(cons[,seq.int(exons[1]),drop=FALSE],prior.params=prior.params)
    pwm2<-PWM(cons[,seq.int(exons[2])+exons[1],drop=FALSE],prior.params=prior.params)
    pwm3<-PWM(cons[,seq.int(exons[3])+sum(exons[1:2]),drop=FALSE],prior.params=prior.params)
    plus_intron1<-sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    plus_intron2<-sprintf('^(G[TC].{%d,%d}AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    minus_intron1<-sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    minus_intron2<-sprintf('^(CT.{%d,%d}[GA]C)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    
    plus<-mclapply(genome, function(x) {
        e1<-suppressWarnings(matchPWM(pwm1,x,min.score))
        e3<-suppressWarnings(matchPWM(pwm3,x,min.score))
        ov<-findOverlapPairs(e1,e3,maxgap=span)
        ov<-ov[start(ov@second)>=end(ov@first)+exons[2]]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first;e3<-ov@second
        intervals<-pgap(e1,e3)
        me<-lapply(seq_along(intervals), function(i) {
            xx<-Views(x,intervals[i])
            m<-suppressWarnings(matchPWM(pwm2,xx,min.score))
            if(length(m)==0) return(IRanges())
            l<-grepl(plus_intron1,flank(m,start(m)-start(xx)),perl=T)
            if (!any(l)) l<-grepl(plus_intron2,flank(m,start(m)-start(xx)),perl=T)
            r<-grepl(plus_intron1,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            if (!any(r)) r<-grepl(plus_intron2,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            m<-m[l&r]
            if(length(m)==0) return(IRanges())
            e2<-m@ranges
            mcols(e2)$NT<-DNAStringSet(paste0(e1[i],m,e3[i]))
            mcols(e2)$AA<-Biostrings::translate(mcols(e2)$NT,if.fuzzy.codon='solve')
            mcols(e2)$region<-IRanges(start(e1[i]),end(e3[i]))
            mcols(e2)$block.starts<-as.matrix(data.frame(
                start(e1)[i],end(e1)[i]+1,
                start(e2),end(e2)+1,start(e3)[i],fix.empty.names=F))
            mcols(e2)$block.sizes<-as.matrix(data.frame(
                exons[1],distance(e1[i],e2),
                exons[2],distance(e2,e3[i]),exons[3],fix.empty.names=F))
            return(e2)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    plus<-unlist(IRangesList(plus[elementNROWS(plus)>0]))
    if (length(plus)==0) plus<-GRanges() else plus<-GRanges(names(plus),plus,'+')
    
    minus<-mclapply(genome, function(x) {
        e1<-suppressWarnings(matchPWM(reverseComplement(pwm1),x,min.score))
        e3<-suppressWarnings(matchPWM(reverseComplement(pwm3),x,min.score))
        ov<-findOverlapPairs(e1,e3,maxgap=span)
        ov<-ov[start(ov@first)>=end(ov@second)+exons[2]]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first;e3<-ov@second
        intervals<-pgap(e1,e3)
        me<-lapply(seq_along(intervals), function(i) {
            xx<-Views(x,intervals[i])
            m<-suppressWarnings(matchPWM(reverseComplement(pwm2),xx,min.score))
            if(length(m)==0) return(IRanges())
            l<-grepl(minus_intron1,flank(m,start(m)-start(xx)),perl=T)
            if (!any(l)) l<-grepl(minus_intron2,flank(m,start(m)-start(xx)),perl=T)
            r<-grepl(minus_intron1,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            if (!any(r)) r<-grepl(minus_intron2,flank(m,end(xx)-end(m),start=FALSE),perl=TRUE)
            m<-m[l&r]
            if(length(m)==0) return(IRanges())
            e2<-m@ranges
            mcols(e2)$NT<-reverseComplement(DNAStringSet(paste0(e3[i],m,e1[i])))
            mcols(e2)$AA<-Biostrings::translate(mcols(e2)$NT,if.fuzzy.codon='solve')
            mcols(e2)$region<-IRanges(start(e3[i]),end(e1[i]))
            mcols(e2)$block.starts<-as.matrix(data.frame(
                start(e1)[i],end(e2)+1,start(e2),end(e3)[i]+1,start(e3)[i],
                fix.empty.names=FALSE))
            mcols(e2)$block.sizes<-as.matrix(data.frame(
                exons[1],distance(e1[i],e2),exons[2],distance(e2,e3[i]),
                exons[3],fix.empty.names=FALSE))
            return(e2)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    minus<-unlist(IRangesList(minus[elementNROWS(minus)>0]))
    if (length(minus)==0) minus<-GRanges() else minus<-GRanges(names(minus),minus,'-')
    
    res<-suppressWarnings(append(plus,minus))
    res<-res[grep("\\*",res$AA,invert = TRUE)]
    if(length(res)>0) {
        res$score<-round(sapply(res$NT, function(x) PWMscoreStartingAt(pwm,x)),4)
        cluster<-reduce(GRanges(seqnames(res),res$region,strand(res)),with.revmap=TRUE)$revmap
        if(any(elementNROWS(cluster)>1)) {
            res<-res[unlist(cluster[lapply(extractList(res$score,cluster), which.max)])]
        }
        res$isMicroExon<-res$block.sizes[,2]!=0 & res$block.sizes[,4]!=0
    }
    cat ("---Number of sequences found:\t", length(res),'\n')
    unname(res)
}

mapPWM<-function(cons,exons,genome,focus=2, min.score='80%',include.intronLoss=TRUE,
                 span=20000, min.intron=20, max.intron=10000, 
                 prior.params=letterFrequency(genome[[1]],DNA_BASES,as.prob = T),
                 cores=1, check.params=TRUE){
    require(Biostrings)
    require(GenomicRanges)
    n<-length(exons)
    if (n==3) {
        return(mapPWM.0(cons,exons,genome,min.score,include.intronLoss, span, 
                        min.intron, max.intron, prior.params, cores))
    }
    if (check.params) {
        if (n<3) stop("'exons' must be a numeric vector of length >=3")
        if (sum(exons)!=ncol(cons))
            stop ("sum(exons)==ncol(cons) is not TRUE")
        if (class(genome)!='DNAStringSet') stop("'genome' must be a DNAStringSet object")
        if (length(focus)!=1 & focus[1]<=1 & focus[1]>=n)
            stop("'focus' must be a numeric in (0,length(exons))!")
    }
    plus_intron1<-sprintf('^(GT[AG].{%d,%d}[CT]AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    plus_intron2<-sprintf('^(G[TC].{%d,%d}AG)%s$',min.intron,max.intron,
                          ifelse(include.intronLoss,"?",""))
    minus_intron1<-sprintf('^(CT[AG].{%d,%d}[CT]AC)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    minus_intron2<-sprintf('^(CT.{%d,%d}[GA]C)%s$',min.intron,max.intron,
                           ifelse(include.intronLoss,"?",""))
    
    pwm<-PWM(cons,prior.params=prior.params)
    for (i in seq(n)) {
        index<-seq.int(exons[i]) + cumsum(c(0,exons))[i]
        assign(paste0('pwm',i),PWM(cons[,index,drop=FALSE],prior.params=prior.params))
    }
    
    plus<-mclapply(genome, function(x) {
        # find matches for the firest and the last exon
        e1<-suppressWarnings(matchPWM(pwm1,x,min.score))@ranges
        en<-suppressWarnings(matchPWM(get(paste0('pwm',n)),x,min.score))@ranges
        # find regions spanning the first exon and the last exon 
        ov<-findOverlapPairs(e1,en,maxgap=span)
        ov<-ov[start(ov@second)-end(ov@first)>=sum(exons[-c(1,n)])]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first; en<-ov@second
        intervals<-pgap(e1,en)
        mcols(intervals)$id<-seq_along(intervals)
        # screen the intervals containing all the other exons
        for(i in seq(2,n-1)) {
            mi<-suppressWarnings(matchPWM(get(paste0('pwm',i)),x,min.score))@ranges
            ov<-findOverlapPairs(mi,intervals,type='within')
            if (length(ov)==0) return(IRanges())
            intervals<-unique(ov@second)
            assign(paste0('m',i),unique(ov@first))
        }
        e1<-e1[mcols(intervals)$id]
        en<-en[mcols(intervals)$id]
        # check the order of the exons and introns
        me<-lapply(seq_along(intervals), function(i) {
            nt<-x[e1[i]] # fist exon
            block<-data.frame(start(e1[i]))
            ir<-intervals[i]
            for (j in seq(2,n-1)) {
                ov<-findOverlapPairs(get(paste0('m',j)),ir,type='within')
                if (length(ov)==0) return(IRanges())
                l<-grepl(plus_intron1,Views(x,flank(
                    ov@first,start(ov@first)-start(ov@second))),perl=T)
                if (!any(l)) l<-grepl(plus_intron2,Views(x,flank(
                    ov@first,start(ov@first)-start(ov@second))),perl=T)
                if (!any(l)) return(IRanges())
                ov<-ov[l]
                nt<-paste0(nt,Views(x,ov@first))
                block<-data.frame(data.frame(lapply(block, rep, length(ov))),
                                  rep(start(ov@first),each=nrow(block)))
                ir<-pgap(ov@first,rep(en[i],length(ov)))
            }
            r<-grepl(plus_intron1,Views(x,ir),perl=TRUE)
            if (!any(r)) r<-grepl(plus_intron2,Views(x,ir),perl=TRUE)
            if (!any(r)) return(IRanges())
            nt<-DNAStringSet(paste0(nt[r],Views(x,en[i])))
            block<-as.matrix(data.frame(block[r,,drop=F],start(en)[i]))
            dimnames(block)<-NULL
            ir<-IRanges(block[,focus],width=exons[focus],
                        NT=nt,AA=Biostrings::translate(nt,if.fuzzy.codon='solve'))
            mcols(ir)$region<-IRanges(start(e1[i]),end(en[i]))
            mcols(ir)$block.starts<-t(apply(block,1,function(xx) {
                head(sort(c(xx,xx+exons)),-1)
            }))
            mcols(ir)$block.sizes<-t(apply(block,1,function(xx) {
                diff(sort(c(xx,xx+exons)))
            }))
            return(ir)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    plus<-unlist(IRangesList(plus[elementNROWS(plus)>0]))
    if (length(plus)==0) plus<-GRanges() else plus<-GRanges(names(plus),plus,'+')
    
    minus<-mclapply(genome, function(x) {
        # find matches for the firest and the last exon
        e1<-suppressWarnings(matchPWM(reverseComplement(pwm1),x,min.score))@ranges
        en<-suppressWarnings(matchPWM(reverseComplement(
            get(paste0('pwm',n))),x,min.score))@ranges
        # find regions spanning the first exon and the last exon 
        ov<-findOverlapPairs(e1,en,maxgap=span)
        ov<-ov[start(ov@first)-end(ov@second)>=sum(exons[-c(1,n)])]
        if (length(ov)==0) return(IRanges())
        e1<-ov@first; en<-ov@second
        intervals<-pgap(e1,en)
        mcols(intervals)$id<-seq_along(intervals)
        # screen the intervals containing all the other exons
        for(i in seq(2,n-1)) {
            mi<-suppressWarnings(matchPWM(reverseComplement(
                get(paste0('pwm',i))),x,min.score))@ranges
            ov<-findOverlapPairs(mi,intervals,type='within')
            if (length(ov)==0) return(IRanges())
            intervals<-unique(ov@second)
            assign(paste0('m',i),unique(ov@first))
        }
        e1<-e1[mcols(intervals)$id]
        en<-en[mcols(intervals)$id]
        # check the order of the exons and introns
        me<-lapply(seq_along(intervals), function(i) {
            nt<-x[e1[i]] # fist exon
            block<-data.frame(start(e1[i]))
            ir<-intervals[i]
            for (j in seq(2,n-1)) {
                ov<-findOverlapPairs(get(paste0('m',j)),ir,type='within')
                if (length(ov)==0) return(IRanges())
                r<-grepl(minus_intron1,Views(x,flank(
                    ov@first,end(ov@second)-end(ov@first),start = FALSE)),
                    perl=T)
                if  (!any(r)) r<-grepl(minus_intron2,Views(x,flank(
                    ov@first,end(ov@second)-end(ov@first),start = FALSE)),
                    perl=T)
                if (!any(r)) return(IRanges())
                ov<-ov[r]
                nt<-paste0(Views(x,ov@first),nt)
                block<-data.frame(data.frame(lapply(block, rep, length(ov))),
                                  rep(start(ov@first),each=nrow(block)))
                ir<-pgap(ov@first,rep(en[i],length(ov)))
            }
            l<-grepl(minus_intron1,Views(x,ir),perl=TRUE)
            if (!any(l)) l<-grepl(minus_intron2,Views(x,ir),perl=TRUE)
            if (!any(l)) return(IRanges())
            nt<-reverseComplement(DNAStringSet(paste0(Views(x,en[i]),nt[l])))
            block<-as.matrix(data.frame(block[l,,drop=F],start(en)[i]))
            dimnames(block)<-NULL
            ir<-IRanges(block[,focus],width=exons[focus],NT=nt,
                        AA=Biostrings::translate(nt,if.fuzzy.codon='solve'))
            mcols(ir)$region<-IRanges(start(en[i]),end(e1[i]))
            mcols(ir)$block.starts<-t(apply(block,1,function(xx) {
                sort(c(xx,xx+exons),decreasing = T)[-1]
            }))
            mcols(ir)$block.sizes<-t(apply(block,1,function(xx) {
                0L-diff(sort(c(xx,xx+exons),decreasing = T))
            }))
            return(ir)
        })
        me<-do.call(c,me[elementNROWS(me)>0])
        me
    },mc.cores = cores)
    minus<-unlist(IRangesList(minus[elementNROWS(minus)>0]))
    if (length(minus)==0) minus<-GRanges() else minus<-GRanges(names(minus),minus,'-')
    
    res<-suppressWarnings(append(plus,minus))
    res<-res[grep("\\*",res$AA,invert = TRUE)]
    if(length(res)>0) {
        res$score<-round(sapply(res$NT, function(x) PWMscoreStartingAt(pwm,x)),4)
        cluster<-reduce(GRanges(seqnames(res),res$region,strand(res)),with.revmap=TRUE)$revmap
        if(any(elementNROWS(cluster)>1)) {
            res<-res[unlist(cluster[lapply(extractList(res$score,cluster), which.max)])]
        }
        res$isMicroExon<-res$block.sizes[,2*(focus-1)]!=0 & res$block.sizes[,2*focus]!=0
    }
    cat ("---Number of sequences found:\t", length(res),'\n')
    unname(res)
}