extract.stats <- function(id , path){
    ## path w/ trailing '/'
    ## DiffRatio max and N>9 based on DR4! --> avoid samples that have big peak at high EditDis (U shape)
    ## use the same spec_strain for mapdam and readdis!
    ## MapDamage for last and first base for node with max N in editDis1-4
    ## uniquePerReference for N>9 for node with max N in editDis1-4
    ind <- id
    out <- list()
    for (run in c('default','ancient')){
        ed.dis <- read.table(paste( path,run,'/editDistance/',ind ,'_editDistance.txt', sep="" ),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # kill last column (">5")
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste( path,run,'/damageMismatch/',ind ,'_damageMismatch.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste( path,run,'/readDist/',ind ,'_readDist.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## get species data
            ed.dis.spec <- ed.dis[ grep(spec,rownames(ed.dis)) , ]
            mp.dam.spec <- mp.dam[ grep(spec,rownames(mp.dam)) , ]
            rd.dis.spec <- rd.dis[ grep(spec,rownames(rd.dis)) , ]
            ## get Spearman: loop through each node in table of given spec
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) > 9 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
            mp.dam.spec.sum <- sum(mp.dam.spec[ rowMax ,"C>T_1"] , mp.dam.spec[ rowMax ,"G>A_20"])
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.sum , read.dis.uniq )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.sum , read.dis.uniq )
            }
        }
    }
    out2 <- do.call(rbind,out)
    colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd')
    return(out2)
}

extract.stats2 <- function(id , path){
    ## version2 is adapted to RMAex v04 which changes the names of the readdis files to alignment incl. some headers. readdistable fun also changed --> v3
    ## path w/ trailing '/'
    ## DiffRatio max and N>9 based on DR4! --> avoid samples that have big peak at high EditDis (U shape)
    ## use the same spec_strain for mapdam and readdis!
    ## MapDamage for last and first base for node with max N in editDis1-4
    ## uniquePerReference for N>9 for node with max N in editDis1-4
    ind <- id
    out <- list()
    for (run in c('default','ancient')){
        ed.dis <- read.table(paste( path,run,'/editDistance/',ind ,'_editDistance.txt', sep="" ),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # kill last column (">5")
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste( path,run,'/damageMismatch/',ind ,'_damageMismatch.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste( path,run,'/readDist/',ind ,'_alignmentDist.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## get species data
            ed.dis.spec <- ed.dis[ grep(spec,rownames(ed.dis)) , ]
            mp.dam.spec <- mp.dam[ grep(spec,rownames(mp.dam)) , ]
            rd.dis.spec <- rd.dis[ grep(spec,rownames(rd.dis)) , ]
            ## get Spearman: loop through each node in table of given spec
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) > 9 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
            mp.dam.spec.sum <- sum(mp.dam.spec[ rowMax ,"C>T_1"] , mp.dam.spec[ rowMax ,"G>A_20"])
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.sum , read.dis.uniq )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.sum , read.dis.uniq )
            }
        }
    }
    out2 <- do.call(rbind,out)
    colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd')
    return(out2)
}

extract.stats3 <- function(id , path){
    ## extract.stats3 updates:
    ## - only report exact node/taxon match. After updating the interrogated Node list we would otherwise extract too many strains!
    ## extract.stats2 updates [partially kept for v3]:
    ## version2 is adapted to RMAex v04 which changes the names of the readdis files to alignment incl. some headers. readdistable fun also changed --> v3
    ## path w/ trailing '/'
    ## DiffRatio max and N>9 based on DR4! --> avoid samples that have big peak at high EditDis (U shape)
    ## use the same spec_strain for mapdam and readdis!
    ## MapDamage for last and first base for node with max N in editDis1-4
    ## uniquePerReference for N>9 for node with max N in editDis1-4
    ind <- id
    out <- list()
    for (run in c('default','ancient')){
        ed.dis <- read.table(paste( path,run,'/editDistance/',ind ,'_editDistance.txt', sep="" ),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # kill last column (">5")
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste( path,run,'/damageMismatch/',ind ,'_damageMismatch.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste( path,run,'/readDist/',ind ,'_alignmentDist.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## get species data
            ## ed.dis.spec <- ed.dis[ grep(spec,rownames(ed.dis)) , ]
            ## mp.dam.spec <- mp.dam[ grep(spec,rownames(mp.dam)) , ]
            ## rd.dis.spec <- rd.dis[ grep(spec,rownames(rd.dis)) , ]
            ed.dis.spec <- ed.dis[ spec , ]
            mp.dam.spec <- mp.dam[ spec , ]
            rd.dis.spec <- rd.dis[ spec , ]

            ## get RatioOfDifferences and N (editdistance0-4 and 0-6)
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) > 9 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
            mp.dam.spec.sum <- sum(mp.dam.spec[ rowMax ,"C>T_1"] , mp.dam.spec[ rowMax ,"G>A_20"])
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.sum , read.dis.uniq )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.sum , read.dis.uniq )
            }
        }
    }
    out2 <- do.call(rbind,out)
    colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd')
    return(out2)
}

extract.stats4 <- function(id , path){
    ## extract.stats3 updates:
    ## - only report exact node/taxon match. After updating the interrogated Node list we would otherwise extract too many strains!
    ## extract.stats2 updates [partially kept for v3]:
    ## version2 is adapted to RMAex v04 which changes the names of the readdis files to alignment incl. some headers. readdistable fun also changed --> v3
    ## path w/ trailing '/'
    ## DiffRatio max and N>9 based on DR4! --> avoid samples that have big peak at high EditDis (U shape)
    ## use the same spec_strain for mapdam and readdis!
    ## MapDamage for last and first base for node with max N in editDis1-4
    ## uniquePerReference for N>9 for node with max N in editDis1-4
    ind <- id
    out <- list()
    for (run in c('default','ancient')){
        ed.dis <- read.table(paste( path,run,'/editDistance/',ind ,'_editDistance.txt', sep="" ),header=F,sep="\t",skip=1,row.names=1,comment.char='')[,1:6] # kill last column (">5")
        colnames(ed.dis) <- c('0','1','2','3','4','5') # R does not like reading numeric header's
        mp.dam <- read.table(paste( path,run,'/damageMismatch/',ind ,'_damageMismatch.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        rd.dis <- read.table(paste( path,run,'/readDist/',ind ,'_alignmentDist.txt', sep="" ),header=T,sep="\t",row.names=1,check.names=F,comment.char='')
        if( run == 'ancient' ){ rhoNM <- 2:6 ;keptDiff <- 1:2} else { rhoNM <- 1:6; keptDiff <- 1:3}
        for ( spec in unq.spec ){
            ## get species data
            ## ed.dis.spec <- ed.dis[ grep(spec,rownames(ed.dis)) , ]
            ## mp.dam.spec <- mp.dam[ grep(spec,rownames(mp.dam)) , ]
            ## rd.dis.spec <- rd.dis[ grep(spec,rownames(rd.dis)) , ]
            ed.dis.spec <- ed.dis[ spec , ]
            mp.dam.spec <- mp.dam[ spec , ]
            rd.dis.spec <- rd.dis[ spec , ]

            ## get RatioOfDifferences and N (editdistance0-4 and 0-6)
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            res <- matrix(ncol=5,nrow=nrow(ed.dis.spec)); rownames(res) <- rownames(ed.dis.spec); colnames(res) <- paste(run,c('node','dr6','n6','dr4','n4'),sep=".")
            for (subset in rownames(ed.dis.spec)){
                a <- diff(as.numeric(ed.dis.spec[ subset , rhoNM ]))
                dr6 <- round(sum(abs(a[a<0]))/sum(abs(a)),3)
                b <- a[ keptDiff ] # only diffs 1:3 considered. When ancient only diffs 1:2
                dr4 <- round(sum(abs(b[b<0]))/sum(abs(b)),3)
                res[subset, ] <- c( subset, dr6 , sum(ed.dis.spec[ subset , 1:6 ]) , dr4 , sum(ed.dis.spec[ subset , 1:4 ]))
            }
            ## require minimum of 10 reads present dr4 analysis and pick the one with highest number of reads
            ## NOTE: Could be shorter if indeed always only 1 Node presented!
            rowMax <- which(as.numeric(res[,paste(run,'.n4',sep="")])==max(as.numeric(res[,paste(run,'.n4',sep="")])))[1]
            if( !is.na(rowMax) & as.numeric(res[ rowMax , paste(run,'.n4',sep="") ]) >= 5 ){
                top.dr <- res[ rowMax , ]
            } else {
                top.dr <- rep(NA,5)
            }
            ## extract map.damage:sum(C>T or G>A pos 1 2 (-1 -2 respectively)) for TopScorer@EdDis
            mp.dam.spec.max <- max(mp.dam.spec[ rowMax ,"C>T_1"] , mp.dam.spec[ rowMax ,"C>T_2"] , mp.dam.spec[ rowMax ,"G>A_20"], mp.dam.spec[ rowMax ,"G>A_19"])
            ## extract max readDis:uniquePerReference for TopScorer@EdDis
            read.dis.uniq <- rd.dis.spec[ rowMax ,'uniquePerReference']
            if(length(read.dis.uniq) == 0){ read.dis.uniq <- NA }
            ## write results list
            if( paste(ind,spec,sep="_") %in% names(out) ){
                out[[ paste(ind,spec,sep="_") ]] <- c( out[[ paste(ind,spec,sep="_") ]] , top.dr , mp.dam.spec.max , read.dis.uniq )
            } else {
                out[[ paste(ind,spec,sep="_") ]] <- c( ind , spec , top.dr , mp.dam.spec.max , read.dis.uniq )
            }
        }
    }
    out2 <- do.call(rbind,out)
    colnames(out2) <- c('id','spec','def.node','def.dr6','def.n6','def.dr4','def.n4','def.mapDam','def.rd','anc.node','anc.dr6','anc.n6','anc.dr4','anc.n4','anc.mapDam','anc.rd')
    return(out2)
}


plot.editDis.perID <- function(id,tax,folders){
    ## function writes 8 barplots that has the edit distance and %identity plots for four different RMAex runs (default, ancient, clipped, complexity)
    ## function takes id of individual and of species (strains will be grepped too)
    ## folders should be a vector with 4 entries (each the RMAex output folders) with names
    for (i in 1:length(folders)){
        ## Edit Distance: read data and grep taxon(s) of interest
        res <- read.table(paste(folders[i],'editDistance/',id,'_editDistance.txt',sep=''),header=F,sep="\t",skip=1,comment.char='',row.names=1) # R does not get the numbers as headers
        colnames(res) <- c('0','1','2','3','4','5','>5')
        res <- res[ grep(tax,rownames(res)) , ]
        mc1 <- c(colorRampPalette(c("lightgreen","darkgreen"))(nrow(res)))
        ## %Identity: read data and grep taxon(s) of interest
        res2 <- read.table(paste(folders[i],'percentIdentity/',id,'_percentIdentity.txt',sep=''),comment.char='',header=F,sep="\t",skip=1,row.names = 1) #auto-assigns spec/strains as row names.
        colnames(res2) <- c('77.5-82.5%','82.5-87.5%','87.5-92.5%','92.5-97.5%','>97.5%') # that might be made dynamic in a future version, when trailing tabs are gone!
        res2 <- res2[ grep(tax,rownames(res2)) , ]
        mc2 <- c(colorRampPalette(c("lightcyan","midnightblue"))(nrow(res2)))
        ## plot
        barplot(as.matrix(res),beside=T,col=mc1,main=paste(names(folders)[i],id,tax),xlab="Edit Distance",ylab="Read count")
        legend("top",legend=paste(rownames(res),rowSums(res)),fill=mc1, cex = 0.6, bty = "n")
        barplot(as.matrix(res2),beside=T,col=mc2,main=paste(names(folders)[i],id,tax),xlab="Percent Identity",ylab="Read count")
        legend("top",legend=paste(rownames(res2),rowSums(res2)),fill=mc2, cex = 0.6, bty = "n")
    }
} 

plot.editDis.cov <- function(id,tax,folders){
    ## function writes barplot EditDistance and Coverage Distribution
    ## function takes id of node
    swtch <- FALSE
    for (i in 1:length(folders)){
        ## Edit Distance: read data and grep taxon(s) of interest
        res <- read.table(paste(folders[i],'editDistance/',id,'_editDistance.txt',sep=''),header=F,sep="\t",skip=1,comment.char='',row.names=1) # R does not get the numbers as headers
        rownames(res) <- chartr("><","..",rownames(res)) # Unusual character fix by JFY, based on plot_summary_rmaex_v05
	colnames(res) <- c('0','1','2','3','4','5','>5')
        res <- res[ tax , ]
        mc1 <- c(colorRampPalette(c("lightgreen","darkgreen"))(nrow(res))) # remnant from multiple spec
        ## Coverage: read data and grep taxon(s) of interest
        res2 <- read.table(paste(folders[i],'readDist/',id,'_coverageHist.txt',sep=''),comment.char='',header=F,sep="\t",skip=1,row.names = 1) #auto-assigns spec/strains as row names.
        rownames(res2) <- chartr("><","..",rownames(res2)) # Unusal character fix by JFY, based on plot_summary_rmaex_v05
	colnames(res2) <- c('Reference','0','1','2','3','4','5','6','7','8','9','10','>10') # that might be made dynamic in a future version, when trailing tabs are gone!
        res2 <- res2[ tax , ]
        res2[ res2==0 ] <- 1 # turn 0to1 (for log10 plotting)
        mc2 <- c(colorRampPalette(c("lightcyan","midnightblue"))(nrow(res2))) # remnant from multiple spec
        ## plot
        barplot(as.matrix(res),beside=T,col=mc1,main=paste(names(folders)[i],id,tax),xlab="edit distance",ylab="read count")
        legend("top",legend=paste('Node:',tax,sum(res)),fill=mc1, cex = 0.6, bty = "n")
        if (swtch==FALSE){
            barplot(as.matrix(log10(res2[-1])),beside=T,col=mc2,main=paste(names(folders)[i],id,tax),xlab="RMAex default: coverage distribution",ylab="bases covered (log10)")
            legend("top",legend=paste('Top Ref:',res2[,'Reference']),fill=mc2, cex = 0.6, bty = "n")
            swtch <- TRUE
        }
    }
} 

plot.mapDamage <- function(id,tax,folder){
    ## Plot damage pattern observed w/ RMAex
    ## Note all species and strain nodes hit are collapsed here for simplicity. Pattern might be idfferential affected by different nodes!
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ## dam[ dam[,21]=='null' , 21] <- NA
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    Nstrains <- nrow(dam)
    dam <- colSums(dam,na.rm=T)
    ## check if mapDamage was calculated at all (apparenlty sometimes everything is 0)
    ## We plot the cumulative damage pattern over all species, normalized by the number of strains that have been detected. The numbers provided by RMAex are fullly corrected per strain!
    if(dam[idxN]){C2T <- dam[1:10]/Nstrains} else { return(paste(id,tax,'no reads considered for mapDamage calculation.')) }
    nonC2T <- dam[21:30]/Nstrains # nonC2T are all substitutions that are non-C2T (which are 11 options, considering also directionality)
    if(dam[idxN]){G2A <- dam[11:20]/Nstrains} else { return(paste(id,tax,'no reads considered for mapDamage calculation.')) }
    nonG2A <- dam[31:40]/Nstrains # nonG2A are all substitutions that are non-G2A (which are 11 options, considering also directionality)
    ## ## plot C2T
    ## plot( nonC2T ,col="grey",ylim=c(0,(max(C2T)+0.1)), type="l",ylab="C2T rate",main=paste("Damage Pattern based on ",dam[idxN]," reads",sep=""))
    ## lines(C2T, col='red',lwd=2)
    ## ## plot G2A
    ## if(dam[idxN]){G2A <- dam[11:20]/dam[idxN]} else { return(paste(id,tax,'no reads considered for mapDamage calculation.')) }
    ## nonG2A <- (dam[31:40]/11)/dam[idxN] # nonG2A are all substitutions that are non-G2A (which are 11 options, considering also directionality)
    ## plot( nonG2A ,col="grey",ylim=c(0,(max(G2A)+0.1)), type="l",ylab="G2A rate",main=paste("Damage Pattern based on ",dam[idxN]," reads",sep=""))
    ## lines(G2A, col='blue',lwd=2)

    ## Plot C>T and G>A in one plot:
    ymax <- max(c(C2T,G2A))+(max(c(C2T,G2A)) * 0.2)
    plot( nonC2T ,col="grey",ylim=c(0,ymax), xlim=c(1,20),type="l",lwd=2,ylab="C2T/G2A rate",xlab="Read Position",xaxt='n',main=paste("Damage Pattern based on ",dam[idxN]," reads",sep=""))
    axis(1, at=seq(1,20,2), labels=c(seq(1,10,2),seq(-10,-1,2)))
    lines(C2T, col='red',lwd=2)
    lines(c(rep("",10),nonG2A),col='grey',lwd=2)
    lines(c(rep("",10),G2A), col='blue',lwd=2)
}


plot.mapDamage2 <- function(id,tax,folder){
    ## Plot damage pattern observed w/ RMAex
    ## Note all species and strain nodes hit are collapsed here for simplicity. Pattern might be idfferential affected by different nodes!
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    dam <- dam[ grep(tax,rownames(dam)) , ]
    maxMatch <- which( dam[, dim(dam)[2] ] == max(dam[ , dim(dam)[2] ]) )[1]
    maxY <- (max(dam[, -dim(dam) ]) * 1.1)
    ## check if mapDamage was calculated at all (apparenlty sometimes everything is 0)
    ## We plot the cumulative damage pattern over all species, normalized by the number of strains that have been detected. The numbers provided by RMAex are fullly corrected per strain!
    ## Plot C>T and G>A in one plot:   
    plot("",col="grey",ylim=c(0,maxY), xlim=c(1,20),type="l",lwd=2,ylab="C2T/G2A rate",xlab="Read Position",xaxt='n')
    for (i in nrow(dam)){
        if(i != maxMatch){
            lines(x=1:10,dam[i,21:30],col='grey')
            lines(x=11:20,dam[i,31:40],col='grey')
            lines(x=1:10,dam[i,1:10],col='pink')
            lines(x=11:20,dam[i,11:20],col='lightblue')
        }
    }
    lines(x=1:10,dam[maxMatch,21:30],col='grey',lwd=1.5)
    lines(x=11:20,dam[maxMatch,31:40],col='grey',lwd=1.5)
    lines(x=1:10,y=dam[maxMatch,1:10],col='red',lwd=1.5)
    lines(x=11:20,dam[maxMatch,11:20],col='blue',lwd=1.5)
    axis(1, at=seq(1,20,2), labels=c(seq(1,10,2),seq(-10,-1,2)))
    legend("top",legend=paste(rownames(dam),dam[, dim(dam)[2] ] ), cex = 0.6, bty = "n")
}


plot.readDisTable <- function(id,tax,folder,runsum){
    ## Plot valuble info about sample-species pair
    ## NOTE: First row in table is species node (however, print name of best match strain!)
    library(gridBase)
    library(gridExtra)
    ## read run sammary
    runsum <- runsum[ grep(tax,rownames(runsum)) , id ]
    ## read stuff from readDis file
    rd <- read.table(paste(folder,'readDist/',id,'_readDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F)
    rd <- rd[ grep(tax,rownames(rd)) , ]
    topNode <- names(runsum)[ which( runsum == max(runsum))[1] ]
    ## read mapDamage stuff
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F)
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    data <- c(topNode,runsum[topNode],rd[topNode,'Reference'],rd[topNode,'TotalReadsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],dam[topNode,'C>T_1'],dam[topNode,'G>A_20'])
    data <- cbind( c('Node','NrunSum','rdNodeRef','rdTotReads','nonDup','readDis','C>T_1','G>A_-1'), data)
    colnames(data)=NULL; rownames(data)=NULL
    plot.new()
    grid.table(data, vp=baseViewports()$figure)
}

plot.readDisTable2 <- function(id,tax,folder){
    ## Plot valuble info about sample-species pair
    ## NOTE: First row in table is species node (however, print name of best match strain!)
    library(gridBase)
    library(gridExtra)
    ## read stuff from readDis file
    rd <- read.table(paste(folder,'readDist/',id,'_readDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    rd <- rd[ grep(tax,rownames(rd)) , ]
    topNode <- rownames(rd)[ which( rd[,'nonDuplicatesonReference'] == max(rd[,'nonDuplicatesonReference']) ) ][1]
    ## read mapDamage stuff
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    data <- c(topNode,rd[topNode,'Reference'],rd[topNode,'TotalReadsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],dam[topNode,'C>T_1'],dam[topNode,'G>A_20'])
    data <- cbind( c('Node','rdNodeRef','rdTotReads','nonDup','readDis','C>T_1','G>A_-1'), data)
    colnames(data)=NULL; rownames(data)=NULL
    plot.new()
    grid.table(data, vp=baseViewports()$figure)
}

plot.readDisTable3 <- function(id,tax,folder){
    ## v3: because Ron changed file names for readDis (alignmentDis) and added nonStacked
    ## Plot valuble info about sample-species pair
    ## NOTE: First row in table is species node (however, print name of best match strain!)
    library(gridBase)
    library(gridExtra)
    ## read stuff from readDis file
    rd <- read.table(paste(folder,'readDist/',id,'_alignmentDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    rd <- rd[ grep(tax,rownames(rd)) , ]
    topNode <- rownames(rd)[ which( rd[,'nonDuplicatesonReference'] == max(rd[,'nonDuplicatesonReference']) ) ][1]
    ## read mapDamage stuff
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    data <- c(topNode,rd[topNode,'Reference'],rd[topNode,'TotalAlignmentsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],rd[topNode,'nonStacked'],dam[topNode,'C>T_1'],dam[topNode,'G>A_20'])
    data <- cbind( c('Node','rdNodeRef','rdTotAlignments','nonDup','readDis','nonStacked','C>T_1','G>A_-1'), data)
    colnames(data)=NULL; rownames(data)=NULL
    plot.new()
    grid.table(data, vp=baseViewports()$figure)
}

plot.readDisTable4 <- function(id,tax,folder){
    ## v4: add length distribution: mean (sd)
    ## v3: because Ron changed file names for readDis (alignmentDis) and added nonStacked
    ## Plot valuble info about sample-species pair
    ## NOTE: First row in table is species node (however, print name of best match strain!)
    library(gridBase)
    library(gridExtra)
    ## read stuff from readDis file
    rd <- read.table(paste(folder,'readDist/',id,'_alignmentDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    rd <- rd[ tax , ]
    topNode <- rownames(rd)
    ## read mapDamage stuff
    dam <- read.table(paste(folder,'damageMismatch/',id,'_damageMismatch.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    idxN <- dim(dam)[2]
    dam <- dam[ grep(tax,rownames(dam)) , ]
    ## length distribution info
    ld <- read.table(paste(folder,'readDist/',id,'_readLengthDist.txt',sep=''),header=T,row.names=1,check.names=F,stringsAsFactors=F,comment.char='')
    ld <- paste( round(ld[ tax , 'Mean' ],0),' (',round(ld[ tax , 'StandardDev' ],3),')',sep="")
    
    data <- c(topNode,rd[topNode,'Reference'],rd[topNode,'TotalAlignmentsOnReference'],rd[topNode,'nonDuplicatesonReference'],rd[topNode,'uniquePerReference'],rd[topNode,'nonStacked'],round(dam[topNode,'C>T_1'],4),round(dam[topNode,'G>A_20'],4),ld)
    data <- cbind( c('Node','rdNodeRef','rdTotAlignments','nonDup','readDis','nonStacked','C>T_1','G>A_-1','mean length (sd)'), data)
    colnames(data)=NULL; rownames(data)=NULL
    plot.new()
    mytheme <- gridExtra::ttheme_default(base_size=8) # https://github.com/baptiste/gridextra/wiki/tableGrob
    grid.table(data, vp=baseViewports()$figure,theme=mytheme)
}
