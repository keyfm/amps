#!/usr/bin/env Rscript
library(getopt); # for args parsing
library(parallel); # sumstat calculation is parallized across libraries

## internal functions
source('/projects1/users/key/anc5h/scripts.backup/functions_anc5h.r')
source('/Users/fmk/Documents/shh/sysEva/scripts/amps/postprocessing.v5.r')

## INFO
## This scripts gathers signatures of species presence at nodes interrogated by RMAextractor.
## Evidence is plotted in a heatmap for all samples and solely the samples with species specific evidence.
## For all sample-species pairs with eidence the profile-signature-pdf is written.

## USAGE
## Input requires outpath of the RMAextractor "def_anc" run and the taxon-of-interest list (node.list) [e.g. RMAex taxon list ] .
## Rscript ./plot_summary_rmaex_vLatest.r -h

## get options, using the spec as defined by the enclosed list.
## we read the options from the default: commandArgs(TRUE).
spec = matrix(c(
    "rmaex.out.fld",  "r" , 1, "character", "RMAextractor output folder (incl. trailing '/').",
    "help"    ,  "h" , 0, "logical", "Print this help.",
    "node.list"   ,  "n" , 1, "character","List (\\n separated) of nodes/taxons to be reported on (e.g. input taxon list used for RMAex)."
), byrow=TRUE, ncol=5);
opt = getopt(spec);
## if help was asked for print a friendly message
## and exit with a non-zero error code
if ( !is.null(opt$help) ) {
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}

## assign args and modify node.vec (tr ' ' '_')
path <- opt$rmaex.out.fld
unq.spec <- unique(gsub(" ","_",scan(file=opt$node,sep="\n",what='character'))) # unique is solely sanity control ;)

## START DATA PROCESSING
all.inds <- colnames(as.matrix(read.table(paste(path,'/default/RunSummary.txt',sep=''),sep="\t",header=T,stringsAsFactors=F,row.names=1,check.names=FALSE,comment.char='')))

### Extract MetaData for all Sample-Species Pairs
out.lists <- mclapply(1:length(all.inds), function(j) extract.stats3( all.inds[j],path ), mc.cores=16 )
data <- do.call(rbind, out.lists)
data <- data.frame(data,stringsAsFactors=F)
data[, c(4:9,11:16) ] = apply(data[ , c(4:9,11:16)], 2, function(x) as.numeric(as.character(x)))


## Step1: DiffRatio0-4: > 0.9
## Step2: Terminal Damage Present
## Step3: DiffRatio1-4: > 0.8

## Step1: DiffRatio0-4: > 0.9
trg1 <- data[ data[,'def.dr4'] >= 0.9 & !is.na(data[,'def.dr4']) , ]

## Step2: Terminal Damage Present
trg2 <- data[ data[,'def.mapDam'] > 0 & !is.na(data[,'def.mapDam']) , ]

## Step3: DiffRatio1-4: > 0.8
trg3 <- data[ data[,'anc.dr4'] > 0.8 & !is.na(data[,'anc.dr4']) , ]

#### Build Matrix for Heatmap
res <- matrix(1L,nrow=length(unq.spec),ncol=length(all.inds),dimnames=list(a=unq.spec,b=all.inds))
for (p in rownames(trg1)){
    if( !p %in% rownames(trg2) ){
        res[ trg1[p,'spec']  , trg1[p,'id'] ] <- 2
    } else if( !p %in% rownames(trg3) ){
        res[ trg2[p,'spec']  , trg2[p,'id'] ] <- 3
    } else {
        res[ trg3[p,'spec']  , trg3[p,'id'] ] <- 4
    }
}


##############
## Plotting
##############
## plot reduced overview heatmap only for samples and tax with evidence
mycol=c('lightgray','yellow','orange','red')
red.res <- res[, colSums(res) > dim(res)[1] , drop = FALSE ]
red.res <- red.res[ rowSums(red.res) > dim(red.res)[2] , , drop = FALSE ]

## NT & organelle results hat <> in species name (equus), caused bug
rownames(red.res)=chartr("><","..",rownames(red.res))

pdf.height <- max(dim(red.res)[1]/2.5,20)
pdf.width <- max(dim(red.res)[2]/10 , 20)
pdf(paste(path,'/heatmap_overview_Wevid.pdf',sep=""),height=pdf.height,width=pdf.width)
par(mar=c(5.1,30.1,25.1,2.1))
image(x=1:ncol(red.res),y=1:nrow(red.res),z=t(red.res),col=mycol,axes=F,ylab="",xlab="",zlim=c(1,4))
axis(side=2,at=1:nrow(red.res),labels=rownames(red.res),las=1,cex.axis=2)
axis(side=3,at=1:ncol(red.res),labels=colnames(red.res),las=2,cex.axis=2,tick=F)
abline(h=1:length(rownames(red.res))+0.5,col='darkgrey') # add horizontal lines for improved vision
abline(v=1:length(colnames(red.res))+0.5,col='darkgrey') # add vertical lines for improved vision
xleg <- ncol(red.res)-(ncol(red.res)*1.35)
yleg <- nrow(red.res)+5
legend(x=xleg,y=yleg, legend=c('edit distance','+damage','+dam. edit dis.'), fill = mycol[-1],xpd=T,cex=3)
dev.off()


########################
###### Candidate Profile PDFs
########################
## plot summary pdf's for candidates if outpath specified
folder.names <- c(paste(path,'default/',sep=""),paste(path,'ancient/',sep=""))

for (spl in colnames(red.res)){
    for (tax in rownames(red.res)){
        if( red.res[tax,spl] > 1 ){
            system(paste('mkdir -p ',path,'pdf_candidate_profiles/',tax,sep='')) #mk pdf output folder
            pdf(paste(path,'pdf_candidate_profiles/',tax,'/stp',red.res[tax,spl]-1,'_',spl,'_',tax,'_summary.pdf',sep=''))
            par(mfrow=c(3,2))
            plot.editDis.cov(spl,tax,folder.names)
            plot.mapDamage2(spl,tax,folder.names[1])
            plot.readDisTable4(spl,tax,folder.names[1])
            dev.off()
        }
    }
}

#########
## Save RData
#########
save.image(paste(path,'analysis.RData',sep=''))

