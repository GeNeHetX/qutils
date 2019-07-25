

# setwd("cccp")
# projectTitle="pdxCCCP76s"
# data=subpdx$exp
# weights=w
# meth.Clustering = c("ward.D2", "complete","average")
# meth.Distance = c("pearson", "euclidean","spearman")
# nbootstrap=1000
# pGenesProbes=0.8
# pSamples = 0.8
# qselectd= exp(seq(log(0.01),log(0.2),length.out=10))*nrow(data)
# Kmax=10
# finalLinkage ="complete"
# ANNOT=NULL



CCCP = function(projectTitle,
    data,
    weights=NULL,
    meth.Clustering = c("ward.D2", "complete","average"),
    meth.Distance = c("pearson", "euclidean","spearman"),
    nbootstrap=1000,
    pGenesProbes=0.8,
    pSamples = 0.8,
    qselectd= exp(seq(log(0.01),log(0.5),length.out=10))*nrow(data),
    Kmax=10,
    finalLinkage ="complete",
    ANNOT=NULL){

    require(ConsensusClusterPlus)
    qselectd= ceiling(qselectd)

    cccpproject = paste(projectTitle)
    dir.create(projectTitle,showWarnings=F)
    setwd(projectTitle)
    on.exit(setwd(".."))
    dat.norm=as.matrix(data)
    psok  <- rownames(dat.norm)
    if(is.null(weights)){
        weights = rowMads(dat.norm)
        names(weights) = rownames(dat.norm)
    }
    pdf("ALLCONSENSUS.pdf")
    on.exit({setwd("..");dev.off()})

    customCCCPrun=mclapply(meth.Distance,function(dis){
        print(dis)
        arun=mclapply(meth.Clustering,function(linkage){
            arun=lapply(  qselectd,function(aq){
                titre=paste(projectTitle,aq,dis,linkage,sep="_")
                print(titre)
                if(aq >0 & aq <= 1) psok <- names(weights)[which(weights >= quantile(weights, probs=aq))]   # aq = quantile
                if(aq > 1 & aq <= nrow(dat.norm)) psok <- rownames(dat.norm)[which(rank(weights) >= (length(weights)-(aq-1)))]  # q = number of pbs
                if(aq==0)psok <- rownames(dat.norm)

                    ccp <- ConsensusClusterPlus(dat.norm[psok, ], maxK =Kmax, reps = nbootstrap,pItem=pSamples,pFeature=pGenesProbes,
                        plot=NULL,clusterAlg="hc",verbose=T,title=titre,
                        distance=dis,                  innerLinkage=linkage,
                        finalLinkage=finalLinkage,corUse="pairwise.complete.obs")
                return(lapply(ccp[2:Kmax],function(x)x$consensusMatrix))
            })
            lapply(1:(Kmax-1),function(k){
                X=Reduce('+',lapply(arun,function(x)x[[k]]))/length(qselectd)
                rownames(X)=colnames(dat.norm)
                colnames(X)=colnames(dat.norm)
                X})

        })
        lapply(1:(Kmax-1),function(k){
            X=Reduce('+',lapply(arun,function(x)x[[k]]))/length(meth.Clustering)
            rownames(X)=colnames(dat.norm)
            colnames(X)=colnames(dat.norm)
            X})
    })
    dev.off()
    on.exit(setwd(".."))

    print("OK")

    cccpConsMat=lapply(1:(Kmax-1),function(k){
        X=Reduce('+',lapply(customCCCPrun,function(x)x[[k]]))/length(meth.Distance)
        rownames(X)=colnames(dat.norm)
        colnames(X)=colnames(dat.norm)
        X})

    pdf("consensusClusterings.pdf")
    lapply(1:length(cccpConsMat),function(i){
        if(is.null(ANNOT)){

            heatmap(cccpConsMat[[i]][colnames(data),colnames(data)],scale="none",col=colorRampPalette(c("white","blue"))(32),main=paste("Consensus clustering k=",i+1))
        }else{
            X=data.frame(cccpConsMat[[i]],stringsAsFactors=F)
            rownames(X)=rownames(cccpConsMat[[i]])
            colnames(X)=colnames(cccpConsMat[[i]])
            cit.heatmap(X,scale="none",heatmapcolors=colorRampPalette(c("white","blue"))(32),
                colclust.hclust=cit.hclust(dx=as.dist(1-cccpConsMat[[i]]),method="complete"),colclust.k=i+1,
                rowclust.hclust=cit.hclust(dx=as.dist(1-cccpConsMat[[i]]),method="complete"),rowclust.k=i+1,
                colannot=ANNOT[colnames(cccpConsMat[[i]]),],
            #          lw = c(10, 5, 1, 20), lh = c(5, 5, 1, 20)
                lw = c(12, 7, 1, 16), lh = c(3, 10, 1, 15),
                test="chisq.test"
                )
        }
    })

    dev.off()

    cccpConsHierarchy = lapply(cccpConsMat,function(x)return(hclust(as.dist(1-x),method=finalLinkage)))
    cccpConsPart = lapply(1:(Kmax-1),function(k)return(cutree(cccpConsHierarchy[[k]],k+1)))

    save(cccpConsPart,file=file.path("cccpConsensusPartition.RData"))
    save(cccpConsHierarchy,file=file.path("cccpConsensusHierarchy.RData"))
    save(cccpConsMat,file=file.path("cccpConsensusCoocurenceMatrix.RData"))

    q=qselectd[length(qselectd)]
    if(q >0 & q <= 1) qkeep <- names(weights)[which(weights >= quantile(weights, probs=q))]   # q = quantile
    if(q > 1 & q <= nrow(dat.norm)) qkeep <- rownames(dat.norm)[which(rank(weights) >= (length(weights)-(q-1)))]  # q = number of pbs
    if(q==0)qkeep <- rownames(dat.norm)


        subdat.norm=as.matrix(dat.norm[qkeep,])

    cccpmetrics_subdata     =.cit.clustMetric(subdat.norm,cccpConsMat,cccpConsHierarchy,cccpConsPart,2:Kmax,plot=T,pdf.name="CCCP.clustMetric_subdata.pdf")
    cccpmetrics_alldata     =.cit.clustMetric(dat.norm,cccpConsMat,cccpConsHierarchy,cccpConsPart,2:Kmax,plot=T,pdf.name="CCCP.clustMetric_alldata.pdf")





    return(list(cccpConsPart=cccpConsPart,cccpConsHierarchy=cccpConsHierarchy,cccpConsMat=cccpConsMat,  cccpmetrics_subdata= cccpmetrics_subdata,    cccpmetrics_alldata=cccpmetrics_alldata))

}



.cit.clustMetric=function (expdata, consensuscooccurencematrix, hierarchy, partitions,
    ks = 2:ifelse(is.list(hierarchy), length(hierarchy) + 1,
        length(unique(partitions))), plot = T, pdf.name = "clustMetric.pdf",
    doGAPstat = F)
{
    require(NMF)
    require(SAGx)
    originalDist = dist(t(expdata))
    if (plot)
        pdf(pdf.name)
    if (is.list(consensuscooccurencematrix) & length(consensuscooccurencematrix) ==
        length(hierarchy) & length(consensuscooccurencematrix) ==
        length(partitions)) {
        message("running on a list of clustering results.")
        rez = data.frame(coph = sapply(hierarchy, function(h) {
            return(cor(as.numeric(cophenetic(h)), as.numeric(originalDist)))
        }), disp = sapply(consensuscooccurencematrix, dispersion),
        silhouettes = sapply(partitions, function(part) {
            p <- as.numeric(as.factor(part))
            names(p) <- names(part)
            sil = silhouette(p, dist = originalDist)
            return(summary(sil)$avg.width)
        }))
        if (doGAPstat)
            rez$gap = sapply(partitions, function(part) {
                try(gap(t(expdata), part[colnames(expdata)])[1])
            })
        if (plot) {
            if (doGAPstat) {
                plot(1:length(rez$gap), rez$gap, xlab = "Cluster number",
                  ylab = "Gap statistics", main = "Gap Statistics",
                  type = "b", axes = F, pch = 16)
                axis(2, las = 2)
                axis(1, at = 1:length(rez$gap), labels = paste("K",
                  2:(length(rez$gap) + 1), sep = ""))
            }
            if (require(ggplot2)) {
                y = melt(as.matrix(rez))
                colnames(y) = c("K", "Metric", "value")
                if (!is.factor(y$K))
                  y$K = factor(paste("K=", y$K + 1, sep = ""),
                    levels = paste("K=", ks, sep = ""))
              print(ggplot(y, aes(x = K, y = value, col = Metric)) +
                  geom_point() + theme_classic())
          }
          else {
            matplot(rez, type = "l", lwd = 3, x = 2:Kmax,
              xlab = "K", ylab = "Metric", ylim = c(min(as.matrix(rez)),
                1.1))
            legend("topright", colnames(rez), col = 1:ncol(rez),
              horiz = T, lty = 1:ncol(rez))
        }
        plot(1:length(rez$coph), rez$coph, xlab = "Cluster number",
            ylab = "Cophenetic correlation", main = "Cophenetic coefficient",
            type = "b", axes = F, pch = 16)
        axis(2, las = 2)
        axis(1, at = 1:length(rez$coph), labels = paste("K",
            2:(length(rez$coph) + 1), sep = ""))
        plot(1:length(rez$disp), rez$disp, xlab = "Cluster number",
            ylab = "Dispersion", main = "Dispersion", type = "b",
            axes = F, pch = 16)
        axis(2, las = 2)
        axis(1, at = 1:length(rez$disp), labels = paste("K",
            2:(length(rez$disp) + 1), sep = ""))
        plot(1:length(rez$silhouettes), rez$silhouettes,
            xlab = "Cluster number", ylab = "Mean silhouette",
            main = "Mean silhouettes", type = "b", axes = F,
            pch = 16)
        axis(2, las = 2)
        axis(1, at = 1:length(rez$silhouettes), labels = paste("K",
            2:(length(rez$silhouettes) + 1), sep = ""))
        N <- length(consensuscooccurencematrix)
        colv <- rainbow(N + 1)
        cdf <- .cit.computeCDF(consensuscooccurencematrix[[1]])
        plot(as.numeric(names(cdf)), cdf, xlab = "Consensus index value",
            ylab = "CDF", col = colv[1], type = "l", lwd = 2.5,
            main = "Cumulative Distribution Function")
        for (i in 2:N) {
            tmp <- .cit.computeCDF(consensuscooccurencematrix[[i]])
            lines(as.numeric(names(tmp)), tmp, col = colv[i],
              lwd = 2.5)
        }
        legend("topleft", legend = paste("K=", 2:(N + 1)),
            fill = colv[1:N], bty = "n")
        x = sapply(partitions, function(part) {
            p <- as.numeric(as.factor(part))
            names(p) <- names(part)
            sil = silhouette(p, dist = originalDist)
            plot(sil)
            return(summary(sil)$avg.width)
        })
    }
}
else {
    sil = silhouette(partitions, dist = originalDist)
    rez = data.frame(coph = cor(cophenetic(hierarchy), originalDist),
        disp = dispersion(consensuscooccurencematrix), silhouettes = summary(sil)$avg.width)
    if (doGAPstat)
        rez$gap = gap = try(gap(t(expdata), partitions[names(expdata)])[1])
    if (plot) {
        cdf <- .cit.computeCDF(consensuscooccurencematrix)
        plot(as.numeric(names(cdf)), cdf, xlab = "Consensus index value",
            ylab = "CDF", col = "red", type = "l", lwd = 2.5,
            main = "Cumulative Distribution Function")
        plot(sil)
    }
}
if (plot)
    dev.off()
rownames(rez) = paste("k=", ks, sep = "")
rez
}


.cit.computeCDF=function (coocurencematrix)
{
    if (nrow(coocurencematrix) != ncol(coocurencematrix))
        stop("A co-classification (s x s) matrix is required")
    consind <- seq(0, 1, 0.05)
    coocurencematrix <- as.matrix(coocurencematrix)
    upmat <- coocurencematrix[upper.tri(coocurencematrix)]
    N <- length(upmat)
    rez <- NULL
    for (ci in consind) rez <- c(rez, length(which(upmat <= ci))/N)
        names(rez) <- consind
    rez
}
