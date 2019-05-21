# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

qfgsea=function(signedcor,pathways,n=1000,thresh=0.01){
  require(fgsea)
  set.seed(1)
  xgsea=fgsea(  pathways,signedcor,nperm=n,nproc=1)
  xgsea=xgsea[which(xgsea$pval<thresh),]
  xgsea=xgsea[order(-abs(xgsea$NES)),]
  xgsea
}

qloadPath=function(){

  data("gs.MSigDB.PathwaysOnly.toSymbol.hs_RN_02032016.dump",package="qutils")
  data("KEGG_2016",package="qutils")
  data("WikiPathways_2016",package="qutils")

  all=c(msigdb,NewKegg,WikiPath)
  rm(msigdb,NewKegg,WikiPath,pos=1)
  all
}

qloadImmu=function(){

  newImmuneGenes = c("CXCL9","CXCL10","CXCL16","IFNG","IL15","PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","DPP4","LILRA3","CCL2","VEGFB","VEGFC","PDGFC","CXCL12","TGFB1","TGFB3","LGALS1","C1QB","C1R","C1S","C3","C3AR1","C5AR1","C7","CFD","CFH","CFI","CXCL13","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M","FOXP3","TNFRSF18","CCR7","SELL","CD44","CD70","CD27","CD163","CD80")

  newImmuneGenesGroups = c(rep("T cells chemotaxis and activation",5),rep("T cells specific inhibition",8),"Myeloid cells chemotactism",rep("Angiogenesis",3),rep("Immunosuppression",4),rep("Complement",10),"Tertiary lymphoid structures",rep("Major Histocompatibility Complex I",7),rep("Regulatory T cells",2),rep("Memory and naive T cells",3),rep("T cell survival",2),rep("Macrophages polarization",2))
  newImmuneGenesGroups
}

#pathways=allpaths[selspaths];stats=nediffprotstat;fgseaRes=neprotgsea;gseaParam=1;colwidths = c(10,3, 0.8, 1.2, 1.2);rename=NULL
qGseaTable=function (pathways, stats, fgseaRes, gseaParam = 1,
colwidths = c(5,3, 0.8, 1.2, 1.2),rename=NULL) {
  require(fgsea)
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  n <- length(statsAdj)
  pathways <- lapply(pathways, function(p) {
    unname(as.vector(na.omit(match(p, names(statsAdj)))))
  })
  if(!is.null(rename)&length(rename)==length(pathways)){
      names(rename)=names(pathways)

  }else{
    rename=names(pathways)
    names(rename)=names(pathways)

  }

  ps <- lapply(names(pathways), function(pn) {
    p <- pathways[[pn]]
    annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
    list(textGrob(rename[pn], just = "right", x = unit(0.95, "npc")),
         ggplot() + geom_segment(aes(x = p, xend = p, y = 0,
                                     yend = statsAdj[p],col=p), size = 0.2) + scale_x_continuous(limits = c(0,
                   length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1,
      1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) +
           scale_colour_gradient2(midpoint=n/2,low = '#ca0020', mid="#f7f7f7",  high = '#0571b0') +
           theme(panel.background = element_blank(), legend.position="none", axis.line = element_blank(),
                 axis.text = element_blank(), axis.ticks = element_blank(),
                 panel.grid = element_blank(), axis.title = element_blank(),
                 plot.margin = rep(unit(0, "null"), 4), panel.spacing = rep(unit(0,
                                                                                 "null"), 4)), textGrob(sprintf("%.2f", annotation$NES)),
         textGrob(sprintf("%.1e", annotation$pval)), textGrob(sprintf("%.1e",
                                                                      annotation$padj)))
  })
  rankPlot <- ggplot() + geom_blank() + scale_x_continuous(limits = c(0,
                                                                      length(statsAdj)), expand = c(0, 0)) + scale_y_continuous(limits = c(-1,
                                                                                                                                           1), expand = c(0, 0)) + xlab(NULL) + ylab(NULL) + theme(panel.background = element_blank(),
                                                                                                                                                                                                   axis.line = element_blank(), axis.text.y = element_blank(),
                                                                                                                                                                                                   axis.ticks.y = element_blank(), panel.grid = element_blank(),
                                                                                                                                                                                                   axis.title = element_blank(), plot.margin = unit(c(0,
                                                                                                                                                                                                                                                      0, 0.5, 0), "npc"), panel.spacing = unit(c(0, 0,
                                                                                                                                                                                                                                                                                                 0, 0), "npc"))
  grobs <- c(lapply(c("Pathway", "Gene ranks", "NES", "pval",
                      "padj"), textGrob), unlist(ps, recursive = FALSE), list(nullGrob(),
                                                                              rankPlot, nullGrob(), nullGrob(), nullGrob()))
  grobsToDraw <- rep(colwidths != 0, length(grobs)/length(colwidths))
  gridExtra::grid.arrange(grobs = grobs[grobsToDraw], ncol = sum(colwidths !=
                                                                   0), widths = colwidths[colwidths != 0])


}






#pathwayName="StromaActiv";stats=x;fgseaRez=puleosec;pathwayDatabase=puleoComps
qGseaPlot=function(pathwayName,stats,fgseaRez,pathwayDatabase,main=NULL,signifN=3,gseaParam=1,doprint=T)
{
  data.table::setkey(fgseaRez,"pathway")

  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathwayDatabase[[pathwayName]], names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL

  print("proc")
  if(abs(max(tops))>abs(min(bottoms))){
    upper=max(tops) +0.15
    lower=0- diff/2
    laby=upper-0.05
  }else{
    lower=min(bottoms) -0.1
    upper= diff/2+0.1
    laby=upper-0.05
  }



  g <- ggplot(toPlot, aes(x = x, y = y))  +
    geom_hline(yintercept = max(tops), colour = "red",
               linetype = "dashed") + geom_hline(yintercept = min(bottoms),
                                                 colour = "red", linetype = "dashed") + geom_hline(yintercept = 0,
                                                                                                   colour = "black") + geom_line(color = "green") + theme_bw() +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x,
                                                               y = -diff/2, xend = x, yend = diff/2), size = 0.2) +
    geom_segment(data = data.frame(x = pathway), mapping = aes(x = x,
                                                               y = -diff/2, xend = x, yend = diff/2,colour=x), size = 0.2) +

    scale_colour_gradient2(midpoint=n/2,low = '#ca0020', mid="#f7f7f7",  high = '#0571b0') +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
    labs(x = "rank", y = "enrichment score")+
    guides(colour=FALSE)+
    ylim(c(lower,upper))+
    xlim(0,n+500)+
    annotate("text",x=n/2,y=laby,label=
               paste("NES:",signif(fgseaRez[pathwayName]$NES,signifN),"   p-value:",signif(fgseaRez[pathwayName]$pval,signifN))
    )

  if(!is.null(main)){
    g=g+ggtitle(main)
  }

  g=g+theme_minimal()
  if(doprint){
    print(g)
    invisible(g)
  }else{
    return(g)
  }

}


