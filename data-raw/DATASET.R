
mcpgenes=readRDS("/Users/remy.nicolle/Workspace/PDAC_meta/dev/mcpgenes.rda")

ChanSengYueSigs=qload("/Users/remy.nicolle/Workspace/PDAC_meta/dev/ChanSengYueSigs.RData")
names(ChanSengYueSigs)[c(2,10,1,6)]=c("BasalA","BasalB","ClassicA","ClassicB")

puleoica=qload("/Users/remy.nicolle/Workspace/PDAC_meta/dev/puleoICA.RData")



turleyCaf=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/dev/TurleySupTab5.xlsx")
turleyCaf=split(turleyCaf[,2],turleyCaf[,1])
names(turleyCaf)=c("hCAF0TGFB","hCAF1early","hCAF2il6lif")




myofibtab=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/CAFscript/CAFsigs/jem_20162024_tables1.xlsx",sheet=2)
myofibtab$padj=as.numeric(myofibtab$padj)
myofibtab$log2FoldChange=as.numeric(myofibtab$log2FoldChange)
icaftab=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/CAFscript/CAFsigs/jem_20162024_tables1.xlsx",sheet=4)
icaftab$padj=as.numeric(icaftab$padj)
icaftab$log2FoldChange=as.numeric(icaftab$log2FoldChange)
myovsicaftab=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/CAFscript/CAFsigs/jem_20162024_tables1.xlsx",sheet=3)
myovsicaftab$padj=as.numeric(myovsicaftab$padj)
myovsicaftab$log2FoldChange=as.numeric(myovsicaftab$log2FoldChange)


myofibDEg=intersect(
  myofibtab$id[which(myofibtab$log2FoldChange< -1 & myofibtab$padj <0.05)],
  myovsicaftab$id[which(myovsicaftab$log2FoldChange< -1 & myovsicaftab$padj <0.05)]
)
icafDEg=intersect(
  icaftab$id[which(icaftab$log2FoldChange> 1 & icaftab$padj <0.05)],
  myovsicaftab$id[which(myovsicaftab$log2FoldChange> 1 & myovsicaftab$padj <0.05)]
)




icafsc=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/CAFscript/CAFsigs/215864_2_supp_5552227_ps91bp.xlsx",sheet=1)[,2]
myosc=openxlsx::read.xlsx("/Users/remy.nicolle/Workspace/PDAC_meta/CAFscript/CAFsigs/215864_2_supp_5552227_ps91bp.xlsx",sheet=2)[,2]





# DATASET=list(
.geneSets=unlist(list(
  list(
    iCAFsc=icafsc,
    myCAFsc=myosc
  ),
  Turley=turleyCaf,
  CSYNotta=ChanSengYueSigs

),recursive=F)



.puleoICAgw=puleoica$gw
colnames(.puleoICAgw)=rownames(puleoica$annot)
.puleoICAgw= t(t(.puleoICAgw)*puleoica$annot$orientation)

.mcpgenes=mcpgenes
usethis::use_data(.geneSets,.puleoICAgw,.mcpgenes,internal=T,overwrite=T)
