# qutils

## install
```R
devtools::install_url(" https://github.com/RemyNicolle/qutils/archive/master.tar.gz")
```

## Pathway

Upload pathways using the exCITingpath

```R
gstat=getUniqueGeneVec(selpICA$cor[,proti$i],protanot[rownames(selpICA$cor),"GENE"],abs(selpICA$cor[,proti$i])))
	protstat[intersect(names(protstat),seli(selcICA$cor,cnai$i,100))]
	protgsea=qfgsea(protstat,allpaths,n=10000)
	
	selupgsea=protgsea[which(protgsea$NES>2.5& protgsea$size>=20),]
	qGseaTable(allpaths[selupgsea$pathway],protstat,selupgsea,gseaParam=0.8,colwidths=c(10,2,1,1,1))
```
