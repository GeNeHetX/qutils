
myUnivariateForest=function(listCoxph,signifnum=3,printArgHelp=F, ...){
	if(printArgHelp){
		print(paste("clip: upper and lower limit of HR to plots (adds arrow if necessary)",
		"xlab for label of plot",
		"graph.pos : column to put plot of HR",
		"hrzl_lines = list('2'=gpar(col='#444444')) : to add a line under first row    ",
		"col = fpColors(box='royalblue',line='darkblue') : to color plots",
		"w",sep="\n"))
	}
	require(forestplot)
	tabletext=list(
		c(NA,names(listCoxph)),
		c("HR (95% CI)",sapply(listCoxph,function(x){
			xx=summary(x);
			paste0(signif(xx$conf.int[1,1],signifnum)," (",
			signif(xx$conf.int[1,3],signifnum),"-",signif(xx$conf.int[1,4],signifnum),")")
	 })),
	 c("P-value (Wald)",sapply(listCoxph,function(x){
		 xx=summary(x); as.character(signif(xx$waldtest[3],signifnum))
	 })),
	 c("P-value (Score)",sapply(listCoxph,function(x){
		xx=summary(x); as.character(signif(xx$sctest[3],signifnum))
	})),
	c("BIC",sapply(listCoxph,function(x){as.character(round(BIC(x),signifnum))}))
	)
	toplotvals=rbind(rep(NA,3),t(sapply(listCoxph,function(x){summary(x)$conf.int[1,c(1,3,4)]})))
	forestplot(tabletext,toplotvals,align=rep("l",5),xlog=T, ...)
}
