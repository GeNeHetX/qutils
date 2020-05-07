
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

fpDrawSimpleCI =function (lower_limit, estimate, upper_limit, size=1, y.offset = 0.5,
	clr.line, clr.marker, lwd, lty = 1, vertices, vertices.height = 0.1,
	...
){
	if (is.na(lower_limit) || is.na(estimate) || is.na(upper_limit))
	return()
	forestplot:::prFpDrawLine(lower_limit = lower_limit, upper_limit = upper_limit,
		clr.line = clr.line, lwd = lwd, lty = lty, y.offset = y.offset,
		vertices = vertices, vertices.height = vertices.height
	)
	box <- convertX(unit(estimate, "native"), "npc", valueOnly = TRUE)
	skipbox <- box < 0 || box > 1
	if (!skipbox) {
		h=unit(0.2, "snpc")
		w=unit(0.01, "snpc")

		grid.rect(x = unit(estimate, "native"), y = y.offset,
		width =w, height = h, gp = gpar(fill = clr.marker,
			col = clr.marker)
		)
	}
}






multivariateForest=function(acoxph,signifnum=3, ...){
	require(forestplot)
	sumcox=summary(acoxph)
	library(forestplot)
	tabletext=list(
		c(NA,rownames(sumcox$conf.int)),
		c("HR (95% CI)",apply(sumcox$conf.int[,c(1,3,4)],1,cipaste)),
		c("P-value (Wald)",signif(sumcox$coefficients[,5], signifnum))
	)
	toplotvals=rbind(rep(NA,3),sumcox$conf.int[,c(1,3,4)])

	forestplot(tabletext,toplotvals,align=rep("l",5),xlog=T,
	is.summary=c(TRUE,rep(FALSE,length(listCoxph))),
	hrzl_lines = gpar(col="#444444"),graph.pos=2,
	fn.ci_norm=fpDrawSimpleCI,...)
}


univariateForest=function(listCoxph,signifnum=3,printArgHelp=F, ...){
	if(printArgHelp){
		print(paste("clip: upper and lower limit of HR to plots (adds arrow if necessary)",
		"xlab for label of plot",
		"graph.pos : column to put plot of HR",
		"hrzl_lines = list('2'=gpar(col='#444444')) : to add a line under first row    ",
		"col = fpColors(box='royalblue',line='darkblue') : to color plots",
		"w",sep="\n"))
	}
	require(forestplot)

	signifnum=3

	library(forestplot)
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
		}))
	)

	toplotvals=rbind(rep(NA,3),t(sapply(listCoxph,function(x){summary(x)$conf.int[1,c(1,3,4)]})))

	forestplot(tabletext,toplotvals,align=rep("l",5),xlog=T,
	is.summary=c(TRUE,rep(FALSE,length(listCoxph))),
	hrzl_lines = gpar(col="#444444"),graph.pos=2,
	fn.ci_norm=fpDrawSimpleCI,...)
}



dualunivariateForest=function(listCoxph1,listCoxph2,signifnum=3,cols=c("green","brown") ,...){
	signifnum=3
	if(any(names(listCoxph1)!=names(listCoxph2)) ){
		stop("not same names of list of cox")
	}


	tabletext=list(
		c(NA,names(listCoxph1)),
		c("HR (95% CI)",sapply(names(listCoxph1),function(x){
			xx1=summary(listCoxph1[[x]]);
			xx2=summary(listCoxph2[[x]]);
			paste0(
				signif(xx1$conf.int[1,1],signifnum)," (",
				signif(xx1$conf.int[1,3],signifnum),"-",signif(xx1$conf.int[1,4],signifnum),")\n",
				signif(xx2$conf.int[1,1],signifnum)," (",
				signif(xx2$conf.int[1,3],signifnum),"-",signif(xx2$conf.int[1,4],signifnum),")"
			)

		})),
		c("P-value (Wald)",sapply(names(listCoxph1),function(x){
			xx1=summary(listCoxph1[[x]]);
			xx2=summary(listCoxph2[[x]]);
			paste0(as.character(signif(xx1$waldtest[3],signifnum)),"\n",
			as.character(signif(xx2$waldtest[3],signifnum)) )
		})),
		c("P-value (Score)",sapply(names(listCoxph1),function(x){
			xx1=summary(listCoxph1[[x]]);
			xx2=summary(listCoxph2[[x]]);
			paste0(as.character(signif(xx1$sctest[3],signifnum)),"\n",
			as.character(signif(xx2$sctest[3],signifnum)) )
		}))
	)

	toplotvals1=rbind(rep(NA,3),t(sapply(listCoxph1,function(x){summary(x)$conf.int[1,c(1,3,4)]})))
	toplotvals2=rbind(rep(NA,3),t(sapply(listCoxph2,function(x){summary(x)$conf.int[1,c(1,3,4)]})))



	forestplot(tabletext,
		mean = cbind(toplotvals1[,1],toplotvals2[,1]),
		lower = cbind(toplotvals1[,2],toplotvals2[,2]),
		upper =  cbind(toplotvals1[,3],toplotvals2[,3]),
		xlog=T,align=rep("l",5),
		is.summary=c(TRUE,rep(FALSE,length(listCoxph1))),
		hrzl_lines = gpar(col="#444444"),graph.pos=2,
		fn.ci_norm=fpDrawSimpleCI,
		col=fpColors(box=cols),...
	)

}


cipaste=function(xv){
	xv=signif(xv,3)
	paste0(xv[1]," (",xv[2],"-",xv[3],")")
}
