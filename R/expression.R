

qNumerote=function(n,l=max(nchar(1:n))) {
  str_pad(1:n,width=l,side="left",pad=0)
}

qNormalize=function(x,type){
  switch(type,
         raw={x},
         gsc={t(scale(t(x)))},
         gc={t(scale(t(x),scale=F))},
         sc={(scale((x),scale=F))},
         ssc={(scale((x)))}
         ,x
  )
}



qICA=function(X,k,maxiter = 10^6,eps = 10^-6){
  resICA=NULL
  try({
    resICA <- JADE::JADE(X, n.comp = k, maxiter = maxiter, eps = eps)
    rownames(resICA$A) <- colnames(X)
    rownames(resICA$S) <- rownames(X)
    colnames(resICA$S)=paste("ICA",qutils::qNumerote(ncol(resICA$S)),sep="")
    colnames(resICA$A)=paste("ICA",qutils::qNumerote(ncol(resICA$A)),sep="")
  })
  resICA
}




qProjICA=function(icarez,dataset,geneNormType="sc",projNormType="raw",removeZeros=F){

  # expg=getUniqueGeneMat(dataset$exp,dataset$probeannot[rownames(dataset$exp),dataset$genecol],rowSds(as.matrix(dataset$exp)))
  # comg = intersect(rownames(icarez$S), rownames(expg))
  scexp = qutils::qNormalize(dat$exp[rownames(icarez$S), ],type=geneNormType)


  invs = MASS::ginv(as.matrix(icarez$S))

  proj=t(qutils::qNormalize(invs %*%  scexp,projNormType))

  # proj = scale(t(t((t(scexp) %*% t(invs)))))

  colnames(proj)=colnames(icarez$S)
  proj

}
