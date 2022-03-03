

qNumerote=function(n,l=max(nchar(1:n))) {
  str_pad(1:n,width=l,side="left",pad=0)
}



loadREDCAP=function(rscript,dirscript){



  system(paste("Rscript ",system.file("processRedcap.R",package="qutils"),dirscript,rscript))

  invisible(qload(file.path(dirscript,gsub("r$|R$","RData",gsub("_R_","_DATA_",rscript)))))
}

qnarm=narm=function(x)x[which(!is.na(x))]


qemptyrm=emptyrm=function(x)x[which(x!="")]



qload=function (filename) {
  if (file.exists(filename))
    return(eval(parse(text = load(filename))))
  cat(paste("error - function qload : file ", filename,
    " doesn't exist\n"))
  NULL
}


# qfilt=function(x,mincounts=3,minsamples=0.05*ncol(x)){
#   x[which(rowSums(x>=mincounts)>=minsamples),]
# }


#' Title Quantile normalization
#'
#' @param x new expression dataset to be normalied
#' @param ref reference dataset to be used as reference for quantile noralisation
#'
#' @return quntile normalized x
#' @export
#'
#' @examples
quantNorm2ref <- function(x,ref){
  if(nrow(ref)!=nrow(x)){
    stop("both datsets need to have same number of rows")
  }

  xr <- apply(x,2,rank,ties.method="min")
  refs <- data.frame(apply(ref, 2, sort))

  refm=rowMeans(as.matrix(refs))

  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(xr, 2, index_to_mean, my_mean=refm)
  dimnames(df_final) <- dimnames(x)


  return(df_final)
}


#' getUniqueGeneMat: get matrix of unique values (rows) per gene
#'
#' @param m matrix of value to simplify
#' @param g vector of gene names
#' @param w vector of weights, for each redndant gene, the highest weight will be kept
#'
#' @return a matrix m with one value per gene
#' @export
#'
#' @examples
getUniqueGeneMat=function(m,g,w){

  if(!all.equal(nrow(m),length(g),length(w))){
    stop("nrow of m should be equal to lenght of g and w")
  }
  i=order(w,decreasing=T)

  oki=i[which(!duplicated(g[i]) & !g[i]%in% c("---"," ","",NA))]

  okm=m[oki,]
  rownames(okm)=g[oki]
  okm
}

getUniqueGeneVec=function(m,g,w){
  if(!all.equal(length(m),length(g),length(w))){
    stop("length of m should be equal to lenght of g and w")
  }

  i=order(w,decreasing=T)

  oki=i[which(!duplicated(g[i]) & !g[i]%in% c("---"," ","",NA))]

  okm=m[oki]
  names(okm)=g[oki]
  okm
}
