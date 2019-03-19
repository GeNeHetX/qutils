qnarm=function(x)x[which(!is.na(x))]
qgivename=function(x,y){names(x)=y;x}
qemptyrm=function(x)x[which(x!="")]

qload=function (filename) {
  if (file.exists(filename))
    return(eval(parse(text = load(filename))))
  cat(paste("error - function cit.load : file ", filename,
            " doesn't exist\n"))
  NULL
}


qfilt=function(x,mincounts=3,minsamples=0.05*ncol(x)){
  x[which(rowSums(x>=mincounts)>=minsamples),]
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
