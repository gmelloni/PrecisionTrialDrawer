# Trick to fool magrittr in checks
if(getRversion() >= "2.15.1"){
  utils::globalVariables(c(".", "%>%"))
}

# Implementation of a special for a "not in"
'%notin%' <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}

# Implementation of two specials for comparison of vectors with NA
# These functions differ from "==" and "!=" because they can also compare NA
# NA==NA is TRUE and "foo"==NA is FALSE (generally they both return NA)
"%eq%" <- function(x , y) {
    out <- x==y
    nacomp <- is.na(x)==is.na(y)
    out[is.na(out)] <- nacomp[is.na(out)]
    return(out)
}
"%noteq%" <- Negate("%eq%")

#change all factor columns to character
.changeFactor <- function(df) {
  for( i in colnames(df)){
    if(is.factor(df[,i]))
      df[ , i] <- as.character(df[,i])
  }
  return(df)
}

# divide a vector in chunks of the specified length
.subsetter <- function(x , len) {
  subs <- seq(1,length(x),len)
  xBlocks <- lapply(2:length(subs) , function(i) {
                    start <- subs[i-1]
                    end <- subs[i]-1
                    return(x[start:end])
                    })
  xBlocks[[length(xBlocks)+1]] <- x[subs[length(subs)]:length(x)]
  return(xBlocks)
}

# fast trimmer of trailing spaces
.myTrimmer <- function (object, ...) 
{
   s <- sub("^[\t\n\f\r ]*", "", as.character(object))
   s <- sub("[\t\n\f\r ]*$", "", s)
   s
}

# like merge but mantaining the order of the first dataframe
.mergeOrder <- function(df1 , df2 , ...){
  if("INDEXVAR" %in% colnames(df1) | 
     "INDEXVAR" %in% colnames(df2)){
    stop("There should not be any column of df1 with the name INDEXVAR")
  }
  df1$INDEXVAR <- seq_len(nrow(df1))
  out <- merge(x=df1 , y=df2 , ...)
  out <- out[ order(out$INDEXVAR , na.last=TRUE) , ]
  out$INDEXVAR <- NULL
  return(out)
}

# find all duplicated elements in a vector
.superdup <- function(x) {
  duplicated(x) | duplicated(x , fromLast=TRUE)
}

# Substitute all NA of a vector with the specified value, default is ""
.noNA <- function(x , subs="") {
  ifelse(is.na(x) , subs , x)
}

# MFROW autosetter
.mfrow <- function(nplots, ncolPlot=FALSE) {
  
  # If we have specified the number of columns we want
  if(ncolPlot){
    # calculate the number of rows we need
    nrowPlot <- ceiling(nplots/ncolPlot)
    return(c(nrowPlot, ncolPlot))
  } 
  
  if(nplots <= 3) {
    return(c(1, nplots))
  }
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    abs(x - round(x)) < tol
  }
  sqrt_nplots <- sqrt(nplots)
  if(is.wholenumber(sqrt_nplots)){
    return(c(sqrt_nplots , sqrt_nplots))
  } else {
    if(round((sqrt_nplots))^2 > nplots){
      return(c(round(sqrt_nplots) , round(sqrt_nplots)))
    }else{
      return(c(round(sqrt_nplots) , round(sqrt_nplots) + 1))
    }
  }
}

# Split a numeric vector in chunks where the difference between
# the element and the next one is 1
# e.g. if v=c(1,2,3,7,8,10) it returns list("1"=c(1,2,3) , "2"=c(7,8) , "3"=10)
.splitter <- function(v){
  if(length(v)==1){
    return(list("1"=v))
  }
  mydiff <- diff(v)
  mysplit <- c(1 , rep(NA , length(mydiff)))
  for( i in seq_len(length(mydiff))) {
    if(mydiff[i]==1){
      mysplit[i+1] <- mysplit[i]
    } else {
      mysplit[i+1] <- mysplit[i] + 1
    }
  }
  split(v , mysplit)
}

# Plot an empty plot with a text message in the center
.emptyPlotter <- function(message){
      par(mar = c(0,0,0,0))
      plot(c(0, 1), c(0, 1), ann = FALSE, bty = 'n'
           , type = 'n', xaxt = 'n', yaxt = 'n')
      text(x = 0.5, y = 0.5, message
        ,cex = 3, col = "black")
}

# Merge a list of dataframe by a common key
.mergeThemAll <- function(biglist , by=NULL , ...){
  if(length(biglist)==0)
    stop("List of length 0")
  if(length(biglist)==1)
    return(biglist[[1]])
  if(is.null(by)) {
    if(length(biglist)==2)
      return(merge(biglist[[1]] , biglist[[2]]) , ...)
    out <- merge(biglist[[1]] , biglist[[2]] , ...)
    for( i in 3:length(biglist)){
      out <- merge(out , biglist[[i]])
    }
    return(out)
  } else {
    if(length(biglist)==2)
      return(merge(biglist[[1]] , biglist[[2]] , by=by , ...))
    out <- merge(biglist[[1]] , biglist[[2]] , by=by , ...)
    for( i in 3:length(biglist)){
      out <- merge(out , biglist[[i]] , by=by , ...)
    }
    return(out)
  }
}

# Stolen from package plyr
.mapvalues <- function (x, from, to, warn_missing = FALSE) {
  if (length(from) != length(to)) {
    stop("`from` and `to` vectors are not the same length.")
  }
  if (!is.atomic(x)) {
    stop("`x` must be an atomic vector.")
  }
  if (is.factor(x)) {
    levels(x) <- .mapvalues(levels(x), from, to, warn_missing)
    return(x)
  }
  mapidx <- match(x, from)
  mapidxNA <- is.na(mapidx)
  x[!mapidxNA] <- to[mapidx[!mapidxNA]]
  x
}
