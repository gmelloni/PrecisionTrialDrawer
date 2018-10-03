.checkGene_to_geneID <- function(genes, myUni, myAlias) {
  # This function checks the gene ids provided by the user
  # transforming eventual aliases to official HugoSymbols,
  # removing duplicated itmes and returning the corresponding EntrezID
  myAliasUnmapped <- myAlias[ myAlias$MappedByLowMACA=="no" , ]
  myAlias <- myAlias[ myAlias$MappedByLowMACA=="yes" , ]
  # make all symbols uppercase in order to avoid
  # ambiguities
  genes <- toupper(genes)
  # good genes are considered when provided as EntrezID or HugoSymbols
  good_genes <- c(as.character(myUni$Gene_Symbol), 
                  as.character(myUni$Entrez))
  if( all(genes %in% good_genes) ) {
    # message("All Gene Symbols correct!")
    # collect the annotations of id provided as gene symbols
    isGeneSymbol <- myUni$Gene_Symbol %in% genes
    Official_gs <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
    # collect the annotations of id provided as EntrezID
    isEntrez <- as.character(myUni$Entrez) %in% genes
    Official_entrez <- myUni[isEntrez, c("Gene_Symbol", "Entrez")]
    # merge the two annotations
    Official <- rbind(Official_gs, Official_entrez)
    Official$Alias <- rep(NA, nrow(Official))
    # check for duplicated itmes, in case there are make a
    # warning with the duplicated items
    if( any(duplicated(Official$Entrez)) ) 
    {
      warning("Either there were duplicated Gene Symbols or Entrez IDs 
              or you put a Gene Symbol along with its Entrez ID:"
              , immediate.=TRUE)
      print(Official[duplicated(Official$Entrez), ])
    }
    # remove duplicated items, drop the factors' levels 
    # and return
    Official <- unique(Official[, c("Gene_Symbol", "Entrez")])
    Official <- droplevels(Official)
    return(Official)
  } else {
    # collect the annotations of id provided as gene symbols
    isGeneSymbol <- myUni$Gene_Symbol %in% genes
    Official_gs <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
    # collect the annotations of id provided as EntrezID
    isEntrez <- as.character(myUni$Entrez) %in% genes
    Official_entrez <- myUni[isEntrez, c("Gene_Symbol", "Entrez") ]
    # merge the two annotations
    Official <- rbind(Official_gs, Official_entrez)
    Official$Alias <- rep(NA, nrow(Official))
    # find genes who were not provided as either EntrezID or HugoSymbol
    notOfficial <- setdiff(genes, good_genes)
    # try to assign them through alias
    bad_alias <- setdiff(notOfficial, myAlias$Alias)
    if( length(bad_alias)==0 ) 
    {
      notOfficial <- lapply(notOfficial, function(x) {
        gs <- unique(myAlias[myAlias$Alias==x,"Official_Gene_Symbol"])
        isGeneSymbol <- myUni$Gene_Symbol %in% gs
        official <- myUni[isGeneSymbol, c("Gene_Symbol", "Entrez")]
        official$Alias <- rep(x, nrow(official))
        return(official)
      })
      if( all(vapply(notOfficial, nrow , numeric(1))==1)) {
        out <- do.call("rbind", notOfficial)
        message("These Genes were reverted to their official Gene Symbol:")
        print(out)
        out <- rbind(Official, out)
        if( any( duplicated(out$Entrez) ) ) {
          warning("There were duplicated Gene Symbols or Entrez IDs 
                  or you put a Gene Symbol along with its Entrez ID:"
                  , immediate.=TRUE)
          print(out[duplicated(out$Entrez), ])
        }
        out <- unique(out[, c("Gene_Symbol", "Entrez")])
        out <- droplevels(out)
        return(out)
        # return(droplevels( unique(out[, c(1, 2)]) ) )
      } else {
        message("There is an ambiguity with some aliases:")
        bad_alias_2 <- lengths(notOfficial)!=1
        bad_alias_3 <- do.call("rbind", notOfficial[bad_alias_2])
        print( bad_alias_3 )
        message("Choose a correct Gene Symbol and start over :(")
        return( bad_alias_3 )
      }
    } else {
      wrongGenes <- !( bad_alias %in% 
                         c(myAliasUnmapped$Official_Gene_Symbol ,
                                        myAliasUnmapped$Alias) )
      if( any(wrongGenes) ) {
        wrongGenes <- bad_alias[wrongGenes]
        message("There are invalid Gene Symbol or Entrez IDs:")
        print(wrongGenes)
        stop("Check manually and start over :(")                
      }
      unmappedGenes <- ( bad_alias %in% 
                           c(myAliasUnmapped$Official_Gene_Symbol ,
                                          myAliasUnmapped$Alias) )
      if( any(unmappedGenes) ) {
        unmappedGenes <- bad_alias[unmappedGenes]
        message("There are valid genes that have not been mapped by LowMACA:")
        print(unmappedGenes)
        stop("We are sorry, remove these genes and start over :(")
      }
    }
    }
}


.clustalOAlign <- function(genesData)
{
  seq_names <- rownames(genesData)
  len <- nchar(as.character(genesData[1, 'AMINO_SEQ']))
  aln <- data.frame(Gene_Symbol=rep(genesData[, 'Gene_Symbol'], len)
                    , Protein=rep(genesData[, 'UNIPROT'])
                    , Entrez=rep(genesData[, 'Entrez'])
                    , Align=seq_len(len)
                    , Ref=seq_len(len)
  )
  see_aln <- NULL
  score <- NA
  # Fasta length
  aln <- data.frame(Gene_Symbol=rep(genesData[ , 'Gene_Symbol'] , len)
                    , domainID=rep(genesData[ , 'Pfam_ID'])
                    , Entrez=rep(genesData[ , 'Entrez'])
                    , Envelope_Start=rep(genesData[ , 'Envelope_Start'])
                    , Envelope_End=rep(genesData[ , 'Envelope_End'] , len)
                    , Align=seq_len(len)
                    , Ref=seq_len(len)
  )
  score <- NA
  return(list(ALIGNMENT=aln, SCORE=score, CLUSTAL=see_aln))
}

#-------------------------------------------
# TAKE repos and accept mutations from PTD
#-------------------------------------------
.getLocalGeneMutations <- function(myGenes=myGenes
    ,localData=NULL  
    ,mutation_type=c("missense", "all", "truncating" , "silent") 
    ,NoSilent=TRUE 
    ,Exonic=TRUE
    ,tumor_type="all_tumors"
)
{
  if( is.null(localData) ) stop('no local file provided')
  ## all mutatios from local data
  mut <- localData
  #Delete all non-SNVs mutation and all non-TCGA MutationType
  mut <- mut[ !is.na(mut$Mutation_Type) , ]
  bad_mut_types <- c("Fusion" , "COMPLEX_INDEL" 
                     , "vIII deletion" , "Splice_Site_SNP" , "Indel")
  mut <- mut[ !(mut$Mutation_Type %in% bad_mut_types) , ]
  mut <- mut[ !(mut$Amino_Acid_Change=="MUTATED") , ]
  ## filter: genes
  chosenGenes <- myGenes$Gene_Symbol
  mut <- mut[mut$Gene_Symbol %in% chosenGenes, ]
  ## filter: mutation type
  mutation_user_choiche <- mutation_type[1]
  chosenMutations <- unique(mut$Mutation_Type)
  ## check flags
  if(mutation_user_choiche=="silent") NoSilent=FALSE
  if(NoSilent)
    chosenMutations <- chosenMutations[chosenMutations != "Silent"]
  if(Exonic) {
    notTransc <- c("3'UTR", "3'Flank", "5'UTR", "5'Flank"
                   ,"IGR1", "IGR", "Intron", "RNA", "Targeted_Region")
    chosenMutations <- chosenMutations[!chosenMutations %in% notTransc]
  }
  ## check user choiche
  if( mutation_user_choiche=="missense" ) {
    chosenMutations <- chosenMutations[chosenMutations %in% 
                                         c("Missense_Mutation"
                                           , "In_Frame_Del"
                                           , "In_Frame_Ins")]
  } else if ( mutation_user_choiche=="silent" ) {
    chosenMutations <- chosenMutations[chosenMutations == 'Silent']
  } else if( mutation_user_choiche=="truncating" ) {
    truncating <- c("Frame_Shift_Del"
                    ,"Nonsense_Mutation"
                    ,"Translation_Start_Site"
                    ,"Frame_Shift_Ins"
                    ,"Nonstop_Mutation"
                    ,"Splice_Site"
                    ,"Indel"
    )
    chosenMutations <- chosenMutations[chosenMutations %in% truncating]
  }
  mut <- mut[mut$Mutation_Type %in% chosenMutations, ]
  ## filter: tumor types
  chosenTumors <- unique(mut$Tumor_Type)
  if( tumor_type[1]!="all_tumors" )
    chosenTumors <- chosenTumors[chosenTumors %in% tumor_type]
  mut <- mut[mut$Tumor_Type %in% chosenTumors, ]
  mut <- unique(mut)
  return(list( Mutations=mut , AbsFreq=NA ))
}

.alnWeights <- function(aln)
{
  aln_agg <- aln$ALIGNMENT
  aln_agg$pos_existance <- ifelse(is.na(aln_agg$Ref) , 0 , 1)
  aln_agg2 <- aggregate(pos_existance~Align , data=aln_agg 
                        , FUN=sum , simplify=TRUE)
  return(aln_agg2$pos_existance/sum(aln_agg2$pos_existance))
}

.MAD <- function(x) {
  med <- median(x , na.rm=TRUE)
  MAD <- median( abs(x - med) )
  return(1.4826 * MAD)
}

######################
###### calculate entropy
############################

.shannon <- function(q) {
  diff <- diff(q$x)[1]
  p <- q$y[q$y != 0]
  shan <- -sum(p*log(p))*diff
  return(shan)
}

.profileDensity <- function(profile, bw=NULL)
{
  nPos <- length(profile)
  positions <- which(as.logical(profile))
  positions <- rep(positions, times=profile[positions])
  if( is.null(bw) ) {
    d <- density(positions, from=1, to=nPos, n=nPos)
  } else {
    if( bw==0 ) {
      d <- list(x=seq_len(nPos), y=profile, bw=0)
    } else {
      d <- density(positions, bw=bw, from=1, to=nPos, n=nPos)
    }
  }
  # normalize before output
  d$y <- d$y/sum(d$y)
  return(d)
}

.profileEntropy <- function(profile, bw=NULL, norm=TRUE
                            , model=NULL, weights=NULL, ...) 
{
  d <- .profileDensity(profile, bw=bw, ...)
  ent <- .shannon(d)
  if( is.null(bw) ) bw <- d$bw
  if( norm ) {
    if( !is.null(model) ) {
      unif <- model(sum(profile))
    } else {
      len <- length(profile)
      nmut <- sum(profile)
      unif <- .sampleUnifEntropyL(len, nmut, bw=bw, weights=weights)
    }
    mean <- unif[[3]]-unif[[1]]
    var <- unif[[2]]^2
    ## check: if variance is 0, put the
    ## profile pvalue to zero
    if( var==0 ) {
      pval <- 1
    } else {
      shape <- mean^2/var
      scale <- var/mean
      pval <- pgamma(unif[[3]]-ent, shape=shape, scale=scale, lower.tail=FALSE)            
    }
    return(log10(pval))
  } else {
    return(ent)
  }
}

.sampleUnifEntropyL.old <- function(len, nmut, bw, nboot=1000, weights=NULL, 
                                    center=median, variability=.MAD)
{
  if(is.null(weights)) weights <- rep(1/len , len)
  boots <- vapply(seq_len(nboot), function(i) {
    d <- density(sample(seq_len(len), nmut, replace=TRUE , prob=weights)
                 , bw=bw, from=1, to=len, n=len)
    .shannon(d)
  } , numeric(1))
  return(list(center=center(boots), variability=variability(boots), max=max(boots)))
}

.sampleUnifEntropyL <- function(len, nmut, bw, nboot=1000, weights=NULL, 
                                center=median, variability=.MAD)
{
  if(is.null(weights)) weights <- rep(1/len , len)
  boots <- vapply(seq_len(nboot), function(i) {
    positions <- sample(seq_len(len), nmut, replace=TRUE , prob=weights)
    tab <- table(positions)
    profile <- rep(0, len)
    profile[as.numeric(names(tab))] <- tab
    .profileEntropy(profile, bw=bw, norm=FALSE)
  } , numeric(1))
  return(list(center=center(boots), variability=variability(boots), max=max(boots)))
}

.makeUniformModel <- function(mat, bw, nboot=1000, plotOUT=TRUE, 
                              weights=NULL, center=median, variability=.MAD, parallelize=FALSE ) 
{
  if( parallelize ) {
    applyfun <- mclapply
  } else applyfun <- lapply
  geneLen <- ncol(mat)
  if( is.null(weights) ) weights <- rep(1/geneLen , geneLen)
  minNMut <- floor(sum(mat)/10)*10 #round to the upper ten
  minNMut <- ifelse(minNMut==0 , 1 , minNMut)
  maxNMut <- ceiling(sum(mat)/10)*10 #round to the upper ten
  maxNMut <- ifelse(maxNMut==0 , 1 , maxNMut)
  nMutInt <- unique(c(minNMut, maxNMut))
  if(length(nMutInt)==1) nMutInt <- c(nMutInt , nMutInt+1)
  outReal <- applyfun(nMutInt, function(i) 
    .sampleUnifEntropyL(geneLen, i, bw=bw, nboot=nboot , weights=weights,
                        center=center, variability=variability))
  outReal <- do.call('cbind',outReal)
  polynomialModel <- function(x, par) {
    vapply( x , function(x_i){
        out <- vapply(seq_len(length(par)), function(i) {
            if(is.na(par[i])) 0 else x_i^(i-1) * par[i]
        } , numeric(1))
        sum(out)
    } , numeric(1))
  }
  pn.optim.aic <- function( tpts , experiment, variance=NULL ) {
    if( length(experiment) < 2 ) return(NA)
    polyOrderChisq <- function(i) {
      model <- lm( experiment~poly( tpts , i , raw=TRUE ) )
      return(list(par=model$coeff, value=AIC(model)))
    }
    lapply(seq_len(min(30,length(tpts)-1)), polyOrderChisq) %>% 
      do.call("cbind" , .)
  }
  pnout <- pn.optim.aic(nMutInt, unlist(outReal[1,]), 1)
  degree <- min(which.min(unlist(pnout[2,])))
  par.mean <- pnout[1,degree]$par
  model.mean <- function(mut) polynomialModel(mut, par.mean)
  pnout <- pn.optim.aic(nMutInt, unlist(outReal[2,]), 1)
  degree <- min(which.min(unlist(pnout[2,])))
  par.sd <- pnout[1,degree]$par
  model.sd <- function(mut) polynomialModel(mut, par.sd)
  pnout <- pn.optim.aic(nMutInt, unlist(outReal[3,]), 1)
  degree <- min(which.min(unlist(pnout[2,])))
  par.max <- pnout[1,degree]$par
  model.max <- function(mut) polynomialModel(mut, par.max)
  modelUnif <- function(nmut) 
    list(mean=model.mean(nmut), sd=model.sd(nmut), max=model.max(nmut))
  if( plotOUT ) {
    par(mfrow=c(1,3))
    plot(nMutInt, unlist(outReal[1,]), xlab='n of mutations', ylab='',
         main='entropy center measure')
    lines(nMutInt, model.mean(nMutInt), col='red', lwd=3)
    plot(nMutInt, unlist(outReal[2,]), xlab='n of mutations', ylab='',
         main='entropy variability measure')
    lines(nMutInt, model.sd(nMutInt), col='red', lwd=3)     
    plot(nMutInt, unlist(outReal[3,]), xlab='n of mutations', ylab='',
         main='max entropy measure')
    lines(nMutInt, model.max(nMutInt), col='red', lwd=3)     
  }
  return(modelUnif)
}

.makeNullProfile <- function(mat, bw, nboot=1000, plotOUT=TRUE, 
                             weights=NULL, center=median, variability=.MAD) 
{
  geneLen <- ncol(mat)
  if( is.null(weights) ) weights <- rep(1/geneLen , geneLen)
  nMut <- sum(mat)
  if( bw==0 ) {
    # to prevent having center and variability 
    # that both equals 0
    center <- mean
    variability <- sd
  }
  boots <- lapply(seq_len(nboot), function(i) {
    # density(sample(1:geneLen, nMut, replace=TRUE , prob=weights), bw=bw, from=1, to=geneLen, n=geneLen)
    positions <- sample(seq_len(geneLen), nMut, replace=TRUE , prob=weights)
    t <- table(positions)
    profile <- rep(0, geneLen)
    profile[as.numeric(names(t))] <- t
    .profileDensity(profile, bw=bw)
  })
  nullDensities <- lapply(boots, '[[', 'y') %>% do.call("cbind" , .)
  # calulate parameters for the gamma distribution
  mu <- apply(nullDensities, 1, center)
  s <- apply(nullDensities, 1, variability)
  s2 <- s^2
  # apply gamma distribution to find thresholds
  upperThreshold <- qgamma(.95, shape=mu^2/s2, scale=s2/mu)
  lowerThreshold <- qgamma(.05, shape=mu^2/s2, scale=s2/mu)
  
  # pvalue of every aa
  d <- .profileDensity(colSums(mat), bw=bw) #, from=1, to=geneLen, n=geneLen)
  pvals <- pgamma(d$y, shape=mu^2/s2, scale=s2/mu, lower.tail=FALSE)
  
  
  return(data.frame(
    mean=mu, lTsh=lowerThreshold, uTsh=upperThreshold, profile=d$y, pvalue=pvals#, qvalue=qvals
  ))
}
