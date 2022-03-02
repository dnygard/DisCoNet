# functions and script for NRC project
# reads in network, and tests network (node) similarity across edge inclusion thresholds
#SPECIFIC VERSION FOR CCLE Dataset
#RUNS IN PARALLEL, should run 10 cases at a time to see if that works well enough locally 


library("Rfast")
library(Matrix)
library(SparseM)
library(pracma)
library(igraph)
library(ggplot2)
library(dgof)
library(plyr)
library(doParallel)
registerDoParallel()
# deltacon implementation from https://gist.github.com/bxshi/b9ce0bbabc458447d1d5
# Port of original DeltaCon to R
# Baoxu(Dash) Shi
# Data Sciense Group
# University of Notre Dame

output_time <- function(debug, tim, s) {
  if(debug) {
    print(paste(s, "user time:", tim[1], "system time:", tim[2], "elapsed time:", tim[3]))
  }
}

inverse_lbp <- function(graph, nnodes, priors=NULL, times = 10, debug = FALSE) {
  .MAX_POWER = 10
  .p = 0.51
  # Sparse identity matrix
  I = NULL
  tim <- system.time(
    {
      I <- .sparseDiagonal(nnodes, x = 1, shape = "t")
    })
  output_time(debug, tim, "Create sparse identity matrix")
  
  # Sparse degree-diagonal matrix, D[i,i] = sum(graph[i,])
  x = NULL
  tim <- system.time(
    {
      x <- rowSums(graph, sparseResult = TRUE)
    })
  output_time(debug, tim, "Calculate degree of each node")
  
  D = NULL
  tim <- system.time(
    {
      D <- sparseMatrix(c(1:nnodes), c(1:nnodes), x = x, dims = c(nnodes, nnodes))  
    })
  output_time(debug, tim, "Create degree diagonal matrix")
  
  # Compute about-half homophily factor to guarantee covergence
  c1 = sum(D) + 2
  c2 = sum(D^2) - 1
  h_h = sqrt((-c1 + sqrt(c1^2 + 4 * c2)) / (8 * c2))
  
  # Compute constant ah and ch
  ah = 4 * h_h^2 / (1 - 4 * h_h^2)
  ch = 2 * h_h / (1 - 4 * h_h^2)
  
  # Invert matrices M1 and M2
  M = NULL
  tim <-system.time({
    M = ch * graph  - ah * D
  })
  output_time(debug, tim, "Initialize Invert matrix")
  
  
  # Calculate inverse of M
  if (is.null(priors)) {
    inv_ = I
    mat_ = M
    pow = 1
    tim <- system.time({
      while(max(mat_) > 1e-09 && pow < .MAX_POWER) {
        inv_ = inv_ + mat_
        mat_ = mat_ %*% M
        pow = pow + 1
      }
    })
    output_time(debug, tim, "Invert of matrix")
    return(inv_)
  } else {
    final_mat <- NULL
    tim <- system.time({
      for(i in c(1:times)) {
        inv_ = matrix(priors[, i], nnodes, 1)
        mat_ = matrix(priors[, i], nnodes, 1)
        pow = 1
        while(max(mat_) > 1e-09 && pow < .MAX_POWER) {
          mat_ = M %*% mat_
          inv_ = inv_ + mat_
          pow = pow + 1
        }
        if (i == 1) {
          final_mat <- matrix(inv_)
        } else {
          final_mat <- cbind(final_mat, matrix(inv_))
        }
      }
    })
    output_time(debug, tim, "Approximate invert of matrix")
    return(final_mat)
  }
  
}

init_priors_percent <- function(percent, nnodes) {
  .MAX_POWER = 10
  .p = 0.51
  times <- ceiling(1 / percent)
  init_nodes <- floor(percent * nnodes)
  
  rand_vector <- runif(nnodes) 
  
  rand_mat <- repmat(matrix(rand_vector, nnodes, 1), 1, times)
  
  for(i in c(1:times)) {
    rand_mat[rand_mat[,i] >= (i-1)*percent & rand_mat[,i] < i *percent, i] <- 1
    rand_mat[rand_mat[,i] != 1, i] <- 0
  }
  
  return(rand_mat)
}

delta_con <- function(graph1, graph2, nnodes,
                      method = "naive", percent = 0.1, debug = FALSE, symmetrical = TRUE) {
  .MAX_POWER = 10
  .p = 0.51
  if(ncol(graph1)!=2 || ncol(graph2)!=2) {
    print("Input file should be data.frame with two cols(src dst).")
    return(0);
  }
  colnames(graph1) <- c("src", "dst")
  colnames(graph2) <- c("src", "dst")
  
  # Construct sparse adjacent matrix from edge list
  node_vector <- NULL
  tim <- system.time({
    g1 <- sparseMatrix(graph1$src, graph1$dst, x=1, dims=c(nnodes, nnodes))
    g2 <- sparseMatrix(graph2$src, graph2$dst, x=1, dims=c(nnodes, nnodes))
  })
  output_time(debug, tim, "Construct sparse adjacent matrix")
  
  # Change directed graph to undirected
  if(symmetrical) {
    tim <- system.time({
      g1 <- g1 + t(g1)
      g2 <- g2 + t(g2)
    })
    output_time(debug, tim, "Change to directed graph")
  }
  
  if (method == "fast") {
    # Number of groups to be initialized
    ngroups <- ceiling(1 / percent)
    
    repetitions <- 10
    t_all <- rep(0, repetitions)
    sim <- rep(0, repetitions)
    for(i in c(1:repetitions)) {
      priors <- NULL
      tim <- system.time({
        priors <- init_priors_percent(percent, nnodes)
      })
      output_time(debug, tim, "Calculate priors")
      
      inv1 <- inverse_lbp(g1, nnodes, priors, ngroups, debug = debug) * (.p - 0.5)
      inv2 <- inverse_lbp(g2, nnodes, priors, ngroups, debug = debug) * (.p - 0.5)
      sim[i] <- 1 / 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    }
    delta_con <- mean(sim)
    return(delta_con)
  } else {
    # Naive FaBP
    inv1 <- inverse_lbp(g1, nnodes, debug = debug) * (.p - 0.5)
    inv2 <- inverse_lbp(g2, nnodes, debug = debug) * (.p - 0.5)
    
    # Compute DeltaCon similarity score
    delta_con <- 1 / (1 + sqrt(sum( (sqrt(inv1) - sqrt(inv2))^2 )))
    return(delta_con)
  }
}
# end of deltacon implementation from https://gist.github.com/bxshi/b9ce0bbabc458447d1d5


# Edge similarity based on jaccard index editted from COGENT paper

checkMatrixList <- function(A){
  if (class(A)!="list" | length(A)!=2)
    stop("Input should be a list of length 2.")
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (nrow(A1)!=ncol(A1)| nrow(A2)!=ncol(A2))
    stop('Adjacency matrix not square.')
  if(is.null(colnames(A1))!=is.null(colnames(A2)))
    stop('Label mismatch: one matrix missing node labels.')
  return(TRUE)
}


getEdgeSimilarity <- function(A, align=FALSE, reduce=TRUE){
  check <- checkMatrixList(A)
  if (align){
    A <- alignMatrices(A)
    if (is.null(A)){
      warning('Alignment failed.')
      return(NULL)
    }
  }
  A1 <- A[[1]]; A2 <- A[[2]]; rm(A)
  if (ncol(A1)!=ncol(A2))
    stop('Different network sizes. Consider setting align=TRUE.')
  intersectionMatrix <- pmin(A1, A2)
  intersectionDegrees <- rowSums(intersectionMatrix)
  unionMatrix <- pmax(A1, A2)
  unionDegrees <- rowSums(unionMatrix)
  isolatedNodesIdx <- which(unionDegrees==0)
  if (!reduce)
    isolatedNodesIdx <- c()
  nodeLabels <- colnames(A1)
  nodeCount <- ncol(A1)
  if (length(isolatedNodesIdx)!=0){
    intersectionMatrix <- intersectionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    intersectionDegrees <- intersectionDegrees[-isolatedNodesIdx]
    unionMatrix <- unionMatrix[-isolatedNodesIdx, -isolatedNodesIdx]
    unionDegrees <- unionDegrees[-isolatedNodesIdx]
    nodeLabels <- nodeLabels[-isolatedNodesIdx]
    nodeCount <- nodeCount-length(isolatedNodesIdx)
  }
  localSimilarity <- intersectionDegrees/unionDegrees
  names(localSimilarity) <- nodeLabels
  globalSimilarity <- sum(intersectionDegrees)/sum(unionDegrees)
  # return(list(
  #   "nodeCount"=nodeCount,
  #   "globalSimilarity"=globalSimilarity,
  #   "localSimilarity"=localSimilarity
  # ))
  return(globalSimilarity)  # TODO remove unneeded computation from this function
}
# end of edge similarity from COGENT


# normalization function to be applied to each row
nrmlz <- function(myVector){
  sqvect <- as.numeric(lapply(myVector, function(x) x*x)) #get the squares of every element in vector
  rootsumsq<- sqrt(sum(sqvect)) # sum them and square root to get normalization divisor (determinant?)
  normvect <- lapply(myVector, function(x) x/rootsumsq) # normalize vector by dividing elements by above
  return(normvect)
}

# replaces NaN with 0 in DF
is.nan.data.frame <- function(x){
  do.call(cbind, lapply(x, is.nan))
}

# creates distance correlation matrix for given dataframe
# current implementation took 1h15min to execute
getDistCorrMatrix <- function(myDF){
  dcormat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  pmat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  itercombs <- expand.grid(range(1, ncol(myDF)), range(1, ncol(myDF)))
  
  # parallel construction of correlation matrix (broken)
  dcormat <-
    foreach(i = 1:ncol(myDF), .combine='cbind') %:%
    foreach(j = 1:ncol(myDF), .combine='c') %dopar% {
      dcor(myDF[,i], myDF[,j])$dcor
    }
  
  pmat <-
    foreach(i = 1:ncol(myDF), .combine='cbind') %:%
    foreach(j = 1:ncol(myDF), .combine='c') %dopar% {
      as.numeric(dcor.ttest(as.matrix(myDF[,i]), as.matrix(myDF[,j]))["p-value"])
    }
  pmat[is.nan.data.frame(pmat)] <- 0  # bc for some reason this keeps putting NaNs on diagonal
  
  #  # in case using foreach dopar sucks 
  #   pb <- txtProgressBar(min = 0, max = ncol(myDF), style = 3)
  #   for(i in 1:ncol(myDF)){
  #     for(j in 1:ncol(myDF)){
  #       thisdcor <- dcor(myDF[,i], myDF[,j])
  #       dcormat[i,j] <- thisdcor$dcor
  #     }
  #     setTxtProgressBar(pb, i)
  #   }
  #   
  #   for(i in 1:ncol(myDF)){
  #     for(j in 1:ncol(myDF)){
  #       thispmat <- dcor.ttest(myDF[,i], myDF[,j])
  #       pmat[i,j] <- thispmat["p-value"]
  #     }
  #     setTxtProgressBar(pb, i)
  #   }
  
  dcormat <- data.frame(dcormat)
  pmat <- data.frame(pmat)
  colnames(dcormat) <- colnames(myDF)
  rownames(dcormat) <- colnames(myDF)
  colnames(pmat) <- colnames(myDF)
  rownames(pmat) <- colnames(myDF)
  
  return(list(dcormat, pmat))
}

# creates distance correlation matrix for given dataframe
getPearsonCorrMatrix <- function(myDF){
  spcormat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  pmat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  pb <- txtProgressBar(min = 0, max = ncol(myDF), style = 3)
  
  for(i in 1:ncol(myDF)){
    for(j in 1:ncol(myDF)){
      spcormat[i,j] <- abs(cor(myDF[,i], myDF[,j], method = "pearson"))
    }
    setTxtProgressBar(pb, i)
  }
  
  for(i in 1:ncol(myDF)){
    for(j in 1:ncol(myDF)){
      pmat[i,j] <- cor.test(myDF[,i], myDF[,j], method = "pearson")[["p.value"]]
    }
    setTxtProgressBar(pb, i)
  }
  
  spcormat <- data.frame(spcormat)
  pmat <- data.frame(pmat)
  colnames(spcormat) <- colnames(myDF)
  rownames(spcormat) <- colnames(myDF)
  colnames(pmat) <- colnames(myDF)
  rownames(pmat) <- colnames(myDF)
  
  return(list(spcormat, pmat))
}

# creates distance correlation matrix for given dataframe
getSpearmanCorrMatrix <- function(myDF){
  spcormat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  pmat <- matrix(, ncol = ncol(myDF), nrow = ncol(myDF))
  pb <- txtProgressBar(min = 0, max = ncol(myDF), style = 3)
  
  for(i in 1:ncol(myDF)){
    for(j in 1:ncol(myDF)){
      spcormat[i,j] <- abs(cor(myDF[,i], myDF[,j], method = "spearman"))
    }
    setTxtProgressBar(pb, i)
  }
  
  for(i in 1:ncol(myDF)){
    for(j in 1:ncol(myDF)){
      pmat[i,j] <- cor.test(myDF[,i], myDF[,j], method = "spearman")[["p.value"]]
    }
    setTxtProgressBar(pb, i)
  }
  
  spcormat <- data.frame(spcormat)
  pmat <- data.frame(pmat)
  colnames(spcormat) <- colnames(myDF)
  rownames(spcormat) <- colnames(myDF)
  colnames(pmat) <- colnames(myDF)
  rownames(pmat) <- colnames(myDF)
  
  return(list(spcormat, pmat))
}

# constructs adjacency matrix from csv file with rows=samples, cols=features
makeCorrGraph <- function(filename, corr="distance", split=T, set="ALL", transpose=F){
  ogdat <- read.csv(filename, header=T) # read in network
  
  # take only rows from specified set
  if(set != "ALL"){
    ogdat <- ogdat[ogdat$DX_bl==set,]
  }
  
  if(!is.numeric(typeof(as.vector(ogdat[,1])))){  # if first col is not numeric, remove it
    while(typeof(as.vector(ogdat[,1])) != "double"){
      ogdat <- ogdat[,-1]
    }
  }
  
  normdat <- ogdat  # normalize data
  for(i in 1:nrow(ogdat)){
    normrow <- as.numeric(as.vector(ogdat[i,]))
    normrow <- nrmlz(normrow)
    normdat[i,] <- normrow
  }
  
  # if transpose is true, transpose the matrix
  if(transpose == T){
    normdat <- as.data.frame(t(as.matrix(normdat)))
  }
  
  # if split is TRUE, split normalized matrix into two equal matrices before
  if(split==T){
    if(nrow(normdat) %% 2 != 0) normdat <- normdat[-1,] # if odd number of rows, cut one so can split evenly
    splitindex <- sample(nrow(normdat), nrow(normdat)/2) # randomly pick half of the rows
    dat1 <- normdat[splitindex,]
    dat2 <- normdat[-splitindex,]
    
    if(corr == "distance") {
      mat1s <- getDistCorrMatrix(dat1)
      cormat1 <- mat1s[[1]]
      pmat1 <- mat1s[[2]]
    }
    else if(corr == "pearson") {
      mat1s <- getPearsonCorrMatrix(dat1)
      cormat1 <- mat1s[[1]]
      pmat1 <- mat1s[[2]]
    }
    else {
      mat1s <- getSpearmanCorrMatrix(dat1)
      cormat1 <- mat1s[[1]]
      pmat1 <- mat1s[[2]]
    }
    
    if(corr == "distance") {
      mat2s <- getDistCorrMatrix(dat2)
      cormat2 <- mat2s[[1]]
      pmat2 <- mat2s[[2]]
    }
    else if(corr == "pearson") {
      mat2s <- getPearsonCorrMatrix(dat2)
      cormat2 <- mat2s[[1]]
      pmat2 <- mat2s[[2]]
    }
    else {
      mat2s <- getSpearmanCorrMatrix(dat2)
      cormat2 <- mat2s[[1]]
      pmat2 <- mat2s[[2]]
    }
  } 
  
  else{
    if(corr == "distance") {
      mats <- getDistCorrMatrix(normdat)
      cormat <- mats[[1]]
      pmat <- mats[[2]]
    }
    else if(corr == "pearson") {
      mats <- getPearsonCorrMatrix(normdat)
      cormat <- mats[[1]]
      pmat <- mats[[2]]
    }
    else {
      mats <- getSpearmanCorrMatrix(normdat)
      cormat <- mats[[1]]
      pmat <- mats[[2]]
    }
  }
  
  if(split==T) return(list(data.frame(cormat1), data.frame(cormat2), data.frame(pmat1), data.frame(pmat2)))
  else return(list(data.frame(cormat), data.frame(pmat)))
}

# accepts an adjacency matrix as created by makeCorrGraph and trims all edges below threshold
# edits in place but if that becomes a problem I'll change so it edits a copy
trimEdgesByThreshold <- function(adjmatrix, threshold){
  for(i in 1:nrow(adjmatrix)){
    for(j in 1:ncol(adjmatrix)){
      if(as.double(adjmatrix[i,j]) < as.double(threshold)){
        adjmatrix[i,j] <- 0
      }
    }
  }
  return(adjmatrix)
}

trimEdgesByPval <- function(adjmatrix, pmatrix, threshold=0.05){
  for(i in 1:nrow(adjmatrix)){
    for(j in 1:ncol(adjmatrix)){
      if(as.double(pmatrix[i,j]) > as.double(threshold)){
        adjmatrix[i,j] <- 0
      }
    }
  }
  return(adjmatrix)
}

# gets pagerank for all vertices and saves a histogram of the values
doPageRank <- function(edgelist, plotname="PRplot.png"){
  g <- graph_from_data_frame(edgelist, directed = FALSE, vertices = NULL)
  pr <- page_rank(g)
  prframe <- as.data.frame(pr['vector'])
  
  # uncomment if you want individual PR plots
  # bins <- seq(from=0, to=1, by=0.05)
  # gg <- ggplot(prframe, aes(x=pr['vector'])) + geom_bar(stat="bin")
  # ggsave(plotname, plot=gg)
  # create and save histogram of pagerank values
  
  return(prframe)
}

# TODO FIX THIS, PR1 AND PR2 CAN BE OF DIFFERING LENGTHS SO YOU CAN'T DO PEARSON ON THEM. USE UNION?
# FOR NOW RETURN NA
# ACTUALLY NOW IT DOES DISTCORR WHICH DOESN'T NEED EQUAL CLASS SIZES
prPearson <- function(pr1, pr2){
  # pcor <- cor(pr1, pr2, method = "pearson")
  # pval <- cor.test(pr1, pr2, method = "pearson")[["p.value"]]
  if(isempty(pr1['vector'])|isempty(pr2['vector'])){
    pcor <- NA
    pval <- NA
  } else {
    pcor <- dcor(as.matrix(pr1), as.matrix(pr2))$dcor
    pval <- as.numeric(dcor.ttest(as.matrix(pr1), as.matrix(pr2))["p-value"])
  }
  return(c(pcor, pval))
}