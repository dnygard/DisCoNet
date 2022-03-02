library(ggplot2)
library(wrapr)
library(dplyr)
library(stringr)
library(tidyr)

# load functions
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

plotRateOfChange <- function(filename, seqmin=0, seqmax=1, seqstep=0.1, method="distance", pthresh=0.05){

# init vars from arguments
myfile <- filename 
thresholds <- seq(from=as.double(seqmin), to=as.double(seqmax), by=as.double(seqstep))
cormethods <- method

# init vars that don't change
resultsdf <- data.frame(Method="dummy", Boot=-1, Set="dummy", Threshold=0.99, Deltacon=0.99, Jaccard=0.99, PRdistcorrM1=0.99, PRdistcorrM2=0.99, PRdistcorrM1pval=0.99, PRdistcorrM2pval=0.99,
                        PRksM1=0.99, PRksM2=0.99, PRksM1pval=0.99, PRksM2pval=0.99, stringsAsFactors=F)
colnames(resultsdf) <- c("Method", "Boot", "Set" , "Threshold", "Deltacon", "Jaccard", "PRdistcorrM1", "PRdistcorrM2", "PRdistcorrM1pval", "PRdistcorrM2pval", "PRksM1", "PRksM2", "PRksM1pval", "PRksM2pval")


for(cor in cormethods){
  # make network
  matrix <- makeCorrGraph(myfile, corr=cor, split=F, set="ALL")
  m <- matrix[[1]] 
  p <- matrix[[2]]
  
  # create nameslist
  for(i in seq(seqmin, seqmax, seqstep)){
    if(firststep){
      stepname = paste(str(firststep), "-", str(i))
      firststep = i
      
      if(nameslist) c(nameslist, stepname)
      else nameslist = stepname
    }
    else firststep = i
  }
  # nameslist <- c("0.0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5", "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1.0")
  
  # create df of edges, pr, transitivity, and degree stats for each threshold
  atthreshdf <- data.frame(Threshold=double(), numedges=integer(), meandegree=double(), mediandegree=double(), meanpr=double(), medianpr=double(), transitivity=double)
  names(atthreshdf) <- c("Threshold", "numedges", "meandegree", "mediandegree", "meanpr", "medianpr", "transitivity")
  for(threshold in seq(seqmin, seqmax, seqstep)){
    # trim matrix
    trimmedmat <- trimEdgesByThreshold(m, threshold)
    trimmedmat <- trimEdgesByPval(trimmedmat, p, threshold=pthresh)
    # if trimming makes the matrix empty, replace with identity
    if(all(sapply(trimmedmat, function(x) x==0))){
      trimmedmat <- as.data.frame(diag(ncol(m)))
      names(trimmedmat) <- names(m)
      rownames(trimmedmat) <- rownames(m)
    }
    # calculate scores for both matrices at threshold
    medges <- data.frame(which(trimmedmat!=0, arr.ind=T)) #convert to edgelist for DConn and PageRank
    colnames(medges) <- c("src", "dst")
    
    # get mean (median?) PR
    pr <- doPageRank(medges)
    meanpr <- mean(as.numeric(as.character(unlist(pr["vector"]))))
    medianpr <- median(as.numeric(as.character(unlist(pr["vector"]))))
    
    # get mean degree(median?)
    g <- graph_from_data_frame(medges, directed = FALSE, vertices = NULL)  # convert to igraph graph
    degreedistrn <- as.vector(degree(g))
    meandegree <- mean(degreedistrn)
    mediandegree <- median(degreedistrn)
    
    # get number of edges
    numedges <- nrow(medges)
    # get transitivity
    thistransitivity <- transitivity(g, type="undirected")
    
    # get PR vector
    pr <- doPageRank(medges)
    prlist[[length(prlist)+1]] <- pr['vector']
    
    # get degree vector
    g <- graph_from_data_frame(medges, directed = FALSE, vertices = NULL)  # convert to igraph graph
    degreedistrn <- as.vector(degree_distribution(g))
    deglist[[length(deglist)+1]] <- degreedistrn
    
    atthreshdf <- rbind(atthreshdf, c(threshold, numedges, meandegree, mediandegree, meanpr, medianpr, thistransitivity))
  }
  
  names(atthreshdf) <- c("Threshold", "numedges", "meandegree", "mediandegree", "meanpr", "medianpr", "transitivity")
  
  # make threshold comparison df
  comparedf <- data.frame(Threshold=character(), numedges=double(), meandegree=double(), mediandegree=double(), meanpr=double(), medianpr=double(), transitivity=double(), stringsAsFactors = F)
  names(comparedf) <- c("Threshold", "numedges", "meandegree", "mediandegree", "meanpr", "medianpr", "transitivity")
  for(i in 1:10){
    casename <- paste(as.character(atthreshdf$Threshold[i]), "-", as.character(atthreshdf$Threshold[i+1]))
    # get edge number ratio
    enumratio <- atthreshdf$numedges[i]/atthreshdf$numedges[i+1]
    
    # get mean degree ratio
    meandegratio <- atthreshdf$meandegree[i]/atthreshdf$meandegree[i+1]
    mediandegratio <- atthreshdf$mediandegree[i]/atthreshdf$mediandegree[i+1]
    
    # get mean page rank ratio
    meanprratio <- atthreshdf$meanpr[i]/atthreshdf$meanpr[i+1]
    medianprratio <- atthreshdf$medianpr[i]/atthreshdf$medianpr[i+1]
    
    # get transitivity ratio
    transitivityratio <- atthreshdf$transitivity[i]/atthreshdf$transitivity[i+1]
    
    
    print(c(casename, enumratio, meandegratio, mediandegratio, meanprratio, medianprratio, transitivityratio))
    comparedf <- rbind(comparedf, data.frame(casename, enumratio, meandegratio, mediandegratio, meanprratio, medianprratio, transitivityratio, stringsAsFactors = F))
  }
  names(comparedf) <- c("Threshold", "numedges", "meandegree", "mediandegree", "meanpr", "medianpr", "transitivity")
  
### plot KS figs
  # initialize empty dfs for eventually plugging into ggplot
  degfigdf <- data.frame(name=character(), degree=integer(), freq1 = double(), freq2 = double())
  prfigdf <- data.frame(name=character(), bin=character(), count1 = integer(), count2 = integer())
  
  # get KS and make df for dists
  for(i in 1:length(nameslist)){
    tempdegdf <- data.frame(name=character(), degree=integer(), freq1 = double(), freq2 = double())
    tempprdf <- data.frame(name=character(), bin=double(), count1 = integer(), count2 = integer())
    prdistrn1 <- prlist[[i]]
    prdistrn2 <- prlist[[i+1]]
    degdistrn1 <- deglist[[i]]
    degdistrn2 <- deglist[[i+1]]
    
    # check if vectors identical for ks
    if(identical(as.vector(prdistrn1), as.vector(prdistrn2))) prks <- NA
    else prks <- ks.test(as.matrix(prdistrn1), as.matrix(prdistrn2))
    if(identical(as.vector(degdistrn1), as.vector(degdistrn2))) degks <- NA
    else degks <- ks.test(as.matrix(degdistrn1), as.matrix(degdistrn2))
    
    if(!is.na(prks$statistic)){  # most confusing way to write an equality lol (is not NA = TRUE)
      prname <- paste(nameslist[[i]], "KS:", as.character(signif(prks$statistic, 2)), "(p-value: ", as.character(signif(prks$p.value, 2)), ")")
    } else {
      prname <- paste(nameslist[[i]], "KS: NA (PR vectors identical)")
    }
    if(!is.na(degks$statistic)){  # it's not NOT opposite day
      degname <- paste(nameslist[[i]], "KS:", as.character(signif(degks$statistic, 2)), "(p-value: ", as.character(signif(degks$p.value, 2)), ")")
    } else {
      degname <- paste(nameslist[[i]], "KS: NA (Degree vectors identical)")
    }
    
    # degree distribution df
    for(j in 1:max(length(degdistrn1), length(degdistrn2))){
      if(j<=length(degdistrn1)){deg1 <- degdistrn1[[j]]}
      else{deg1 <- NA}
      if(j<=length(degdistrn2)){deg2 <- degdistrn2[[j]]}
      else{deg2 <- NA}
      tempdegdf <- rbind(tempdegdf, list(name=degname, degree=j, freq1=deg1, freq2=deg2))
    }
    
    # pr distribution df
    prdistrn1 <- as.numeric(as.character(unlist(prdistrn1)))
    prdistrn2 <- as.numeric(as.character(unlist(prdistrn2)))
    prhist1 <- hist(prdistrn1, breaks=seq(0,1,0.0001))
    prhist2 <- hist(prdistrn2, breaks=seq(0,1,0.0001))
    
    # fill df with distribution and relevant data in long format
    for(j in 1:length(prhist1$breaks)){
      tempprdf <- rbind(tempprdf, list(name=prname, bin = prhist1$breaks[j], count1=prhist1$count[j], count2=prhist2$count[j]))
    }
    
    # get bind rows to big df
    degfigdf <- rbind(degfigdf, tempdegdf)
    prfigdf <- rbind(prfigdf, tempprdf)
  }
  
  # plot deg KS fig
  ggtitle <- paste("KS of degree distribution for", cor, "networks with threshold difference of ", as.character(seqstep)) 
  ggdeg <- ggplot(data = degfigdf) + 
    labs(title=ggtitle) +
    geom_bar(aes(x=degree, y=freq1, fill="blue", alpha=0.2), stat="identity") +
    geom_bar(aes(x=degree, y=freq2, fill="red", alpha=0.2), stat="identity") +
    scale_y_continuous(name="Frequency", breaks=seq(0,0.05,0.01), limits=c(0,0.05)) +
    facet_wrap(~ name, ncol=5, nrow=2) +
    theme_classic()
  
  plotname <- paste("KSvDegree", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot PR KS fig
  ggtitle <- paste("KS of PR distribution for", cor, "networks with threshold difference of 0.1 ", as.character(seqstep)) 
  ggpr <- ggplot(data = prfigdf) + 
    labs(title=ggtitle) +
    geom_area(aes(x=bin, y=count1, fill="former \ndistribution", alpha=0.2), stat="identity") +
    geom_area(aes(x=bin, y=count2, fill="latter \ndistribution", alpha=0.2), stat="identity") +
    scale_x_continuous(name="PR", breaks=seq(0,0.005,0.0025), limits=c(0,0.005)) +
    facet_wrap(~ name, ncol=5, nrow=2) +
    theme_classic()
  
  plotname <- paste("KSvPR", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggpr)
  dev.off()
}
  

  ### plot rate of change figs
  # plot edge number rate of change fig
  ggtitle <- paste("Ratio of edge number between", cor, "metabolomics networks with threshold difference of ", as.character(seqstep)) 
  ggdeg <- ggplot(data = comparedf) + 
    labs(title=ggtitle, x="Thresholds being compared", y="Edge Number at 1st Threshold / Edge Number at 2nd Threshold") +
    geom_line(aes(x=Threshold, y=numedges, group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsEdgenumRateofChange", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot deg Rate of Change fig
  ggtitle <- paste("Ratio of average degree between", cor, "metabolomics networks with threshold difference of ", as.character(seqstep)) 
  ggdeg <- ggplot(data = comparedf) + 
    labs(title=ggtitle, x="Thresholds being compared", y="Avg. Degree at 1st Threshold / Avg. Degree at 2nd Threshold") +
    geom_line(aes(x=Threshold, y=meandegree, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=mediandegree, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsDegreeRateofChange", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot pr Rate of Change fig
  ggtitle <- paste("Ratio of average PageRank between", cor, "metabolomics networks with threshold difference ", as.character(seqstep)) 
  ggdeg <- ggplot(data = comparedf) + 
    labs(title=ggtitle, x="Thresholds being compared", y="Avg. PageRank at 1st Threshold / Avg. PageRank at 2nd Threshold") +
    geom_line(aes(x=Threshold, y=meanpr, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=medianpr, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsPageRankRateofChange", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot transitivity rate of change fig
  ggtitle <- paste("Ratio of transitivity between", cor, "metabolomics networks with threshold difference of ", as.character(seqstep)) 
  ggdeg <- ggplot(data = comparedf) + 
    labs(title=ggtitle, x="Thresholds being compared", y="Transitivity at 1st Threshold / Transitivity at 2nd Threshold") +
    geom_line(aes(x=Threshold, y=meanpr, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=medianpr, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsTransitivityRateofChange", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  
  ### plot non rate of change figs
  # plot edge number rate of change fig
  ggtitle <- paste("Edge number in", cor, "metabolomics networks at various thresholds") 
  ggdeg <- ggplot(data = atthreshdf) + 
    labs(title=ggtitle, x="Correlation Threshold", y="Number of Edges") +
    geom_line(aes(x=Threshold, y=numedges, group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsEdgenum", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot deg Rate of Change fig
  ggtitle <- paste("Average degree of", cor, "metabolomics networks at various thresholds") 
  ggdeg <- ggplot(data = atthreshdf) + 
    labs(title=ggtitle, x="Correlation Threshold", y="Avg. Degree") +
    geom_line(aes(x=Threshold, y=meandegree, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=mediandegree, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsDegree", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot pr Rate of Change fig
  ggtitle <- paste("Average PageRank of", cor, "metabolomics networks at various thresholds") 
  ggdeg <- ggplot(data = atthreshdf) + 
    labs(title=ggtitle, x="Correlation Threshold", y="Avg. PageRank") +
    geom_line(aes(x=Threshold, y=meanpr, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=medianpr, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsPageRank", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
  # plot transitivity rate of change fig
  ggtitle <- paste("Transitivity of", cor, "metabolomics networks at various thresholds") 
  ggdeg <- ggplot(data = comparedf) + 
    labs(title=ggtitle, x="Correlation Threshold", y="Transitivity") +
    geom_line(aes(x=Threshold, y=meanpr, colour="Mean", group=1)) +
    geom_line(aes(x=Threshold, y=medianpr, colour="Median", group=1)) +
    theme_classic()
  plotname <- paste("metabolomicsTransitivity", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(ggdeg)
  dev.off()
  
}


