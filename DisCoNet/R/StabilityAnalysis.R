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

stabilityAnalysis <- function(filename, seqmin=0, seqmax=1, seqstep=0.1, method="distance", pthresh=0.05, boots=25){
myfile = filename
cormethods = method
thresholds = seq(from=as.double(seqmin), to=as.double(seqmax), by=as.double(seqstep))

resultsdf <- data.frame(Method="dummy", By="byDummy", Set="dummy", Threshold=0.99, Deltacon=0.99, Jaccard=0.99, PRpearsonM1=0.99, PRpearsonM2=0.99, PRpearsonM1pval=0.99, PRpearsonM2pval=0.99,
                        PRksM1=0.99, PRksM2=0.99, PRksM1pval=0.99, PRksM2pval=0.99, stringsAsFactors=F)
colnames(resultsdf) <- c("Method", "Boot", "Set" , "Threshold", "Deltacon", "Jaccard")
set="ALL"

#TODO: more for loops! we want one for [bySample and byFeature] and one for [all, AD, CN, LMCI], and maybe one for bootstrapping?
resultsdf <- foreach(boot=1:boots, .combine=rbind) %:%
  foreach(cor=cormethods, .combine=rbind) %:%
    foreach(threshold=thresholds, .combine=rbind) %dopar% {
      # calculate adj matrix for both sets
      matrices <- makeCorrGraph(myfile, corr=cor, split=T, set=set)
      m1 <- matrices[[1]]
      m2 <- matrices[[2]]
      p1 <- matrices[[3]]
      p2 <- matrices[[4]]

        # at each threshold, remove edges with weights below threshold
        trimmedmat1 <- trimEdgesByThreshold(m1, threshold)
        trimmedmat2 <- trimEdgesByThreshold(m2, threshold)
        trimmedmat1 <- trimEdgesByPval(trimmedmat1, p1, threshold=pthresh)
        trimmedmat2 <- trimEdgesByPval(trimmedmat2, p2, threshold=pthresh)
        
        # if trimming makes the matrix empty, replace with identity
        if(all(sapply(trimmedmat1, function(x) x==0))){
          trimmedmat1 <- as.data.frame(diag(ncol(m1)))
          names(trimmedmat1) <- names(m1)
          rownames(trimmedmat1) <- rownames(m1)
        }
        
        if(all(sapply(trimmedmat2, function(x) x==0))){
          trimmedmat2 <- as.data.frame(diag(ncol(m2)))
          names(trimmedmat2) <- names(m2)
          rownames(trimmedmat2) <- rownames(m2)
        }
        
        # calculate scores for both matrices at threshold
        m1edges <- data.frame(which(trimmedmat1!=0, arr.ind=T)) #convert to edgelist for DConn and PageRank
        colnames(m1edges) <- c("src", "dst")
        m2edges <- data.frame(which(trimmedmat2!=0, arr.ind=T)) 
        colnames(m2edges) <- c("src", "dst")
        
        # check that edgelists aren't empty before doing dcon
        if(empty(m1edges)|empty(m2edges)){
          dcon <- NA
          jaccard <- NA
        } else {
          dcon <- delta_con(m1edges, m2edges, ncol(trimmedmat1))
          jaccard <- getEdgeSimilarity(list(trimmedmat1, trimmedmat2))
        }
        
        # output both matrices and score table to file or print
        basename <- paste(cor, as.character(threshold), as.character(set), "_", as.character(boot), sep="")
        
        write(paste(as.character(boot), set, threshold, "\n"), file="StabilityBatchProgress.txt", append=T)
        
        return(c(cor, boot, set, threshold, dcon, jaccard))
      }
colnames(resultsdf) <- c("Method", "Boot", "Set", "Threshold", "Deltacon", "Jaccard")
fname <- "NRCthresholdResultsCCLE25_2.csv"
write.csv(resultsdf, fname) # writes to file with dummy row removed

### do plotting
resultsfile <- fname
df <- read.csv(resultsfile)

### make deltacon figure ###
for(cormethod in c("distance", "pearson", "spearman")){
  myTitle <- paste("Mean DeltaCon between", str_to_title(cormethod), " Metabolomics Correlation Networks at Various Edge Thresholds")
  myXLab <- paste(str_to_title(cormethod), "Correlation Threshold" )

  newdf <- df %>%
    filter(Method == cormethod) %>%
    group_by(Threshold, Set) %>%
    summarize(mean = mean(Deltacon), sd = sd(Deltacon))

  beedf <- filter(df, Method == cormethod)

  # TODO fix bees
  gg1 <- ggplot(data = newdf) + geom_point(aes(x=Threshold, y=mean)) +
    geom_errorbar(aes(x=Threshold, y=mean, ymin=mean-sd, ymax=mean+sd)) +
    scale_x_continuous(name=myXLab, breaks=seq(0,1,0.1)) +
    scale_y_continuous(name="DeltaCon", breaks=seq(0,1,0.1), limits=c(0,1.1)) +
    labs(title=myTitle) +
    geom_abline(intercept=0, slope=0) +
    geom_quasirandom(data=beedf, aes(x=Threshold, y=Deltacon, shape=".", alpha=0.01)) +
    facet_wrap(~ Set, ncol=4) +
    theme_classic()
  plotname <- paste("metabolomicsDeltacon", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(gg1)
  dev.off()
}


### make jaccard figure ###
for(cormethod in c("distance", "pearson", "spearman")){
  myTitle <- paste("Mean Jaccard Index between", str_to_title(cormethod), "Metabolomics Correlation Networks at Various Edge Thresholds")
  myXLab <- paste(str_to_title(cormethod), "Correlation Threshold" )

  newdf <- df %>%
    filter(Method == cormethod) %>%
    group_by(Set, Threshold) %>%
    summarize(mean = mean(Jaccard), sd = sd(Jaccard))

  beedf <- filter(df, Method == cormethod)

  # TODO fix bees
  gg1 <- ggplot(data = newdf) + geom_point(aes(x=Threshold, y=mean)) +
    geom_errorbar(aes(x=Threshold, y=mean, ymin=mean-sd, ymax=mean+sd)) +
    scale_x_continuous(name=myXLab, breaks=seq(0,1,0.1)) +
    scale_y_continuous(name="Jaccard Index", breaks=seq(0,1,0.1), limits=c(0,1.1)) +
    labs(title=myTitle) +
    geom_quasirandom(data=beedf, aes(x=Threshold, y=Jaccard, shape=".", alpha=0.01)) +
    geom_abline(intercept=0, slope=0) +
    facet_wrap(~ Set, ncol=4) +
    theme_classic()
  plotname <- paste("metabolomicsJaccard", cor, ".png", sep="")
  png(plotname, width = 1600, height = 900)
  plot(gg1)
  dev.off()
}
}





