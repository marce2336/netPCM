 #---------For plotting gene expression, metabolite concentration, etc...-------------#


                                                                                         ### You may even skip everything before...
normToMedian <- function (X){                                                            ### this is the part that will make the heatmap
   row_medians = apply( X, 1, median, na.rm = T)
   names(row_medians) = NULL
   X1 = log2(X/row_medians)
   return(X1)
}
                                               ### want to visualize (DATA, NORM_DATA, MED_DATA)
plotHeatMap <- function (X){

    X[X==0] <- NA 
    
    ind.r <- apply(X, 1, var) == 0
    X <- X[!ind.r,]

    hc.c <- hclust(dist(t(X)), "single")

    

    X = normToMedian(X)

    library(gplots)
    breaks <- seq(-2,2, length = 61)

    X.matrix = as.matrix(X)   

    heatmap.2(X.matrix, col = redblue(60) , trace = "none", scale = "none", breaks = breaks, margins = c(10, 10),
     Colv=as.dendrogram(hc.c), Rowv='none', dendrogram = 'column', cexRow = 0.5, cexCol = 0.8
     #Colv=FALSE, Rowv =FALSE, dendrogram = 'none', cexRow = 0.8, cexCol = 0.8,
     #Colv=FALSE, Rowv=as.dendrogram(hc.r), dendrogram = 'row', cexRow = 0.7, cexCol = 0.7,
    )
          
}

plotHeatMap2 <- function (X , cellnote){

    #hc.r <- hclust(dist(X), "average")
    #hc.c <- hclust(dist(t(X)), "average")

    #X = normToMedian(X)

    library(gplots)
    breaks <- seq(-2,2, length = 61)

    X.matrix = as.matrix(X)
    #X1 = round(DATA1, digits = 2)


    
    heatmap.2(X.matrix, col = redblue(60) , trace = "none", scale = "none", breaks = breaks,
      margins = c(10, 10),
      cellnote = cellnote, notecex=0.8,
    #colsep = c(6, 12, 17, 22, 27, 32, 37, 43, 48, 53, 59, 64, 70, 75),
    #sepcolor = "black", sepwidth = 0.1,
    #colsep = c(14,33,41,47,61,71,82,92,97), sepcolor = c("black","black"), sepwidth = c(0.15,0.15),
     #rowsep = c(14,33,41,47,61,71,82,92,97), #sepcolor = "black", sepwidth = 0.1,
     #labRow = DATA0[,1],
     #labCol = labCol,
     #key = FALSE,

   
    #       notecol="black",
    #       na.color=par("bg"),

     Colv=FALSE, Rowv =FALSE, dendrogram = 'none', cexRow = 0.5, cexCol = 0.5,)
     #Colv=FALSE, Rowv=as.dendrogram(hc.r), dendrogram = 'row', cexRow = 0.7, cexCol = 0.7,)
     #Colv=as.dendrogram(hc.c), Rowv='none', dendrogram = 'column', cexRow = 0.5, cexCol = 0.5,)

}

plotHeatMap3 <- function (X, colv){

   # X[X==0] <- NA 
    
    ind.r <- apply(X, 1, var) == 0
    X <- X[!ind.r,]

    hc.c <- hclust(dist(t(X)), "single")

    

    X = normToMedian(X)

    library(gplots)
    breaks <- seq(-2,2, length = 61)

    X.matrix = as.matrix(X)   

    heatmap.2(X.matrix, col = redblue(60) , trace = "none", scale = "none", breaks = breaks, margins = c(10, 10),
     Colv= colv, Rowv='none', dendrogram = 'column', cexRow = 0.4, cexCol = 0.4
     #Colv=FALSE, Rowv =FALSE, dendrogram = 'none', cexRow = 0.5, cexCol = 0.5,)
     #Colv=FALSE, Rowv= as.dendrogram(hc.r) , dendrogram = 'row', cexRow = 0.7, cexCol = 0.7,
    )

}



#---------------------------Normalize your data to the internal standard ----------------#

#normToStd <- function (X){

#    m = median(X["PC 34:0" ,], na.rm = T)
#    norm.factor = X["PC 34:0",]/m
#    norm.factor = as.numeric(norm.factor)
#    names(norm.factor) = colnames(X) 
                                                          
#    NORM_DATA = t(apply(X, 1, function(x) x/norm.factor)) 
    
#    return(NORM_DATA)
#}


#----------------Normalize to ion count excluding compounds with high variance--------------#

normToTIC2 <- function (X){

    tags = grep("TAG",row.names(X))
    no_tags = setdiff(1:nrow(X), tags)
    DATA_NO_TAGS = X[no_tags,]


#    sdev = apply(X, 1, sd, na.rm=T)
#    sorted.sdev = sort(sdev, decreasing = T)
#    changing.more = sorted.sdev[1:78]
#    changing.less = setdiff(row.names(X), names(changing.more))
#    DATA_CHANGING_LESS = X[changing.less,]
#    sum4 = apply(DATA_CHANGING_LESS, 2, sum, na.rm=T)
#    NORM_4 = t(apply(X, 1, function(x) x/sum4 * 10000))
#

    row_medians = apply( DATA_NO_TAGS, 1, median, na.rm = T)
    sdev = apply(DATA_NO_TAGS, 1, sd, na.rm=T)
    c.of.v = sdev/row_medians
    sorted.c.of.v = sort(c.of.v, decreasing = T)
    changing.more = sorted.c.of.v[1:round(nrow(DATA_NO_TAGS)/2)]
    changing.less = setdiff(row.names(DATA_NO_TAGS), names(changing.more))
    DATA_CHANGING_LESS = DATA_NO_TAGS[changing.less,]
    sum4 = apply(DATA_CHANGING_LESS, 2, sum, na.rm=T)
    NORM_4 = t(apply(X, 1, function(x) x/sum4 * 10000))
    return(NORM_4)

}

#-------------------------------- Average replicates-----------------------------------------#

getMedData <- function (X){                                                                    #### This section is useful in case you don't want
                                                                                               #### to visualize each replicate separate and you
                                                                                               #### prefer to visualize the median instead
    reps =   colnames(X)
    
#    reps1 = sub("_Pi_", "_ctrlPi_", reps, perl =T )
#    reps1 = sub("_Pi-_", "_lowPi_", reps1, perl =T )
#    reps1 = sub("_Pi\\+_", "_highPi_", reps1, perl =T)

    reps1 = gsub("\\db?_", "", reps, perl =T )                                                #### This depends on the name format of your replicates

    tpts = unique(reps1)                                           

    MED_X = X[,1:length(tpts)]
    colnames(MED_X) = tpts

    for(i in 1:nrow(X)){

        for(j in 1:length(tpts)){
            indexes = grep( tpts[j], reps1, perl =T)
            MED_X[i,j] = median(X[i,indexes], na.rm = T)

        }

    }

  return(MED_X)
}



#-----------------------------------------------------------------------------------------#

plotCors <- function (X){

    #hc.r <- hclust(dist(X), "complete")
    #hc.c <- hclust(dist(t(X)), "complete")
    library(gplots)
    breaks <- seq(-1,1, length = 61)
    X.matrix = as.matrix(X)
    heatmap.2(X.matrix, col = redblue(60) , trace = "none", scale = "none", breaks = breaks,
     margins = c(5, 5),
    #heatmap.2(X.matrix, col = colorpanel(60, "deepskyblue3","gold","firebrick1") , trace = "none", scale = "none", breaks = breaks,
    #colsep = c(6, 12, 17, 23, 28, 33, 38, 44, 50, 55, 61, 66, 72, 77, 83),# sepcolor = c("black","black"), sepwidth = c(0.15,0.15),
     #rowsep = c(14,33,41,47,61,71,82,92,97), 
     #sepcolor = "black", sepwidth = 0.1,
     #labRow = labRow,
     #labCol = labCol,
     
     #### Para mantener el orden existente
     
     Colv=FALSE, Rowv =FALSE, dendrogram = 'none', cexRow = 0.5, cexCol = 0.5,)
     
     #Colv=FALSE, Rowv=as.dendrogram(hc.r), dendrogram = 'row', cexRow = 0.4, cexCol = 0.4,)
     
     ###Para ordenar segun dendrograma de las correlaciones .....................................................
     
     #Colv=as.dendrogram(hc.c), Rowv=FALSE, dendrogram = 'column', cexRow = 0.4, cexCol = 0.4,)

}


getDistanceMatrix<- function (X){
 # X = NORM_DATA3[,X]
  DMATRIX <- cor(X, method=c("spearman"), use = "pairwise.complete.obs")
  return(DMATRIX)
}



#---------------------------------------------------------------------------------------#


    
anovaAndTtest <- function (X){

    samp.c1 = colnames(X)
    samp.c2 = sub("_.+","", samp.c1, perl = T )
    samp.c3 = sub(".+_","", samp.c1, perl = T )

    samp = cbind(samp.c1, samp.c2, samp.c3)

    colnames(samp)[1] = ""
    colnames(samp)[2] = "sample"
    colnames(samp)[3] = "replicate"

    groups <- as.factor(samp[,2])         ##Creates the array "groups" which contains the column with repeated variables
    group.length <- length(levels(groups))    ##Records the number of variables that repeat

        # P-values
    anova.2 <- function(x,y) anova(lm(x ~ y))$Pr[1]
    p.values <- apply(X,1, anova.2, groups)

    #t.test
    ncol = group.length * (group.length-1) / 2
    p.t.test <- matrix(nrow = nrow(X), ncol = ncol )
    colnames(p.t.test) = 1:ncol

    t.test.2 <- function(x, y) {
        z <- try(t.test(x, y)$p.value, silent=TRUE)
        if(is.numeric(z))
            return(z)
        else
            return(NA)
    }

    c = 1
    for(i in 1:(group.length-1)) {

        for(j in (i+1):group.length){

            idx.i <- which(groups == unique(groups)[i])
            idx.j <- which(groups == unique(groups)[j])
            p.t.test[,c] <- apply(X, 1, function(x) { t.test.2(x[idx.i], x[idx.j]) })
            col.name = as.character(paste(unique(groups)[i],unique(groups)[j], sep = "_") )
            colnames(p.t.test)[c] = col.name 
            c = c+1

        }

    }

    return(data.frame(anova.p = p.values, p.t.test))

}  

plotBarPlots <- function (X, Classes){

#  NORM_DATA["PC 34:0",] = rep(NA, ncol(NORM_DATA))
  X[ X== Inf] = NA
  X = mergeIsomers(X)
 reps = colnames(X)
  lipid.names = row.names(X)
  #glycerolipids =  c("PC","PI","PE","PS","PG","SQDG","MGDG","DGDG")
  glycerolipids =  c("PC","PI","PE","PS","PG","SQDG","MGDG","DGDG","DAG")
  lysos =  c("Lyso-PC","Lyso-PE","Lyso-MGDG","Lyso-DGDG","MAG")
  pos.mode =  c("PC","PE","MGDG","DGDG") #  "DAG",
  neg.mode =  c("PG","PS","PI","SQDG") #  "DAG",
  
  c.34  = c("34:0","34:1","34:2","34:3","34:4","34:5","34:6")
      
  all.species = c("32:0","32:1","32:2","32:3","32:5","32:4","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "34:6","36:1","36:2","36:3","36:4", "36:5","36:6","36:7","38:0","38:1", "38:2","38:3", "38:4", "38:5", "38:6","40:0", 
        "40:1","40:2", "40:3","40:4","40:5","42:1","42:2","42:3","42:4")
        
  smaller.than.40 = c("32:0","32:1","32:2","32:3","32:5","32:4","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "34:6","36:1","36:2","36:3","36:4", "36:5","36:6","36:7","38:0","38:1", "38:2","38:3", "38:4", "38:5", "38:6")
         

        

  low.abundant =  c("32:0","32:1","32:2","32:3","32:5","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "36:1","38:2","38:3","38:4","38:5","38:6")
       
            
  tag = c(  "48:0","48:1","48:2","48:3","48:4","48:5","48:6","48:7","48:8","48:9",
                  "50:0","50:1","50:2","50:3","50:4","50:5","50:6","50:7","50:8","50:9",
                  "52:0","52:1","52:2","52:3","52:4","52:5","52:6","52:7","52:8","52:9",
                  "54:0","54:1","54:2","54:3","54:4","54:5","54:6","54:7","54:8","54:9",
                  "56:0","56:1","56:2","56:3","56:4","56:5","56:6","56:7","56:8","56:9",
                  "58:0","58:1","58:2","58:3","58:4","58:5","58:6","58:7","58:8","58:9",
                  "60:0","60:1","60:2","60:3","60:4","60:5","60:6","60:7","60:8","60:9"
              
  ) 
      
  Lysos = c("16:0","16:1","16:2","16:3","18:0","18:1","18:2","18:3","20:1","20:1","22:0","22:1")
  
  if(Classes == "glycerolipids"){
      classes =  glycerolipids
      names = all.species
      par(mfrow=c(3,3), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "pos.mode"){
      classes =  pos.mode
      names = smaller.than.40
      par(mfrow=c(2,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "neg.mode"){
      classes =  neg.mode
      names = all.species
      par(mfrow=c(2,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "lysos"){
      classes =  lysos
      names = Lysos
      par(mfrow=c(3,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "tag"){
      classes =  "TAG"
      names = tag

  }
  
  library(gplots)

  for( j in 1:length(classes)){
  
      #lipid.x.inds = grep (classes[j], lipid.names )
      lipid.x.inds = grep (paste("^",classes[j], sep = ""), lipid.names, perl = T )
      X_DATA = X[lipid.x.inds, , drop = FALSE]
      x = rep(NA, length(lipid.x.inds))
      X_DATA = cbind(X_DATA, x)
      row.names(X_DATA) = row.names(X)[lipid.x.inds]

      #means = c()
      means = rep(NA, length(names))
      ci.l = rep(NA, length(names))
      ci.u = rep(NA, length(names))

      names(means) = names
      names(ci.l) = names
      names(ci.u) = names

      row.names = row.names(X_DATA)
      row.names = sub(" $", "", row.names, perl = T )
      row.names = sub(".+ ", "", row.names, perl = T )
      row.names(X_DATA) = row.names

#      intersect(row.names, names)
      for(i in 1:nrow(X_DATA)){
          row.vector = as.numeric(X_DATA[i,])

          means[row.names(X_DATA)[i]] = mean(row.vector, na.rm=T)
          sdev = sd(row.vector, na.rm=T)
          ci.l[row.names(X_DATA)[i]] =  mean(row.vector, na.rm=T) - sdev
          ci.u[row.names(X_DATA)[i]] =  mean(row.vector, na.rm=T) + sdev
      }


      col.vector = sub(".+:", "", names, perl = T )
      #col.vector = sub(":", "", row.names, perl = T )
      #col.vector = 12*(7 - (as.numeric(col.vector)))
      col.vector = 10*(10 - (as.numeric(col.vector)))
      col.vector = paste("gray", col.vector, sep="")




      barplot2(means,  width=1,  las= 1,cex.names = 0.9, main=classes[j],
          #col = col.vector,
          col = "gray",
          #xlim= c(0,size),
          ylim= c(0,max(ci.u, na.rm = T)),
          #horiz = T,
          names.arg = names, las = 3,
          plot.ci = T, ci.l = ci.l, ci.u = ci.u
      )

     # dev.off()

  }
}   

mergeIsomers <- function (NORM_DATA){

      compounds = gsub(" +\\(\\d\\)", "", row.names(NORM_DATA), perl =T )
      isomers = which(duplicated(compounds))
      no.isomers = setdiff(1:nrow(NORM_DATA), isomers)
      NO_ISOMERS = NORM_DATA[no.isomers,]
      row.names(NO_ISOMERS) = gsub(" +\\(\\d\\)", "", row.names(NO_ISOMERS), perl =T )
      compounds.with.isomers = unique(compounds[isomers])
      for(i in 1:length(compounds.with.isomers)){
           inds = grep(compounds.with.isomers[i], row.names(NORM_DATA))
           COMPOUND = NORM_DATA[inds,]
           compound.sum = apply( COMPOUND, 2, sum, na.rm = T)
           NO_ISOMERS[compounds.with.isomers[i],] = compound.sum
      }
      
      return(NO_ISOMERS)
}

mergeAndPlotBarPlots <- function (X, Y, Classes){

#  X = mergeIsomers(X)
#  Y = mergeIsomers(Y)
    
  X[ X== Inf] = NA
  Y[ Y== Inf] = NA

  reps = colnames(X)
  lipid.names = row.names(X)
  #glycerolipids =  c("PC","PI","PE","PS","PG","SQDG","MGDG","DGDG")
  glycerolipids =  c("PC","PI","PE","PS","PG","SQDG","MGDG","DGDG","DAG")
  lysos =  c("Lyso-PC","Lyso-PE","Lyso-MGDG","Lyso-DGDG","MAG")
  pos.mode =  c("PC","PE","MGDG","DGDG") #  "DAG",
  neg.mode =  c("PG","PS","PI","SQDG") #  "DAG",
  
  c.34  = c("34:1","34:2","34:3","34:4","34:5","34:6")
      
  all.species = c("32:0","32:1","32:2","32:3","32:5","32:4","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "34:6","36:1","36:2","36:3","36:4", "36:5","36:6","36:7","38:0","38:1", "38:2","38:3", "38:4", "38:5", "38:6","40:0", 
        "40:1","40:2", "40:3","40:4","40:5","42:1","42:2","42:3","42:4")
        
  smaller.than.40 = c("32:0","32:1","32:2","32:3","32:5","32:4","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "34:6","36:1","36:2","36:3","36:4", "36:5","36:6","36:7","38:0","38:1", "38:2","38:3", "38:4", "38:5", "38:6")
         

        

  low.abundant =  c("32:0","32:1","32:2","32:3","32:5","32:6","34:0","34:1","34:2","34:3","34:4","34:5",
        "36:1","38:2","38:3","38:4","38:5","38:6")
       
            
  tag = c(  "48:0","48:1","48:2","48:3","48:4","48:5","48:6","48:7","48:8","48:9",
                  "50:0","50:1","50:2","50:3","50:4","50:5","50:6","50:7","50:8","50:9",
                  "52:0","52:1","52:2","52:3","52:4","52:5","52:6","52:7","52:8","52:9",
                  "54:0","54:1","54:2","54:3","54:4","54:5","54:6","54:7","54:8","54:9",
                  "56:0","56:1","56:2","56:3","56:4","56:5","56:6","56:7","56:8","56:9",
                  "58:0","58:1","58:2","58:3","58:4","58:5","58:6","58:7","58:8","58:9",
                  "60:0","60:1","60:2","60:3","60:4","60:5","60:6","60:7","60:8","60:9"
              
  ) 
      
  Lysos = c("16:0","16:1","16:2","16:3","18:0","18:1","18:2","18:3","20:1","20:1","22:0","22:1")
  
  if(Classes == "special"){
      classes =  c("PC", "DGDG", "MGDG" )
      names = c.34
      par(mfrow=c(1,3), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  if(Classes == "glycerolipids"){
      classes =  glycerolipids
      names = all.species
      par(mfrow=c(3,3), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "pos.mode"){
      classes =  pos.mode
      names = smaller.than.40
      par(mfrow=c(2,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "neg.mode"){
      classes =  neg.mode
      names = all.species
      par(mfrow=c(2,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "lysos"){
      classes =  lysos
      names = Lysos
      par(mfrow=c(3,2), mar=c(3.1, 2.5, 2.5, 2.1))
  
  }
  
  if(Classes == "tag"){
      classes =  "TAG"
      names = tag

  }
  
  library(gplots)

  for( j in 1:length(classes)){
  
      lipid.x.inds = grep (paste("^",classes[j], sep = ""), lipid.names, perl = T )
      #X_DATA = COND_DATA[lipid.x.inds,]
      X_DATA = X[lipid.x.inds, , drop = FALSE]
      x = rep(NA, length(lipid.x.inds))
      X_DATA = cbind(X_DATA, x)
      row.names(X_DATA) = row.names(X)[lipid.x.inds]
      Y_DATA = Y[lipid.x.inds, , drop = FALSE]
      y = rep(NA, length(lipid.x.inds))
      Y_DATA = cbind(Y_DATA, y)
      row.names(Y_DATA) = row.names(Y)[lipid.x.inds]

      #means = c()

      row.names = row.names(X_DATA)
      row.names = sub(" +$", "", row.names, perl = T )
      row.names = sub(".+ ", "", row.names, perl = T )
      row.names(X_DATA) = row.names
      
      X_DATA = X_DATA[intersect(row.names(X_DATA),names),]
            
      means = rep(NA, length(names))
      ci.l = rep(NA, length(names))
      ci.u = rep(NA, length(names))
      
      names(means) = names
      names(ci.l) = names
      names(ci.u) = names



      #names = intersect(names, row.names)
      for(i in 1:nrow(X_DATA)){

          row.vector = as.numeric(X_DATA[i,])

          #means[names[i]] = mean(row.vector, na.rm=T)
          #sdev = sd(row.vector, na.rm=T)
          #ci.l[names[i]] = means[names[i]] - sdev
          #ci.u[names[i]] = means[names[i]] + sdev
          
          
          means[row.names(X_DATA)[i]] = summary(row.vector)["Median"]
          ci.l[row.names(X_DATA)[i]] = summary(row.vector)["1st Qu."]
          ci.u[row.names(X_DATA)[i]] = summary(row.vector)["3rd Qu."]
          
#          means[names[i]] = median(row.vector, na.rm = T) 
#          ci.l[names[i]] =  min(row.vector, na.rm = T)
#          ci.u[names[i]] =  max(row.vector, na.rm = T)
    
#          means[names[i]] = summary(row.vector)["Median"]
#          ci.l[names[i]] =  summary(row.vector)["1st Qu."]
#          ci.u[names[i]] =  summary(row.vector)["3rd Qu."]
    
 
      }


      

      
      row.names(Y_DATA) = row.names
      Y_DATA = Y_DATA[intersect(row.names(Y_DATA),names),]
      
      means1 = rep(NA, length(names))
      ci.l1 = rep(NA, length(names))
      ci.u1 = rep(NA, length(names))
      
      names(means1) = names
      names(ci.l1) = names
      names(ci.u1) = names
      
  
      for(i in 1:nrow(Y_DATA)){

          row.vector = as.numeric(Y_DATA[i,])

        #  means1[names[i]] = mean(row.vector, na.rm=T)
        #  sdev = sd(row.vector, na.rm=T)
        #  ci.l1[names[i]] = means1[names[i]] - sdev
        #  ci.u1[names[i]] = means1[names[i]] + sdev
          
          
          means1[row.names(Y_DATA)[i]] = summary(row.vector)["Median"]
          ci.l1[row.names(Y_DATA)[i]] = summary(row.vector)["1st Qu."]
          ci.u1[row.names(Y_DATA)[i]] = summary(row.vector)["3rd Qu."]
          
#          means1[names[i]] = median(row.vector, na.rm = T) 
#          ci.l1[names[i]] =  min(row.vector, na.rm = T)
#          ci.u1[names[i]] =  max(row.vector, na.rm = T)
          
#          means1[names[i]] = summary(row.vector)["Median"]
#          ci.l1[names[i]] =  summary(row.vector)["1st Qu."]
#          ci.u1[names[i]] = summary(row.vector)["3rd Qu."]
          
    
          
      }

      both.means = rep(NA, 2*(length(means)))
      both.ci.l =  rep(NA, 2*(length(ci.l)))
      both.ci.u =  rep(NA, 2*(length(ci.u)))
      both.names = rep("", 2*(length(names)))

      y.inds = 2*(1:length(means))
      x.inds = y.inds -1

      both.means[x.inds] = means
      both.means[y.inds] = means1

      both.ci.l[x.inds] = ci.l
      both.ci.l[y.inds] = ci.l1

      both.ci.u[x.inds] = ci.u
      both.ci.u[y.inds] = ci.u1

      both.names[y.inds] = names

      both.means[both.means == Inf] = NA 
      both.ci.l[both.ci.l == Inf] = NA
      both.ci.u[both.ci.u == Inf] = NA
      barplot2(both.means,  width=1,  las= 1,cex.names = 0.9, main=classes[j],
          #col = col.vector,
          col = rep(c("white","gray"), 28) ,
          #xlim= c(0,size),
          #ylim= c(0,size),
          #horiz = T,
          names.arg = both.names,
          las = 3,
          plot.ci = T, ci.l = both.ci.l, ci.u = both.ci.u
      )


     # dev.off()

  }
}


##############################################################################################
# Here is where the normalization starts wo_WT
##############################################################################################
#--------------------------------------------------------------------------------------------#
# This part of the code was used to merge all the raw data sets into one file: DO NOT RUN AGAIN!!!!

RAW_DATA_POS_1 =  read.delim("raw_intensities_121204_pos.txt", row.names = 1 )
RAW_DATA_POS_2 =  read.delim("raw_intensities_130122_pos.txt", row.names = 1 )
RAW_DATA_POS_3 =  read.delim("raw_intensities_130314_pos.txt", row.names = 1 )

myb = grep("X1_|X2_|X3_|X4_|X5_|X6_|X7_|X9_|X10_|X11_|Tray_1|Tray_2|Tray_3",colnames(RAW_DATA_POS_2))
not.myb = setdiff(1:ncol(RAW_DATA_POS_2), myb)
RAW_DATA_POS_2 = RAW_DATA_POS_2[,not.myb]

#NORM_DATA_2 = normToTIC2(RAW_DATA_POS_2)
#write.table(NORM_DATA, file='norm_intensities_130122.txt', sep="\t", quote = F)

labels =  read.delim("sample_description.txt" , row.names = NULL)

row.int = intersect(row.names(RAW_DATA_POS_1), row.names(RAW_DATA_POS_2) )
row.int = intersect( row.int , row.names(RAW_DATA_POS_3))

RAW_DATA = cbind(RAW_DATA_POS_1[row.int,], RAW_DATA_POS_2[row.int,], RAW_DATA_POS_3[row.int,])

colnames(RAW_DATA) = sub("X","", colnames(RAW_DATA) )
colnames(RAW_DATA) = sub("_pos","", colnames(RAW_DATA) )
colnames(RAW_DATA) = sub("tray","T", colnames(RAW_DATA) )  
colnames(RAW_DATA) = sub("t|Tray","T", colnames(RAW_DATA) )
colnames(RAW_DATA) = sub("Col_0","_Col-0", colnames(RAW_DATA) )
colnames(RAW_DATA) = sub("Col0","_Col-0", colnames(RAW_DATA) )
colnames(RAW_DATA) = gsub("_+","_", colnames(RAW_DATA) )
colnames(RAW_DATA) = gsub("T_","T", colnames(RAW_DATA) )

setdiff(colnames(RAW_DATA), labels$Sample.name)
setdiff(labels$Sample.name, colnames(RAW_DATA))
RAW_DATA1 =  RAW_DATA[, intersect(labels$Sample.name,colnames(RAW_DATA) )]

write.table(RAW_DATA1, file='raw_intensities.txt', sep="\t", quote = F)


##############################################################################################
#--------------------------------------------------------------------------------------------#

# Here comes the normalization to TIC:


setwd("C:/Users/sandr/OneDrive - Universidad de Antioquia/Publications/Paper 5/KO_lines/DataNormalization/Data_processing")
#RAW_DATA_POS =  read.delim("raw_intensities_All.txt", row.names = 1 )
#RAW_DATA_NEG =  read.delim("raw_intensities_neg.txt", row.names = 1 )
#RAW_DATA_POS =  read.delim("raw_intensities_pos.txt", row.names = 1 )


# Import .csv file with data:
# Sheet names:   data_processed_minValue.csv    |    data_processed_PPCA.csv
RAW_DATA_POS <- read.csv("data_processed_PPCA.csv", head = TRUE, sep=";", row.names = 1)


# Delete row with label information and set column names:
RAW_DATA_POS <- RAW_DATA_POS[-1,]


#colnames(RAW_DATA_POS) = sub("_pos","",colnames(RAW_DATA_POS),)
#colnames(RAW_DATA_NEG) = sub("_neg","",colnames(RAW_DATA_NEG),)

#int = intersect(colnames(RAW_DATA_POS), colnames(RAW_DATA_NEG))
#RAW_DATA = rbind(RAW_DATA_POS[,int], RAW_DATA_NEG[,int])

# Carry-out normalization to TIC:
RAW_DATA <- RAW_DATA_POS
NORM_DATA = normToTIC2(RAW_DATA)
write.table(NORM_DATA, file='norm_to_TIC.txt', sep="\t", quote = F)



##############################################################################################
##############################################################################################
# Here are other functions to probably plot the results but they are not working properly!!!!!
#--------------------------------------------------------------------------------------------#


MED_DATA = getMedData(NORM_DATA)
MED_DATA1 = mergeIsomers(MED_DATA)
#-----------------------------------------------------------------------------------------#
wt = grep("Col", colnames(NORM_DATA))
mut = setdiff(1:ncol(NORM_DATA), wt)

MUTANTS = NORM_DATA[,mut]
MUTANTS_MED = getMedData(MUTANTS)

genotype = sub("_.+","", colnames(MUTANTS_MED))

labels1 =  read.delim("mutant_labels_and_descriptions~.txt" , row.names = NULL)
sample.name = labels1$Gene
names(sample.name) = labels1$Tube.number
setdiff(genotype, labels1$Tube.number)
new.colnames = sample.name[genotype]
order = intersect(labels1$Tube.number, genotype)
new.colnames1 = new.colnames[order]

MUTANTS_MED1 = MUTANTS_MED[,order]
colnames(MUTANTS_MED1) = new.colnames1

MUTANTS_MED2 = mergeIsomers(MUTANTS_MED1)
plotHeatMap(MUTANTS_MED2)

BIG_CHANGES = locateBigChanges(MUTANTS_MED2)
plotHeatMap4(BIG_CHANGES)

CLASS_TOTALS = getClassTotal(MUTANTS_MED2)
CLASS_TOTALS1 = CLASS_TOTALS[,order(CLASS_TOTALS["TAG",])]
plotHeatMap(CLASS_TOTALS1)

wt = grep("Col", colnames(NORM_DATA))
WT = NORM_DATA[,wt]

#-------------------------------------------------------------------------------#

mut1.number = which(labels1$Gene == "NMT3") 
mut1.number1 = paste("_", mut1.number,"_", sep = "")
mut1.inds = grep(mut1.number1, paste("_", colnames(NORM_DATA), sep = "") )
mut1 = NORM_DATA[,mut1.inds]

mut2.number = which(labels1$Gene == "LPD2") 
mut2.number1 = paste("_", mut2.number,"_", sep = "")
mut2.inds = grep(mut2.number1, paste("_", colnames(NORM_DATA), sep = "") )
mut2 = NORM_DATA[,mut2.inds]

mergeAndPlotBarPlots(mut1, mut2, "tag")

mut = NORM_DATA[,grep("_69_", paste("_",colnames(NORM_DATA), sep = ""))]
wt = NORM_DATA[,grep("T10_", colnames(NORM_DATA), perl = T)]

srs = NORM_DATA[,grep("_114_", paste("_",colnames(NORM_DATA), sep = ""))]
wt = NORM_DATA[,grep("T15_", colnames(NORM_DATA), perl = T)]

mergeAndPlotBarPlots(wt, srs, "tag")




#--------------------------------------------------------------------------------#

row.names(labels) = labels$Sample.name
labels = labels[intersect(labels$Sample.name,colnames(RAW_DATA)),]

NORM_BY_M_DAY = normToFactor(RAW_DATA1, "Measurement.day")
NORM_BY_TRAY = normToFactor(RAW_DATA1, "Tray")

NORM_DATA = normToTIC2(NORM_BY_TRAY)

NORM_DATA = normToTIC2(RAW_DATA1)

write.table(NORM_DATA, file='norm_intensities~.txt', sep="\t", quote = F)

col.sum = apply(NORM_DATA, 2, sum, na.rm=T)
valid = setdiff(colnames(NORM_DATA),names(which(col.sum == 0)))
NORM_DATA = NORM_DATA[,valid]

MED_DATA =getMedData(NORM_DATA)

wt = grep("Col", colnames(NORM_DATA))
mut = setdiff(1:ncol(NORM_DATA), wt)

MUTANTS = NORM_DATA[,mut]
MUTANTS_MED = getMedData(MUTANTS)

genotype = sub("_.+","", colnames(MUTANTS_MED))

labels1 =  read.delim("mutant_labels_and_descriptions~.txt" , row.names = NULL)
sample.name = labels1$Gene
names(sample.name) = labels1$Tube.number
setdiff(genotype, labels1$Tube.number)
new.colnames = sample.name[genotype]
order = intersect(labels1$Tube.number, genotype)
new.colnames1 = new.colnames[order]

MUTANTS_MED1 = MUTANTS_MED[,order]
colnames(MUTANTS_MED1) = new.colnames1

jpg("ads_heatmap.jpg")
MUTANTS_MED2 = mergeIsomers(MUTANTS_MED1)
plotHeatMap(MUTANTS_MED2)
dev.off()

pdf("ads_heatmap.pdf")
MUTANTS_MED2 = mergeIsomers(MUTANTS_MED1)
plotHeatMap(MUTANTS_MED2)
dev.off()


myb = grep("myb", colnames(MUTANTS_MED2))
no.myb = setdiff(1:ncol(MUTANTS_MED2), myb)
MUTANTS_MED3 = MUTANTS_MED2[,no.myb]
plotHeatMap(MUTANTS_MED3)

DELTAS = locateBigChanges(MUTANTS_MED2)
plotHeatMap4(DELTAS)
