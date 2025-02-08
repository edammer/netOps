netORA <- function(inputDat1,
                   inputDat2,
                   net1,
                   net2,
                   species1 = "hsapiens",
                   species2 = "hsapiens",
                   outfiletitle = "Network1_vs_Network2",
                   heatmapTitle = "Network1 Modules vs Network2 Modules",
                   net1.prefix = "Net1",
                   net2.prefix = "Net2",
                   oneOrTwoTailed = 2,
                   orderByRelatedness = TRUE,
                   columnReorderOfAutoReorderedOutput = c(),
                   rootdir = getwd()) {

  # Load required libraries
  cat(paste0("- Loading WGCNA package ...\n"))
  require(WGCNA, quietly=TRUE)
  
  ###############################################
  ## Helper functions for Fisher’s Exact Test ##
  ###############################################
  ORoneTailed <- function(q, k, m, t) {
    fisher.out <- fisher.test(matrix(c(q, k - q, m - q, t - m - k + q), 2, 2),
                                conf.int = TRUE,
                                alternative = "greater")
    OR <- fisher.out$estimate
    pval <- fisher.out$p.value
    upCI <- fisher.out$conf.int[1]
    downCI <- fisher.out$conf.int[2]
    
    output <- c(OR, pval, upCI, downCI)
    names(output) <- c("OR", "Fisher p", "-95%CI", "+95%CI")
    return(output)
  }
  
  ORtwoTailed <- function(q, k, m, t) {
    fisher.out <- fisher.test(matrix(c(q, k - q, m - q, t - m - k + q), 2, 2),
                                conf.int = TRUE)
    OR <- fisher.out$estimate
    pval <- fisher.out$p.value
    upCI <- fisher.out$conf.int[1]
    downCI <- fisher.out$conf.int[2]
    
    output <- c(OR, pval, upCI, downCI)
    names(output) <- c("OR", "Fisher p", "-95%CI", "+95%CI")
    return(output)
  }
  
  ORA.core <- function(testpath, refpath, testbackground, refbackground) {
    q <- length(intersect(testpath, refpath))          ## overlap count
    k <- length(intersect(refpath, testbackground))      ## size of reference list
    m <- length(intersect(testpath, refbackground))      ## size of input (module) list
    t <- length(intersect(testbackground, refbackground))## total background
    
    if (oneOrTwoTailed == 1) {
      empvals <- ORoneTailed(q, k, m, t)
    } else {
      empvals <- ORtwoTailed(q, k, m, t)
    }
    
    tmpnames <- names(empvals)
    empvals <- as.character(c(empvals, q, k, m, t, 100 * signif(q/k, 3)))
    names(empvals) <- c(tmpnames,
                        "Overlap",
                        "Reference List",
                        "Input List",
                        "Background",
                        "% List Overlap")
    return(empvals)
  }
  
  ###############################################
  ## Prepare gene info from the input datasets ##
  ###############################################
  geneInfo.net1 <- as.data.frame(inputDat1)
  geneInfo.net1$Symbol <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(inputDat1)), "[|]")))[, 1]
  geneInfo.net1$UniqueID<-as.character(rownames(inputDat1))
  
  geneInfo.net2 <- as.data.frame(inputDat2)
  geneInfo.net2$Symbol <- as.data.frame(do.call(rbind, strsplit(as.character(rownames(inputDat2)), "[|]")))[, 1]
  geneInfo.net2$UniqueID<-as.character(rownames(inputDat2))
  
  ## Build a gene list from the first dataset
  genelist1 <- data.frame(unlist(geneInfo.net1$Symbol), ncol = 1)
  colnames(genelist1) <- "species1"
  genelist1 <- data.frame(genelist1[which(!genelist1[, 1] == ""), "species1"], ncol = 1)
  colnames(genelist1) <- "species1"
  genelist1 <- data.frame(genelist1$species1[!grepl("'", genelist1$species1)], ncol = 1)
  colnames(genelist1) <- c("species1", "species2")
  
  ###############################################
  ## Convert gene symbols if species differ  ##
  ###############################################
  if(!species1==species2) {
    require(biomaRt, quietly=TRUE)
    if(species1=="hsapiens") {
      cat(paste0("- Converting ",species2," to ",species1," for lists in Second Network ... \n"))

      genelist2 <- data.frame(unlist(geneInfo.net2$Symbol), ncol = 1)
      colnames(genelist2) <- "species2"
      genelist2 <- data.frame(genelist2[which(!genelist2[, 1] == ""), "species2"], ncol = 1)
      colnames(genelist2) <- "species2"
      genelist2 <- data.frame(genelist2$species2[!grepl("'", genelist2$species2)], ncol = 1)
      colnames(genelist2) <- c("species2", "species1")

      conv.fromSpecies=species2
      conv.toSpecies=species1
    } else {
      cat(paste0("- Converting ",species1," to ",species2," for lists in First Network ... \n"))

      conv.fromSpecies=species1
      conv.toSpecies=species2
    }
  
    this.heatmapTitle=paste0(heatmapTitle," in ",species1," homologs")
    #library(biomaRt)
  
    #human = useEnsembl("genes", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021
    #mouse = useEnsembl("genes", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org") 
  
    this.species = useEnsembl("genes", dataset=paste0(conv.toSpecies,"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")  #ver=105 equivalent to dec2021  ; old code: #useMart("ensembl",dataset=paste0(species1,"_gene_ensembl"))
    other.species = useEnsembl("genes", dataset=paste0(conv.fromSpecies,"_gene_ensembl"), host="https://dec2021.archive.ensembl.org")   # old code: #useMart("ensembl",dataset=paste0(species2,"_gene_ensembl"))
	
    #category species to other species conversion (first column is other species)
    if(species1=="hsapiens") {
#      genelist.conv<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=genelist2$species2, mart=other.species, attributesL="external_gene_name",martL = this.species)  # timeout
      all.ids<-getBM(attributes=c("ensembl_gene_id"), mart=other.species)
      genelist.conv<-getLDS(attributes=c("external_gene_name","ensembl_gene_id"), filters="ensembl_gene_id", values=all.ids, mart=other.species, attributesL="hgnc_symbol",martL = this.species)
      genelist2$convertedSymbol <- genelist.conv[match(genelist2$species2,genelist.conv$Gene.name),"HGNC.symbol"]

      geneInfo.net2$Symbol.original<-geneInfo.net2$Symbol
      geneInfo.net2$Symbol<-genelist2$convertedSymbol[match(geneInfo.net2$Symbol.original,genelist2$species2)]
    } else {
      if(species2=="hsapiens") {
#        genelist.conv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=genelist1$species1, mart=other.species, attributesL="hgnc_symbol",martL = this.species)  # timeout
        all.ids<-getBM(attributes=c("ensembl_gene_id"), mart=other.species)
        genelist.conv<-getLDS(attributes=c("external_gene_name","ensembl_gene_id"), filters="ensembl_gene_id", values=all.ids, mart=other.species, attributesL="hgnc_symbol",martL = this.species)
        genelist1$convertedSymbol <- genelist.conv[match(genelist1$species1,genelist.conv$Gene.name),"HGNC.symbol"]
      } else {
#        genelist.conv<-getLDS(attributes=c("external_gene_name"), filters="external_gene_name", values=genelist1$species1, mart=other.species, attributesL="external_gene_name",martL = this.species)  # timeout
        all.ids<-getBM(attributes=c("ensembl_gene_id"), mart=other.species)
        genelist.conv<-getLDS(attributes=c("external_gene_name","ensembl_gene_id"), filters="ensembl_gene_id", values=all.ids, mart=other.species, attributesL="external_gene_name",martL = this.species)
        genelist1$convertedSymbol <- genelist.conv[match(genelist1$species1,genelist.conv$Gene.name),"Gene.name.1"]
      }
      geneInfo.net1$Symbol.original<-geneInfo.net1$Symbol
      geneInfo.net1$Symbol<-genelist1$convertedSymbol[match(geneInfo.net1$Symbol.original,genelist1$species1)]
    }
  }
	
  #Reduce available data to module members with non-missing symbols
  geneInfo.net1$net.colors<-net1$colors
  geneInfo.net2$net.colors<-net2$colors
  
  net1.keep.idx=intersect(which(!is.na(geneInfo.net1$Symbol)),which(geneInfo.net1$Symbol!=""))
  net2.keep.idx=intersect(which(!is.na(geneInfo.net2$Symbol)),which(geneInfo.net2$Symbol!=""))
  geneInfo.net1<-geneInfo.net1[net1.keep.idx,]
  geneInfo.net2<-geneInfo.net2[net2.keep.idx,]
    
  net1$colors<-geneInfo.net1$net.colors
  net2$colors<-geneInfo.net2$net.colors
	
  #Define backgrounds for each network
  background.net1 <- na.omit(as.character(unique(geneInfo.net1$Symbol)))
  background.net2 <- na.omit(as.character(unique(geneInfo.net2$Symbol)))

  
  ##############################################################
  ## Loop over two types of tests ("FDR" and "p") and produce ##
  ## heatmaps (including a “reordered” version)               ##
  ##############################################################
  for (test in c("FDR", "p")) {
    
    ## Set module ordering – either according to relatedness (default) or by size‐rank.
    if (orderByRelatedness) {
      if (length(which(gsub("ME", "", colnames(net1$MEs)) == "grey")) > 0) {
        net1.colorOrder <- gsub("ME", "", colnames(net1$MEs))[-which(gsub("ME", "", colnames(net1$MEs)) == "grey")]
      } else {
        net1.colorOrder <- gsub("ME", "", colnames(net1$MEs))
      }
      uniquemodcolors.net1 <- net1.colorOrder
      if (length(which(gsub("ME", "", colnames(net2$MEs)) == "grey")) > 0) {
        net2.colorOrder <- gsub("ME", "", colnames(net2$MEs))[-which(gsub("ME", "", colnames(net2$MEs)) == "grey")]
      } else {
        net2.colorOrder <- gsub("ME", "", colnames(net2$MEs))
      }
      uniquemodcolors.net2 <- net2.colorOrder
      
      nModules <- max(length(net1.colorOrder), length(net2.colorOrder))
      orderedModules <- data.frame(Mnum = seq(1:nModules), Color = labels2colors(1:nModules))
      Emory.Order <- as.numeric(orderedModules$Mnum[match(uniquemodcolors.net2, orderedModules$Color)])
      BLSA.Order <- as.numeric(orderedModules$Mnum[match(uniquemodcolors.net1, orderedModules$Color)])
    } else {
      net1.colorCount <- length(unique(net1$colors)) - length(which(unique(net1$colors) == "grey"))
      net2.colorCount <- length(unique(net2$colors)) - length(which(unique(net2$colors) == "grey"))
      
      Emory.Order <- seq(1, net2.colorCount)
      BLSA.Order <- seq(1, net1.colorCount)
      uniquemodcolors.net1 <- labels2colors(BLSA.Order)
      uniquemodcolors.net2 <- labels2colors(Emory.Order)
    }

    cat(paste0("- [",test,"] Calculating Fisher's exact tests for network overlaps by gene product ... \n"))
    
    ## Initialize matrices to store ORA results.
    ORmat <- matrix(NA, nrow = length(uniquemodcolors.net2), ncol = 1)
    Pmat <- matrix(NA, nrow = length(uniquemodcolors.net2), ncol = 1)
    overlapList <- matrix(NA, ncol = length(uniquemodcolors.net1), nrow = length(uniquemodcolors.net2))
    
    ## Loop through each pair of modules (from network 1 and network 2)
    for (i in 1:length(uniquemodcolors.net1)) {
      thismod <- uniquemodcolors.net1[i]
      thisGene <- geneInfo.net1$Symbol[net1$colors == thismod]
      testpath <- as.character(unique(thisGene))
      oraMat1 <- matrix(NA, ncol = 9, nrow = length(uniquemodcolors.net2))
      
      for (j in 1:length(uniquemodcolors.net2)) {
        thismod1 <- uniquemodcolors.net2[j]
        thisGene1 <- geneInfo.net2$Symbol[net2$colors == thismod1]
        refpath <- as.character(unique(thisGene1))
        testbackground <- background.net1
        refbackground <- background.net2
        
        oraout <- ORA.core(testpath, refpath, testbackground, refbackground)
        oraMat1[j, ] <- oraout
        
        ## Also record the actual overlapping symbols.
        overlapList1 <- intersect(refpath, testpath)
        if (length(overlapList1) > 0)
          overlapList[j, i] <- paste(overlapList1, collapse = ",")
      }
      
      ORmat <- cbind(ORmat, as.numeric(oraMat1[, 1]))
      Pmat <- cbind(Pmat, as.numeric(oraMat1[, 2]))
    }
    
    ORmat.Array <- ORmat[, -1]
    Pmat.Array <- Pmat[, -1]
    FDRmat.Array <- matrix(p.adjust(Pmat.Array, method = "BH"),
                           nrow = nrow(Pmat.Array),
                           ncol = ncol(Pmat.Array))
    
    colnames(overlapList) <- colnames(ORmat.Array) <- colnames(Pmat.Array) <-
      colnames(FDRmat.Array) <- paste(net1.prefix, uniquemodcolors.net1, sep = ".")
    rownames(overlapList) <- rownames(ORmat.Array) <- rownames(Pmat.Array) <-
      rownames(FDRmat.Array) <- paste(net2.prefix, uniquemodcolors.net2, sep = ".")
    
    #########
    ## Choose display matrix based on test type.
    if (test == "p") {
      dispMat <- -log10(Pmat.Array) * sign(log2(ORmat.Array))
      FDRmat.Array <- Pmat.Array
      testtext <- "-log10(P Value)"
      testfileext <- "pval"
    } else {
      outputData <- rbind("FET pValue", Pmat.Array,
                          "FDR (BH) corrected", FDRmat.Array,
                          "Overlap (Symbols)", overlapList)
      write.csv(outputData,
                file = paste0(rootdir, "/ORA-", oneOrTwoTailed, "tailed-", outfiletitle, ".FULL.csv"))
      dispMat <- -log10(FDRmat.Array) * sign(log2(ORmat.Array))
      testtext <- "-log10(FDR, BH)"
      testfileext <- "FDR"
    }
    
    ## Create text matrix for heatmap labels.
    txtMat <- dispMat
    txtMat[FDRmat.Array >= 0.05] <- ""
    txtMat[FDRmat.Array < 0.05 & FDRmat.Array > 0.01] <- "*"
    txtMat[FDRmat.Array < 0.01 & FDRmat.Array > 0.005] <- "**"
    txtMat[FDRmat.Array < 0.005] <- "***"
    
    txtMat1 <- signif(dispMat, 2)
    txtMat1[txtMat1 < 1.5] <- ""
    textMatrix1 <- matrix(paste(txtMat1, "\n", txtMat, sep = ""),
                          ncol = ncol(Pmat.Array),
                          nrow = nrow(Pmat.Array))
    
    ## Set up colors for the heatmap.
    bw <- colorRampPalette(c("#0058CC", "white"))
    wr <- colorRampPalette(c("white", "#CC3300"))
    colvec <- if (oneOrTwoTailed == 1) { 
      wr(100) 
    } else { 
      c(bw(50), wr(50)) 
    }
    zvec <- if (oneOrTwoTailed == 1) { 
      c(0, 10) 
    } else { 
      c(-10, 10) 
    }
    
    ## Plot the primary heatmap.
    pdf(file = paste0(rootdir, "/ORA-", oneOrTwoTailed, "tailed-", outfiletitle,
                      ".FULL.", testfileext, ".pdf"), width = 16, height = 8)
    par(mar = c(8, 12, 3, 3))
    par(mfrow = c(1, 1))
    
    suppressWarnings(
      labeledHeatmap(Matrix = dispMat,
                     colorLabels = TRUE,
                     setStdMargins = FALSE,
                     yLabels = paste("ME", uniquemodcolors.net2, sep = ""),
                     ySymbols = as.vector(paste(net2.prefix, "-M", Emory.Order, sep = "")),
                     xLabels = paste("ME", uniquemodcolors.net1, sep = ""),
                     xSymbols = as.vector(paste(net1.prefix, "-M", BLSA.Order, sep = "")),
                     colors = colvec,
                     textMatrix = textMatrix1,
                     cex.text = 0.55,
                     cex.lab.x = 1,
                     zlim = zvec,
                     main = paste0(heatmapTitle, " ORA/FET Signed ", testtext))
    )
    dev.off()
    
    #############################################
    ## Reorder the matrix to create a “diagonal” ##
    ## positive overlap pattern (if any)         ##
    #############################################
    orderVecEmory <- matrix(NA, ncol = 1, nrow = 1)
    orderVecBLSA <- matrix(NA, ncol = 1, nrow = 1)
    BetterThanCutoff <- 1  ## adjust this cutoff as desired
    
    for (j in 1:length(uniquemodcolors.net1)) {
      bestMatch <- which(dispMat[, j] == max(dispMat[, j])[1])[1]
      if (dispMat[bestMatch, j] > BetterThanCutoff) {
        orderVecBLSA <- rbind(orderVecBLSA, j)
        truncatedRow <- dispMat[-bestMatch, j]
        if (is.na(orderVecEmory[1, 1])) {
          orderVecEmory <- matrix(bestMatch, ncol = 1)
        } else {
          if (is.na(match(bestMatch, orderVecEmory[, 1])))
            orderVecEmory <- rbind(as.matrix(orderVecEmory, ncol = 1), bestMatch)
        }
        bestMatch1 <- which(dispMat[, j] == max(truncatedRow)[1])[1]
        truncatedRow1 <- truncatedRow[-match(max(truncatedRow), truncatedRow)]
        if (is.na(match(bestMatch1, orderVecEmory)) & abs(dispMat[bestMatch1, j]) > BetterThanCutoff) {
          orderVecEmory <- rbind(as.matrix(orderVecEmory, ncol = 1), bestMatch1)
          bestMatch2 <- which(dispMat[, j] == max(truncatedRow1)[1])
          truncatedRow2 <- truncatedRow1[-which(truncatedRow1 == bestMatch2[1])]
          if (is.na(match(bestMatch2[1], orderVecEmory)) & abs(dispMat[bestMatch2[1], j]) > BetterThanCutoff) {
            orderVecEmory <- rbind(as.matrix(orderVecEmory, ncol = 1), bestMatch2[1])
          }
        }
      }
    }
    
    if (!dim(orderVecBLSA)[1] == 1) {  ## if any modules pass the cutoff
      Emory.ReOrder <- as.vector(orderVecEmory)
      BLSA.ReOrder <- as.vector(orderVecBLSA[-1, 1])
      
      ## Add back any modules that did not pass the cutoff.
      BLSA.ReOrder <- c(BLSA.ReOrder,
                        setdiff(seq(1, (length(table(net1$colors)) -
                                           length(which(unique(net1$colors) == "grey")))),
                                BLSA.ReOrder))
      Emory.ReOrder <- c(Emory.ReOrder,
                         setdiff(seq(1, (length(table(net2$colors)) -
                                           length(which(unique(net2$colors) == "grey")))),
                                 Emory.ReOrder))
      
      ## Optionally re-order columns manually if the length matches.
      if (length(columnReorderOfAutoReorderedOutput) == length(BLSA.ReOrder))
        BLSA.ReOrder <- as.numeric(as.character(BLSA.ReOrder)[columnReorderOfAutoReorderedOutput])
      
      uniquemodcolors.ordered.net2 <- labels2colors(Emory.Order[Emory.ReOrder])
      uniquemodcolors.ordered.net1 <- labels2colors(BLSA.Order[BLSA.ReOrder])
      dispMat.ordered <- dispMat[Emory.ReOrder, BLSA.ReOrder]
      FDRmat.Array.ordered <- FDRmat.Array[Emory.ReOrder, BLSA.ReOrder]
      
      txtMat <- dispMat.ordered
      txtMat[FDRmat.Array.ordered >= 0.05] <- ""
      txtMat[FDRmat.Array.ordered < 0.05 & FDRmat.Array.ordered > 0.01] <- "*"
      txtMat[FDRmat.Array.ordered < 0.01 & FDRmat.Array.ordered > 0.005] <- "**"
      txtMat[FDRmat.Array.ordered < 0.005] <- "***"
      
      txtMat1 <- signif(dispMat.ordered, 2)
      txtMat1[txtMat1 < 1.5] <- ""
      textMatrix.ordered <- matrix(paste(txtMat1, "\n", txtMat, sep = ""),
                                   ncol = length(BLSA.ReOrder),
                                   nrow = length(Emory.ReOrder))
      
      dispMat.ordered.untransposed <- t(dispMat.ordered)
      textMatrix.ordered.untransposed <- t(textMatrix.ordered)
      Emory.ReOrder2 <- Emory.ReOrder
      
      pdf(file = paste0(rootdir, "/ORA-", oneOrTwoTailed, "tailed-", outfiletitle,
                        ".reordered.GTcutoff", BetterThanCutoff, ".", testfileext, ".pdf"),
          width = 16, height = 8)
      par(mar = c(8, 12, 3, 3))
      par(mfrow = c(1, 1))
      
      suppressWarnings(
        labeledHeatmap(Matrix = t(dispMat.ordered.untransposed),
                       colorLabels = TRUE,
                       setStdMargins = FALSE,
                       yLabels = paste("ME", uniquemodcolors.ordered.net2, sep = ""),
                       ySymbols = as.vector(paste(net2.prefix, "-M", Emory.Order[Emory.ReOrder2], sep = "")),
                       xLabels = paste("ME", uniquemodcolors.ordered.net1, sep = ""),
                       xSymbols = as.vector(paste(net1.prefix, "-M", BLSA.Order[BLSA.ReOrder], sep = "")),
                       colors = colvec,
                       textMatrix = t(textMatrix.ordered.untransposed),
                       cex.text = 0.55,
                       cex.lab.x = 1,
                       zlim = zvec,
                       main = paste0(heatmapTitle, " (Reordered) ORA/FET Signed ", testtext))
      )
      suppressWarnings(
        labeledHeatmap(Matrix = t(dispMat.ordered),
                       colorLabels = TRUE,
                       setStdMargins = FALSE,
                       xLabels = paste("ME", uniquemodcolors.ordered.net2, sep = ""),
                       xSymbols = as.vector(paste(net2.prefix, "-M", Emory.Order[Emory.ReOrder2], sep = "")),
                       yLabels = paste("ME", uniquemodcolors.ordered.net1, sep = ""),
                       ySymbols = as.vector(paste(net1.prefix, "-M", BLSA.Order[BLSA.ReOrder], sep = "")),
                       colors = colvec,
                       textMatrix = t(textMatrix.ordered),
                       cex.text = 0.55,
                       cex.lab.x = 1,
                       zlim = zvec,
                       main = paste0(heatmapTitle, " (Reordered, Transposed) ORA/FET Signed ", testtext))
      )
      dev.off()
    }
  } # end test loop
}
