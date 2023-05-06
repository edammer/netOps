buildIgraphs <- function(dummyVar="", env=.GlobalEnv) {

  if (!exists("symbols2Add")) symbols2Add=c()
  if (!exists("GOIlist")) GOIlist=c()
  if (!exists("recalcMEs")) recalcMEs=TRUE
  if (!exists("outFilePrefix")) { outFilePrefix="" } else { outFilePrefix=paste0(outFilePrefix,".") }
  if (!exists("outFileSuffix") & exists("FileBaseName")) { outFileSuffix=FileBaseName } else { if(!exists("outFileSuffix")) outFileSuffix="Unnamed_Net" }
  if (!exists("showTOMhubs")) showTOMhubs=FALSE
  if (showTOMhubs & (!exists("keepTopEdges") & !exists("keepTopPercentEdges"))) { cat(" - showTOMhubs=TRUE without specifying edges to keep/use. Using 150.\n"); keepTopEdges=150; }
  if (!exists("keepTopEdges") & !exists("keepTopPercentEdges"))  { cat(" - # or % of edges to show in iGraph not specified. Using 700.\n"); keepTopEdges=700; }
  if (suppressWarnings(is.na(as.integer(keepTopEdges))))  { cat(" - keepTopEdges is not coercible to an integer. Using 150.\n"); keepTopEdges=150; }
  if (exists("keepTopPercentEdges")) if(suppressWarnings(is.na(as.numeric(keepTopPercentEdges))))  { cat(" - keepTopPercentEdges is not coercible to numeric value. Using 15.\n"); keepTopPercentEdges=15; }
  if (!exists("vertexsize")) { cat(" - vertexsize not specified. Using 16 for large nodes.\n"); vertexsize=16; }
  if (showTOMhubs & !is.numeric(power)) { cat(" - showTOMhubs=TRUE, but we need a power to calculate TOM. Trying power=8.\n"); power=8; }
  if (!exists("PPIedges")) { cat(" - not using PPIedges from BioGRID.\n"); PPIedges=FALSE; } else {
    if(!exists("myHumanBioGrid.tsvFile")) myHumanBioGrid.tsvFile="nonexistent.file"
    if(!file.exists(myHumanBioGrid.tsvFile))  {
      if (interactive()) {
        suppressPackageStartupMessages(require(rvest,quietly=TRUE))
        linksText<-html_nodes(read_html("https://downloads.thebiogrid.org/BioGRID"), xpath="//a")
        bioGRIDcurrentRelease.link <- html_attr(linksText[which(grepl("[Cc]urrent.[Rr]elease",linksText))], "href")
        # bioGRIDcurrentRelease.link <- bioGRIDcurrentRelease.link[grepl("^[A-z].*\\/$",species.links)]
        linksText2<-html_nodes(read_html(bioGRIDcurrentRelease.link), xpath="//a")
        bioGRIDorganism.mitabPage.link <- html_attr(linksText2[which(grepl("[Dd]ownload/.*BIOGRID-ORGANISM-.*.mitab.zip",linksText2))], "href")
  
        cat(" - myHumanBioGrid.tsvFile (Human BioGRID TSV file) not specified or not found: ", myHumanBioGrid.tsvFile,"\n\n")
        print(data.frame(File=bioGRIDorganism.mitabPage.link))
        input.idx <- readline(paste0("[INTERACTIVE]\nChoose a BioGRID .ZIP file download with organism-specific binary PPI files inside\n(only current ver shown from https://thebiogrid.org/) [1-",length(bioGRIDorganism.mitabPage.link),"]: "))
        input.idx <- as.integer(input.idx)
        
        get.url=bioGRIDorganism.mitabPage.link[input.idx]
        full.dl.file=gsub("^.*/(BIOGRID-.*)$","\\1",get.url) #gmt.links[grepl("*\\_GO\\_AllPathways\\_with\\_GO\\_iea\\_.*\\.[Gg][Mm][Tt]",gmt.links)]
        
        BioGRIDtargetPath=getwd()  #gsub("\\/\\/","/", gsub("(.*\\/).*$","\\1",GMTdatabaseFile) )
        if(file.exists(file.path(BioGRIDtargetPath,full.dl.file))) {
          cat(paste0(" - Found that the full current .ZIP file online matches a file name you already have:\n  ",full.dl.file," [skipping download]\n"))
          BioGRIDzip=paste0(BioGRIDtargetPath,full.dl.file)
        } else {
          cat("Current Release .ZIP file online:  ",get.url,"\n")
          cat("Download this file to folder:  ",BioGRIDtargetPath,"\n")
          input.dlYN <- readline("[Y/n]?")
          if(input.dlYN == "Y" | input.dlYN == "y" | input.dlYN == "") {
            suppressPackageStartupMessages(require(curl,quietly=TRUE))
            if (!dir.exists(BioGRIDtargetPath)) dir.create(BioGRIDtargetPath)
            curr.dir<-getwd()
            setwd(BioGRIDtargetPath)
            cat("Downloading BioGRID .zip (multi)organism mitab file specified above...\n")
            curl_download(url=get.url, destfile=full.dl.file, quiet = TRUE, mode = "wb")
            setwd(curr.dir)
          }
        }
        cat("Unzipping downloaded file: ", paste0(BioGRIDtargetPath,full.dl.file),"...\n")
        BioGRIDzip=paste0(BioGRIDtargetPath,"/",full.dl.file)
        zipFilenames<-unzip(full.dl.file,list=TRUE)
        human.bioGRIDfile=zipFilenames$Name[which(grepl("BIOGRID-ORGANISM-Homo_sapiens-.*.mitab.txt",zipFilenames$Name))[1]]
        cat("..specifically, ",human.bioGRIDfile,"\n")
        if (!file.exists(human.bioGRIDfile)) unzip(full.dl.file,files=human.bioGRIDfile)
  
        cat("\nProcessing entrezIDs to gene symbols for binary PPI interactor list.\n")
        table<-read.delim(file=human.bioGRIDfile,header=TRUE,sep="\t")
        table.entrez<-apply(table[,1:2],2,function(x) { as.data.frame(do.call(rbind,strsplit(x,"[:]")))[,2] })
        colnames(table.entrez)<-c("From","To")
        
        suppressPackageStartupMessages(require(biomaRt,quietly=TRUE))
        human= useEnsembl("genes", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
        entrez.uniqueVec<-unique(c( table.entrez[which(!duplicated(table.entrez[,1])),1], table.entrez[which(!duplicated(table.entrez[,2])),2] ))
        #length(entrez.uniqueVec)
        # [1] 27886

        Ez2Symbol<-getBM(attributes=c("hgnc_symbol","entrezgene_id"), filters="entrezgene_id", values=entrez.uniqueVec, mart=human, uniqueRows=TRUE)
        Ez2Symbol$entrezgene_id<-as.character(Ez2Symbol$entrezgene_id)
  
        suppressPackageStartupMessages(require(doParallel,quietly=TRUE))
        if (!exists("parallelThreads")) { cat (" - parallelThreads not set. Trying with 4.\n"); parallelThreads=4; }
        clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
        registerDoParallel(clusterLocal)
  
        # reduced.table.Entrez<- t(parApply(cl=clusterLocal, table.entrez,1, sort, chunk.size=500))
        parallel::clusterExport(cl=clusterLocal, varlist=c("table.entrez"),envir=environment())
        parts <- splitIndices(nrow(table.entrez), length(clusterLocal))
        table.entrezParts <- lapply(parts, function(i) table.entrez[i,,drop=FALSE])
#        thisList<-parallel::clusterApply(cl=clusterLocal, table.entrezParts, function(x) { t(apply(x,1, sort)) })
#        reduced.table.Entrez<- do.call(rbind, thisList)
        reduced.table.Entrez <- foreach(i=1:length(table.entrezParts), .combine=rbind) %dopar% { t(apply(table.entrezParts[[i]],1,sort)) }

        #dim(reduced.table.Entrez)
        # [1] 1147670      2  #February v4.4.219 BioGRID
        reduced.table.Entrez<-unique(reduced.table.Entrez)
        #dim(reduced.table.Entrez)
        #[1] 862661      2  # as done
        #[1] 913581      2  #without pre-sort of columns 1 and 2 on every row.
  
        # Note: parApply call from within function does not parallelize.
        parallel::clusterExport(cl=clusterLocal, varlist=c("Ez2Symbol"),envir=environment())
        #table.HGNC<-parApply(cl=clusterLocal, reduced.table.Entrez, 2, function(x) { Ez2Symbol$hgnc_symbol[match(x,Ez2Symbol$entrezgene_id)] } )
        parts <- splitIndices(nrow(reduced.table.Entrez), length(clusterLocal))
        reduced.table.EntrezParts <- lapply(parts, function(i) reduced.table.Entrez[i,,drop=FALSE])
        #table.HGNC2 <- do.call(rbind, parallel::clusterApply(cl=clusterLocal, reduced.table.EntrezParts, function(x) { t(apply(x,1,function(x2) { Ez2Symbol$hgnc_symbol[match(x,Ez2Symbol$entrezgene_id)] })) }) )
        table.HGNC <- foreach(i=1:length(reduced.table.EntrezParts), .combine=rbind) %dopar% { apply(reduced.table.EntrezParts[[i]],2,function(x) Ez2Symbol$hgnc_symbol[match(x,Ez2Symbol$entrezgene_id)] ) }
        
#        parallel::clusterExport(cl=clusterLocal, varlist=c("table.HGNC"),envir=environment())
#        # Note: parApply call from within function does not parallelize.
##        removeRows<-unlist( parSapply(cl=clusterLocal, 1:nrow(table.HGNC), function(x) { if(length(which(is.na(table.HGNC[x,])))>0) x }, chunk.size=500) )
#        removeRows<-foreach(i=1:nrow(table.HGNC), .combine=c) %dopar% {
#          if(length(which(is.na(table.HGNC[i,])))>0) i
#        }
        keepIdx=which((as.numeric(is.na(table.HGNC[,1]))+as.numeric(is.na(table.HGNC[,2])))==0)
        
        #nrow(table.HGNC)-length(keepIdx)  # (remove rows count)
        #[1] 65816
  
        table.HGNC.noNA<-table.HGNC[keepIdx,]
        #dim(table.HGNC.noNA)
        # [1] 796845      2
  
        colnames(table.HGNC.noNA)<-c("From","To")
  
        myHumanBioGrid.tsvFile=gsub(".mitab.txt",".mitab.HUMANsimple.txt",human.bioGRIDfile)
        cat("Writing file for human binary protein-protein interactions (PPIs): ",myHumanBioGrid.tsvFile,"\n")
        write.table(table.HGNC.noNA, file=myHumanBioGrid.tsvFile,quote=FALSE, sep="\t", row.names=FALSE,col.names=TRUE)
      } else { stop(paste0("This is not an interactive session and required binary PPI file specified in variable myHumanBioGrid.tsvFile not found.\nThis file is created via processing after selecting a human organism-specific mitab.txt in a .ZIP file download from https://thebiogrid.org .\n\n")) }
    }
  }
  
  if (!exists("showAllPPIs")) showAllPPIs=FALSE
  if (!exists("boldEdgeColor")) boldEdgeColor="#483D8B66"   #darkslateblue, "66" makes the color transparent, about 30-40%; only used if PPIedges=TRUE
  if (!exists("netSpecies")) { cat(" - Variable species not specified. Assuming human symbols are in your network.\n   You can specify netSpecies='mouse' for conversion of BioGRID PPIs to mouse for PPI mapping.\n"); netSpecies="human"; }

  suppressPackageStartupMessages(require(igraph,quietly=TRUE))
  suppressPackageStartupMessages(require(RColorBrewer,quietly=TRUE))
  suppressPackageStartupMessages(require(WGCNA,quietly=TRUE))
  suppressPackageStartupMessages(require(Cairo,quietly=TRUE))
  if(PPIedges | showTOMhubs) suppressPackageStartupMessages(require(doParallel,quietly=TRUE))
  
  if(recalcMEs) {
    cat(" - Recalculating MEs from cleanDat and net$colors module assignments...\n")
    # (re)Calculate MEs (no grey)
    MEs<-tmpMEs<-data.frame()
    MEList = suppressMessages(moduleEigengenes(t(cleanDat), colors = net$colors))
    MEs = orderMEs(MEList$eigengenes)
    colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
    tmpMEs <- MEs #net$MEs
    colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
    MEs[,"grey"] <- NULL
    tmpMEs[,"MEgrey"] <- NULL
  } else {
    if(!exists("MEs")) stop(paste0("MEs not in memory but recalcMEs=FALSE was specified.\nPlease provide MEs eigenprotein or eigengene values as a data frame called MEs.\n\n"))
    colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
    tmpMEs <- MEs #net$MEs
    colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
    MEs[,"grey"] <- NULL
    tmpMEs[,"MEgrey"] <- NULL
  }  
  
  #(re)Establish M# for module colors table
  nModules<-length(table(net$colors))-1
  modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
  orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
  
  
  #load interactome [full bioGRID 4.4 (03/21/2023)]
  if(PPIedges) bioGrid <- read.table(file=myHumanBioGrid.tsvFile,sep="\t",header=TRUE)
  if (PPIedges==FALSE) { bioGrid<-data.frame(From=c(1),To=c(1)) }  # skip bolding edges for PPIs from BioGrid
  #------
  
  if (netSpecies=="mouse" & PPIedges) {
    cat(" - Starting biomaRt conversion of Symbols to Human for PPI mapping.\n")
    suppressPackageStartupMessages(require(biomaRt,quietly=TRUE))
    human = useEnsembl("genes", dataset="hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
    mouse = useEnsembl("genes", dataset="mmusculus_gene_ensembl",host="https://dec2021.archive.ensembl.org")
    
    ## If converting bioGrid from one species to another
    genelist.convert<-getLDS(attributes=c("hgnc_symbol"), filters="hgnc_symbol", values=unique(c(as.vector(bioGrid[,"From"]),as.vector(bioGrid[,"To"]))), mart=human, attributesL=c("mgi_symbol","external_gene_name"),martL = mouse)
  
    bioGrid.mouse.convert<-cbind(genelist.convert[match(bioGrid[,1],genelist.convert[,1]),3],genelist.convert[match(bioGrid[,2],genelist.convert[,1]),3])
    bioGrid.human<-bioGrid
    bioGrid<-bioGrid.mouse.convert
    #dim(bioGrid)
    bioGrid<-as.data.frame(bioGrid[-which(is.na(bioGrid)),])
    colnames(bioGrid)<-colnames(bioGrid.human)
    #dim(bioGrid)
  }
  
  
  ## Get sigend kME values
  kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
  # KME.vis=signedKME(cleanDat, kMEdat,corFnc="bicor",outputColumnName = "");
  KME.vis=kMEdat[,] # -ncol(kMEdat) #minus grey column if calculating signedKME on the fly
  annot=as.data.frame(do.call("rbind",strsplit(as.character(rownames(cleanDat)),"[|]")))
  # Rescue bad IDs with secondary ID (ENSGid or UniprotID)
  annot$V1[is.na(annot$V1)]<-annot$V2[is.na(annot$V1)]
  annot$V1[annot$V1==""]<-annot$V2[annot$V1==""]
  annot$V1[as.character(annot$V1)=="0"]<-annot$V2[as.character(annot$V1)=="0"]
  
  KME.vis$Symbol=annot$V1
  
  geneInfo.vis=as.data.frame(cbind(KME.vis$Symbol,net$colors, KME.vis))  
  
  geneInfo.vis=geneInfo.vis[,-ncol(geneInfo.vis)] # check if last column is Ensembl gene id
  colnames(geneInfo.vis)[1]= "Symbol"
  colnames(geneInfo.vis)[2]= "Module.Color"
  colnames(geneInfo.vis)<-gsub("kMEME","kME",colnames(geneInfo.vis))
  
  if(showTOMhubs) {
    cat(" - Starting adjacency calculation...[This could take some time and cause the R session to appear nonresponsive.]\n")
    if(exists("clusterLocal")) tryCatch(suppressWarnings(stopCluster(clusterLocal)),error=function(err) { })
    clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
    registerDoParallel(clusterLocal)
    enableWGCNAThreads(nThreads=parallelThreads)
    # Allows us to get the top connected genes in the module
    softPower = power
    adjacency = adjacency(t(cleanDat), power=softPower, type="signed",corFnc="bicor")
    cat(" - Starting TOM calculation...[This could take some time, as well.]\n")
    TOM = TOMsimilarity(adjacency)
  
    TOM.matrix = as.matrix(TOM);
  }
  
  uniquemodcolors = gsub("kME","",gsub("kMEME","",colnames(kMEdat[,]))); #-ncol(kMEdat) last column if KMEdat.vis was calculated on the fly... moduleColors
  # (OR SELECT MODULES INSTEAD OF ALL)
  # uniquemodcolors = c("violet")
  
  cat(paste0(" - Beginning iGraph builds ",if(showTOMhubs) { "with" } else { "without" }," TOM hubs broken out by top coexpression edges.\n   ",length(uniquemodcolors)," module iGraphs to be built ",if(PPIedges) { paste0("with ",nrow(bioGrid)) } else { "without" }," BioGRID PPI edges to map.\n\n"))
  
  if(PPIedges) {
    # require(doParallel,quietly=TRUE)
    if(exists("clusterLocal")) tryCatch(suppressWarnings(stopCluster(clusterLocal)),error=function(err) { })
    clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
    registerDoParallel(clusterLocal)
    clusterExport(cl=clusterLocal, varlist=c("bioGrid"),envir=environment())
  }
  
  listboldEdges<-list()
  firstPDFdone=FALSE
  
  
  for (CAIRO in c(TRUE,FALSE)) {
    if (CAIRO) {
      CairoPDF(file=paste0("./",outFilePrefix,"iGraph_Modules-",FileBaseName,"-CAIRO-",if(!recalcMEs) { "noMErecalc" }, if(showTOMhubs) {"withTOMhubs"}, ".pdf"),width= if(showTOMhubs) {32} else {16}, height=12)
      } else {
        cat("\n - Replotting above iGraphs in non-Cairo version of PDF output...\n")
        pdf(paste0("./",outFilePrefix,"iGraph_Modules-",FileBaseName,"-nonCAIRO-",if(!recalcMEs) { "noMErecalc" }, if(showTOMhubs) {"withTOMhubs"}, ".pdf"),height=9,width= if(showTOMhubs) {20} else {10})
      }
    
    if(showTOMhubs) par(mfrow=c(1,2))
    
    # for (i in 1:length(sigmodcolors))  {
      # mod=sigmodcolors[i];  
      # numgenesingraph = 50;
      # numconnections2keep = 1500;
    
    for (mod in uniquemodcolors)  {
      #mod="darkslateblue"
      numgenesingraph = 100;
      numconnections2keep = keepTopEdges;
      if(exists("keepTopPercentEdges")) percentconnections2keep=keepTopPercentEdges;
      # cat('module:',mod,'\n');
      geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="NA",]
      geneInfo.vis=geneInfo.vis[geneInfo.vis$Symbol!="",]
    
      colind = which(colnames(geneInfo.vis)== paste("kME",mod,sep=""));
      rowind = which(geneInfo.vis[,2]==mod);
      # cat(' ',length(rowind),'probes in module\n');
      submatrix = geneInfo.vis[rowind,];
      orderind = order(submatrix[,colind],decreasing=TRUE);
      if (length(rowind) < numgenesingraph) {
        numgenesingraph = length(rowind);
        numconnections2keep = numgenesingraph/2 * (numgenesingraph/6 - 1); #added /2 and /6 9/14/2015
      }
      if (length(rowind)<10) { innercircleNum=length(rowind) } else { innercircleNum=10 }
      cat(paste0(orderedModules[which(orderedModules[,"Color"]==mod),"Mnum"]," ",mod," module"),'- Calculating network graph(s), using top',paste0('[',numgenesingraph,'/',length(rowind),']'),'probes & showing hub connections using the top',round(numconnections2keep,0),if(showTOMhubs) {'coexpression edges of TOM.\n'} else {'pseudoedges (showTOMhubs=/=TRUE).'});
      submatrix = submatrix[orderind[1:numgenesingraph],];
      #Identify the columns in the TOM that correspond to these hub probes
      matchind = match(submatrix$Symbol,KME.vis$Symbol)  #annot$Symbol);   #*** KME.vis$Symbol would be appropriate second term to match to, but only inner hub edges drawn then
      
    #x*  orderind = order(reducedTOM,decreasing=TRUE);
    #x*  connections2keep = orderind[1:numconnections2keep];
    #x*  reducedTOM = matrix(0,nrow(reducedTOM),ncol(reducedTOM));
    #x*  reducedTOM[connections2keep] = 1;
    
      if(showTOMhubs) {
        reducedTOM = TOM.matrix[matchind,matchind];
        # keep up to numconnections2keep
        reducedTOM.diag0<-reducedTOM
        diag(reducedTOM.diag0)<-0
        if(!exists("keepTopPercentEdges")) {
          # keep up to numconnections2keep
          orderedTOMvals = sort(reducedTOM.diag0,decreasing=TRUE);
          reducedTOM[which(reducedTOM<orderedTOMvals[numconnections2keep])]<-0
        } else {
          # keep top percentconnections2keep
          connectionsNotZero<-reducedTOM[reducedTOM!=0]
          reducedTOM[reducedTOM<quantile(connectionsNotZero,(1-percentconnections2keep/100))]<-0
        }
      } else { # adjacency and TOM calculations were not performed -- draw generic connections
        reducedTOM = matrix(0,length(matchind),length(matchind))
        connections2keep.idx = 1:numconnections2keep
        reducedTOM[connections2keep.idx] = 1;
      }
      
      g0 <- graph.adjacency(as.matrix(reducedTOM[1:innercircleNum,1:innercircleNum]),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMata <- layout.circle(g0)
    
      if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) {   
        g0 <- graph.adjacency(as.matrix(reducedTOM[11:ncol(reducedTOM),11:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMatb <- layout.circle(g0)
    
        g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
        layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75)
      } else {
        if (ncol(reducedTOM) > 10) { #****
    
          g0 <- graph.adjacency(as.matrix(reducedTOM[11:50,11:50]),mode="undirected",weighted=TRUE,diag=FALSE)
          layoutMatb <- layout.circle(g0)
        
          g0 <- graph.adjacency(as.matrix(reducedTOM[51:ncol(reducedTOM),51:ncol(reducedTOM)]),mode="undirected",weighted=TRUE,diag=FALSE)
          layoutMatc <- layout.circle(g0)
          g1 <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
          layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc)
        } else {
          g1 <- g0
          layoutMat <- layoutMata
        } #****
      }
     
    #PLOT DONE BELOW WITH ADDED PHYSICAL INTERACTION EDGES (OR USE THIS SINGLE LINE FOR NO PHYSICAL INT HIGHLIGHTING)
    #par(mfrow=c(1,2))
    #  plot(g1,edge.color="grey",vertex.color=mod,vertex.label=as.character(rbind(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=submatrix[,colind]^2*vertexsize,main=paste(mod,"module"))
    
     if(showTOMhubs) {
       # prepare second plot on same page, with option enabled
       g.alt <- graph.adjacency(as.matrix(reducedTOM),mode="undirected",weighted=TRUE,diag=FALSE)
       layoutMat.alt <- layout.fruchterman.reingold(g.alt)
       isolated = which(degree(g.alt)==0)
       g.alt = delete.vertices(g.alt, isolated)
       layoutMat.alt = layoutMat.alt[-isolated,]
       labels.alt=as.character(rbind(submatrix$Symbol[-isolated]))
       degrees.alt=degree(g.alt)
  
     ## Plot is called below, after first triple-circular layout plot
     #  plot(g.alt,edge.color="grey", vertex.color=mod, vertex.label=sapply(1:length(labels.alt),function(x) as.expression(bquote(.(labels.alt[x])^.(degrees.alt[x]) )) ),
     #       vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat.alt,vertex.size=submatrix[,colind]^2*vertexsize,
     #       main=as.expression(bquote(.(paste0(mod," module hubs connected by top ", sum(degree(g.alt)), " TOM edges: HUB"))^.("degree"))) )
     }
    
    
    #Add symbols2Add
    iter=as.numeric(0)
    
    if (length(symbols2Add)==4) {
      for (symbol2Add in symbols2Add) {
        iter=iter+1
        if (iter==1) { position=cbind(-0.7071068,-0.7071068) } 
        if (iter==2) { position=rbind(position,cbind(0.7071068,-0.7071068)) }
        if (iter==3) { position=rbind(position,cbind(-0.7071068,0.7071068)) } 
        if (iter==4) { position=rbind(position,cbind(0.7071068,0.7071068)) }
        # ADD APP to graph g3 (to be used if APP not already in graph)
        g2 <- vertex(1) #,size=40,color="green", label=symbol2Add
        # layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,cbind(-0.7071068,-0.7071068)))
        g3 <- g1+g2
        V(g3)$name <- c(1:(numgenesingraph+1))
        #WORKS# plot(g3,edge.color="grey",vertex.color=c(rep(mod,nrow(layoutMat)-1),"green"),vertex.label=as.character(c(submatrix$Symbol,symbol2Add)),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat,vertex.size=c(submatrix[,colind]^2*vertexsize,15),main=paste(mod,"module"))
    
        g1<-g3
        numgenesingraph=numgenesingraph+1
      }
    }
    
    #x***moved out of loop
    if (length(symbols2Add)==4) {  if (ncol(reducedTOM) < 51) { layoutMat <- rbind(layoutMata*0.25,rbind(layoutMatb*0.75,position)) } else { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, rbind(layoutMatc,position*1.33)) }
                                 } else { if (ncol(reducedTOM) < 51 & ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75) } else { if(ncol(reducedTOM) > 10) { layoutMat <- rbind(layoutMata*0.25,layoutMatb*0.75, layoutMatc) } } #do not add 4 node positions if symbols2Add is empty (length 0)
                                         }
    symbolList<-c(as.character(submatrix$Symbol),symbols2Add)
    if (length(symbols2Add)==4) { vertexColorVec<-c(rep(mod,numgenesingraph-4),"green","darkgreen","steelblue","darkslateblue"); vertexSizeMat<-c(submatrix[,colind]^2*vertexsize,rep(15,4)); } else { vertexColorVec<-rep(mod,numgenesingraph); vertexSizeMat<-submatrix[,colind]^2*vertexsize; }
    
    
    
    ## FIND EDGES OVERLAPPING WITH BIOGRID PAIRS
    if(!firstPDFdone) listboldEdges[[mod]]<-matrix(ncol=2,nrow=0)
    
    noPPIfound.thisMod=FALSE # (re)set this flag

    if(nrow(bioGrid)>1) {
      ## non-parallel (slow) version of below collection of PPI (bold) edges with nested parSapply, sapply...
      #for(i in 1:numgenesingraph) {
      # for(j in i:numgenesingraph) {
      #  if(!(length(which(bioGrid$From==symbolList[i]&bioGrid$To==symbolList[j]))+length(which(bioGrid$From==symbolList[j]&bioGrid$To==symbolList[i])))==0) { listboldEdges <- rbind(listboldEdges,c(i,j)) }
      # } }
    
    
      clusterExport(cl=clusterLocal, varlist=c("numgenesingraph","symbolList"),envir=environment())
      if(!firstPDFdone) {
        #  listboldEdges[[mod]] <- do.call(rbind, lapply(as.list(
        #                 parSapply(cl=clusterLocal, 1:numgenesingraph, function(i) { edges<-sapply(i:numgenesingraph,function(j) {
        #                                                                                        if(!(length(which(bioGrid$From==symbolList[i]&bioGrid$To==symbolList[j]))+length(which(bioGrid$From==symbolList[j]&bioGrid$To==symbolList[i])))==0) { c(i,j) }  }
        #                                                                                          )
        #                                                                             if(is.list(edges)) { do.call(rbind, edges) } else { do.call(rbind, list(edges)) }  },
        #                           chunk.size=2)  ), function(x) if(!is.null(x)) { if(ncol(x)==1) { t(x) } else { x } })  # this lapply transposes self-edges that may appear as 2-row, 1 column data inside list elements, so that the outer rbind always works.
        #                                                 )
        ## Parallelizes within function                                                 
        if(!firstPDFdone) holderList <- foreach(i = 1:numgenesingraph) %:% foreach(j=i:numgenesingraph) %dopar% {
          edges <- if(!(length(which(bioGrid$From==symbolList[i]&bioGrid$To==symbolList[j]))+length(which(bioGrid$From==symbolList[j]&bioGrid$To==symbolList[i])))==0) { c(i,j) }
          if(is.list(edges)) { do.call(rbind, edges) } else { do.call(rbind, list(edges)) }
        }
        holderList2<-list()
        for (x2 in 1:length(holderList)) holderList2[[x2]] <- do.call(rbind, lapply(holderList[[x2]], function(x3) { if(!is.null(x3)) { if(ncol(x3)==1) { t(x3) } else { x3 }; }}))
        listboldEdges[[mod]] <- do.call(rbind, lapply(holderList2, function(x3) { if(!is.null(x3)) { if(ncol(x3)==1) { t(x3) } else { x3 }; }}))
      }
      
      
      ##Remove self-loops
      if (PPIedges) {
        if (length(nrow(listboldEdges[[mod]]))>0) {
          if (nrow(listboldEdges[[mod]])>0) { 
            for(i in nrow(listboldEdges[[mod]]):1) {
              if (i<=nrow(listboldEdges[[mod]])) {
                if(listboldEdges[[mod]][i,1]==listboldEdges[[mod]][i,2]) {  #identify self-loop
                  listboldEdges[[mod]] <-listboldEdges[[mod]][-i,]
                  #if (!length(nrow(listboldEdges[[mod]]))>0) { 
                  #  cat('NO PPI FOUND FOR THIS MODULE.\n')
                  #  break #HANDLE NO BOLD EDGES CASE
                  #}
                }
              }
            }
          }
        } else {
          cat('NO PPI FOUND FOR THIS MODULE.\n')
          noPPIfound.thisMod=TRUE  # HANDLE NO BOLD EDGES CASE
        }
      }
    }
    
    if (is.vector(listboldEdges[[mod]])) { listboldEdges[[mod]]<-t(data.frame(listboldEdges[[mod]])) }
    
    newgraph <- g1 %>%
    #  delete_edges(listboldEdges[[mod]]) %>%
      set_edge_attr("color",value="lightgrey") %>%
      set_edge_attr("width",value=1) %>%
      set_edge_attr("curved",value=0) 
    
    
    if(!noPPIfound.thisMod) {
      if (length(symbols2Add)==4) {
        if (is.vector(listboldEdges[[mod]])) { 
          cat('ONLY 1 SINGLE INTERACTION FOUND.\n')
          newgraph <- newgraph %>%
          add_edges(c(listboldEdges[[mod]][1],listboldEdges[[mod]][2]),color="steelblue",width=2,curved=0)
        } else {
          if (dim(listboldEdges[[mod]])[1]>0) { 
            for(k in 1:nrow(listboldEdges[[mod]])) {
              if((length(which(submatrix$Symbol==symbols2Add[1]))==0)&(listboldEdges[[mod]][k,1]==(numgenesingraph-3)|listboldEdges[[mod]][k,2]==(numgenesingraph-3))) { 
                newgraph <- newgraph %>%
                add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color="#33BB33",width=3,curved=0)
              }
            }
            for(k in 1:nrow(listboldEdges[[mod]])) {
              if((length(which(submatrix$Symbol==symbols2Add[2]))==0)&(listboldEdges[[mod]][k,1]==(numgenesingraph-2)|listboldEdges[[mod]][k,2]==(numgenesingraph-2))) { 
                newgraph <- newgraph %>%
                add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color="#338833",width=3,curved=0)
              }
            }
            for(k in 1:nrow(listboldEdges[[mod]])) {
              if((length(which(submatrix$Symbol==symbols2Add[3]))==0)&(listboldEdges[[mod]][k,1]==(numgenesingraph-1)|listboldEdges[[mod]][k,2]==(numgenesingraph-1))) { 
                newgraph <- newgraph %>%
                add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color="steelblue",width=3,curved=0)
              }
            }
            for(k in 1:nrow(listboldEdges[[mod]])) {
              if((length(which(submatrix$Symbol==symbols2Add[4]))==0)&(listboldEdges[[mod]][k,1]==numgenesingraph|listboldEdges[[mod]][k,2]==numgenesingraph)) { 
                newgraph <- newgraph %>%
                add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color="darkslateblue",width=3,curved=0)
              }
            }
      
            if (showAllPPIs) { 
              for(k in 1:nrow(listboldEdges[[mod]])) {
                if(length(which(listboldEdges[[mod]][k,] %in% c(numgenesingraph-3,numgenesingraph-2,numgenesingraph-1,numgenesingraph)))==0) {
                  newgraph <- newgraph %>%
                  add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color=boldEdgeColor,width=3,curved=0)
                }
              }
            }
          }
        }
      } else {
        if (dim(listboldEdges[[mod]])[1]>0) { 
          for(k in 1:nrow(listboldEdges[[mod]])) {
            newgraph <- newgraph %>%
            add_edges(c(listboldEdges[[mod]][k,1],listboldEdges[[mod]][k,2]),color=boldEdgeColor,width=3,curved=0)
          }
        }
      }
    }
    
     ## make dark nodes contrast vs. their black text by adding transparency (checked luminance of colors 1:100)
     if(mod=="black") { vertexColorVec[vertexColorVec=="black"]<- "#444444BB" }
     if(mod=="midnightblue") { vertexColorVec[vertexColorVec=="midnightblue"]<- "#191970BB" }
     if(mod=="darkmagenta") { vertexColorVec[vertexColorVec=="darkmagenta"]<- "#8B008BBB" }
     if(mod=="brown4") { vertexColorVec[vertexColorVec=="brown4"]<- "#8B2323BB" }
     if(mod=="magenta4") { vertexColorVec[vertexColorVec=="magenta4"]<- "#8B008BBB" }
    
     highlightNodes<-match(GOIlist,symbolList) #ADmagmaList$Gene.Symbol,symbolList
     highlightNodes<-sort(na.omit(highlightNodes))
     if(length(highlightNodes)>0) { if(mod=="yellow") {vertexColorVec[highlightNodes]<-"cyan" } else {vertexColorVec[highlightNodes]<-"yellow" } }
    
     plot(newgraph,vertex.color=vertexColorVec,vertex.label=as.character(symbolList),vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout=layoutMat,vertex.size=vertexSizeMat,main=paste0(orderedModules[which(orderedModules[,"Color"]==mod),"Mnum"]," ",mod," module"))
    
     if(showTOMhubs) {
       if(ncol(reducedTOM)>10) {
         vertexColorVec.alt<-vertexColorVec[1:(length(vertexColorVec)-length(symbols2Add))]
         if(length(isolated)>0) vertexColorVec.alt<-vertexColorVec.alt[-isolated]
         plot(g.alt,edge.color="grey", vertex.color=vertexColorVec.alt, vertex.label=sapply(1:length(labels.alt),function(x) as.expression(bquote(.(labels.alt[x])^.(degrees.alt[x]) )) ),
              vertex.label.cex=1.1,vertex.label.dist=0.45,vertex.label.degree=-pi/4,vertex.label.color="black",layout= layoutMat.alt,vertex.size=submatrix[,colind]^2*vertexsize, font.main=2,
              main=as.expression(bquote(.(paste0(orderedModules[which(orderedModules[,"Color"]==mod),"Mnum"]," ",mod," module hubs connected by top ", sum(degree(g.alt)), " TOM edges: HUB"))^.("degree"))) )
       } else {
         plot.new()
         mtext("Hub connectivity plot skipped for module with < 11 nodes.",3)
       } 
     }
    }
    dev.off();
    firstPDFdone=TRUE
  
  } # for (CAIRO in c(TRUE,FALSE)) loop... <repeat>
} # close buildIgraphs()
