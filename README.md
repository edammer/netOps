# netOps
Network Operations for Systems Biology

<b>buildIgraphs</b> - Build network graphs for your modules from the Seyfried Systems Biology pipeline.

Additional network visualization routine functions will be added to this repository later.

Simple Wrapper for buildIgraphs():
(All variables set in global environment, but defaults are tried if none are specified).
```
unzip("sampleInput.zip")  # Download from this repo; contains net.csv and cleanDat.csv as loaded below, for sample input

rootdir="e:/OneDrive/SystemsBioPipeline/"
setwd(rootdir)

cleanDat=read.csv("cleanDat.csv",header=TRUE,row.names=1)
net<-as.list(read.csv("net.csv",header=TRUE))


outFileSuffix="MEGATMT_NatNeuro2022-44mods+PPIs"
parallelThreads=31       # needed if PPIedges=TRUE; set to # of threads on your computer
PPIedges=TRUE            # total time to map ~800000 bioGRID binary PPIs to the 44 modules is ~40 minutes.
                         # use FALSE for a faster output.
showTOMhubs=TRUE         # calculates TOM from cleanDat, so that hubs can be found and shown on a second graph for each module.
                         # the second graph uses Fruchterman-Reingold layout and only nodes with a TOM value calculated by WGCNA
                         # that is equal to or greater than the Nth top such edge. Default N=150. See full parameter list to change.

source("buildIgraphs.R")
buildIgraphs()
```
## Other Notes
Requires R packages: igraph, curl, rvest, biomaRt, doparallel, RColorBrewer, WGCNA, Cairo.

BioGRID PPIs will be downloaded and processed to a simple 2-column human gene SYMBOL format tab-seperated value .txt file from the current BioGRID mitab-formatted release on <a href="https://thebiogrid.org/">https://thebiogrid.org/</a>.

Full parameter list for buildIgraphs():
```
############################################################################################
# iGRAPHs (Multiple Toggle Options, e.g. BioGRID interactome overlap) // CONNECTIVITY PLOT #
############################################################################################
showTOMhubs=FALSE        # default is FALSE if not set. Calculates adjacency and TOM for exact determination of top coexpression-based network edges;
                         # a second plot on each page will show these hubs in isolation.
power=8                  # network power beta for adjacency calculation; this is needed if showTOMhubs=TRUE
keepTopEdges=150         # default is 150 if not set.  This many edges are plotted as grey connections between nodes
keepTopPercentEdges=NULL # 5 to 20, typically, for pulling hubs together in second same-page iGraphs with Fruchterman-Reingold layout; overrides keepTopEdges setting

parallelThreads=4        # needed if PPIedges=TRUE; set to # of threads (-1) on your computer
#####################
PPIedges=TRUE            # Draw bold edges for PPIs in your network modules. TRUE will be somewhat slower, and if few parallelThreads specified below...
myHumanBioGrid.tsvFile = "nonexistent.file"
                         # "BIOGRID-ORGANISM-Homo_sapiens-4.4.219.mitab.HUMANsimple.txt" #v4.4.219 downloaded 3/21/2023; ONLY NEEDED IF PPIedges=TRUE
species="human"          # current option "mouse" will convert bioGRID to mouse symbols before drawing PPI edges in your network.
boldEdgeColor=paste0(gplots::col2hex("darkslateblue"),"66")
                         #"66" makes the color transparent, about 30-40%; only used if PPIedges=TRUE
symbols2Add=c()          # If non-blank, must specify exactly 4 gene symbols. Interactions to these extra nodes (not necessarily in the modules)
                         # will be drawn for every module. These nodes appear in the 4 corners of each plot.
showAllPPIs=FALSE        # if there are 4 symbols2Add, unless this flag is true, only PPIs to the 4 corner nodes will be drawn.
recalcMEs=TRUE           # Default. Usually TRUE; unless MEs are from another cleanDat and stored in variable 'MEs'.
#####################
GOIlist<- c()            # Genes of interest to highlight-could be all RNA binding proteins, or a list of MAGMA-significant genes for a disease...
                         # ...they get highlighted yellow, or cyan if the module is yellow.
vertexsize=16            # 8 for regular, 16 for large balls
#####################
outFilePrefix="4",         # typically "4", or step # in pipeline being run
outFileSuffix=FileBaseName # A description of the project, used as a filename suffix
############################################################################################
```
