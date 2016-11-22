# ========== Cab_seq_tree Version 1.0 ============ #
# Anthony Bowen
# anthony.bowen@med.einstein.yu.edu
# Laboratory of Arturo Casadevall
# Albert Einstein College of Medicine, Bronx NY
# Last Modified 10/17/16; updated to work with IgMotif output results and previously analyzed PDB database

# This script will show positional seq variation in antibody sequences

### Setup
initial.dir<-getwd() # store the current directory
temp.dir<-  #"B:/Tony/Documents/R Working Directory/Temp" #2500K
            file.path("~","R Working Directory","Temp") #Macbook Pro
pdb.dir<-   #"B:/Tony/Documents/Casadevall Lab/Non-redundant PDB files" #2500K
            file.path("~","Documents","Casadevall Lab","Non-redundant PDB files") #Macbook Pro
src.dir<-   #"B:/Tony/Google Drive/R Programming/IgMotif" #2500K
            file.path("~","Google Drive","R Programming","IgMotif") #Macbook Pro

setwd(src.dir) # change to the src directory
options(scipen=10) # number of digits before scientific notation is used

#load needed packages
library("data.table")
library("ggplot2")
library("xlsx")
library("seqinr")
library("bio3d")
library("Biostrings")
library("dendextend")
library("circlize")
library("RColorBrewer")

### Change variables
imgt_append<- FALSE #whether to include imgt list of antibodies
paired<- FALSE # whether to remove unpaired chains from the analysis (1 H and 1 L per PDB ID)
refscore.min<- -350 # minimum ref gene alignment score for counting H and L chains, not used
write.fasta<- FALSE # whether to write new fasta files of raw sequences
max.position<-120 # for alignment figures
pdb_list_mode<-"antibodies" # "models", "antibodies", "cabs", "hydrolases"
mean_temp_mode<-"antibodies" #"models", "antibodies", "hydrolases", "serprot", "cysprot", "thrprot", "metprot", "aspprot", "glycosylase"
rms_method<-"rmsd or DRMS" # "rmsd", "DRMS", "rmsd or DRMS", "rmsd and DRMS" for distance measure cutoff
rms_cutoff<-1 #rms cutoff for determining matches
meanT<- #choose one mean template from below by uncommenting
      "4hdi-1L-98L-26L" #3E5 G3 1
      #"4hdi-1L-98L-28L" #3E5 G3 2
      #"2h1p-1L-98L-26L" #2H1 G1 1
      #"2h1p-1L-98L-28L" #2H1 G1 2
      #"18b7_G1_2H1Pmodel-1L-98L-26L" #18B7 G1 model 1
      #"18b7_G1_2H1Pmodel-1L-98L-28L" #18B7 G1 model 2
      #"18b7_G3_4HDImodel-1L-98L-26L" #18B7 G3 model 1
      #"18b7_G3_4HDImodel-1L-98L-28L" #18B7 G3 model 2
      #"1utn-102A-57A-195A" #1UTN- bovine trypsin: A- H57 S195 D102
      #"1st2-32A-64A-221A" #1ST2- bacillus subtilisin: A- H64 S221 D32
      #"1bqi-175A-159A-25A" #1BQI- papain triad 1: A- N175 H159 C25, must allow Asn
      #"1bqi-158A-159A-25A" #91BQI- papain triad 2: A- H159 C25 D158
      #"1lay-157A-63A-132A" #1LAY- CMV protease triad: A- H157 H63 S132

#=============================BEGIN Define custom functions=============================#
circos_Ldend<-function(dend,seq_dt,color_codes,filename="Ldend.pdf",width=12,height=12){
      dend<-dend
      seq_dt<-seq_dt
      dend<-dendextend::set(dend,"labels",seq_dt[chain=="L",name][order.dendrogram(dend)])
      dend<-color_branches(dend, k=4, col=brewer.pal(4,"Set1"))
      circle_max<-length(labels(dend))
      pdf(file=file.path("figs",filename), width=width,height=height)
      par(mar = c(1, 1, 1, 1))
      
      circos.par(track.height=0.08,track.margin=c(0.001,0),cell.padding=c(0,1,0,1))
      circos.initialize("foo",xlim=c(0,circle_max))
      
      
      ### plot color bar tracks
      # cab track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="L",cab])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      # cab text labels
      offset<-1.25*circle_max/64
      circos.rect(9*circle_max/16+offset,0.2,8.2*circle_max/16-offset,0.8, track.index = 1,
                  col=rgb(192/255,192/255,192/255,0.8), border = NA)
      circos.text(9*circle_max/16,.5,track.index=1,"Yes Catalytic",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["TRUE"], cex=1)
      circos.text(8.2*circle_max/16,.5,track.index=1,"? Catalytic",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["FALSE"], cex=1)
      
      # match track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="L",match])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      # match track labels
      offset<-1.25*circle_max/64
      circos.rect(9*circle_max/16+offset,0.2,8.2*circle_max/16-offset,0.8, track.index = 2,
                  col=rgb(192/255,192/255,192/255,0.8), border = NA)
      circos.text(9*circle_max/16,.5,track.index=2,"Matching",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["TRUE"], cex=1)
      circos.text(8.2*circle_max/16,.5,track.index=2,"Not Matching",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["FALSE"], cex=1)
      
      #family track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="L",family])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      #family track labels
      offset<-1*circle_max/64
      circos.rect(11.5*circle_max/16+offset,0.2,11*circle_max/16-offset,0.8, track.index = 3,
                  col=rgb(216/255,216/255,216/255,0.8), border = NA)
      circos.text((11.5*circle_max/16),.5,track.index=3,"IgKV",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["IgKV"], cex=1)
      circos.text((11*circle_max)/16,.5,track.index=3,"IgLV",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["IgLV"], cex=1)
      
      # species track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="L",species])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      #species track labels
      offset<-1*circle_max/64
      circos.rect(13.5*circle_max/16+offset,0.2,10.5*circle_max/16-offset,0.8, track.index = 4,
                  col=rgb(1,1,1,0.8), border = NA)
      circos.text((13.5*circle_max/16),.5,track.index=4,"Human",facing="bending.inside", niceFacing=FALSE,
                  col=color_codes["human"], cex=1)
      circos.text((13*circle_max)/16,.5,track.index=4,"Mouse",facing="bending.inside", niceFacing=FALSE,
                  col=color_codes["mouse"], cex=1)
      circos.text((12.5*circle_max)/16,.5,track.index=4,"Macaque",facing="bending.inside", niceFacing=FALSE,
                  col=color_codes["rhesus monkey"], cex=1)
      circos.text((12*circle_max)/16,.5,track.index=4,"Rat",facing="bending.inside", niceFacing=FALSE,
                  col=color_codes["rat"], cex=1)
      circos.text((11.5*circle_max)/16,.5,track.index=4,"Rabbit",facing="bending.inside", niceFacing=FALSE,
                  col=color_codes["rabbit"], cex=1)
      circos.text((11*circle_max)/16,.5,track.index=4,"Alpaca",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["alpaca"], cex=1)
      circos.text((10.5*circle_max)/16,.5,track.index=4,"Pig",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["pig"], cex=1)
      
      
      # plot labels track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.text(1:circle_max-0.5, rep(0, circle_max), seq_dt[name %in% labels(dend),refgene],
                        #col = color_codes[as.character(seq_dt[chain=="L",refgene])[order.dendrogram(dend)]],
                        col = "gray75",
                        cex=0.1, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      }, bg.border = NA, track.height = 0.022, track.margin=c(0,0.002))
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.text(1:circle_max-0.5, rep(0, circle_max), labels(dend), col = "black", cex=0.1,
                        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      }, bg.border = NA, track.height = 0.01, track.margin=c(0.002,0))
      
      # plot dendrogram
      max_height = attr(dend, "height")
      circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
            circos.dendrogram(dend, max_height = max_height)
      }, track.height = 0.5, bg.border = NA)
      
      circos.par(RESET=TRUE)
      circos.clear()
      dev.off()
}
circos_Hdend<-function(dend,seq_dt,color_codes,filename="Hdend.pdf",width=12,height=12){
      dend<-dend
      seq_dt<-seq_dt
      dend<-dendextend::set(dend,"labels",seq_dt[chain=="H",name][order.dendrogram(dend)])
      dend<-color_branches(dend, k=4, col=brewer.pal(4,"Set1"))
      circle_max<-length(labels(dend))
      pdf(file=file.path("figs",filename), width=width,height=height)
      par(mar = c(1, 1, 1, 1))
      
      circos.par(track.height=0.08,track.margin=c(0.001,0),cell.padding=c(0,1,0,1))
      circos.initialize("foo",xlim=c(0,circle_max))
      
      ### plot color bar tracks
      # cab track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="H",cab])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      # cab text labels
      offset<-1.25*circle_max/64
      circos.rect(13.5*circle_max/16+offset,0.2,12.75*circle_max/16-offset,0.8, track.index = 1,
                  col=rgb(192/255,192/255,192/255,0.8), border = NA)
      circos.text(13.5*circle_max/16,.5,track.index=1,"Yes Catalytic",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["TRUE"], cex=1)
      circos.text(12.75*circle_max/16,.5,track.index=1,"? Catalytic",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["FALSE"], cex=1)
      
      # match track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="H",match])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      # match track labels
      offset<-1.25*circle_max/64
      circos.rect(13.5*circle_max/16+offset,0.2,12.75*circle_max/16-offset,0.8, track.index = 2,
                  col=rgb(192/255,192/255,192/255,0.8), border = NA)
      circos.text(13.5*circle_max/16,.5,track.index=2,"Matching",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["TRUE"], cex=1)
      circos.text(12.75*circle_max/16,.5,track.index=2,"Not Matching",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["FALSE"], cex=1)
      
      # species track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.rect(1:circle_max-1, rep(0, circle_max), 1:circle_max, rep(1,circle_max), 
                        col = color_codes[as.character(seq_dt[chain=="H",species])[order.dendrogram(dend)]],
                        border = NA)
      }, bg.border = NA)
      
      #species track labels
      offset<-1*circle_max/64
      circos.rect(4*circle_max/16+offset,0.2,1*circle_max/16-offset,0.8, track.index = 3,
                  col=rgb(1,1,1,0.8), border = NA)
      circos.text((3.5*circle_max/16),.5,track.index=3,"Human",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["human"], cex=1)
      circos.text((4*circle_max)/16,.5,track.index=3,"Mouse",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["mouse"], cex=1)
      circos.text((3*circle_max)/16,.5,track.index=3,"Macaque",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["rhesus monkey"], cex=1)
      circos.text((2.5*circle_max)/16,.5,track.index=3,"Rat",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["rat"], cex=1)
      circos.text((2*circle_max)/16,.5,track.index=3,"Rabbit",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["rabbit"], cex=1)
      circos.text((1.5*circle_max)/16,.5,track.index=3,"Alpaca",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["alpaca"], cex=1)
      circos.text((1*circle_max)/16,.5,track.index=3,"Pig",facing="bending.inside", niceFacing=TRUE,
                  col=color_codes["pig"], cex=1)
      
      # plot labels track
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.text(1:circle_max-0.5, rep(0, circle_max), seq_dt[name %in% labels(dend),refgene],
                        #col = color_codes[as.character(seq_dt[chain=="L",refgene])[order.dendrogram(dend)]],
                        col = "gray75",
                        cex=0.1, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      }, bg.border = NA, track.height = 0.022, track.margin=c(0,0.002))
      circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
            circos.text(1:circle_max-0.5, rep(0, circle_max), labels(dend), col = "black", cex=0.1,
                        facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
      }, bg.border = NA, track.height = 0.01, track.margin=c(0.002,0))
      
      # plot dendrogram
      max_height = attr(dend, "height")
      circos.track(ylim = c(0, max_height), panel.fun = function(x, y) {
            circos.dendrogram(dend, max_height = max_height)
      }, track.height = 0.6, bg.border = NA)
      
      circos.par(RESET=TRUE)
      circos.clear()
      dev.off()
}
#==============================Begin Script====================================#

#Read in IMGT reference sequences and create table on consensus sequences
imgt_files<-list.files(file.path(src.dir,"IMGT Gene-DB AA V gaps"),pattern=".fasta")
imgt_con_dt<-data.table(file=imgt_files) #consensus sequence data table
imgt_consq<-character()
imgt_species<-character()
imgt_family<-character()
imgt_allseq_list<-list() #total sequence list
for (file in imgt_files){
      species<-strsplit(gsub(".fasta","",file),"_")[[1]][2] #mouse, human, rat, etc...
      family<-strsplit(gsub(".fasta","",file),"_")[[1]][1] #IgHV, IgKV, IgLV
      sequences<-bio3d::read.fasta(file.path(src.dir,"IMGT Gene-DB AA V gaps",file))
      
      #total sequence data table
      seq_vector<-apply(sequences$ali,1,function(x){paste0(x,collapse="")})
      seq_dt<-data.table(name=names(seq_vector))
      seq_dt[,seq:=seq_vector]
      seq_dt[,species:=rep(species,length(seq_vector))]
      seq_dt[,family:=rep(family,length(seq_vector))]
      imgt_allseq_list[[file]]<-seq_dt
      
      #consensus sequence data table
      consq<-paste0(consensus(sequences)$seq,collapse="")
      imgt_consq<-c(imgt_consq,consq)
      imgt_family<-c(imgt_family,family)
      imgt_species<-c(imgt_species,species)
}
imgt_con_dt[,seq:=imgt_consq]
imgt_con_dt[,family:=imgt_family]
imgt_con_dt[,species:=imgt_species]
imgt_allseq<-rbindlist(imgt_allseq_list) #merge data tables from list into a single table with all sequences

#read in IgMotif results and PDB database files
pdb_dt<-readRDS(file=file.path("output","pdb_dt.RDS"))# open pdb_dt file
meanT_run<-readRDS(file=file.path("mean temps",paste0(meanT,"_runMean_ABNalign_",pdb_list_mode,"_",rms_method,".rds")))
#meanT_gen<-readRDS(file=file.path("mean temps",paste0(meanT,"_genMean_ABNalign_",mean_temp_mode,"_",rms_method,".rds")))

#Define lists of known catalytic antibodies or other groups by PDB ID
pdbids<-data.table(read.xlsx2("PDBid_lists.xlsx", sheetName = "ids"))
cabs<-unique(tolower(as.character(na.omit(pdbids$cabs_chk))))
cabs<-cabs[nchar(cabs)>0]
matches<-unique(meanT_run[mean_rms<=rms_cutoff]$pdb) #matches defined by rms_cutoff from mean template
imgt_abs<-data.table(read.xlsx2("PDBid_lists.xlsx", sheetName = "IMGT_pdb_Igs",stringsAsFactors=FALSE,
                               colClasses=c(rep("character",7),"numeric","character")))
imgt_abs[,PDB.release.date:=as.Date(imgt_abs$PDB.release.date,"%d-%B-%y")]
sw20_abs<-tolower(as.character(list.files(path=file.path(pdb.dir,"SW20- all Abs"),pattern=".pdb")))
all_abs<-unique(meanT_run$pdb) #all antibodies to analyze
rm(meanT_run)

#Clean up seq_dt table
seq_dt<-pdb_dt[PDBID %in% all_abs & len > 0] # remove any empty and non-Ab sequences 
rm(pdb_dt)
seq_dt[,match:=PDBID %in% matches]
seq_dt[,PDBID_count:=sapply(seq_dt[,PDBID], function(x){length(seq_dt[PDBID == x,PDBID])})] # count chains per PDBID
if(paired){seq_dt<-seq_dt[PDBID_count==2,]} #remove chains that are not paired (1 H and 1 L per PDBID)
seq_dt[,refgene:=sapply(strsplit(seq_dt[,refgene],"\\|"),function(x){paste(x[2],sep="|")})] #remove species from refgene label
setkey(seq_dt,name)

#Write fasta files of all Ab, Cab, and Matching antibody sequences to output directory
if(write.fasta){
      seqinr::write.fasta(as.list(seq_dt[chain=="L",seq]),names=seq_dt[chain=="L",name],
                  file.out=file.path("output","Lseqs_all.fasta"), as.string=TRUE) #write L chain fasta file
      seqinr::write.fasta(as.list(seq_dt[chain=="H",seq]),names=seq_dt[chain=="H",name],
                  file.out=file.path("output","Hseqs_all.fasta"), as.string=TRUE) #write H chain fasta file
      seqinr::write.fasta(as.list(seq_dt[,seq]),names=seq_dt[,name],
                  file.out=file.path("output","seqs_all.fasta"), as.string=TRUE) #write combined fasta file
      
      seqinr::write.fasta(as.list(seq_dt[chain=="L" & cab,seq]),names=seq_dt[chain=="L" & cab,name],
                  file.out=file.path("output","Lseqs_cabs.fasta"), as.string=TRUE) #write L chain fasta file
      seqinr::write.fasta(as.list(seq_dt[chain=="H" & cab,seq]),names=seq_dt[chain=="H" & cab,name],
                  file.out=file.path("output","Hseqs_cabs.fasta"), as.string=TRUE) #write H chain fasta file

      seqinr::write.fasta(as.list(seq_dt[chain=="L" & match,seq]),names=seq_dt[chain=="L" & match,name],
                  file.out=file.path("output","Lseqs_matches.fasta"), as.string=TRUE) #write L chain fasta file
      seqinr::write.fasta(as.list(seq_dt[chain=="H" & match,seq]),names=seq_dt[chain=="H" & match,name],
                  file.out=file.path("output","Hseqs_matches.fasta"), as.string=TRUE) #write H chain fasta file
      set.seed(1212)
      randL<-sample(seq_dt[chain=="L" & !cab,name],length(seq_dt[chain=="L" & cab, name])) #list of random control L seqs
      set.seed(1987)
      randH<-sample(seq_dt[chain=="H" & !cab,name],length(seq_dt[chain=="H" & cab, name])) #list of random control H seqs
      seqinr::write.fasta(as.list(seq_dt[name %in% randL,seq]),names=seq_dt[name %in% randL,name],
                          file.out=file.path("output","Lseqs_rand.fasta"), as.string=TRUE) #write L chain fasta file
      seqinr::write.fasta(as.list(seq_dt[name %in% randH,seq]),names=seq_dt[name %in% randH,name],
                          file.out=file.path("output","Hseqs_rand.fasta"), as.string=TRUE) #write H chain fasta file
}

### must do clustal omega alignment and download results in fasta & trees in Newick format
### msa package and ape package
#Lalign<-msaClustalOmega(seq_dt[len>0 & chain=="L",seq], type="protein")
#Halign<-msaClustalOmega(seq_dt[len>0 & chain=="H",seq], type="protein")
#Lalign<-read.fasta("Lalign.fasta",seqtype="AA",seqonly=FALSE,as.string=TRUE)
#Halign<-read.fasta("Halign.fasta",seqtype="AA",seqonly=FALSE,as.string=TRUE)
#Htree<-read.tree("Htree.txt")
#Ltree<-read.tree("Ltree.txt")

### Create alignment tables
#align_dt<-data.table("name"=c(attr(Lalign,"name"),attr(Halign,"name")))
#align_dt[,align:=c(as.character(Lalign),as.character(Halign))]
#seq_dt<-merge(align_dt,seq_dt,by=c("name"),all=TRUE)
#Lmat<-do.call(rbind,strsplit(seq_dt[len>0 & chain=="L",align],""))
#Hmat<-do.call(rbind,strsplit(seq_dt[len>0 & chain=="H",align],""))
#Lx<-seq(1,dim(Lmat)[2])
#Hx<-seq(1,dim(Hmat)[2])

###Do clustering analysis BIO3D package###
#Total_align<-bio3d::read.fasta(paste0("align_",sequences,".fasta"))


### Create VL phylogeny graph
distances<-Biostrings::stringDist(seq_dt[chain == "L",seq], method = "levenshtein")
#distances<-Biostrings::stringDist(seq_dt[chain == "L",seq], method = "substitutionMatrix",
#                                  substitutionMatrix="BLOSUM80", gapOpening=15,
#                                  gapExtension=2)

#kclust<-kmeans(distances,2)
hcluster<-hclust(distances, method = "ward.D2")
dend<-as.dendrogram(hcluster)

#Assign colors to category labels using RColorBrewer
species_cols<-brewer.pal(length(unique(seq_dt[,species])),"Set2")
names(species_cols)<-unique(seq_dt[,species])
family_cols<-brewer.pal(length(unique(seq_dt[,family])),"RdYlGn")
names(family_cols)<-unique(seq_dt[,family])
TF_cols<-brewer.pal(length(as.character(unique(seq_dt[,match]))),"Blues")[c(1,3)]
names(TF_cols)<-as.character(unique(seq_dt[,match]))
refgene_cols<-rep_len(brewer.pal(12,"Set3"),length(unique(seq_dt[chain=="L",refgene])))
names(refgene_cols)<-unique(seq_dt[chain=="L",refgene])
color_codes<-c(species_cols,family_cols,TF_cols,refgene_cols)

# Circos dendrogram plots
circos_Ldend(dend,seq_dt,color_codes,filename="Ldend_all.pdf")

## Create VH phylogey plot
distances<-Biostrings::stringDist(seq_dt[chain == "H",seq], method = "levenshtein")

#kclust<-kmeans(distances,2)
hcluster<-hclust(distances, method = "ward.D2")
dend<-as.dendrogram(hcluster)

#Assign colors to category labels using RColorBrewer
species_cols<-brewer.pal(length(unique(seq_dt[,species])),"Set2")
names(species_cols)<-unique(seq_dt[,species])
family_cols<-brewer.pal(length(unique(seq_dt[,family])),"RdYlGn")
TF_cols<-brewer.pal(length(as.character(unique(seq_dt[,match]))),"Blues")[c(1,3)]
names(TF_cols)<-as.character(unique(seq_dt[,match]))
refgene_cols<-rep_len(brewer.pal(12,"Set3"),length(unique(seq_dt[chain=="L",refgene])))
names(refgene_cols)<-unique(seq_dt[chain=="L",refgene])
color_codes<-c(species_cols,family_cols,TF_cols,refgene_cols)

# Circos dendrogram plots
circos_Hdend(dend,seq_dt,color_codes,filename="Hdend_all.pdf")

###====================Make graphs of positional sequence entropy=======================###
Lalign_cabs<-bio3d::read.fasta(file.path("alignments",paste0("Lalign_cabs.fasta")))
Lalign_rand<-bio3d::read.fasta(file.path("alignments",paste0("Lalign_rand.fasta")))
sequences<-c("rand", "cabs") #cabs or rand

Lalign_cabs$ali
sum(Lalign_cabs$ali[,2]=="D")
sum(Lalign_cabs$ali[,2]=="D")/nrow(Lalign_cabs$ali)
sum(Lalign_cabs$ali[,27]=="S")
sum(Lalign_cabs$ali[,27]=="S")/nrow(Lalign_cabs$ali)
sum(Lalign_cabs$ali[,100]=="H")
sum(Lalign_cabs$ali[,100]=="H")/nrow(Lalign_cabs$ali)

Lalign_rand$ali
sum(Lalign_rand$ali[,116]=="D")
sum(Lalign_rand$ali[,116]=="D")/nrow(Lalign_rand$ali)
sum(Lalign_rand$ali[,141]=="S")
sum(Lalign_rand$ali[,141]=="S")/nrow(Lalign_rand$ali)
sum(Lalign_rand$ali[,214]=="H")
sum(Lalign_rand$ali[,214]=="H")/nrow(Lalign_rand$ali)

# Plot entropy and residue frequencies (excluding positions >=60 percent gaps)
for (set in sequences){
      Lalign<-eval(parse(text=paste0("Lalign_",set)))
      h   <- entropy(Lalign)
      con <- bio3d::consensus(Lalign)
      names(h$H)=con$seq
      
      ## find IMGT start by finding location of first conserved Cys and setting to position 23
      cys23<-which((h$H<0.15) & (names(h$H)=="C"))[1]
      if(cys23>22){start<-cys23-22}else{start<-1}
      
      ## bar colors for catalytic triad
      bar_cols<-rep("grey",max.position)
      bar_cols[c(1,26,99)]<-"firebrick"
      
      h.sub <- entropy(Lalign$ali[,start:(start+max.position-1)])# Entropy for sub-alignment (positions start to max) 
      H <- h.sub$H.norm
      H[ apply(h.sub$freq[21:22,],2,sum)>=0.6 ] = 0
      
      
      ### BIG GRAPH
      big.graph<-paste0("VL_",set,".pdf")
      pdf(file=file.path("figs",big.graph),width=22,height=10)
      col <- mono.colors(32)
      aa  <- rev(rownames(h.sub$freq))
      oldpar <- par(no.readonly=TRUE)
      layout(matrix(c(1,2),2,1,byrow = TRUE), widths = 7, 
             heights = c(2, 8), respect = FALSE)
      
      # Plot 1: entropy
      par(mar = c(0, 4, 2, 2))
      barplot(H, border="white", ylab = "Entropy", space=0, xlim=c(4.2,max.position-4.2), yaxt="n",
              col=bar_cols)
      axis(side=2, at=c(0.2,0.4, 0.6, 0.8))
      axis(side=3, at=(seq(5,max.position-1,by=5)-0.5), labels=seq(5,max.position-1,by=5))
      box()
      
      # Plot2: residue frequencies 
      par(mar = c(5, 4, 0, 2))
      image(x=1:ncol(con$freq[,start:(start+max.position-1)]),
            y=1:nrow(con$freq[,start:(start+max.position-1)]),
            z=as.matrix(rev(as.data.frame(t(con$freq[,start:(start+max.position-1)])))),
            col=col, yaxt="n", xaxt="n",
            xlab="Alignment Position", ylab="Residue Type")
      axis(side=1, at=seq(5,max.position-1,by=5))
      axis(side=2, at=c(1:22), labels=aa)
      axis(side=3, at=c(1:max.position), labels =con$seq[start:(start+max.position-1)],cex.axis=.5)
      axis(side=4, at=c(1:22), labels=aa)
      grid(max.position, length(aa))
      box()
      
      for(i in start:(start+max.position-1)) {
            text(i-start+1, which(aa==con$seq[i]),con$seq[i],col="white", cex=0.75)
      }
      #abline(h=c(2.5, 5.5, 5.5, 3.5, 7.5, 9.5,
      #           12.5, 14.5, 16.5, 19.5), col="gray") # grouping of AAs in table
      par(oldpar)
      dev.off()
      
      ### SMALL GRAPH
      small.graph<-paste0("VLsm_",set,".pdf")
      # Plot entropy and residue frequencies (excluding positions >=60 percent gaps)
      pdf(file=file.path("figs",small.graph),width=22,height=6)

      col <- mono.colors(32)
      aa  <- rev(rownames(h.sub$freq))
      oldpar <- par(no.readonly=TRUE)
      
      # Plot 1: entropy
      par(mar = c(2, 4, 2, 2))
      barplot(H, border="white", ylab = "Entropy",space=0, xlim=c(4.2,max.position-4.2), yaxt="n",
              col=bar_cols)
      axis(side=1, at=c(1:max.position)-0.5, labels =con$seq[start:(start+max.position-1)],cex.axis=.5)
      axis(side=2, at=c(0.2,0.4, 0.6, 0.8))
      axis(side=3, at=(seq(5,max.position-1,by=5)-0.5),
           labels=seq(5,max.position-1,by=5))
      box()
      par(oldpar)
      dev.off()
}

setwd(initial.dir)
