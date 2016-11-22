# ========================= PDB Templates Version 1.3 ============================ #
# Anthony Bowen
# anthony.bowen@med.einstein.yu.edu
# Laboratory of Arturo Casadevall
# Albert Einstein College of Medicine, Bronx NY
# Last Modified 10/15/16; -faster rot.lsq, fit.xyz, and rmsd functions,
#                        -added hydrolase set,
#                        -output the residue names of triads,
#                        -fixed error when only one acid, base, or nucleophile is present in a
#                         given protein
#                        -fixed bug that did not properly remove duplicate sequences
#                        -added plot outputs of 3D coordinates of seed and mean templates

# Future additions:
#     -classify canonical CDR groups
#     -compare sequence alignment to structural template alignment
#     -option for whether or not to include backbone atoms and protons in SCg and SCm point calculations
#     -use GPU to speed up calculations
#     -glycine cannot be used in templates yet becuse it has no CB atoms
#     -transition to mmCIF file types for PDB structures
#     -implement database for file storage and better memory usage
#     -allow asparagine, glutamine residues, OCS hetatm for CYS (i.e. in 9PAP) or other modifications

# This script will take lists of PDB files, generate templates, superimpose them, and analyse them

#load needed packages
library("data.table")
library("bio3d")
library("lubridate")
library("xlsx")
library("XML")
library("Biostrings")
library("plyr")
library("dplyr")
library("dbscan")
library("ggplot2")
library("parallel")
library("rgl")
#library("gpuR")
#library("R.utils")

##################################### Setup #####################################
version<-"1.3, 10/15/2016"
set.seed(12121987) #set seed for reproducibility
ostype<-Sys.info()["sysname"]
if(!(ostype %in% c("Windows","Darwin"))){stop("Please specify directories for unrecognized OS")}
#use.gpu<-FALSE #whether to use available gpu for matrix computations
use.cpu.cores<-FALSE #whether to parallelize with CPU cores (not supported in Windows)
initial.dir<-getwd() # store the current directory
if(ostype=="Windows"){temp.dir<-"B:/Tony/Documents/R Working Directory/IgMotif_temp"}
if(ostype=="Darwin"){temp.dir<-file.path("~","R Working Directory","IgMotif_temp")}
if(ostype=="Windows"){pdb_db<-"B:/Tony/Documents/R Working Directory/IgMotif_temp/PDB_DB/PDB_files"} #2500K
if(ostype=="Darwin"){pdb_db<-file.path(temp.dir,"PDB_DB")} #Macbook Pro
if(ostype=="Windows"){pdb.dir<-"B:/Tony/Documents/Casadevall Lab/Non-redundant PDB files"} #2500K
if(ostype=="Darwin"){pdb.dir<-file.path("~","Documents","Casadevall Lab","Non-redundant PDB files")} #Macbook Pro
if(ostype=="Windows"){src.dir<-"B:/Tony/Google Drive/R Programming/IgMotif"} #2500K
if(ostype=="Darwin"){src.dir<-file.path("~","Google Drive","R Programming","IgMotif")} #Macbook Pro
setwd(src.dir) # change to the src directory
options(scipen=10) # number of digits before scientific notation is used

##################################### Script options #####################################
template_gen<-FALSE #whether to generate templates from PDB files
pdb_list_mode<-"antibodies" # "models", "antibodies", "cabs", "hydrolases"
temp_min<-0 #minimum distance for template residues
temp_max<-11 #maximum distance for template residues
base_atms<-c("HISND1","HISNE2","LYSNZ","ARGNH1","ARGNH2","ARGNE") #possible bases for creating templates
acid_atms<-c("HISND1","HISNE2","ASPOD1","ASPOD2","GLUOE1","GLUOE2") #possible acids for creating templates
nuc_atms<-c("THROG1","CYSSG","SEROG","TYROH") #possible nucleophiles for creating templates
pointnames<-c("CA","CB","O","C","N","SCa","SCg","SCm") #points to use for each template *not currently used, update in future version to choose specific template points
amino_acids<-c("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")
Lscore<- -350 #minimum alignment score to classify an L chain for antibodies
Hscore<- -350 #minimum alignment score to classify an H chain for antibodies
rem_duplicates<-TRUE #whether to remove duplicate sequences from analysis
temp_res_cutoff<-NA #resolution cutoff for all template structures, NA for no cutoff
ab_HL_required<-FALSE #require antibodies to have both H and L chain or not

gen_mean_temp<-FALSE #whether to generate a mean template based on the seed and abn table
mean_temp_mode<-"serprot" #"models", "antibodies", "hydrolases", "serprot", "cysprot", "thrprot", "metprot", "aspprot", "glycosylase"
rms_method<-"rmsd or DRMS" # "rmsd", "DRMS", "rmsd or DRMS", "rmsd and DRMS" for distance measure cutoff
avg_rms_cutoff<-1 #RMS difference cutoff for averaging templates
avg_res_cutoff<-2.5 #resolution cutoff for averaging templates
seed<- #choose one seed template from below by uncommenting
#"4hdi-1L-98L-26L" #3E5 G3 1
#"4hdi-1L-98L-28L" #3E5 G3 2
#"2h1p-1L-98L-26L" #2H1 G1 1
#"2h1p-1L-98L-28L" #2H1 G1 2
#"18b7_G1_2H1Pmodel-1L-98L-26L" #18B7 G1 model 1
#"18b7_G1_2H1Pmodel-1L-98L-28L" #18B7 G1 model 2
#"18b7_G3_4HDImodel-1L-98L-26L" #18B7 G3 model 1
#"18b7_G3_4HDImodel-1L-98L-28L" #18B7 G3 model 2
#"1utn-102A-57A-195A" # 1UTN- bovine trypsin: A- D102 H57 S195
"1st2-32A-64A-221A" # 1ST2- bacillus subtilisin: A- D32 H64 S221
#"1bqi-175A-159A-25A" #1BQI- papain triad 1: A- N175 H159 C25, must allow Asn
#"1bqi-158A-159A-25A" #1BQI- papain triad 2: A- D158 H159 C25
#"1lay-157A-63A-132A" #1LAY- CMV protease triad: A- H157 H63 S132

run_mean_temp<-TRUE #whether to open a saved mean template and run against the abn table
meanT<- #choose one mean template from below by uncommenting
#"4hdi-1L-98L-26L" #3E5 G3 1
#"4hdi-1L-98L-28L" #3E5 G3 2
#"2h1p-1L-98L-26L" #2H1 G1 1
#"2h1p-1L-98L-28L" #2H1 G1 2
#"18b7_G1_2H1Pmodel-1L-98L-26L" #18B7 G1 model 1
#"18b7_G1_2H1Pmodel-1L-98L-28L" #18B7 G1 model 2
#"18b7_G3_4HDImodel-1L-98L-26L" #18B7 G3 model 1
#"18b7_G3_4HDImodel-1L-98L-28L" #18B7 G3 model 2
#"1utn-102A-57A-195A" #1UTN- bovine trypsin: A- H57 S195 D102
"1st2-32A-64A-221A" #1ST2- bacillus subtilisin: A- H64 S221 D32
#"1bqi-175A-159A-25A" #1BQI- papain triad 1: A- N175 H159 C25, must allow Asn
#"1bqi-158A-159A-25A" #91BQI- papain triad 2: A- H159 C25 D158
#"1lay-157A-63A-132A" #1LAY- CMV protease triad: A- H157 H63 S132
template_chunks<-500000 #break template table into chunks if longer than this

enriched_temps<-FALSE #whether to find most enriched templates in a given group compared to ref group
seed_group<-"model" # temp, pdb, model, cabs; group to use templates as seeds
seed_list<-c("18b7_G3_4HDImodel","18b7_G1_2H1Pmodel") #add custom list of templates or pdb ids here
mean_group<-"antibodies" #cabs, antibodies; group to calculate mean templates from
enriched_group<-"cabs" # list, cabs; group to look for template enrichment
enriched_list<-c() #add custom list of PDB IDs here to find enriched templates
reference_group<-"antibodies" # list, antibodies; group to use as referece
reference_list<-c() #add custom list of PDB IDs here as reference
find_mean_templates<-TRUE #whether to create a mean for each template in the enriched group

##################################### Custom functions #####################################
#sum x+(x-1)+(x-2)... down to the number n
desc_sum<-function(x,n=0){
      sum<-0
      while(x>=n){
            sum<-sum+x
            x<-x-1
      }
      return(sum)
}
#euclidean distance between two vectors or adjacent rows of two matrices
euc_dist<-function(x,y){
      if(is.vector(x) & is.vector(y)){
            return(sqrt(sum((x-y)^2)))
      }
      if(is.matrix(x) & is.matrix(y)){
            return(sqrt(rowSums((x-y)^2)))
      }
      warning("euc_dist does not recognize x and y as vectors or matrices")
}
#template_table to create abn_dt for a given set of acid base nucleophile matrices and the coordinate table
template_table<-function(coords,base_mat,acid_mat,nuc_mat){
      
      #get atom coordinates and check residue numbering
      ba_dist<-apply(base_mat,1,function(row){
            row<-data.table(row)
            mat<-cbind(row,t(acid_mat))
            mat<-t(mat)
            dist(mat)[1:(length(acid_mat[,x]))]
      })
      if(is.vector(ba_dist)){
            ba_dist<-t(matrix(ba_dist))
      }
      bn_dist<-apply(base_mat,1,function(row){
            row<-data.table(row)
            mat<-cbind(row,t(nuc_mat))
            mat<-t(mat)
            dist(mat)[1:(length(nuc_mat[,x]))]
      })
      if(is.vector(bn_dist)){
            bn_dist<-t(matrix(bn_dist))
      }
      
      #ID base-acid and base-nucleophile pairs between temp_min and temp_max distances
      # then combine ba and bn pairs into threesomes
      ba_logical<- ba_dist>=temp_min & ba_dist<=temp_max
      bn_logical<- bn_dist>=temp_min & bn_dist<=temp_max
      ba_index<-which(ba_logical==TRUE,arr.ind=TRUE)
      bn_index<-which(bn_logical==TRUE,arr.ind=TRUE)
      ba_dt<-data.table(ba_index)
      setnames(ba_dt,c("acid_i","base_i"))
      bn_dt<-data.table(bn_index)
      setnames(bn_dt,c("nuc_i","base_i"))
      abn_dt<-rbindlist(apply(ba_dt,1,function(x){
            abn<-bn_dt[base_i==x[2]]
            abn[,acid_i:=x[1]]
      }))
      
      #find residues of abn indices and get info from the coord table. remove duplicate triads
      abn_dt[,nuc_atomkey:=coords[nuc==TRUE,atomkey][nuc_i]]
      abn_dt[,base_atomkey:=coords[base==TRUE,atomkey][base_i]]
      abn_dt[,acid_atomkey:=coords[acid==TRUE,atomkey][acid_i]]
      abn_dt[,nuc_reskey:=coords[nuc==TRUE,reskey][nuc_i]]
      abn_dt[,base_reskey:=coords[base==TRUE,reskey][base_i]]
      abn_dt[,acid_reskey:=coords[acid==TRUE,reskey][acid_i]]
      abn_dt[,triadkey:=paste(coords[acid==TRUE,reskey][acid_i],coords[base==TRUE,reskey][base_i],
                              coords[nuc==TRUE,reskey][nuc_i],sep="-")]
      setkey(abn_dt,triadkey)
      abn_dt<-unique(abn_dt)
      
      #rename nuc_i, base_i, and acid_i columns with residue id
      abn_dt[,nuc_i:=coords[nuc==TRUE,resid][nuc_i]]
      abn_dt[,base_i:=coords[base==TRUE,resid][base_i]]
      abn_dt[,acid_i:=coords[acid==TRUE,resid][acid_i]]
      
      #label coordinate table with acitve atoms and residue masses, add pairwise distances
      atom_mass<-c(12.01,14.01,32.07,16.00,1.01,1.01)
      names(atom_mass)<-c("C","N","S","O","H","D")
      active_anchors<-c(GLY='CA',ALA='CB',SER='OG',THR='OG1',CYS='SG',VAL='CB',LEU='CG',ILE='CD1',MET='SD',
                        PRO='CG',PHE='CZ',TYR=' OH ',TRP='NE1',ASP='CG',GLU='CD',ASN='ND2',GLN='NE2',
                        HIS='ND1',LYS='NZ',ARG='NE')
      coords[,mass:=atom_mass[substr(elety,1,1)]]
      coords[,sc_active:=(elety==active_anchors[resid])]
      
      #add pairwise distances to abn_dt
      #abn_dt[,ba_d:=euc_dist(coords[atomkey==base_atomkey,c("x","y","z"),with=FALSE],
      #                       coords[atomkey==acid_atomkey,c("x","y","z"),with=FALSE])]
      
      ##get nucleophile coordinates
      setkey(abn_dt,nuc_reskey) #set abn key for getting nucleophile coordinates
      ###get C alpha CA_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & elety == "CA",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCA_n")
      setnames(abn_dt,"y","yCA_n")
      setnames(abn_dt,"z","zCA_n")
      ###get C beta CB_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & elety == "CB",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCB_n")
      setnames(abn_dt,"y","yCB_n")
      setnames(abn_dt,"z","zCB_n")
      ###get O carbonyl O_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & elety == "O",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xO_n")
      setnames(abn_dt,"y","yO_n")
      setnames(abn_dt,"z","zO_n")
      ###get C carbonyl C_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & elety == "C",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xC_n")
      setnames(abn_dt,"y","yC_n")
      setnames(abn_dt,"z","zC_n")
      ###get N backbone N_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & elety == "N",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xN_n")
      setnames(abn_dt,"y","yN_n")
      setnames(abn_dt,"z","zN_n")
      ###get side-chain active atom SCa_n coordinates
      y<-coords[reskey %in% abn_dt$nuc_reskey & sc_active==TRUE,c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xSCa_n")
      setnames(abn_dt,"y","ySCa_n")
      setnames(abn_dt,"z","zSCa_n")
      ###get side-chain geometric mean SCg_n and weighted mean SCm_n coordinates
      y<-summarise(group_by(coords,reskey),mean(x),mean(y),mean(z),
                   sum(x*mass)/sum(mass),sum(y*mass)/sum(mass),sum(z*mass)/sum(mass))
      names(y)<-c("reskey","xSCg_n","ySCg_n","zSCg_n","xSCm_n","ySCm_n","zSCm_n")
      y<-data.table(y)
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="nuc_reskey", by.y="reskey") #do the merge
      
      ##get base coordinates
      setkey(abn_dt,base_reskey) #set abn key for getting base coordinates
      ###get C alpha CA_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & elety == "CA",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCA_b")
      setnames(abn_dt,"y","yCA_b")
      setnames(abn_dt,"z","zCA_b")
      ###get C beta CB_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & elety == "CB",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCB_b")
      setnames(abn_dt,"y","yCB_b")
      setnames(abn_dt,"z","zCB_b")
      ###get O carbonyl O_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & elety == "O",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xO_b")
      setnames(abn_dt,"y","yO_b")
      setnames(abn_dt,"z","zO_b")
      ###get C carbonyl C_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & elety == "C",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xC_b")
      setnames(abn_dt,"y","yC_b")
      setnames(abn_dt,"z","zC_b")
      ###get N backbone N_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & elety == "N",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xN_b")
      setnames(abn_dt,"y","yN_b")
      setnames(abn_dt,"z","zN_b")
      ###get side-chain active atom SCa_b coordinates
      y<-coords[reskey %in% abn_dt$base_reskey & sc_active==TRUE,c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xSCa_b")
      setnames(abn_dt,"y","ySCa_b")
      setnames(abn_dt,"z","zSCa_b")
      ###get side-chain geometric mean SCg_b and weighted mean SCm_b coordinates
      y<-summarise(group_by(coords,reskey),mean(x),mean(y),mean(z),
                   sum(x*mass)/sum(mass),sum(y*mass)/sum(mass),sum(z*mass)/sum(mass))
      names(y)<-c("reskey","xSCg_b","ySCg_b","zSCg_b","xSCm_b","ySCm_b","zSCm_b")
      y<-data.table(y)
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="base_reskey", by.y="reskey") #do the merge
      
      ##get acid coordinates
      setkey(abn_dt,acid_reskey) #set abn key for getting acid coordinates
      ###get C alpha CA_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & elety == "CA",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCA_a")
      setnames(abn_dt,"y","yCA_a")
      setnames(abn_dt,"z","zCA_a")
      ###get C beta CB_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & elety == "CB",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xCB_a")
      setnames(abn_dt,"y","yCB_a")
      setnames(abn_dt,"z","zCB_a")
      ###get O carbonyl O_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & elety == "O",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xO_a")
      setnames(abn_dt,"y","yO_a")
      setnames(abn_dt,"z","zO_a")
      ###get C carbonyl C_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & elety == "C",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xC_a")
      setnames(abn_dt,"y","yC_a")
      setnames(abn_dt,"z","zC_a")
      ###get N backbone N_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & elety == "N",c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xN_a")
      setnames(abn_dt,"y","yN_a")
      setnames(abn_dt,"z","zN_a")
      ###get side-chain active atom SCa_a coordinates
      y<-coords[reskey %in% abn_dt$acid_reskey & sc_active==TRUE,c("reskey","x","y","z"), with = FALSE] #merge
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      setnames(abn_dt,"x","xSCa_a")
      setnames(abn_dt,"y","ySCa_a")
      setnames(abn_dt,"z","zSCa_a")
      ###get side-chain geometric mean SCg_a and weighted mean SCm_a coordinates
      y<-summarise(group_by(coords,reskey),mean(x),mean(y),mean(z),
                   sum(x*mass)/sum(mass),sum(y*mass)/sum(mass),sum(z*mass)/sum(mass))
      names(y)<-c("reskey","xSCg_a","ySCg_a","zSCg_a","xSCm_a","ySCm_a","zSCm_a")
      y<-data.table(y)
      setkey(y,reskey) #set key of coordinates to merge
      abn_dt<-merge(abn_dt, y, by.x="acid_reskey", by.y="reskey") #do the merge
      abn_dt<-unique(abn_dt,by="triadkey") #remove duplicate traids by the reskey column
      setkey(abn_dt,triadkey)
      
      # Generate point distances for dRMS alignments from template coordinates
      coord_mat<-as.matrix(abn_dt[,grep("^x|^y|^z",names(abn_dt)),with=FALSE])
      ca_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="CA"]]
      cb_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="CB"]]
      o_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="O"]]
      c_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="C"]]
      n_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="N"]]
      sca_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="SCa"]]
      scg_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="SCg"]]
      scm_mat<-coord_mat[,colnames(coord_mat)[substring(sub("_.+","",colnames(coord_mat),perl=TRUE),2)=="SCm"]]
      
      dca<-cbind(euc_dist(ca_mat[,1:3],ca_mat[,4:6]),euc_dist(ca_mat[,1:3],ca_mat[,7:9]),euc_dist(ca_mat[,4:6],ca_mat[,7:9]))
      dcb<-cbind(euc_dist(cb_mat[,1:3],cb_mat[,4:6]),euc_dist(cb_mat[,1:3],cb_mat[,7:9]),euc_dist(cb_mat[,4:6],cb_mat[,7:9]))
      do<-cbind(euc_dist(o_mat[,1:3],o_mat[,4:6]),euc_dist(o_mat[,1:3],o_mat[,7:9]),euc_dist(o_mat[,4:6],o_mat[,7:9]))
      dc<-cbind(euc_dist(c_mat[,1:3],c_mat[,4:6]),euc_dist(c_mat[,1:3],c_mat[,7:9]),euc_dist(c_mat[,4:6],c_mat[,7:9]))
      dn<-cbind(euc_dist(n_mat[,1:3],n_mat[,4:6]),euc_dist(n_mat[,1:3],n_mat[,7:9]),euc_dist(n_mat[,4:6],n_mat[,7:9]))
      dsca<-cbind(euc_dist(sca_mat[,1:3],sca_mat[,4:6]),euc_dist(sca_mat[,1:3],sca_mat[,7:9]),euc_dist(sca_mat[,4:6],sca_mat[,7:9]))
      dscg<-cbind(euc_dist(scg_mat[,1:3],scg_mat[,4:6]),euc_dist(scg_mat[,1:3],scg_mat[,7:9]),euc_dist(scg_mat[,4:6],scg_mat[,7:9]))
      dscm<-cbind(euc_dist(scm_mat[,1:3],scm_mat[,4:6]),euc_dist(scm_mat[,1:3],scm_mat[,7:9]),euc_dist(scm_mat[,4:6],scm_mat[,7:9]))
      
      d_coords<-data.table(dca,dcb,do,dc,dn,dsca,dscg,dscm)
      setnames(d_coords,c(1:ncol(d_coords)),c("dCA_nb","dCA_na","dCA_ba",
                                              "dCB_nb","dCB_na","dCB_ba",
                                              "dCO_nb","dO_na","dO_ba",
                                              "dC_nb","dC_na","dC_ba",
                                              "dN_nb","dN_na","dN_ba",
                                              "dSCa_nb","dSCa_na","dSCa_ba",
                                              "dSCg_nb","dSCg_na","dSCg_ba",
                                              "dSCm_nb","dSCm_na","dSCm_ba"))
      d_coords[,triadkey:=abn_dt$triadkey]
      setkey(d_coords,triadkey)
      abn_dt<-merge(abn_dt,d_coords,by="triadkey")
      
      #remove triads where the same residue is included more than once
      abn_dt<-abn_dt[(dSCa_nb != 0) & (dSCa_na != 0) & (dSCa_ba != 0),]
      
      return(abn_dt)
}
#multiplot to organize and plot multiple ggplot2 objects together
# thanks to Winston Chang's Cookbook for R, http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
      library(grid)
      
      # Make a list from the ... arguments and plotlist
      plots <- c(list(...), plotlist)
      
      numPlots = length(plots)
      
      # If layout is NULL, then use 'cols' to determine layout
      if (is.null(layout)) {
            # Make the panel
            # ncol: Number of columns of plots
            # nrow: Number of rows needed, calculated from # of cols
            layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                             ncol = cols, nrow = ceiling(numPlots/cols))
      }
      
      if (numPlots==1) {
            print(plots[[1]])
            
      } else {
            # Set up the page
            grid.newpage()
            pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
            
            # Make each plot, in the correct location
            for (i in 1:numPlots) {
                  # Get the i,j matrix positions of the regions that contain this subplot
                  matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
                  
                  print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
            }
      }
}
#z.prop gives the z stastic of a test between two proportions x1/n1 and x2/n2
z.prop = function(x1,x2,n1,n2){
      numerator = (x1/n1) - (x2/n2)
      p.common = (x1+x2) / (n1+n2)
      denominator = sqrt(p.common * (1-p.common) * (1/n1 + 1/n2))
      z.prop.ris = numerator / denominator
      return(z.prop.ris)
}
#rmsdF is a faster replacement for getting the rmsd between a vector and a matrix
rmsdF<-function(a,b,ncore=1,nseg.scale=1,fit=FALSE){
      ncore <- setup.ncore(ncore)
      if (ncore > 1) {
            R_NCELL_LIMIT_CORE = 268435448
            R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE
            if (nseg.scale < 1) {
                  warning("nseg.scale should be 1 or a larger integer\n")
                  nseg.scale = 1
            }
      }
      if (is.vector(b) || nrow(b) == 1) {
            if (fit) {b <- fit.xyzF(fixed = a, mobile = b)}
            irmsd <- sqrt(sum((a-b)^2)/(length(a)/3))
            return(round(irmsd, 3))
      }
      else{
            if(fit){b<-fit.xyzF(fixed=a,mobile=b,ncore=ncore,nseg.scale=nseg.scale)}
            if(ncore > 1){
                  RLIMIT = R_NCELL_LIMIT
                  nDataSeg = floor((nrow(b) - 1)/RLIMIT) + 1
                  nDataSeg = floor(nDataSeg * nseg.scale)
                  lenSeg = floor(nrow(b)/nDataSeg)
                  irmsd = vector("list", nDataSeg)
                  for (i in 1:nDataSeg) {
                        istart = (i - 1) * lenSeg + 1
                        iend = if (i < nDataSeg) 
                              i * lenSeg
                        else nrow(b)
                        irmsd[[i]] <- mclapply(istart:iend, function(j) sqrt(sum((b[j,] - a[])^2)/(length(a[])/3)), 
                                               mc.preschedule = TRUE)
                  }
                  irmsd <- unlist(irmsd)
            }
            else {
                  irmsd <- sqrt(apply((apply(b, 1, "-", a)^2), 2, sum)/(length(a)/3))
            }
            return(round(irmsd, 3))
      }
}
#fit.xyzF is a faster replacement for fitting a vector to a matrix
fit.xyzF<-function(fixed,mobile,ncore=1,nseg.scale = 1){
      ncore <- setup.ncore(ncore)
      if(is.vector(mobile)){
            if (!is.numeric(mobile)) 
                  stop("input 'mobile' should be numeric")
            if (any(is.na(fixed)) || any(is.na(mobile))) {
                  stop(" NA elements selected for fitting (check indices)")
            }
            fit <- rot.lsqF(xx = mobile, yy = fixed)
            return(fit)
      }
      if (ncore > 1) {
            R_NCELL_LIMIT_CORE = 268435448
            R_NCELL_LIMIT = ncore * R_NCELL_LIMIT_CORE
            if (nseg.scale < 1) {
                  warning("nseg.scale should be 1 or a larger integer\n")
                  nseg.scale = 1
            }
            RLIMIT = floor(R_NCELL_LIMIT/ncol(mobile))
            nDataSeg = floor((nrow(mobile) - 1)/RLIMIT) + 1
            nDataSeg = floor(nDataSeg * nseg.scale)
            lenSeg = floor(nrow(mobile)/nDataSeg)
            fit = vector("list", nDataSeg)
            for (i in 1:nDataSeg) {
                  istart = (i - 1) * lenSeg + 1
                  iend = if (i < nDataSeg) 
                        i * lenSeg
                  else nrow(mobile)
                  fit[[i]] <- mclapply(istart:iend, function(j) rot.lsqF(xx = mobile[j,], yy = fixed), mc.preschedule = TRUE)
            }
            fit <- matrix(unlist(fit), ncol = ncol(mobile), byrow = TRUE)
      }
      else{fit <- t(apply(mobile,1,rot.lsqF,yy=fixed))}
      return(fit)
}
#rot.lsqF is a faster replacement function for the Bio3D package function rot.lsq
rot.lsqF<-function(xx, yy, xfit = rep(TRUE, length(xx)), yfit = xfit, verbose = FALSE){
      xx <- matrix(xx, nrow = 3)
      x <- matrix(xx[xfit], nrow = 3)
      y <- matrix(yy[yfit], nrow = 3)
      if (length(x) != length(y)) 
            stop("dimension mismatch in x and y")
      xbar <- rowMeans(x)
      ybar <- rowMeans(y)
      xx <- xx-xbar
      x <- x-xbar
      y <- y-ybar
      R <- y %*% t(x)
      RR <- t(R) %*% R
      prj <- eigen(RR)
      prj$values[prj$values < 0 & prj$values >= -1e-12] <- 1e-12
      A <- prj$vectors
      b <- A[, 1]
      c <- A[, 2]
      A[1, 3] <- (b[2] * c[3]) - (b[3] * c[2])
      A[2, 3] <- (b[3] * c[1]) - (b[1] * c[3])
      A[3, 3] <- (b[1] * c[2]) - (b[2] * c[1])
      B <- R %*% A
      #B <- sweep(B, 2, sqrt(prj$values), "/")
      B<-t(t(B)/sqrt(prj$values)) #to avoid using sweep above
      b <- B[, 1]
      c <- B[, 2]
      B[1, 3] <- (b[2] * c[3]) - (b[3] * c[2])
      B[2, 3] <- (b[3] * c[1]) - (b[1] * c[3])
      B[3, 3] <- (b[1] * c[2]) - (b[2] * c[1])
      U <- B %*% t(A)
      xx <- U %*% xx
      if (verbose) {
            x <- U %*% x
            frmsd <- sqrt(sum((x - y)^2)/dim(y)[2])
            cat("#rmsd= ", round(frmsd, 6), "\n")
      }
      xx <- xx+ybar
      as.vector(xx)
}

##################################### Script sections #####################################
if(use.cpu.cores==TRUE){
      cpu_cores<-max(1,detectCores()-1) #number of CPU cores to use for Parallel package
      if(Sys.info()["sysname"]=="Windows"){cpu_cores<-1;print("Windows OS detected, using 1 CPU core")}
}else{cpu_cores<-1}

#designate PDB list directory
if (pdb_list_mode=="models"){pdb.dir<-file.path(pdb.dir,"models")}else{pdb.dir<-file.path(pdb_db,"PDB_files")}
#if (pdb_list_mode=="antibodies" | pdb_list_mode=="cabs"){pdb.dir<-file.path(pdb_db,"PDB_files")}
pdb_dt<-readRDS(file=file.path("output","pdb_dt.RDS"))# open pdb_dt file

# 1. generate templates if selected and save the file
if (template_gen){
      # define list of pdb files
      pdb_list<-sub(".pdb","",list.files(pdb.dir)[grep(".pdb",list.files(pdb.dir))])
      if(rem_duplicates){
            Lchains<-pdb_dt[chain=="L" & score>=Lscore & duplicate_seq==FALSE]$PDBID
            Hchains<-pdb_dt[chain=="H" & score>=Hscore & duplicate_seq==FALSE]$PDBID
            enzymes<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & hydrolase == TRUE]$PDBID)
            glycosylase<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & glycosylase == TRUE]$PDBID)
            serprot<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & serprot == TRUE]$PDBID)
            cysprot<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & cysprot == TRUE]$PDBID)
            thrprot<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & thrprot == TRUE]$PDBID)
            metprot<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & metprot == TRUE]$PDBID)
            aspprot<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & aspprot == TRUE]$PDBID)
      }else{
            Lchains<-pdb_dt[chain=="L" & score>=Lscore]$PDBID
            Hchains<-pdb_dt[chain=="H" & score>=Hscore]$PDBID
            enzymes<-unique(pdb_dt[len>0 & hydrolase == TRUE]$PDBID)
            glycosylase<-unique(pdb_dt[len>0 & glycosylase == TRUE]$PDBID)
            serprot<-unique(pdb_dt[len>0 & serprot == TRUE]$PDBID)
            cysprot<-unique(pdb_dt[len>0 & cysprot == TRUE]$PDBID)
            thrprot<-unique(pdb_dt[len>0 & thrprot == TRUE]$PDBID)
            metprot<-unique(pdb_dt[len>0 & metprot == TRUE]$PDBID)
            aspprot<-unique(pdb_dt[len>0 & aspprot == TRUE]$PDBID)
      }
      if(ab_HL_required){abs<-intersect(Lchains,Hchains)}else{abs<-union(Lchains,Hchains)}
      enzymes<-enzymes[!(enzymes %in% abs)] #remove enzymes with antibodies
      glycosylase<-glycosylase[!(glycosylase %in% abs)]
      serprot<-serprot[!(serprot %in% abs)]
      cysprot<-cysprot[!(cysprot %in% abs)]
      thrprot<-thrprot[!(thrprot %in% abs)]
      aspprot<-aspprot[!(aspprot %in% abs)]
      metprot<-metprot[!(metprot %in% abs)]
      
      if(pdb_list_mode=="antibodies"){pdb_list<-abs}
      if(pdb_list_mode=="cabs"){
            if(rem_duplicates){
                  Lchains<-pdb_dt[chain=="L" & score>=Lscore & duplicate_seq==FALSE & cab==TRUE]$PDBID
                  Hchains<-pdb_dt[chain=="H" & score>=Hscore & duplicate_seq==FALSE & cab==TRUE]$PDBID
            }else{
                  Lchains<-pdb_dt[chain=="L" & score>=Lscore & cab==TRUE]$PDBID
                  Hchains<-pdb_dt[chain=="H" & score>=Hscore & cab==TRUE]$PDBID
            }
            if(ab_HL_required){pdb_list<-intersect(Lchains,Hchains)}else{pdb_list<-union(Lchains,Hchains)}
      }
      if(pdb_list_mode=="hydrolases"){pdb_list<-unique(enzymes)}
      if(!is.na(temp_res_cutoff)){pdb_list[pdb_list %in% pdb_dt[resolution<=temp_res_cutoff]$PDBID]}
      
      # create empty abn data table to hold triads for all files
      abn<-data.table(matrix(nrow=0,ncol=108))
      setnames(abn,c( "key","pdb","triadkey","acid_reskey","base_reskey","nuc_reskey", "nuc_i","base_i","acid_i", 
                      "nuc_atomkey","base_atomkey","acid_atomkey","xCA_n","yCA_n","zCA_n",
                      "xCB_n","yCB_n","zCB_n","xO_n","yO_n","zO_n","xC_n","yC_n","zC_n",
                      "xN_n","yN_n","zN_n","xSCa_n","ySCa_n","zSCa_n","xSCg_n","ySCg_n","zSCg_n",
                      "xSCm_n","ySCm_n","zSCm_n","xCA_b","yCA_b","zCA_b","xCB_b","yCB_b","zCB_b",
                      "xO_b","yO_b","zO_b","xC_b","yC_b","zC_b","xN_b","yN_b","zN_b",
                      "xSCa_b","ySCa_b","zSCa_b","xSCg_b","ySCg_b","zSCg_b","xSCm_b","ySCm_b","zSCm_b",
                      "xCA_a","yCA_a","zCA_a","xCB_a","yCB_a","zCB_a","xO_a", "yO_a","zO_a",
                      "xC_a","yC_a","zC_a","xN_a","yN_a","zN_a","xSCa_a","ySCa_a","zSCa_a",
                      "xSCg_a","ySCg_a","zSCg_a","xSCm_a","ySCm_a","zSCm_a","dCA_nb","dCA_na","dCA_ba",
                      "dCB_nb","dCB_na","dCB_ba","dCO_nb","dO_na","dO_ba","dC_nb","dC_na","dC_ba",
                      "dN_nb","dN_na","dN_ba","dSCa_nb","dSCa_na","dSCa_ba","dSCg_nb","dSCg_na","dSCg_ba",
                      "dSCm_nb","dSCm_na","dSCm_ba"))
      abn[,key:=as.character(key)]
      setkey(abn,key)
      empty.files<-character()
      empty.atoms<-character()
      no.files<-character()
      no.viable.chains<-character()
      incomplete.templates<-character()
      error.templates<-character()
      
      # open each file and generate 3-residue template table from the structure
      pdb_list<-paste0(pdb_list,".pdb")
      for (file in pdb_list){
            path<-file.path(pdb_db,file)
            
            # read in file and do checks
            if (file.exists(path)){
                  if (file.size(path)>0){
                        txt<-read.csv(path,colClasses="character",row.names=NULL,header=FALSE)[[1]]
                        sumtest<-sum(sapply(txt,function(x){grepl("^ATOM ",x)}))
                        if(sumtest<2){empty.atoms<-c(empty.atoms,file)}
                  }
                  else{sumtest<-0; empty.files<-c(empty.files,file)}
                  if((sumtest>1) ){
                        data<-bio3d::read.pdb(path)
                  }else{warning(paste0(file," cannot be read due to sumtest ",sumtest));next}
            }else{warning(paste0(path," does not exist"));no.files<-c(no.files,file);next}
            
            #get info about chains, sequences, and remove duplicates
            chains<-unique(as.factor(names(data$seqres[data$seqres %in% amino_acids])))
            seqs<-as.factor(data$seqres)
            table<-data.table(chain=chains)
            table[,seq:=sapply(chain,function(x){paste0(aa321(seqs[names(seqs)==x]),collapse="")})]
            table[,len:=sapply(seq,nchar)] #add column for chain length
            #if using antibodies from pdb_dt, get the right chains
            if(pdb_list_mode=="antibodies" |pdb_list_mode=="cabs"){
                  sequences<-pdb_dt[PDBID==sub(".pdb","",file) & len!=0]
                  table[,antibody:=seq %in% sequences$seq]
                  table<-table[antibody==TRUE]
            }
            table[,duplicate:=duplicated(seq)] #label chains with same sequence
            table<-table[!(seq %in% grep("X",table$seq,value=TRUE))] #remove sequences with "X"
            table<-table[duplicate==FALSE,] #remove duplicate chains
            chains<-as.character(table$chain)
            
            # check residue numbering and chain number
            coords<-data.table(data$atom[atom.select(data,"protein",chain=chains)$atom,])
            coords<-coords[type=="ATOM"] #select only ATOM rows, not HETATM
            if(nrow(coords)==0){
                  no.viable.chains<-c(no.viable.chains,file)
                  warning(paste0(file," has no viable chains. Skipping..."))
                  next
            }
            atom_counts<-table(coords$elety)
            backbone_atoms<-c("CA","C","N","O")
            if(!all(atom_counts["CA"] == atom_counts[names(atom_counts) %in% backbone_atoms])){
                  incomplete.templates<-c(incomplete.templates,file)
                  warning(paste0(file," has incomplete backbone atoms"))
                  next
            }
            alt_conf_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],17,17)) #alternate conformation list
            chain_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],21,22)) #chain ID list
            resno_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],23,29)) #residue num list
            
            atomtxt<-txt[grepl("^ATOM ",txt)]
            atomtxt<-atomtxt[which((alt_conf_txt==""|alt_conf_txt=="A") & chain_txt %in% chains)]
            chain_txt<-gsub("\\s","",substring(atomtxt,21,22)) #chain ID list
            resno_txt<-gsub("\\s","",substring(atomtxt,23,29)) #residue num list
            txt_reskey<-paste0(resno_txt,chain_txt)
            coords[,reskey:=paste0(resno,chain)]
            
            if (sum(nchar(table$seq)) != length(unique(coords$reskey))){ #are there enough unique residue keys from the bio3D read?
                  if(length(resno_txt)>=nrow(coords)){ #replace reskey if text parse is long enough
                        setnames(coords,"resno","resno_num")
                        coords[,resno:=resno_txt[1:nrow(coords)]]
                        coords[,reskey:=txt_reskey[1:nrow(coords)]]
                        if(length(resno_txt)>nrow(coords)){warning(paste0(file,
                                          " has extra residue keys. check for multiple models."))}
                  }else{
                        print(paste0("ERROR: cannot find unique residue keys for ",file))
                        error.templates<-c(error.templates,file)
                        next
                  }
            }
            coords[,atomkey:=paste0(resno,chain,eleno)]
            setkey(coords,atomkey)
            
            # label acid base or nucleophile atoms and compute pairwise distances
            coords[,abnkey:=paste0(resid,elety)]
            coords[,base:=abnkey %in% base_atms]
            coords[,acid:=abnkey %in% acid_atms]
            coords[,nuc:=abnkey %in% nuc_atms]
            base_mat<-coords[base==TRUE,c("x","y","z"),with=FALSE]
            acid_mat<-coords[acid==TRUE,c("x","y","z"),with=FALSE]
            nuc_mat<-coords[nuc==TRUE,c("x","y","z"),with=FALSE]
            if(nrow(base_mat)==0|nrow(acid_mat)==0|nrow(nuc_mat)==0){
                  incomplete.templates<-c(incomplete.templates,file)
                  warning(paste0(file," has no ABN atoms"))
                  next
            }
            
            #generate triad table from matrices of acid, base, nucleophile atoms and coordinate table
            abn_dt<-template_table(coords,acid_mat=acid_mat,base_mat=base_mat,nuc_mat=nuc_mat)
            
            #add PDB ID and unique triad key to table. then merge with the total abn table
            abn_dt[,pdb:=sub(".pdb",'',file)]
            abn_dt[,key:=paste(pdb,triadkey,sep="-")]
            setcolorder(abn_dt,names(abn))
            setkey(abn_dt,key)
            
            abn<-rbind(abn,abn_dt)
      }
      
      attr(abn,"script_version")<-version
      attr(abn,"date")<-today()
      attr(abn,"pdb_list_mode")<-pdb_list_mode
      attr(abn,"pdb_list")<-pdb_list
      attr(abn,"min_dist")<-temp_min
      attr(abn,"max_dist")<-temp_max
      attr(abn,"base_atoms")<-base_atms
      attr(abn,"acid_atoms")<-acid_atms
      attr(abn,"nucleophile_atoms")<-nuc_atms
      #attr(abn,"atom_masses")<-atom_mass #atom_mass is now in the template function
      #attr(abn,"sidechain_active_atoms")<-active_anchors #active_anchors now in template function
      attr(abn,"remove_duplicate_seqs")<-rem_duplicates
      attr(abn,"template_resolution_cutoff")<-temp_res_cutoff
      attr(abn,"antibody_L_min_score")<-Lscore
      attr(abn,"antibody_H_min_score")<-Hscore
      attr(abn,"incomplete_templates")<-incomplete.templates
      attr(abn,"no_viable_chains")<-no.viable.chains
      attr(abn,"error_templates")<-error.templates
      attr(abn,"empty_files")<-empty.files
      attr(abn,"empty_atoms")<-empty.atoms
      attr(abn,"no_files")<-no.files
      saveRDS(abn,file.path("output",paste0("abn_triads_",pdb_list_mode,".rds")))
}

# 2. generate a mean template from the specificed seed template and the abn triad table
#    seed template PDB file must be in seed folder
if (gen_mean_temp){
      #abn<-abn[1:100000] #truncate abn table for testing purposes
      
      if(mean_temp_mode %in% c("serprot","cysprot","thrprot","metprot","aspprot","glycosylase")){
            if(rem_duplicates){
                  Lchains<-pdb_dt[chain=="L" & score>=Lscore & duplicate_seq==FALSE]$PDBID
                  Hchains<-pdb_dt[chain=="H" & score>=Hscore & duplicate_seq==FALSE]$PDBID
                  enzymes<-unique(pdb_dt[len>0 & duplicate_seq==FALSE & get(mean_temp_mode) == TRUE]$PDBID)
            }else{
                  Lchains<-pdb_dt[chain=="L" & score>=Lscore]$PDBID
                  Hchains<-pdb_dt[chain=="H" & score>=Hscore]$PDBID
                  enzymes<-unique(pdb_dt[len>0 & get(mean_temp_mode) == TRUE]$PDBID)
                  
            }
            enzymes<-enzymes[!(enzymes %in% unique(c(Lchains,Hchains)))]
            pdb_list<-unique(enzymes)
            abn<-readRDS(file.path("output",paste0("abn_triads_hydrolases.rds")))
            abn<-abn[pdb %in% pdb_list]
      }else{
            abn<-readRDS(file.path("output",paste0("abn_triads_",mean_temp_mode,".rds")))
            pdb_list<-attr(abn,"pdb_list")
      }
      
      #get seed template if it's not already in abn table
      if (!(seed %in% abn$key)){
            
            seed_pdb<-file.path("mean temps","seeds",paste0(strsplit(seed,"-")[[1]][1],".pdb"))
            seed_data<-bio3d::read.pdb(seed_pdb)
            txt<-read.csv(seed_pdb,colClasses="character",row.names=NULL,header=FALSE)[[1]]
            seed_chains<-unique(sapply(strsplit(seed,"-")[[1]][2:4],function(x){substring(x,nchar(x))}))
            
            #get sequences
            chains<-unique(as.factor(names(seed_data$seqres)))
            seqs<-as.factor(seed_data$seqres)
            table<-data.table(chain=chains)
            table[,seq:=sapply(chain,function(x){paste0(aa321(seqs[names(seqs)==x]),collapse="")})]
            table[,len:=sapply(seq,nchar)] #add column for chain length
            table[,duplicate:=duplicated(seq)] #label chains with same sequence
            table<-table[chain %in% seed_chains]
            
            #check numbering for seed pdb file and get abn coordinates
            coords<-data.table(seed_data$atom[atom.select(seed_data,"protein",chain=seed_chains)$atom,])
            if(nrow(coords)==0){
                  warning(paste0(seed," has no viable chains."))
            }
            if(sum(c("CA","C","N","CB","O") %in% unique(coords$elety))!=5){
                  warning(paste0(seed," has incomplete backbone atoms"))
            }
            alt_conf_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],17,17)) #alternate conformation list
            chain_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],21,22)) #chain ID list
            resno_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],23,29)) #residue num list
            chain_txt<-chain_txt[alt_conf_txt==""|alt_conf_txt=="A"] #select only 1 conformation
            resno_txt<-resno_txt[alt_conf_txt==""|alt_conf_txt=="A"] #select only 1 conformation
            resno_txt<-resno_txt[chain_txt %in% seed_chains] #remove extra chains
            coords[,reskey:=paste0(resno,chain)]
            if (sum(nchar(table$seq)) != length(unique(coords$reskey))){ #are there enough unique residue keys from the bio3D read?
                  if(length(resno_txt)>=nrow(coords)){ #replace reskey if text parse is long enough
                        setnames(coords,"resno","resno_num")
                        coords[,resno:=resno_txt[1:nrow(coords)]]
                        coords[,reskey:=paste0(resno,chain)]
                        if(length(resno_txt)>nrow(coords)){warning(paste0(seed,
                                                                          " has extra residue keys. check for multiple models."))}
                  }else{print(paste0("ERROR: cannot find unique residue keys for ",seed))}
            }
            coords[,atomkey:=paste0(resno,chain,eleno)]
            setkey(coords,atomkey)
            
            # label acid base or nucleophile atoms and compute pairwise distances
            coords[,abnkey:=paste0(resid,elety)]
            coords[,base:=abnkey %in% base_atms]
            coords[,acid:=abnkey %in% acid_atms]
            coords[,nuc:=abnkey %in% nuc_atms]
            base_mat<-coords[base==TRUE,c("x","y","z"),with=FALSE]
            acid_mat<-coords[acid==TRUE,c("x","y","z"),with=FALSE]
            nuc_mat<-coords[nuc==TRUE,c("x","y","z"),with=FALSE]
            if(nrow(base_mat)==0|nrow(acid_mat)==0|nrow(nuc_mat)==0){
                  warning(paste0(seed," has no ABN atoms"))
            }
            
            #generate triad table for seed
            seed_abn<-template_table(coords,acid_mat=acid_mat,base_mat=base_mat,nuc_mat=nuc_mat)
            #add PDB ID and unique triad key to table. then merge with the total abn table
            seed_abn[,pdb:=strsplit(seed,"-")[[1]][1]]
            seed_abn[,key:=paste(pdb,triadkey,sep="-")]
            setcolorder(seed_abn,names(abn))
            setkey(seed_abn,key)
      }else{seed_abn<-abn[key==seed]}
      
      xyz_indices<-grep("^x|^y|^z",names(seed_abn))
      dRMS_indices<-grep("^d",names(seed_abn))
      coord_indices<-grep("^x|^y|^z|^d",names(seed_abn))
      if(!(seed %in% seed_abn$key)){stop(paste0("Error generating seed template for ",seed))}
      seed_xyz<-unlist(seed_abn[key==seed,xyz_indices,with=FALSE]) #original seed
      seed_dists<-unlist(seed_abn[key==seed,dRMS_indices,with=FALSE]) #original seed
      seed_coords<-c(seed_xyz,seed_dists)
      
      
      #remove templates for PDBs that don't meet resolution requirement
      if(!is.na(avg_res_cutoff)){abn<-abn[pdb %in% pdb_dt[resolution<=avg_res_cutoff]$PDBID]}
      
      #profvis({ #profiling
      averaged<-character() #triads already averaged
      averaged_iters<-list() #averaging iteration for averaged templates
      iterations<-numeric()
      find_mean<-TRUE
      iteration<-1
      
      #create alignment table to store superimposed coordinates
      setkey(abn,key)
      align_coords<-data.table(t(c(seed,seed_coords)))
      setnames(align_coords,c("key",names(seed_coords)))
      align_coords<-rbind(align_coords, abn[key!=seed,c("key",names(abn)[coord_indices]),with=FALSE])
      align_coords[,rmsd:=as.numeric(NaN)]
      align_coords[,DRMS:=as.numeric(NaN)]
      align_coords[key==seed,rmsd:=0]
      align_coords[key==seed,DRMS:=0]
      setcolorder(align_coords,c("key","rmsd","DRMS",names(seed_coords)))
      setkey(align_coords,key)
      
      while(find_mean){
            #do the superposition
            print(paste0("Starting mean template superposition iteration #",iteration))
            
            #apply the superposition and rmsd functions to all rows in the table
#            align_coords<-align_coords[,{
#
#                  row<-as.numeric(mget(names(seed_xyz)))
#                  dRMS<-as.numeric(mget(names(seed_dists)))
#                  seed_xyz<-as.numeric(seed_xyz)
#                  seed_dists<-as.numeric(seed_dists)
#                  
#                  t_align<-fit.xyz(seed_xyz,row,fixed.inds=c(1:length(seed_xyz)),
#                                   mobile.inds=c(1:length(row)),ncore=cpu_cores)
#                  t_dRMS<-rmsd(seed_dists,dRMS,ncore=cpu_cores)
#                 t_rmsd<-rmsd(seed_xyz,t_align,ncore=cpu_cores)
#                  as.list(c(t_rmsd,t_dRMS,t_align,dRMS))
#            },by=key]
            
            mat1<-as.matrix(align_coords[,names(align_coords)[grep("^x|^y|^z",names(align_coords))],with=FALSE])
            mat2<-as.matrix(align_coords[,names(align_coords)[grep("^d",names(align_coords))],with=FALSE])
            mat1<-matrix(as.numeric(mat1),ncol=ncol(mat1))
            mat2<-matrix(as.numeric(mat2),ncol=ncol(mat2))
            mseed_xyz<-as.numeric(seed_xyz)
            seed_dists<-as.numeric(seed_dists)
            t_align<-fit.xyzF(seed_xyz,mat1,ncore=cpu_cores)
            t_dRMS<-rmsdF(seed_dists,mat2,ncore=cpu_cores)
            t_rmsd<-rmsdF(seed_xyz,t_align,ncore=cpu_cores)
            t_rmsd[is.nan(t_rmsd)]<-0 #set NaN superpositions to 0 (identical coordinates)
            align_coords2<-data.table(align_coords[,"key",with=FALSE],t_rmsd,t_dRMS,t_align,mat2)
            align_coords<-align_coords2
            rm(align_coords2)
            
            setnames(align_coords,c(1:ncol(align_coords)),c("key","rmsd","DRMS",names(seed_coords)))
            setkey(align_coords,key)
            
            #Find the indices to average
            if(rms_method=="rmsd"){hits<-which(align_coords$rmsd<=avg_rms_cutoff)}
            if(rms_method=="DRMS"){hits<-which(align_coords$DRMS<=avg_rms_cutoff)}
            if(rms_method=="rmsd or DRMS"){hits<-which(align_coords$rmsd<=avg_rms_cutoff | align_coords$DRMS<=avg_rms_cutoff)}
            if(rms_method=="rmsd and DRMS"){hits<-which(align_coords$rmsd<=avg_rms_cutoff & align_coords$DRMS<=avg_rms_cutoff)}
            
            iterations[iteration]<-nrow(align_coords[hits])
            if(iteration==1){
                  if(iterations[iteration]>1){newmean=TRUE}else{newmean=FALSE}
                  if(!(seed %in% abn$key)){
                        abn[,seed_rmsd:=align_coords[key!=seed,rmsd]]
                        abn[,seed_DRMS:=align_coords[key!=seed,DRMS]]
                  }else{
                        abn[,seed_rmsd:=align_coords$rmsd]
                        abn[,seed_DRMS:=align_coords$DRMS]
                  }
            }else{if(iterations[iteration]>iterations[iteration-1]){newmean<-TRUE}else{newmean=FALSE}}
            
            #record new rmsd and DRMS values for this iteration 
            if(!(seed %in% abn$key)){
                  abn[,mean_rmsd:=align_coords[key!=seed,rmsd]]
                  abn[,mean_DRMS:=align_coords[key!=seed,DRMS]]
                  
            }else{
                  abn[,mean_rmsd:=align_coords[,rmsd]]
                  abn[,mean_DRMS:=align_coords[,DRMS]]
            }
            
            #average templates that are below the cutoff, record averaged template keys
            if(newmean){
                  mean_coords<-apply(align_coords[hits,grep("^x|^y|^z|^d",names(align_coords)),with=FALSE],
                                     2,function(column){mean(as.numeric(column))})
                  seed_coords<-mean_coords
                  seed_xyz<-seed_coords[grep("^x|^y|^z",names(seed_coords))]
                  seed_dists<-seed_coords[grep("^d",names(seed_coords))]
                  averaged<-as.character(align_coords[hits]$key) #averaged template
                  averaged_iters[[iteration]]<-averaged
                  iteration<-iteration+1
            }else{find_mean<-FALSE}
      }
      print(paste("Mean template generated with iterations ",paste0(iterations,collapse=", ")))
      #}) #profiling
      
      #save aligned coordinates to abn table
      if(!(seed %in% abn$key)){
            abn[,names(seed_coords):=align_coords[key!=seed,names(seed_coords),with=FALSE]]      
      }else{abn[,names(seed_coords):=align_coords[,names(seed_coords),with=FALSE]]}
      if(rms_method=="rmsd"){abn[,mean_rms:=mean_rmsd]}
      if(rms_method=="DRMS"){abn[,mean_rms:=mean_DRMS]}
      if(rms_method=="rmsd or DRMS"){abn[,mean_rms:=pmin(mean_rmsd,mean_DRMS)]}
      if(rms_method=="rmsd and DRMS"){abn[,mean_rms:=pmax(mean_rmsd,mean_DRMS)]}
      if(rms_method=="rmsd"){abn[,seed_rms:=seed_rmsd]}
      if(rms_method=="DRMS"){abn[,seed_rms:=seed_DRMS]}
      if(rms_method=="rmsd or DRMS"){abn[,seed_rms:=pmin(seed_rmsd,seed_DRMS)]}
      if(rms_method=="rmsd and DRMS"){abn[,seed_rms:=pmax(seed_rmsd,seed_DRMS)]}
      
      #Do clustering of mean_rmsd, mean_DRMS, seed_rmsd, and seed_DRMS points to identify groups
      #rms_cluster<-kmeans(abn[,c("mean_rmsd","mean_DRMS"),with=FALSE],2,nstart=25) #k-means
      rms_cluster<-dbscan(abn[,c("mean_rmsd","mean_DRMS"),with=FALSE],eps=0.2,minPts=5) #density-based
      seed_rms_cluster<-dbscan(abn[,c("seed_rmsd","seed_DRMS"),with=FALSE],eps=0.2,minPts=5)
      abn[,rms_cluster:=rms_cluster$cluster]
      abn[,seed_rms_cluster:=seed_rms_cluster$cluster]
      
      attr(abn,"script_version")<-version
      attr(abn,"seed")<-seed
      attr(abn,"date")<-today()
      attr(abn,"averaged")<-averaged
      attr(abn,"RMS_method")<-rms_method
      attr(abn,"rms_cutoff")<-avg_rms_cutoff
      attr(abn,"iterations")<-iterations
      attr(abn,"averaged")<-averaged_iters
      attr(abn,"avg_resolution_cutoff")<-avg_res_cutoff
      
      attr(mean_coords,"script_version")<-version
      attr(mean_coords,"seed")<-seed
      attr(mean_coords,"date")<-today()
      attr(mean_coords,"averaged")<-averaged
      attr(mean_coords,"iterations")<-iterations
      attr(mean_coords,"RMS_method")<-rms_method
      attr(mean_coords,"avg_rms_cutoff")<-avg_rms_cutoff
      attr(mean_coords,"avg_resolution_cutoff")<-avg_res_cutoff
      
      attr(mean_coords,"template_pdb_list_mode")<-attr(abn,"pdb_list_mode")
      attr(mean_coords,"template_pdb_list")<-attr(abn,"pdb_list")
      attr(mean_coords,"template_min_dist")<-attr(abn,"min_dist")
      attr(mean_coords,"template_max_dist")<-attr(abn,"max_dist")
      attr(mean_coords,"template_base_atoms")<-attr(abn,"base_atoms")
      attr(mean_coords,"template_acid_atoms")<-attr(abn,"acid_atoms")
      attr(mean_coords,"template_nucleophile_atoms")<-attr(abn,"nucleophile_atoms")
      attr(mean_coords,"template_atom_masses")<-attr(abn,"atom_masses")
      attr(mean_coords,"template_sidechain_active_atoms")<-attr(abn,"sidechain_active_atoms")
      attr(mean_coords,"template_remove_duplicate_seqs")<-attr(abn,"remove_duplicate_seqs")
      attr(mean_coords,"template_resolution_cutoff")<-attr(abn,"template_resolution_cutoff")
      attr(mean_coords,"template_antibody_L_min_score")<-attr(abn,"antibody_L_min_score")
      attr(mean_coords,"template_antibody_H_min_score")<-attr(abn,"antibody_H_min_score")
      saveRDS(mean_coords,file.path("mean temps",paste0(seed,"_",mean_temp_mode,"_",rms_method,".rds")))
      saveRDS(abn,file.path("mean temps",paste0(seed,"_genMean_ABNalign_",mean_temp_mode,"_",rms_method,".rds")))
      rm(abn)
}

# 3. run opened mean template against the opened abn triad table
if (run_mean_temp){
      print(paste0("Running mean template against the full dataset"))
      abn<-readRDS(file.path("output",paste0("abn_triads_",pdb_list_mode,".rds")))
      pdb_list<-attr(abn,"pdb_list")
      mean_coords<-readRDS(file.path("mean temps",paste0(meanT,"_",mean_temp_mode,"_",rms_method,".rds")))
      
      meanT<-paste0("mean_",meanT)
      coord_indices<-grep("^x|^y|^z|^d",names(abn))
      xyz_indices<-grep("^x|^y|^z",names(abn))
      dist_indices<-grep("^d",names(abn))
      mean_dists<-mean_coords[grep("^d",names(mean_coords))]
      mean_xyz<-mean_coords[grep("^x|^y|^z",names(mean_coords))]
      
      #break abn table into chunks if it is too long
      if(nrow(abn)>template_chunks){
            chunk_max<-ceiling(nrow(abn)/template_chunks)
            print(paste0("Breaking template table into ",chunk_max," chunks for run."))
      }else{chunk_max<-1}
      chunk_num<-1
      abn2<-abn[0,]
      abn2[,mean_rmsd:=numeric()]
      abn2[,mean_DRMS:=numeric()]
      abn2[,mean_rms:=numeric()]
      
      
      #loop through chunks
      while(chunk_num<=chunk_max){
            print(paste0("Running chunk number ",chunk_num,"."))
            abn<-readRDS(file.path("output",paste0("abn_triads_",pdb_list_mode,".rds")))
            if(chunk_max>1){
                  if(chunk_num<chunk_max){
                        abn<-abn[(((chunk_num-1)*template_chunks)+1):(chunk_num*template_chunks)]
                  }else if (chunk_num==chunk_max){
                        abn<-abn[(((chunk_num-1)*template_chunks)+1):nrow(abn)]
                  }
            }
            
            #create alignment table to store superimposed coordinates
            setkey(abn,key)
            align_coords<-data.table(t(c(meanT,mean_coords)))
            setnames(align_coords,c("key",names(mean_coords)))
            align_coords<-rbind(align_coords, abn[,c("key",names(abn)[coord_indices]),with=FALSE])
            align_coords[,rmsd:=as.numeric(NaN)]
            align_coords[key==meanT,rmsd:=0]
            align_coords[,DRMS:=as.numeric(NaN)]
            align_coords[key==meanT,DRMS:=0]
            setcolorder(align_coords,c("key","rmsd","DRMS",names(mean_coords)))
            setkey(align_coords,key)
            
            #apply the superposition and rmsd functions to all rows in the table
            #      align_coords<-align_coords[,{
            
            #define matrices whether use.gpu==TRUE
            #            row<-as.numeric(mget(names(mean_xyz)))
            #            dRMS<-as.numeric(mget(names(mean_dists)))
            #            meanmat<-as.numeric(mean_xyz)
            #            meand_mat<-as.numeric(mean_dists)
            #            meanc_mat<-as.numeric(mean_xyz)
            #            t_align<-fit.xyz(meanmat,row,fixed.inds=c(1:length(mean_xyz)),
            #                             mobile.inds=c(1:length(row)),ncore=cpu_cores)
            #            t_dRMS<-rmsd(meand_mat,dRMS,ncore=cpu_cores)
            #            t_rmsd<-rmsd(meanc_mat,t_align,ncore=cpu_cores)
            #            as.list(c(t_rmsd,t_dRMS,t_align,dRMS))
            #      },by=key]
            
            mat1<-as.matrix(align_coords[,names(align_coords)[grep("^x|^y|^z",names(align_coords))],with=FALSE])
            mat2<-as.matrix(align_coords[,names(align_coords)[grep("^d",names(align_coords))],with=FALSE])
            mat1<-matrix(as.numeric(mat1),ncol=ncol(mat1))
            mat2<-matrix(as.numeric(mat2),ncol=ncol(mat2))
            mean_xyz<-as.numeric(mean_xyz)
            mean_dists<-as.numeric(mean_dists)
            t_align<-fit.xyzF(mean_xyz,mat1,ncore=cpu_cores)
            rm(mat1)
            t_dRMS<-rmsdF(mean_dists,mat2,ncore=cpu_cores)
            t_rmsd<-rmsdF(mean_xyz,t_align,ncore=cpu_cores)
            t_rmsd[is.nan(t_rmsd)]<-0 #set NaN superpositions to 0 (identical coordinates)
            align_coords2<-data.table(align_coords[,"key",with=FALSE],t_rmsd,t_dRMS,t_align,mat2)
            rm(mat2)
            rm(t_align)
            align_coords<-align_coords2
            rm(align_coords2)
            
            setnames(align_coords,c(1:ncol(align_coords)),c("key","rmsd","DRMS",names(mean_coords)))
            setkey(align_coords,key)
            
            abn[,mean_rmsd:=align_coords[key!=meanT,rmsd]]
            abn[,mean_DRMS:=align_coords[key!=meanT,DRMS]]
            if(rms_method=="rmsd"){abn[,mean_rms:=align_coords[key!=meanT,rmsd]]}
            if(rms_method=="DRMS"){abn[,mean_rms:=align_coords[key!=meanT,DRMS]]}
            if(rms_method=="rmsd or DRMS"){abn[,mean_rms:=pmin(align_coords[key!=meanT]$rmsd,align_coords[key!=meanT]$DRMS)]}
            if(rms_method=="rmsd and DRMS"){abn[,mean_rms:=pmax(align_coords[key!=meanT]$rmsd,align_coords[key!=meanT]$DRMS)]}
            abn[,names(mean_coords):=align_coords[key!=meanT,names(mean_coords),with=FALSE]] #save aligned coordinates
            rm(align_coords)
            
            #join new chunk
            abn2<-rbind(abn2,abn)
            chunk_num<-chunk_num+1
      }
      
      abn<-abn2
      rm(abn2)
      saveRDS(abn,file.path("mean temps",paste0(seed,"_runMean_ABNalign_",pdb_list_mode,"_",rms_method,".rds")))
      
      #Do clustering of mean_rmsd, mean_DRMS points to identify groups
      cluster_mat<-as.matrix(abn[,c("mean_rmsd","mean_DRMS"),with=FALSE])
      rm(abn)
      #rms_cluster<-kmeans(cluster_mat,2,nstart=25)
      tryCatch({rms_cluster<<-dbscan(cluster_mat,eps=0.2,minPts=5)},error=function(x){
                        rms_cluster<<-1
                        print("Error in clustering. All points set to cluster 1.")
                  })
      rm(cluster_mat)
      abn<-readRDS(file.path("mean temps",paste0(seed,"_runMean_ABNalign_",pdb_list_mode,"_",rms_method,".rds")))
      if(is.list(rms_cluster)){
            abn[,rms_cluster:=rms_cluster$cluster]
      }else{abn[,rms_cluster:=rms_cluster]}
      
      attr(abn,"script_version")<-version
      attr(abn,"date")<-today()
      attr(abn,"mean_temp_mode")<-mean_temp_mode
      attr(abn,"mean_template")<-meanT
      attr(abn,"RMS_method")<-rms_method
      attr(abn,"meanT_attributes")<-attributes(mean_coords)
      saveRDS(abn,file.path("mean temps",paste0(seed,"_runMean_ABNalign_",pdb_list_mode,"_",rms_method,".rds")))
      
      # read template alignment files for mean generation and the mean run
      meanT_run<-abn
      rm(abn)
      #meanT_run<-readRDS(file=file.path("mean temps",paste0(seed,"_runMean_ABNalign_",pdb_list_mode,"_",rms_method,".rds")))
      meanT_gen<-readRDS(file=file.path("mean temps",paste0(seed,"_genMean_ABNalign_",mean_temp_mode,"_",rms_method,".rds")))
      
      #make 3D plots of seed and mean template points
      seedT<-meanT_gen[key==seed]
      cols<-c(rep("darkred",8),rep("darkblue",8),rep("gold",8))
      xind<-grep("^x",names(seedT))
      yind<-grep("^y",names(seedT))
      zind<-grep("^z",names(seedT))
      plot3d(x=seedT[,xind,with=FALSE],y=seedT[,yind,with=FALSE],z=seedT[,zind,with=FALSE],
             xlab="x",ylab="y",zlab="z",col=cols,box=FALSE,size=8,axes=TRUE,expand=1.03)
      rgl.postscript(filename=file.path("mean temps",paste0(seed,"_seed_points.pdf")),fmt="pdf")
      xind<-grep("^x",names(mean_coords))
      yind<-grep("^y",names(mean_coords))
      zind<-grep("^z",names(mean_coords))
      plot3d(x=mean_coords[xind],y=mean_coords[yind],z=mean_coords[zind],
             xlab="x",ylab="y",zlab="z",col=cols,box=FALSE,size=8,axes=TRUE,expand=1.03)
      rgl.postscript(filename=file.path("mean temps",paste0(seed,"_mean_",mean_temp_mode,"_points.pdf")),fmt="pdf")
      
      # plot distribution graphs and save to output folder
      pdf(file=file.path("mean temps",paste0(
            meanT,"_",mean_temp_mode,"_run_",pdb_list_mode,"_dists_",rms_method,".pdf")),width=6,height=4)
      par(mar=c(4,3,1,0),ps=8)
      par(mfrow=c(3,2))
      if(nrow(meanT_run[mean_rms<=avg_rms_cutoff])==0){zoom_y<-50}else{
            zoom_y<-round_any(max(table(cut(
                  meanT_run[mean_rms<=avg_rms_cutoff]$mean_rms,breaks=(avg_rms_cutoff/0.05))))*2,
                  10,f=ceiling)
      }
      
      hist(as.numeric(meanT_gen[,seed_rms]),
           breaks=seq(0,ceiling(max(as.numeric(meanT_gen[,seed_rms]))),by=0.05),
           main="Mean Template Generation Distribution",
           ylab="",
           xlab="",
           border=FALSE,
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(as.numeric(meanT_gen[,seed_rms])))/0.05)
                                   -(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the seed (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      hist(as.numeric(meanT_gen[,seed_rms]),
           xlim=c(0,(2*avg_rms_cutoff)),
           ylim=c(0,zoom_y),
           breaks=seq(0,ceiling(max(as.numeric(meanT_gen[,seed_rms]))),by=0.05),
           main="Averaged Templates",
           xlab="",
           ylab="",
           border="white",
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(as.numeric(meanT_gen[,seed_rms])))/0.05)
                                   -(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the seed (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      
      hist(meanT_gen$mean_rms,
           breaks=seq(0,ceiling(max(meanT_gen$mean_rms)),by=0.05),
           main="Mean Template Generation Distribution",
           ylab="",
           xlab="",
           border=FALSE,
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(meanT_gen$mean_rms))/0.05)-(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the mean (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      hist(meanT_gen$mean_rms,
           xlim=c(0,(2*avg_rms_cutoff)),
           ylim=c(0,zoom_y),
           breaks=seq(0,ceiling(max(meanT_gen$mean_rms)),by=0.05),
           main="Averaged Templates",
           xlab="",
           ylab="",
           border="white",
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(meanT_gen$mean_rms))/0.05)-(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the mean (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      
      hist(meanT_run$mean_rms,
           breaks=seq(0,ceiling(max(meanT_run$mean_rms)),by=0.05),
           main="Mean Template Run Distribution",
           xlab="",
           ylab="",
           border=FALSE,
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(meanT_run$mean_rms))/0.05)-(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the mean (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      hist(meanT_run$mean_rms,
           xlim=c(0,(2*avg_rms_cutoff)),
           ylim=c(0,zoom_y),
           breaks=seq(0,ceiling(max(meanT_run$mean_rms)),by=0.05),
           main="Superposition hits <= RMS cutoff",
           xlab="",
           ylab="",
           border="white",
           col=c(rep("cadetblue4",times=(avg_rms_cutoff/0.05)),
                 rep("gray",times=((ceiling(max(meanT_run$mean_rms))/0.05)-(avg_rms_cutoff/0.05)))))
      mtext(side=1,text="Template RMS from the mean (angstroms)",line=2,cex=(1))
      mtext(side=2,text="Number of templates",line=2,cex=(1))
      dev.off()
      
      
      run_hits<-(nrow(meanT_run[mean_rms<avg_rms_cutoff])>0)
      gen_hits<-(nrow(meanT_gen[mean_rms<avg_rms_cutoff])>0)
            
      #plot RMSD vs DRMS graphs in ggplot2
      p1 <- ggplot(meanT_gen, aes(x = seed_DRMS, y = seed_rmsd)) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="rmsd & DRMS") +
            ggtitle("Mean Template Generation: seed")
      if(gen_hits){
            p1<-p1+geom_point(aes(color=factor((seed_rmsd<=avg_rms_cutoff | seed_DRMS<=avg_rms_cutoff),
                  labels=c(paste0("Both >",avg_rms_cutoff),paste0("Either <=",avg_rms_cutoff)))),alpha=0.5)
      }
      
      p2 <- ggplot(meanT_gen, aes(x = seed_DRMS, y = seed_rmsd)) +
            geom_point(aes(color=factor(seed_rms_cluster)),alpha=0.5) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="RMS Cluster") +
            ggtitle("Mean Template Generation Clusters: seed")
      
      p3 <- ggplot(meanT_gen, aes(x = mean_DRMS, y = mean_rmsd)) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="rmsd & DRMS") +
            ggtitle("Mean Template Generation: mean")
      if(gen_hits){
            p3<-p3+geom_point(aes(color=factor((mean_rmsd<=avg_rms_cutoff | mean_DRMS<=avg_rms_cutoff),
                  labels=c(paste0("Both >",avg_rms_cutoff),paste0("Either <=",avg_rms_cutoff)))),alpha=0.5)
      }
      
      p4 <- ggplot(meanT_gen, aes(x = mean_DRMS, y = mean_rmsd)) +
            geom_point(aes(color=factor(rms_cluster)),alpha=0.5) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="RMS Cluster") +
            ggtitle("Mean Template Generation Clusters: mean")
      
      p5 <- ggplot(meanT_run, aes(x = mean_DRMS, y = mean_rmsd)) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="rmsd & DRMS") +
            ggtitle("Mean Template Run")
      if(run_hits){
            p5<-p5+geom_point(aes(color=factor((mean_rmsd<=avg_rms_cutoff | mean_DRMS<=avg_rms_cutoff),
                labels=c(paste0("Both >",avg_rms_cutoff),paste0("Either <=",avg_rms_cutoff)))),alpha=0.5)
      }
      
      p6 <- ggplot(meanT_run, aes(x = mean_DRMS, y = mean_rmsd)) +
            geom_point(aes(color=factor(rms_cluster)),alpha=0.5) +
            #geom_abline(slope=1,intercept=0,linetype="dashed",size=1) +
            geom_hline(yintercept=avg_rms_cutoff,linetype="dashed",size=1) +
            geom_vline(xintercept=avg_rms_cutoff,linetype="dashed",size=1)+
            geom_smooth(method="lm") +
            labs(x="DRMS",y="RMSD",color="RMS Cluster") +
            ggtitle("Mean Template Run Clusters")
      
      #write file and plot the objects
      jpeg(file=file.path("mean temps",paste0(meanT,"_",mean_temp_mode,"_run_",pdb_list_mode,"_RMSplots_",rms_method,".jpg")),
          units="px",width=4000,height=3200,pointsize=100,quality=100)
      multiplot(p1,p3,p5,p2,p4,p6,cols=2)
      dev.off()
      rm(p1,p2,p3,p4,p5,p6)
      rm(meanT_gen)
      
      ### create output table and save to file
      output_data<-meanT_run[,c("key","pdb","nuc_i","base_i","acid_i","mean_rmsd","mean_DRMS","rms_cluster"),with=FALSE]
      rm(meanT_run)
      setkey(output_data,key)
      output_data[,min_rms:=pmin(mean_rmsd, mean_DRMS)]
      output_data[,max_rms:=pmax(mean_rmsd, mean_DRMS)]
      output_data[,averaged:=(key %in% attr(mean_coords,"averaged"))]
      output_data<-merge(output_data,pdb_dt[(chain=="L")&(PDBID %in% output_data$pdb),
                        c("PDBID","resolution","cab","match","pdb_ab","imgt_ab","hydrolase",
                          "glycosylase","serprot","cysprot","thrprot","metprot","aspprot"),with=FALSE],
                        by.x="pdb",by.y="PDBID")
      if(rms_method=="rmsd"){setorder(output_data,mean_rmsd);output_data<-output_data[mean_rmsd<=(2*avg_rms_cutoff)]}
      if(rms_method=="DRMS"){setorder(output_data,mean_DRMS);output_data<-output_data[mean_DRMS<=(2*avg_rms_cutoff)]}
      if(rms_method=="rmsd or DRMS"){setorder(output_data,min_rms);output_data<-output_data[min_rms<=(2*avg_rms_cutoff)]}
      if(rms_method=="rmsd and DRMS"){setorder(output_data,max_rms);output_data<-output_data[max_rms<=(2*avg_rms_cutoff)]}
      write.table(output_data, file=file.path("mean temps",paste0(meanT,"_",mean_temp_mode,"_run_",
                                    pdb_list_mode,"_output_",rms_method,".tab")),
                                    sep="\t",row.names=FALSE,quote=FALSE)
      #rm(output_data)
}

# 4. run all templates from enriched_group against themselves and a reference group to find enriched templates
if (enriched_temps){
      
      #open abn file with all templates in it & remove duplicates if desired
      abn<-readRDS(file.path("output",paste0("abn_triads_",reference_group,".rds")))
      if(rem_duplicates){abn<-abn[pdb %in% pdb_dt[duplicate_seq==FALSE]$PDBID]}
      
      #get list of templates to use as seeds for generating means, if desired
      if(seed_group=="cabs"){seed_list<-unique(pdb_dt[cab==TRUE,]$PDBID)}
      seed_dt<-abn[(pdb %in% seed_list),]
      
      #get seed template if it's a model
      if (seed_group=="model"){
            
            for(seed in seed_list){
                  
                  seed_pdb<-file.path("mean temps","seeds",paste0(seed,".pdb"))
                  seed_data<-bio3d::read.pdb(seed_pdb)
                  txt<-read.csv(seed_pdb,colClasses="character",row.names=NULL,header=FALSE)[[1]]
                  
                  #get sequences
                  chains<-unique(as.character(names(seed_data$seqres)))
                  seqs<-as.factor(seed_data$seqres)
                  table<-data.table(chain=chains)
                  table[,seq:=sapply(chain,function(x){paste0(aa321(seqs[names(seqs)==x]),collapse="")})]
                  table[,len:=sapply(seq,nchar)] #add column for chain length
                  table[,duplicate:=duplicated(seq)] #label chains with same sequence
                  table<-table[duplicate==FALSE] #remove duplicate sequences
                  
                  #check numbering for seed pdb file and get abn coordinates
                  coords<-data.table(seed_data$atom[atom.select(seed_data,"protein",chain=table$chain)$atom,])
                  if(nrow(coords)==0){
                        warning(paste0(seed," has no viable chains."))
                  }
                  if(sum(c("CA","C","N","CB","O") %in% unique(coords$elety))!=5){
                        warning(paste0(seed," has incomplete backbone atoms"))
                  }
                  alt_conf_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],17,17)) #alternate conformation list
                  chain_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],21,22)) #chain ID list
                  resno_txt<-gsub("\\s","",substring(txt[grepl("^ATOM ",txt)],23,29)) #residue num list
                  chain_txt<-chain_txt[alt_conf_txt==""|alt_conf_txt=="A"] #select only 1 conformation
                  resno_txt<-resno_txt[alt_conf_txt==""|alt_conf_txt=="A"] #select only 1 conformation
                  resno_txt<-resno_txt[chain_txt %in% table$chain] #remove extra chains
                  coords[,reskey:=paste0(resno,chain)]
                  if (sum(nchar(table$seq)) != length(unique(coords$reskey))){ #are there enough unique residue keys from the bio3D read?
                        if(length(resno_txt)>=nrow(coords)){ #replace reskey if text parse is long enough
                              setnames(coords,"resno","resno_num")
                              coords[,resno:=resno_txt[1:nrow(coords)]]
                              coords[,reskey:=paste0(resno,chain)]
                              if(length(resno_txt)>nrow(coords)){warning(paste0(seed,
                                                                                " has extra residue keys. check for multiple models."))}
                        }else{print(paste0("ERROR: cannot find unique residue keys for ",seed))}
                  }
                  coords[,atomkey:=paste0(resno,chain,eleno)]
                  setkey(coords,atomkey)
                  
                  # label acid base or nucleophile atoms and compute pairwise distances
                  coords[,abnkey:=paste0(resid,elety)]
                  coords[,base:=abnkey %in% base_atms]
                  coords[,acid:=abnkey %in% acid_atms]
                  coords[,nuc:=abnkey %in% nuc_atms]
                  base_mat<-coords[base==TRUE,c("x","y","z"),with=FALSE]
                  acid_mat<-coords[acid==TRUE,c("x","y","z"),with=FALSE]
                  nuc_mat<-coords[nuc==TRUE,c("x","y","z"),with=FALSE]
                  if(nrow(base_mat)==0|nrow(acid_mat)==0|nrow(nuc_mat)==0){
                        warning(paste0(seed," has no ABN atoms"))
                  }
                  
                  #generate triad table for seed
                  seed_abn<-template_table(coords,acid_mat=acid_mat,base_mat=base_mat,nuc_mat=nuc_mat)
                  #add PDB ID and unique triad key to table. then merge with the total abn table
                  seed_abn[,pdb:=strsplit(seed,"-")[[1]][1]]
                  seed_abn[,key:=paste(pdb,triadkey,sep="-")]
                  setcolorder(seed_abn,names(abn))
                  seed_dt<-rbind(seed_dt,seed_abn)
            }
            
      }
      setkey(seed_dt,key)
      
      #get list of templates for generating means from seeds
      if(mean_group=="cabs"){mean_list<-unique(pdb_dt[cab==TRUE,]$PDBID)}
      if(mean_group=="antibodies"){mean_list<-unique(abn$pdb)}
      mean_dt<-abn[(pdb %in% mean_list[mean_list %in% pdb_dt[resolution<=avg_res_cutoff]$PDBID]),]
      mean_dt<-merge(seed_dt,mean_dt,by=names(seed_dt),all=TRUE) #add seed templates if not present
      
      #get coordinate indices
      xyz_indices<-grep("^x|^y|^z",names(seed_dt))
      dRMS_indices<-grep("^d",names(seed_dt))
      coord_indices<-grep("^x|^y|^z|^d",names(seed_dt))
      
      #mean_dt<-mean_dt[key %in% seed_dt$key] #shorten mean table for testing purposes
      #seed_dt<-seed_dt[1:50] #shorten seed table for testing purposes
      
      #generate a mean template for each seed and put in a data table, if desired
      if(find_mean_templates){
            if(seed_group=="pdb"|seed_group=="model"){title<-paste0(seed_group,seed_list,collapse=",")}else{title<-seed_group}
            if(!file.exists(file.path("mean temps",paste0(title,"_meanTempTable_",rms_method,".rds")))){
                  row<-0
                  seed_dt<-seed_dt[,{
                        row<-row+1
                        seed_coords<-as.numeric(mget(names(seed_dt[,coord_indices,with=FALSE])))
                        seed_xyz<-seed_coords[1:length(xyz_indices)]
                        seed_dists<-seed_coords[(length(xyz_indices)+1):length(coord_indices)]
                        print(paste0("Generating mean template for ",key,": ",row," of ",nrow(seed_dt)))
                        iteration<-1
                        iterations<-c(0)
                        find_mean<-TRUE
                        
                        while(find_mean){
                              print(paste0("Mean template superposition #",iteration))
                              
                              #create gpuMatrix objects if use.gpu==TRUE
                              dRMSmat<-as.matrix(mean_dt[,dRMS_indices,with=FALSE])
                              rmsdmat<-as.matrix(mean_dt[,xyz_indices,with=FALSE])

                              #Calc distanceRMS and rmsd vectors
                              dRMSvec<-rmsdF(seed_dists,dRMSmat,fit=FALSE,ncore=cpu_cores)
                              rmsdvec<-rmsdF(seed_xyz,rmsdmat,fit=TRUE,ncore=cpu_cores)
                              rmsdvec[is.nan(rmsdvec)]<-0 #set NaN superpositions to 0 (identical coordinates)
                              rmsd_fit<-which(rmsdvec!=0) #indices of templates to superimpose
                              rmsd_nonfit<-which(rmsdvec==0) #templates that will not be superimposed
                        
                              #Find the indices to average
                              if(rms_method=="rmsd"){hits<-which(rmsdvec<=avg_rms_cutoff)}
                              if(rms_method=="DRMS"){hits<-which(dRMSvec<=avg_rms_cutoff)}
                              if(rms_method=="rmsd or DRMS"){hits<-which(rmsdvec<=avg_rms_cutoff | dRMSvec<=avg_rms_cutoff)}
                              if(rms_method=="rmsd and DRMS"){hits<-which(rmsdvec<=avg_rms_cutoff & dRMSvec<=avg_rms_cutoff)}
                              iterations[iteration]<-length(hits)
                              
                              #Calculate an average and repeat
                              if(iterations[iteration]>1){
                                    if(iteration>1){
                                          if(!(iterations[iteration]>iterations[iteration-1])){find_mean<-FALSE;next}
                                    }
                                    seed_dists<-colMeans(mean_dt[hits,dRMS_indices,with=FALSE])
                                    meanmat<-as.matrix(mean_dt[rmsd_fit[rmsd_fit %in% hits],xyz_indices,with=FALSE])
                                    seed_xyz<-fit.xyzF(seed_xyz,meanmat,ncore=cpu_cores)
                                    seed_xyz<-colMeans(rbind(mean_dt[rmsd_nonfit,xyz_indices,with=FALSE],data.table(seed_xyz),use.names=FALSE))
                                    seed_coords<-c(seed_xyz,seed_dists)
                                    iteration<-iteration+1
                              }else{find_mean<-FALSE}
                        }
                        print(paste("Mean template generated with iterations ",paste0(iterations,collapse=", ")))
                        
                        #send output mean template as list
                        as.list(seed_coords)
                  },by=key]
                  setnames(seed_dt,c("key",names(abn[,coord_indices,with=FALSE])))
                  
                  #add in missing abn columns
                  add<-mean_dt[key %in% seed_dt$key,c(1,which(!(names(abn) %in% names(seed_dt)))),with=FALSE]
                  seed_dt<-merge(add,seed_dt,by="key")
                  setkey(seed_dt,key)
                  
                  attr(seed_dt,"script_version")<-version
                  attr(seed_dt,"date")<-today()
                  attr(seed_dt,"seed_group")<-seed_group
                  attr(seed_dt,"seed_list")<-seed_list
                  attr(seed_dt,"mean_group")<-mean_group
                  attr(seed_dt,"enriched_group")<-enriched_group
                  attr(seed_dt,"enriched_list")<-enriched_list
                  attr(seed_dt,"reference_group")<-reference_group
                  attr(seed_dt,"reference_list")<-reference_list
                  
                  attr(seed_dt,"RMS_method")<-rms_method
                  attr(seed_dt,"avg_rms_cutoff")<-avg_rms_cutoff
                  attr(seed_dt,"avg_resolution_cutoff")<-avg_res_cutoff
                  attr(seed_dt,"remove_duplicate_seqs")<-rem_duplicates
                  
                  saveRDS(seed_dt,file.path("mean temps",paste0(title,"_meanTempTable_",rms_method,".rds")))
            }else{seed_dt<-readRDS(file.path("mean temps",paste0(title,"_meanTempTable_",rms_method,".rds")))}
      }

      #Do a single alignment for all templates in seed_dt with the reference group and calculate hit numbers
      if(enriched_group=="cabs"){enriched_list<-unique(pdb_dt[cab==TRUE,]$PDBID)}
      row<-0
      alignment_dt<-seed_dt[,{
            row<-row+1
            seed_coords<-as.numeric(mget(names(seed_dt[,coord_indices,with=FALSE])))
            seed_xyz<-seed_coords[1:length(xyz_indices)]
            seed_dists<-seed_coords[(length(xyz_indices)+1):length(coord_indices)]
            print(paste0("Generating alignment results for ",key,": ",row," of ",nrow(seed_dt)))
            
            #create gpuMatrix objects if use.gpu==TRUE
            dRMSmat<-as.matrix(abn[,dRMS_indices,with=FALSE])
            rmsdmat<-as.matrix(abn[,xyz_indices,with=FALSE])

            #Calc distanceRMS and rmsd vectors
            dRMSvec<-rmsdF(seed_dists,dRMSmat,fit=FALSE,ncore=cpu_cores)
            rmsdvec<-rmsdF(seed_xyz,rmsdmat,fit=TRUE,ncore=cpu_cores)
            rmsdvec[is.nan(rmsdvec)]<-0 #set NaN superpositions to 0 (identical coordinates)
            
            #Find the indices of hits
            if(rms_method=="rmsd"){hits<-which(rmsdvec<=avg_rms_cutoff)}
            if(rms_method=="DRMS"){hits<-which(dRMSvec<=avg_rms_cutoff)}
            if(rms_method=="rmsd or DRMS"){hits<-which(rmsdvec<=avg_rms_cutoff | dRMSvec<=avg_rms_cutoff)}
            if(rms_method=="rmsd and DRMS"){hits<-which(rmsdvec<=avg_rms_cutoff & dRMSvec<=avg_rms_cutoff)}

            #Calculate enrichment statistics for the seed and reference groups
            ref_hits<-length(unique(abn[hits]$pdb))
            enriched_hits<-length(unique(abn[hits]$pdb[abn[hits]$pdb %in% enriched_list]))
            refN<-length(unique(abn$pdb))
            enrichedN<-length(unique(enriched_list[enriched_list %in% abn$pdb]))

            #send output mean template as list
            as.list(c(ref_hits,refN,enriched_hits,enrichedN))
      },by=key]
      
      #Compare enriched group to reference group
      setnames(alignment_dt,c("key","refHits","refN","enrichedHits","enrichedN"))
      alignment_dt[,enrichedP:=enrichedHits/enrichedN]
      alignment_dt[,refP:=refHits/refN]
      alignment_dt[,Zscore:=z.prop(x1=enrichedHits,x2=(refHits-enrichedHits),n1=enrichedN,n2=(refN-enrichedN))]
      chisquares<-t(apply(alignment_dt[,c("refHits","refN","enrichedHits","enrichedN"),with=FALSE],1,function(x){
            test<-prop.test(x=c(x["enrichedHits"],(x["refHits"]-x["enrichedHits"])),
                            n=c(x["enrichedN"],(x["refN"]-x["enrichedN"])),correct=FALSE)
            test2<-prop.test(x=c(x["enrichedHits"],(x["refHits"]-x["enrichedHits"])),
                            n=c(x["enrichedN"],(x["refN"]-x["enrichedN"])),correct=TRUE)
            c(test$statistic,test2$statistic)
      }))
      alignment_dt[,chisq:=chisquares[,1]]
      alignment_dt[,chisq_lowN:=chisquares[,2]]
      
      #Clean up alignment output table and save to file
      setkey(alignment_dt,key)
      setorder(alignment_dt,-Zscore)
      if(seed_group=="pdb" | seed_group=="model"){title<-paste0(seed_group,seed_list,collapse=",")}else{title<-seed_group}
      title<-paste0(title,"_EN-",enriched_group,"_REF-",reference_group)
      write.table(alignment_dt,
                  file=file.path("mean temps",paste0(title,"_enrichedTemps_",rms_method,".tab")),
                  sep="\t",row.names=FALSE,quote=FALSE)
}
