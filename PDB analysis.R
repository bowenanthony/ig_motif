# ========== PDB Analysis Version 1.0 ============ #
# Anthony Bowen
# anthony.bowen@med.einstein.yu.edu
# Laboratory of Arturo Casadevall
# Albert Einstein College of Medicine, Bronx NY
# Last Modified 9/9/16;  created, modified to automatically update PDB lists with new downloads and
#                         to re-download current PDB lists if desired.
#                         classify all pdb sequences with imgt reference genes to ID antibodies
#                         pull structure resolution from pdb files

# This script will identify all antibodies in the PDB and generate other lists based on other criteria

#load needed packages
library("data.table")
library("bio3d")
library("lubridate")
library("xlsx")
library("XML")
library("Biostrings")
#library("ggplot2")
#library("rhdf5")

### Setup
initial.dir<-getwd() # store the current directory
temp.dir<-  "B:/Tony/Documents/R Working Directory/IgMotif_temp" #2500K
#file.path("~","R Working Directory","IgMotif_temp") #Macbook Pro
pdb_db<-    "B:/Tony/Documents/R Working Directory/IgMotif_temp/PDB_DB" #2500K
#file.path(temp.dir,"PDB_DB") #Macbook Pro
pdb.dir<-   "B:/Tony/Documents/Casadevall Lab/Non-redundant PDB files" #2500K
#file.path("~","Documents","Casadevall Lab","Non-redundant PDB files") #Macbook Pro
src.dir<-   "B:/Tony/Google Drive/R Programming/IgMotif" #2500K
#file.path("~","Google Drive","R Programming","IgMotif") #Macbook Pro
setwd(src.dir) # change to the src directory
options(scipen=10) # number of digits before scientific notation is used

### Script options
update_pdbs<-FALSE #fetch new PDBs and remove obsolete ones from the RCSB, & update tables
update_pdb_labels<-TRUE #whether to update pdb_dt with new Cab or PDB Ab ID labels
Lscore<- -350 #minimum alignment score to classify an L chain
Hscore<- -350 #minimum alignment score to classify an H chain

### Custom Functions###

# IMGT reference comparison
# input: list of PDB IDs for files downloaded in a certain directory
# output: pdb_dt of sequence alignment scores & other IMGT info for protein chains in PDB files
imgt_ref_compare<-function(pdb_list,pdb_dir){
      
      pdb_list<-tolower(pdb_list)
      pdb_dt<-data.table(rep(pdb_list,2))
      setnames(pdb_dt,"V1","PDBID")
      pdb_dt[,chain:=c(rep("L",length(pdb_list)),rep("H",length(pdb_list)))]
      pdb_dt[,seq:=rep("",length(pdb_dt$PDBID))]
      pdb_dt[,family:=rep("",length(pdb_dt$PDBID))]
      pdb_dt[,species:=rep("",length(pdb_dt$PDBID))]
      pdb_dt[,refgene:=rep("",length(pdb_dt$PDBID))]
      pdb_dt[,score:=as.numeric(rep(NA,length(pdb_dt$PDBID)))]
      pdb_dt[,resolution:=as.numeric(rep(NA,length(pdb_dt$PDBID)))]
      empty.files<-character()
      empty.atoms<-character()
      no.files<-character()
      
      for (id in pdb_list){
            file<-file.path(pdb_dir,paste0(id,".pdb"))
            if (file.exists(file)){
                  if (file.size(file)>0){
                        txt<-read.csv(file,colClasses="character",row.names=NULL,header=FALSE)[[1]]
                        sumtest<-sum(sapply(txt,function(x){grepl("^ATOM ",x)}))
                        if(sumtest<2){empty.atoms<-c(empty.atoms,id)}
                  }
                  else{sumtest<-0; empty.files<-c(empty.files,id)}
                  if((sumtest>1) ){
                        pdbs<-bio3d::read.pdb(file)
                  }else{warning(paste0(id," cannot be read due to sumtest ",sumtest));next}
            }else{warning(paste0(file," does not exist"));no.files<-c(no.files,id);next}
            chains<-unique(as.factor(names(pdbs$seqres)))
            seqs<-as.factor(pdbs$seqres)
            table<-data.table(chain=chains)
            table[,seq:=sapply(chain,function(x){paste0(aa321(seqs[names(seqs)==x]),collapse="")})]
            table[,len:=sapply(seq,nchar)] #add column for chain length
            table[,duplicate:=duplicated(seq)] #label chains with same sequence
            table[,short:=len<=min(nchar(imgt_con_dt$seq))-20] #flag sequences that are too short
            table<-table[duplicate==FALSE & short == FALSE,] #remove duplicate and short chains
            
            #find resolution of structures if available
            if(length(as.numeric(substring(txt[grepl("REMARK   2 RESOLUTION.",txt)],25,30)))==1){
                  pdb_dt[PDBID==id,resolution:=
                               as.numeric(substring(txt[grepl("REMARK   2 RESOLUTION.",txt)],25,30))]
            }
            
            #Do pairwise alignment of chain sequences with reference sequences to guess species & family
            aliscores<-numeric()
            alispecies<-character()
            alifamily<-character()
            aliname<-character()
            for (seq in table[,seq]){
                  if (nchar(seq)>max(nchar(imgt_con_dt$seq)+20)){
                        seq<-substr(seq,1,max(nchar(imgt_con_dt$seq)+20))
                  }
                  ali<-pairwiseAlignment(imgt_allseq$seq,seq,type="global",scoreOnly=TRUE)
                  hit_index<-which.max(ali)
                  aliscores<-c(aliscores,max(ali))
                  alifamily<-c(alifamily,imgt_allseq[which.max(ali),family])
                  alispecies<-c(alispecies,imgt_allseq[which.max(ali),species])
                  aliname<-c(aliname,imgt_allseq[which.max(ali),name])
            }
            table[,family:=alifamily]
            table[,species:=alispecies]
            table[,name:=aliname]
            table[,score:=aliscores]
            
            # Assign sequence info from max alignment score of chains with reference seqs
            Lindex<-which.max(table[which(table$family=="IgLV"|table$family=="IgKV"),score])
            Ltable<-table[which(table$family=="IgLV"|table$family=="IgKV"),][Lindex]
            Hindex<-which.max(table[which(table$family=="IgHV"),score])
            Htable<-table[which(table$family=="IgHV"),][Hindex]
            
            # Fill in info to PDB pdb_dt table
            if(nrow(Htable)==1){
                  pdb_dt[PDBID==id & chain=="H",seq:=Htable$seq]
                  pdb_dt[PDBID==id & chain=="H",family:=Htable$family]
                  pdb_dt[PDBID==id & chain=="H",species:=Htable$species]
                  pdb_dt[PDBID==id & chain=="H",refgene:=Htable$name]
                  pdb_dt[PDBID==id & chain=="H",score:=Htable$score]
            }
            if(nrow(Ltable)==1){
                  pdb_dt[PDBID==id & chain=="L",seq:=Ltable$seq]
                  pdb_dt[PDBID==id & chain=="L",family:=Ltable$family]
                  pdb_dt[PDBID==id & chain=="L",species:=Ltable$species]
                  pdb_dt[PDBID==id & chain=="L",refgene:=Ltable$name]
                  pdb_dt[PDBID==id & chain=="L",score:=Ltable$score]
            }
      }
      pdb_dt[,len:=nchar(pdb_dt[,seq])] #add lengths of sequences
      pdb_dt[,name:=paste0(PDBID,"-",chain)] #add unique seq name
      return(pdb_dt)
}


### Script sections

download_errors<-character()# list of PDBs that could not be downloaded
new_downloads<-character()# list of new PDBs that were just downloaded
if(update_pdbs){
      
      ### download new lists of current and obsolete PDB IDs from the RCSB
      
      current_url<-"http://www.rcsb.org/pdb/rest/getCurrent"
      download.file("http://www.rcsb.org/pdb/rest/getCurrent",
                    file.path("PDB downloads",paste0(today(),"_currentPDBs.xml")))
      current_pdb<-xmlTreeParse(file.path("PDB downloads",
                                          paste0(today(),"_currentPDBs.xml")),useInternal=TRUE)
      unlink(file.path("PDB downloads",paste0(today(),"_currentPDBs.xml")))
      current_pdb<-data.table(xpathSApply(current_pdb,"//PDB",xmlAttrs))
      setnames(current_pdb,"V1","PDBID")
      attr(current_pdb,"accession_date")<-now()
      current_pdb[,mod_date:=""]
      
      obsolete_url<-"http://www.rcsb.org/pdb/rest/getObsolete"
      download.file("http://www.rcsb.org/pdb/rest/getObsolete",
                    file.path("PDB downloads",paste0(today(),"_obsoletePDBs.xml")))
      obsolete_pdb<-xmlTreeParse(file.path("PDB downloads",
                                           paste0(today(),"_obsoletePDBs.xml")),useInternal=TRUE)
      unlink(file.path("PDB downloads",paste0(today(),"_obsoletePDBs.xml")))
      obsolete_pdb<-data.table(xpathSApply(obsolete_pdb,"//PDB",xmlAttrs))
      setnames(obsolete_pdb,"V1","PDBID")
      attr(obsolete_pdb,"accession_date")<-now()
      obsolete_pdb[,mod_date:=""]
      saveRDS(obsolete_pdb,file.path(pdb_db,paste0(today(),"_obsoletePDBs.rds")))
      
      ###   Get general info about each PDB including method, resolution, deposit date, entities,
      #     chains, etc...
      # add columns to current_pdb table for data to be added
      
      # split searches up by certain number of PDB ids
      max_pdb_number<-1000
      max_search_number<-ceiling(length(current_pdb$PDBID)/max_pdb_number)
      pdb_i<-1
      for (search in c(1:max_search_number)){
            if(search==max_search_number){max_pdb_number<-length(current_pdb$PDBID)}
            pdb_list<-current_pdb$PDBID[pdb_i:(pdb_i+max_pdb_number-1)]
            mol_d<-paste0("http://www.rcsb.org/pdb/rest/describeMol?structureId=",
                          paste(pdb_list,collapse=","))
            pdb_d<-paste0("http://www.rcsb.org/pdb/rest/describePDB?structureId=",
                          paste(pdb_list,collapse=","))
            
            #download XML files
            download.file(mol_d, file.path(pdb_db,paste0(search,"_molD.xml")))
            download.file(pdb_d, file.path(pdb_db,paste0(search,"_pdbD.xml")))
            
            pdb_i<-pdb_i+max_pdb_number
      }
      
      # parse and add data to current_pdb table, then delete xml files
      pdb_i<-1
      for (search in c(1:max_search_number)){
            if(search==max_search_number){max_pdb_number<-length(current_pdb$PDBID)}
            pdb_list<-current_pdb$PDBID[pdb_i:(pdb_i+max_pdb_number-1)]
            mol_data<-xmlTreeParse(file.path(pdb_db, paste0(search,"_molD.xml")),useInternal=TRUE)
            pdb_data<-xmlTreeParse(file.path(pdb_db, paste0(search,"_pdbD.xml")),useInternal=TRUE)
            
            
            #data.table(xpathSApply(current_pdb,"//PDB",xmlAttrs))
            
            unlink(file.path(pdb_db, paste0(search,"_molD.xml")))
            unlink(file.path(pdb_db, paste0(search,"_pdbD.xml")))
            pdb_i<-pdb_i+max_pdb_number
      }
      
      
      saveRDS(current_pdb,file.path(pdb_db,paste0(today(),"_currentPDBs.rds")))
      
      ### check for most recent PDB list in the folder and compare to the one just downloaded
      
      recent_cdate<-sort(ymd(sub("_currentPDBs.rds","",list.files(pdb_db)
                                 [grep("currentPDBs.rds",list.files(pdb_db))])),
                         decreasing = TRUE)[2]
      recent_odate<-sort(ymd(sub("_obsoletePDBs.rds","",list.files(pdb_db)
                                 [grep("obsoletePDBs.rds",list.files(pdb_db))])),
                         decreasing=TRUE)[2]
      if(!is.na(recent_cdate)){
            recent_current<-readRDS(file.path(pdb_db, paste0(recent_cdate,"_currentPDBs.rds")))
      }else{recent_current<-data.table()}
      if(!is.na(recent_odate)){
            recent_obsolete<-readRDS(file.path(pdb_db, paste0(recent_odate,"_obsoletePDBs.rds")))
      }else{recent_obsolete<-data.table()}
      
      add_current<-current_pdb[!(current_pdb$PDBID %in% recent_current$PDBID)] #added PDBs
      rem_current<-recent_current[!(recent_current$PDBID %in% current_pdb$PDBID)] #removed PDBs
      if (length(add_current$PDBID)>0 | length(rem_current$PDBID)>0){
            print("Additional PDBs to be added:")
            print(add_current)
            print("PDBs to be removed:")
            print(rem_current)
            
      }else{
            unlink(file.path(pdb_db,paste0(today(),"_currentPDBs.rds")))
            current_pdb<-recent_current
      }
      
      add_obsolete<-obsolete_pdb[!(obsolete_pdb$PDBID %in% recent_obsolete$PDBID)] #added PDBs
      rem_obsolete<-recent_obsolete[!(recent_obsolete$PDBID %in% obsolete_pdb$PDBID)] #removed PDBs
      if (length(add_obsolete$PDBID)>0 | length(rem_obsolete$PDBID)>0){
            print("New obsolete PDBs:")
            print(add_obsolete)
            print("No longer obsolete PDBs:")
            print(rem_obsolete)
      }else{
            unlink(file.path(pdb_db,paste0(today(),"_obsoletePDBs.rds")))
            obsolete_pdb<-recent_obsolete
      }
      rm(recent_obsolete)
      rm(recent_current)
      
      ### delete or download PDB files from RCSB as needed to update list
      
      for (id in add_current$PDBID){
            tryCatch({
                  if(!file.exists(file.path(pdb_db,"PDB_files",paste0(id,".pdb")))){
                        download.file(paste0("http://files.rcsb.org/view/",id,".pdb"),
                                      file.path(pdb_db,"PDB_files",paste0(id,".pdb")))
                        current_pdb[PDBID==id,mod_date:=as.character(now())]
                        new_downloads<-c(new_downloads,id)
                  }else{print(paste0("PDB ",id," already downloaded"))}
            }, error = function(e){
                  cat("ERROR ",id,":",conditionMessage(e),"\n")
                  current_pdb[PDBID==id,mod_date:=paste0("ERROR: ",as.character(now()))]
                  download_errors<<-c(download_errors,id)})
      }
      for (id in rem_current$PDBID){
            unlink(file.path(pdb_db,"PDB_files",paste0(id,".pdb")))
            obsolete_pdb[PDBID==id,mod_date:=as.character(now())]
      }
}

if(!update_pdbs){
      recent_cdate<-sort(ymd(sub("_currentPDBs.rds","",list.files(pdb_db)
                                 [grep("currentPDBs.rds",list.files(pdb_db))])),
                         decreasing = TRUE)[1]
      recent_odate<-sort(ymd(sub("_obsoletePDBs.rds","",list.files(pdb_db)
                                 [grep("obsoletePDBs.rds",list.files(pdb_db))])),
                         decreasing=TRUE)[1]
      if(!is.na(recent_cdate)){
            current_pdb<-readRDS(file.path(pdb_db, paste0(recent_cdate,"_currentPDBs.rds")))
      }else{warning("no current PDB list file found")}
      if(!is.na(recent_odate)){
            obsolete_pdb<-readRDS(file.path(pdb_db, paste0(recent_odate,"_obsoletePDBs.rds")))
      }else{warning("no obsolete PDB list file found")}
}

### check that all PDBs are downloaded
for (id in current_pdb$PDBID){
      tryCatch({
            if(!file.exists(file.path(pdb_db,"PDB_files",paste0(id,".pdb")))){
                  download.file(paste0("http://files.rcsb.org/view/",id,".pdb"),
                                file.path(pdb_db,"PDB_files",paste0(id,".pdb")))
                  current_pdb[PDBID==id,mod_date:=as.character(now())]
                  new_downloads<-c(new_downloads,id)
            }else{print(paste0("PDB ",id," already downloaded"))}
      }, error = function(e){cat("ERROR ",id,":",conditionMessage(e),"\n")
            current_pdb[PDBID==id,mod_date:=paste0("ERROR: ",as.character(now()))]
            download_errors<<-c(download_errors,id)})
}

### save updated lists with newly modified dates
saveRDS(current_pdb,file.path(pdb_db,paste0(
      date(attr(current_pdb,"accession_date")),"_currentPDBs.rds")))
saveRDS(obsolete_pdb,file.path(pdb_db,paste0(
      date(attr(obsolete_pdb,"accession_date")),"_obsoletePDBs.rds")))

### delete files that had download errors and filter pdb_list to remove error files
for (id in download_errors){
      unlink(file.path(pdb_db,"PDB_files",paste0(id,".pdb"))) #delete files with errors
}
pdb_list<-current_pdb[!(PDBID %in% download_errors)]


### read in IMGT files and create table
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
      rm(seq_dt)
      
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

#Define lists of known catalytic antibodies or other groups by PDB ID
pdbids<-data.table(read.xlsx2("PDBid_lists.xlsx", sheetName = "ids"))
cabs<-tolower(as.character(na.omit(pdbids$cabs_chk)))
cabs<-cabs[nchar(cabs)>0]
pdb_abs<-tolower(as.character(na.omit(pdbids$abs8_16)))
pdb_abs<-pdb_abs[nchar(pdb_abs)>0]
matches<-tolower(as.character(na.omit(pdbids$X4hdi_mean)))
matches<-matches[nchar(matches)>0]
hydrolase<-tolower(as.character(na.omit(pdbids$hydrolase9.16)))
hydrolase<-hydrolases[nchar(hydrolase)>0]
serprot<-tolower(as.character(na.omit(pdbids$serprot9.16)))
serprot<-serprot[nchar(serprot)>0]
cysprot<-tolower(as.character(na.omit(pdbids$cysprot9.16)))
cysprot<-cysprot[nchar(cysprot)>0]
thrprot<-tolower(as.character(na.omit(pdbids$thrprot9.16)))
thrprot<-thrprot[nchar(thrprot)>0]
aspprot<-tolower(as.character(na.omit(pdbids$aspprot9.16)))
aspprot<-aspprot[nchar(aspprot)>0]
metprot<-tolower(as.character(na.omit(pdbids$metprot9.16)))
metprot<-metprot[nchar(metprot)>0]
glycosylase<-tolower(as.character(na.omit(pdbids$glycosylase9.16)))
glycosylase<-glycosylase[nchar(glycosylase)>0]

imgt_abs<-data.table(read.xlsx2("PDBid_lists.xlsx", sheetName = "IMGT_pdb_Igs",stringsAsFactors=FALSE,
                                colClasses=c(rep("character",7),"numeric","character")))
imgt_abs[,PDB.release.date:=as.Date(imgt_abs$PDB.release.date,"%d-%B-%y")]
#all_abs<-tolower(as.character(na.omit(pdbids$abs7_16)))
#all_abs<-all_abs[nchar(all_abs)>0]

### get sequences from PDB files, and align with IMGT ref genes to classify,
#     if pdb_dt file doesn't exist already
if (!file.exists(file.path("output","pdb_dt.RDS"))){
      pdb_dt<-imgt_ref_compare(pdb_list = pdb_list$PDBID, pdb_dir = file.path(pdb_db,"PDB_files"))
      pdb_dt[,cab:=PDBID %in% cabs]
      pdb_dt[,match:=PDBID %in% matches]
      pdb_dt[,pdb_ab:=pdb_dt$PDBID %in% pdb_abs]
      pdb_dt[,imgt_ab:=PDBID %in% imgt_abs$IMGT.entry.ID]
      
      #make combined seq table & ID duplicates
      Ldt<-pdb_dt[chain=="L",c("PDBID","seq"),with=FALSE]
      Hdt<-pdb_dt[chain=="H",c("PDBID","seq"),with=FALSE]
      setkey(Ldt,PDBID)
      setkey(Hdt,PDBID)
      mergedt<-merge(Ldt, Hdt,by="PDBID")
      mergedt[,LHseq:=paste0(seq.x,seq.y)]
      mergedt[,duplicate_seq:=duplicated(mergedt[,LHseq])] #label any rows that have duplicate H&L sequence
      pdb_dt[,duplicate_seq:=mergedt[PDBID==PDBID,"duplicate_seq",with=FALSE]]
      #pdb_dt[,duplicate_seq:=duplicated(pdb_dt[,seq])] #label any rows that have duplicate chain sequence
      
      #label protease rows
      pdb_dt[,hydrolase:=pdb_dt$PDBID %in% hydrolase]
      pdb_dt[,serprot:=pdb_dt$PDBID %in% serprot]
      pdb_dt[,cysprot:=pdb_dt$PDBID %in% cysprot]
      pdb_dt[,thrprot:=pdb_dt$PDBID %in% thrprot]
      pdb_dt[,aspprot:=pdb_dt$PDBID %in% aspprot]
      pdb_dt[,metprot:=pdb_dt$PDBID %in% metprot]
      pdb_dt[,glycosylase:=pdb_dt$PDBID %in% glycosylase]
      
      saveRDS(pdb_dt,file=file.path("output","pdb_dt.RDS"))
}else{
      pdb_dt<-readRDS(file=file.path("output","pdb_dt.RDS"))# open existing file
      if (length(new_downloads)>0){
            # create table for new downloads
            new_dt<-imgt_ref_compare(pdb_list = new_downloads, pdb_dir = file.path(pdb_db,"PDB_files"))
            new_dt[,cab:=PDBID %in% cabs]
            new_dt[,match:=PDBID %in% matches]
            new_dt[,duplicate_seq:=duplicated(new_dt[,seq])] #label any rows that have duplicate chain sequence
            
            #label protease rows
            new_dt[,hydrolase:=new_dt$PDBID %in% hydrolase]
            new_dt[,serprot:=new_dt$PDBID %in% serprot]
            new_dt[,cysprot:=new_dt$PDBID %in% cysprot]
            new_dt[,thrprot:=new_dt$PDBID %in% thrprot]
            new_dt[,aspprot:=new_dt$PDBID %in% aspprot]
            new_dt[,metprot:=new_dt$PDBID %in% metprot]
            new_dt[,glycosylase:=new_dt$PDBID %in% glycosylase]
            
            #add new download data to existing list
            pdb_dt<-rbind(pdb_dt,new_dt) #combine tables
            
            #make combined seq table & ID duplicates
            Ldt<-pdb_dt[chain=="L",c("PDBID","seq"),with=FALSE]
            Hdt<-pdb_dt[chain=="H",c("PDBID","seq"),with=FALSE]
            setkey(Ldt,PDBID)
            setkey(Hdt,PDBID)
            mergedt<-merge(Ldt, Hdt,by="PDBID")
            mergedt[,LHseq:=paste0(seq.x,seq.y)]
            mergedt[,duplicate_seq:=duplicated(mergedt[,LHseq])] #label any rows that have duplicate H&L sequence
            pdb_dt[,duplicate_seq:=mergedt[PDBID==PDBID,"duplicate_seq",with=FALSE]]
            
            #pdb_dt[,duplicate_seq:=duplicated(pdb_dt[,seq])] #label any rows that have duplicate chain sequence
            
            saveRDS(pdb_dt,file=file.path("output","pdb_dt.RDS"))
      }
      if(update_pdb_labels){
            pdb_dt<-imgt_ref_compare(pdb_list = pdb_list$PDBID, pdb_dir = file.path(pdb_db,"PDB_files"))
            pdb_dt[,cab:=PDBID %in% cabs]
            pdb_dt[,match:=PDBID %in% matches]
            pdb_dt[,pdb_ab:=pdb_dt$PDBID %in% pdb_abs]
            pdb_dt[,imgt_ab:=PDBID %in% imgt_abs$IMGT.entry.ID]
            
            #make combined seq table & ID duplicates
            Ldt<-pdb_dt[chain=="L",c("PDBID","seq"),with=FALSE]
            Hdt<-pdb_dt[chain=="H",c("PDBID","seq"),with=FALSE]
            setkey(Ldt,PDBID)
            setkey(Hdt,PDBID)
            mergedt<-merge(Ldt, Hdt,by="PDBID")
            mergedt[,LHseq:=paste0(seq.x,seq.y)]
            mergedt[,duplicate_seq:=duplicated(mergedt[,LHseq])] #label any rows that have duplicate H&L sequence
            pdb_dt[,duplicate_seq:=mergedt[PDBID==PDBID,"duplicate_seq",with=FALSE]]
            
            #pdb_dt[,duplicate_seq:=duplicated(pdb_dt[,seq])] #label any rows that have duplicate chain sequence
            
            #label protease rows
            pdb_dt[,hydrolase:=pdb_dt$PDBID %in% hydrolase]
            pdb_dt[,serprot:=pdb_dt$PDBID %in% serprot]
            pdb_dt[,cysprot:=pdb_dt$PDBID %in% cysprot]
            pdb_dt[,thrprot:=pdb_dt$PDBID %in% thrprot]
            pdb_dt[,aspprot:=pdb_dt$PDBID %in% aspprot]
            pdb_dt[,metprot:=pdb_dt$PDBID %in% metprot]
            pdb_dt[,glycosylase:=pdb_dt$PDBID %in% glycosylase]
            
            saveRDS(pdb_dt,file=file.path("output","pdb_dt.RDS"))
      }
}
      



########## Misc code

###   Create HDF5 database if it doesn't exist, and append new PDB lists if selected above
#     stopped using HDF5 because datasets cannot be easily deleted or moved after added to file
# pdb_db<-file.path(temp.dir,"fullPDB.h5")
# if (!file.exists(file.path(temp.dir,"fullPDB"))){
#       if(!update_pdbs){warning("no PDB database found, please change 'update_pdbs' to TRUE")}
#       h5createFile(file.path(temp.dir,"fullPDB.h5"))
#       h5createGroup(pdb_db,"PDB_lists")
#       h5createGroup(pdb_db,"PDB_lists/old")
#       h5createGroup(pdb_db,"PDB_files")
#       h5write(current_pdb,pdb_db,name="PDB_lists/current",write.attributes=TRUE)
#       h5write(obsolete_pdb,pdb_db,name="PDB_lists/obsolete",write.attributes=TRUE)
# }else{
#       if(update_pdbs){
#             # move older pdb lists to old grouping and overwrite with the new lists if there are any changes
#             old_current_pdb<-h5read(pdb_db,"PDB_lists/current",read.attributes=TRUE)
#             old_obsolete_pdb<-h5read(pdb_db,"PDB_lists/obsolete",read.attributes=TRUE)
#             add_current<-current_pdb[!(current_pdb$PDBID %in% old_current_pdb$PDBID)] #added PDBs
#             rem_current<-old_current_pdb[!(old_current_pdb$PDBID %in% current_pdb$PDBID)] #removed PDBs
#             if (length(add_current)>0 | length(rem_current)>0){
#                   print("Additional PDBs to be added:")
#                   print(add_current)
#                   print("PDBs to be removed:")
#                   print(rem_current)
#                   h5write(old_current_pdb,pdb_db,
#                           name=paste0("PDB_lists/old/current_",attr(old_current_pdb,"accession_date")),
#                           write.attributes=TRUE)
#                   h5write(current_pdb,pdb_db,name="PDB_lists/current",write.attributes=TRUE)
#             }
#             add_obsolete<-obsolete_pdb[!(obsolete_pdb$PDBID %in% old_obsolete_pdb$PDBID)] #added PDBs
#             rem_obsolete<-old_obsolete_pdb[!(old_obsolete_pdb$PDBID %in% obsolete_pdb$PDBID)] #removed PDBs
#             if (length(add_obsolete)>0 | length(rem_obsolete)>0){
#                   print("New obsolete PDBs:")
#                   print(add_obsolete)
#                   print("PDBs removed from the obsolete list:")
#                   print(rem_obsolete)
#                   h5write(old_obsolete_pdb,pdb_db,
#                           name=paste0("PDB_lists/old/obsolete_",attr(old_current_pdb,"accession_date")),
#                           write.attributes=TRUE)
#                   h5write(obsolete_pdb,pdb_db,name="PDB_lists/obsolete",write.attributes=TRUE)
#             }
#       }else{
#             current_pdb<-h5read(pdb_db,"PDB_lists/current",read.attributes=TRUE)
#             obsolete_pdb<-h5read(pdb_db,"PDB_lists/obsolete",read.attributes=TRUE)
#       }
# }
# 
# ###   Read in PDB files from web, save data to HDF5 database file if they don't already exist
# db<-H5Fopen(pdb_db)
# for (id in current_pdb$PDBID){
#       if(!H5Lexists(db, paste0("PDB_files/",id))){
#             pdb_data<-readLines(con=paste0("http://files.rcsb.org/view/",id,".pdb"))
#             attr(pdb_data,"downloaded")<-as.character(now())
#             db$PDB_files[id]<-pdb_data
#             h5write(pdb_data, file = pdb_db, name = paste0("PDB_files/",id), write.attributes=TRUE)
#       }
# }
# H5Fclose(db)
