########################
##declare the function
#######################
pair.sets.filter<- function(ids, min.seqs){    
        require(stringr)
        if(missing(min.seqs)) min.seqs = 5 
        '%!in%' <- function(x, y) ! ('%in%'(x, y)) #define operator
        vector.is.empty <- function(x) return(length(x) ==0 ) #define function
        for(id in ids){ #loop every transmission pair 
                cluster.name <- id
                print(paste("Processing ", cluster.name, sep=""))
                info.gb <- read.csv(paste(cluster.name, "/", "info.gb.csv", sep=""), stringsAsFactors = F) #read data
                info <- read.csv(paste(cluster.name, "/", "info.csv", sep=""),stringsAsFactors = F)
                cluster <- unique(info$ClusterId) 
                source <- na.omit(info$PatientId[info$States=='source'])[[1]] 
                recipient <- na.omit(info$PatientId[info$States=='recipient'])[[1]] 
                set.id<- read.csv(paste(id, "/sets.csv", sep=""), row.names = NULL, stringsAsFactors = F) #load sets
                set<- which(set.id$length == (max(set.id$length))) #select set with longest genomic region
                if(length(set)>1) set<- max(set) #if there are hits for longest genomic region, pick the first one
                accessions <- unlist(strsplit(set.id$accessions[set], ",")) #retrieve accessions for set
                accessions.source <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id ==source)]))]
                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$patient_id ==recipient)]))]
                s.terminate <- r.terminate <- FALSE #logical vectors to evualate progress
                
                r.dfi <- unique(info.gb$days_from_infection[info.gb$accession_id %in% accessions.recipient])
                r.dfi[r.dfi==""]<- NA
                r.dfi<- unique(r.dfi)
                if(any(is.na(r.dfi) & length(r.dfi)>1)) r.dfi <- r.dfi[!is.na(r.dfi)]
                backup.control <- r.dfi
                backup.control.v <- numeric()
                if(length(r.dfi) > 1) { #size control
                        size.control <- c()
                        for(j in 1:length(r.dfi)){
                                backup.control.v <- c(backup.control.v, length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection==r.dfi[j])]))]))
                                if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection==r.dfi[j])]))]) >= min.seqs) {
                                        size.control <- c(size.control, r.dfi[j])
                                }
                        }
                        if(!is.null(size.control)) r.dfi <- size.control else {
                                r.dfi <- backup.control
                                backup.control.v <- c(backup.control.v, length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection%!in%r.dfi)]))]))
                                if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection%!in%r.dfi)]))]) >= min.seqs){
                                        accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection%!in%r.dfi)]))]
                                } else {
                                        r.dfi <- backup.control[which(backup.control.v == max(backup.control.v))]
                                        if(is.na(r.dfi)) r.dfi <- ""
                                }
                        }
                        
                }
                if(length(r.dfi)>1){ #prioritize numerical over character. pick lower
                        if(vector.is.empty(unique(na.omit(as.numeric(unlist( strsplit(as.character(r.dfi), "[^0-9]+"))))))) r.dfi<- sort(r.dfi)[1] else r.dfi<-min(unique(na.omit(as.numeric(unlist( strsplit(as.character(r.dfi), "[^0-9]+"))))))
                }
                if(!is.na(r.dfi)) {
                        accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_infection==r.dfi)]))]
                        r.terminate <- TRUE
                }
                #checks dfss
                if(!r.terminate){ # :D
                        r.dfs <- unique(info.gb$days_from_seroconversion[info.gb$accession_id %in% accessions.recipient])
                        r.dfs[r.dfs==""]<- NA
                        r.dfs<- unique(r.dfs)
                        if(any(is.na(r.dfs) & length(r.dfs)>1)) r.dfs <- r.dfs[!is.na(r.dfs)]
                        backup.control <- r.dfs
                        backup.control.v <- numeric()
                        if(length(r.dfs) > 1) { #size control
                                size.control <- c()
                                for(j in 1:length(r.dfs)){
                                        backup.control.v <- c(backup.control.v, length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_seroconversion==r.dfs[j])]))]))
                                        if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_seroconversion==r.dfs[j])]))]) >= min.seqs) {
                                                size.control <- c(size.control, r.dfs[j])
                                        }
                                }
                                if(!is.null(size.control)) r.dfs <- size.control else {
                                        r.dfs <- backup.control
                                        backup.control.v <- c(backup.control.v, length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_seroconversion%!in%r.dfs)]))]))
                                        if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_seroconversion%!in%r.dfs)]))]) >= min.seqs){
                                                r.dfs <- ""
                                        } else {
                                                r.dfs <- backup.control[which(backup.control.v == max(backup.control.v))]
                                                if(is.na(r.dfs)) r.dfs <- ""
                                        }
                                }
                                
                                
                        }
                        if(length(r.dfs)>1){ #prioritize numerical over character. pick lower
                                if(vector.is.empty(unique(na.omit(as.numeric(unlist(strsplit(as.character(r.dfs), "[^0-9]+"))))))) r.dfs<- sort(r.dfs)[1] else r.dfs <- r.dfs[which(unique(na.omit(as.numeric(unlist(strsplit(as.character(as.numeric(r.dfs)), "[^0-9]+")))))==min(unique(na.omit(as.numeric(unlist(strsplit(as.character(as.numeric(r.dfs)), "[^0-9]+")))))))]
                        }
                        if(!is.na(r.dfs)) {
                                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$days_from_seroconversion==r.dfs)]))]
                                r.terminate <- TRUE
                        }
                }
                
                #checks fiebig
                if(!r.terminate){
                        r.fiebig <- unique(info.gb$fiebig_stage[info.gb$accession_id %in% accessions.recipient])
                        r.fiebig[r.fiebig==""]<- NA
                        r.fiebig<- unique(r.fiebig)
                        if(any(is.na(r.fiebig) & length(r.fiebig)>1)) r.fiebig <- r.fiebig[!is.na(r.fiebig)]
                        backup.control <- r.fiebig
                        backup.control.v <- numeric()
                        if(length(r.fiebig) > 1) { #size control
                                size.control <- c()
                                for(j in 1:length(r.fiebig)){
                                        backup.control.v <- c(backup.control.v,length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage==r.fiebig[j])]))]))
                                        if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage==r.fiebig[j])]))]) >= min.seqs) {
                                                size.control <- c(size.control, r.fiebig[j])
                                        }
                                }
                                if(!is.null(size.control)) r.fiebig <- size.control else {
                                        r.fiebig <- backup.control
                                        backup.control.v <- c(backup.control.v, length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage%!in%r.fiebig)]))]))
                                        if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage%!in%r.fiebig)]))]) >= min.seqs){
                                                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage%!in%r.fiebig)]))]
                                        } else {
                                                r.fiebig <- backup.control[which(backup.control.v == max(backup.control.v))]
                                                if(is.na(r.fiebig)) r.fiebig <- ""
                                        }
                                }
                        }
                        if(length(r.fiebig)>1){ #prioritize numerical over character. pick lower
                                if(vector.is.empty(unique(na.omit(as.numeric(unlist( strsplit(as.character(r.fiebig), "[^0-9]+"))))))) r.fiebig<- sort(r.fiebig)[1] else r.fiebig<-min(unique(na.omit(as.numeric(unlist( strsplit(as.character(r.fiebig), "[^0-9]+"))))))
                        }
                        if(!is.na(r.fiebig)) {
                                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$fiebig_stage==r.fiebig)]))]
                                r.terminate <- TRUE
                        }
                }
                #check dates samples
                r.sample <- unique(info.gb$sampling.dates[info.gb$accession_id %in% accessions.recipient])
                r.sample[r.sample==""]<- NA
                if(vector.is.empty(r.sample)) r.sample<- NA
                if(any(is.na(r.sample) & length(r.sample)>1)) r.sample <- r.sample[!is.na(r.sample)]
                r.sample<- unique(r.sample)
                if(length(r.sample)==1) {
                        if(!is.na(r.sample)) {
                                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$sampling.dates==r.sample)]))] 
                                r.max <- max(unlist(lapply(lapply(r.sample, function(x) unlist(str_split(x, "-"))), length)))
                        } else r.max <- NA
                } else {
                        r.max <- max(unlist(lapply(lapply(r.sample, function(x) unlist(str_split(x, "-"))), length)))
                        r.sample<-r.sample[which(unlist(lapply(lapply(r.sample, function(x) unlist(str_split(x, "-"))), length)) ==r.max)]
                        if(r.max == 3) r.sample<- r.sample[which(as.Date(r.sample,'%d-%b-%Y') == min(as.Date(r.sample,'%d-%b-%Y')))] else if (r.max == 2) r.sample <- r.sample[which(as.Date(paste(01, "-", r.sample, sep=""), '%d-%b-%Y') == min(as.Date(paste(01, "-", r.sample, sep=""), '%d-%b-%Y')))] else if (r.max == 1) r.sample <- min(as.integer(r.sample))
                        accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$sampling.dates==r.sample)]))]
                }
                

                #if all NA, check years
                if(all(unlist(lapply(c("r.sample",  "r.fiebig", "r.dfs", "r.dfi"), function(x) exists(x))))){
                        if(all(is.na(c(r.sample, r.fiebig, r.dfs, r.dfi)))){
                                r.year <- unique(info.gb$sampling_year[info.gb$accession_id %in% accessions.recipient])
                                r.year[r.year==""]<- NA
                                r.year<- unique(r.year)
                                if(any(is.na(r.year) & length(r.year)>1)) r.year <- r.year[!is.na(r.year)]
                                if(length(r.year) > 1) { #size control
                                        size.control <- c()
                                        for(j in 1:length(r.year)){
                                                if(length(accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$sampling_year==r.year[j])]))]) >= min.seqs) {
                                                        size.control <- c(size.control, r.year[j])
                                                }
                                        }
                                        r.year <- size.control
                                }
                                r.year<-min(as.numeric(as.character(r.year))) 
                                accessions.recipient <- accessions[which(accessions %in% as.character(info.gb$accession_id[which(info.gb$accession_id%in%accessions.recipient & info.gb$sampling_year==r.year)]))]
                        }
                        
                }

                
                accessions <- c(accessions.source, accessions.recipient) 
                accession.min <- min(info.gb$pos_ini[which(info.gb$accession_id%in%accessions)])
                accession.max <- min(info.gb$pos_end[which(info.gb$accession_id%in%accessions)])
                empirical <- data.frame(
                        "cluster_name"=cluster.name, "cluster_id"=cluster, "donor.id"=source, "recipient.id"=recipient,
                        "subtype"=set.id$subtype[[set]],
                        "ref"=set.id$ref[[set]],
                        "counts" = length(accessions), "source.n"=length(accessions.source), "recipient.n"=length(accessions.recipient),
                        "min" = accession.min, "max" = accession.max,
                        "length"= accession.max - accession.min,
                        "source.dfi"=paste(unique(info.gb$days_from_infection[info.gb$accession_id %in% accessions.source]), collapse=" "),
                        "source.dfs"=paste(unique(info.gb$days_from_seroconversion[info.gb$accession_id %in% accessions.source]), collapse=" "),
                        "source.y"=paste(unique(info.gb$sampling_year[info.gb$accession_id %in% accessions.source]), collapse=" "),
                        "source.fiebig"=paste(unique(info.gb$fiebig_stage[info.gb$accession_id %in% accessions.source]), collapse=" "),
                        "source.sampling.dates"=paste(unique(info.gb$sampling.dates[info.gb$accession_id %in% accessions.source]), collapse=" "),
                        "recipient.dfi"=if (length(unique(info.gb$days_from_infection[info.gb$accession_id %in% accessions.recipient]))==1) unique(info.gb$days_from_infection[info.gb$accession_id %in% accessions.recipient]) else unique(as.vector(na.omit(info.gb$days_from_infection[info.gb$accession_id %in% accessions.recipient]))),
                        "recipient.dfs"=if (length(unique(info.gb$days_from_seroconversion[info.gb$accession_id %in% accessions.recipient]))==1) unique(info.gb$days_from_seroconversion[info.gb$accession_id %in% accessions.recipient]) else unique(as.vector(na.omit(info.gb$days_from_seroconversion[info.gb$accession_id %in% accessions.recipient]))),
                        "recipient.y"=if (length(unique(info.gb$sampling_year[info.gb$accession_id %in% accessions.recipient]))==1) unique(info.gb$sampling_year[info.gb$accession_id %in% accessions.recipient]) else unique(as.vector(na.omit(info.gb$sampling_year[info.gb$accession_id %in% accessions.recipient]))),
                        "recipient.fiebig"= if (length(unique(info.gb$fiebig_stage[info.gb$accession_id %in% accessions.recipient]))==1) unique(info.gb$fiebig_stage[info.gb$accession_id %in% accessions.recipient]) else unique(as.vector(na.omit(info.gb$fiebig_stage[info.gb$accession_id %in% accessions.recipient]))),
                        "recipient.sampling.dates"=if (length(unique(info.gb$sampling.dates[info.gb$accession_id %in% accessions.recipient]))==1) unique(info.gb$sampling.dates[info.gb$accession_id %in% accessions.recipient]) else unique(as.vector(na.omit(info.gb$sampling.dates[info.gb$accession_id %in% accessions.recipient]))),
                        "accessions"= paste(accessions, collapse=","), 
                        "accessions.source"= paste(accessions.source, collapse=","), 
                        "accessions.rec"= paste(accessions.recipient, collapse=","),
                        stringsAsFactors = FALSE)
                empirical[] <- lapply(empirical, as.character)
                empirical <<- empirical
                write.csv(empirical, paste(id, "/set.csv", sep=""), row.names = FALSE)
        }
}
#######################
##use the function:
#######################
#setwd("/Multi.regions") #change your working directory if needed

pair.sets.filter(ids = ids.multi, min.seqs=5)
