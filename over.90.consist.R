#create database of consistent posterior results##

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/HPGH/Atkins Project/TransmissionPairs-master/Clusters/Clusters.90.new")

library(plyr)
library(dplyr)

#filter those with consistent results over 90

ann.90.2 <- ann.consist %>% filter(freq >89.9)


#create list of all annoated files which have 90% or over consistent results
ann.90.files <- list.files()


#create file of all ANN logs from pairs which are 90% and over consistent
ann.90 <- ldply(list.files(), read.table, header=TRUE)

#save as csv to use in other scripts
write.csv(ann.90, "ann.90.csv")

#filter out only the consistent results

ann.90.consist <- ann.90 %>% select(par.direction, lineages.recipient, LANLdb_cluster_name) %>% filter (par.direction == "consistent") 

#test with one to see shape
#ann.90.consist1 <- ann.90 %>% select(par.direction, lineages.recipient, LANLdb_cluster_name) %>% filter (par.direction == "consistent") %>%  filter(LANLdb_cluster_name == "MC11-56")

#list the cluster names in the 90 and over group
ann.90.names <- ann.90 %>% select(par.direction,LANLdb_cluster_name) %>% filter (par.direction == "consistent") %>% group_by(LANLdb_cluster_name) %>% tally() 

#add transmission route - this is done manually from reading papers
ann.90.consist.add <- ann.90 %>% select(par.direction, LANLdb_cluster_name, lineages.recipient) %>% filter (par.direction == "consistent") %>% mutate(Transmission = case_when(LANLdb_cluster_name == "CASE10" ~ "BM",        
                                                                                                                                                                               LANLdb_cluster_name == "HN88" ~ "Intra",                                    LANLdb_cluster_name == "M1001-P1024" ~ "Intra",                           LANLdb_cluster_name == "M1007-P1046" ~ "Intra",                           LANLdb_cluster_name == "MB_pair_i20-m20" ~ "BM", 
LANLdb_cluster_name == "MB_pair_i38-m38" ~ "BM",                          LANLdb_cluster_name == "MB_pair_i40-m40" ~ "BM", 
LANLdb_cluster_name == "MC-DU" ~ "IUT", 
LANLdb_cluster_name == "MC-IP2" ~ "Intra or BM", 
LANLdb_cluster_name == "MC-IU1" ~ "IUT",                                  LANLdb_cluster_name == "MC-IU2" ~ "IUT", 
LANLdb_cluster_name == "MC-IU4" ~ "IUT",                                  LANLdb_cluster_name == "MC10-1" ~ "IUT", 
LANLdb_cluster_name == "MC10-4" ~ "Intra", 
LANLdb_cluster_name == "MC10-8" ~ "IUT", 
LANLdb_cluster_name == "MC10-9" ~ "Intra",               LANLdb_cluster_name == "MC1002-1031" ~ "Intra",                      
LANLdb_cluster_name == "MC1084" ~ "Intra or BM",
LANLdb_cluster_name == "MC11-56" ~ "IUT",
LANLdb_cluster_name == "MC11-95" ~ "IUT",
LANLdb_cluster_name == "MC12-0377" ~ "IUT or Intra (not BM)",
LANLdb_cluster_name == "MC12-0779" ~ "Intra",
LANLdb_cluster_name == "MC12-1005" ~ "IUT",
LANLdb_cluster_name == "MC12-1110" ~ "IUT",
LANLdb_cluster_name == "MC12-1224" ~ "IUT",
LANLdb_cluster_name == "MC1288" ~ "IUT",
LANLdb_cluster_name == "MC5-A" ~ "UT",
LANLdb_cluster_name == "Mother-child-pair1" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair10" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair11" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair13" ~ "IUT",
LANLdb_cluster_name == "Mother-child-pair14" ~ "IUT",
LANLdb_cluster_name == "Mother-child-pair15" ~ "IUT",
LANLdb_cluster_name == "Mother-child-pair19" ~ "UT",
LANLdb_cluster_name == "Mother-child-pair3" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair6" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair7" ~ "BM",
LANLdb_cluster_name == "Mother-child-pair9" ~ "BM",)) 
                                                                                                                                                                               
                                                                                                                                                         #save to csv                      
#write.csv(ann.90.consist, "ann.90.consist.csv")
#write.csv(ann.90.names, "ann.90.names.csv")
write.csv(ann.90.consist.add, "ann.90.consist.add.csv")
#some stats

summary(ann.90.consist$lineages.recipient)

