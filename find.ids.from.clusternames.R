#######Transmission pairs have been found from published sources "clusternames", find the IDs from these files. 
#Load required libraries and specify path
#######
require(tidyverse)
mypath <- "."
#######
#List all files
#######
myfiles <- list.files(path = mypath, recursive = FALSE)
#######
#Create a data frame to save the ids
#######
myids <- data.frame('cluster_name'=character(), 
                     'cluster_id'=integer())
#######
#Loop every file
#######
for(i in myfiles){
        print(paste("reading file", i))
        mydf <- read.delim(i)
        cluster_comb <- unique(mydf$cluster_comb)
        for(j in cluster_comb){
                cluster_name <- str_match_all(j, "^[^\\(]+") %>% toString()
                cluster_id <- str_extract_all(j, "\\([^()]+\\)")[[1]]
                cluster_id <- substring(cluster_id, 2, nchar(cluster_id)-1)
                cluster_id <- as.numeric(cluster_id)
                myrow <- data.frame('cluster_name'=cluster_name, 
                                    'cluster_id'=cluster_id)
                myids <- rbind(myrow, myids)
                
        }; rm(cluster_name, cluster_id, myrow)
        
};rm(i, j, mydf, cluster_comb, myfiles)


library("writexl")
no
write_xlsx(myids,"\\Users\\rivkalim\\Library\\Mobile\\Document\\com\\apple\\CloudDocs\\HPGH\\Atkins\\Project\\id.xlsx")
