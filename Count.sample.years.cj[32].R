#Count samples years, find most recent
setwd("/Users/rivkalim/Library/Mobile Documents/com~apple~CloudDocs/HPGH/Atkins Project/TransmissionPairs-master/Clusters") #set working dorectory to where i want to start the loop

#create a data frame to save your data. Inform what kind of data each column stores ()
uniq.years.df<- data.frame('id'=character(),
                        'count'=integer(),
                        'nYears'=integer(),
                        'maxFreq'=numeric(),
                        'aveFreq'=numeric(),
                        'minFreq'=numeric(),
                        'anyNA'=numerical())

for(i in 1:length(idsno129)) {
  setwd(idsno129[i]) #1:length means start at first one is ids and apply loop to all
  #########
  #read the data
  #########
  x <- read.csv("info.gb.csv") #load data
  z <- plyr::count(x$sampling_year) #count. 
  # the '::' tells R to use the function 'count' from the package 'plyr'
  # we use '::' when you dont want to load the full package as you are using
  #only one function from the package
  
  #########
  #get info from z
  #########
  
  #note that plyr::count creates a data frame. The columns names are "x" and "freq"
  #you can confirm this by running:
  #class(z); names(z)
  myNYears <- length(z$x) #extract the number of years
  myMaxFreq <- max(z$freq/nrow(x)) #extract the max freq (divide counts by # observations in x)
  myMinFreq <- min(z$freq/nrow(x)) #extract the min freq (...)
  myAveFreq <- mean(z$freq/nrow(x)) #extract the average freq (...)
  myNa<- sum(is.na(x$sampling_year)) #check is any of the dates is na
  
  #########
  #save info to a dataframe
  #########
  uniq.years.df <- data.frame('id'=i, #in your loop, i is the name of your cluster
                          'count'=nrow(x), # observations in x
                          'nYears'=myNYears,
                          'maxFreq'=myMaxFreq,
                          'aveFreq'=myMinFreq,
                          'minFreq'=myAveFreq,
                          'anyNA'=myNa)
  
  #########
  #Update uniq.years.df by adding uniq.years to it (row binding)
  #########
  uniq.years.df <- rbind(uniq.years.df, uniq.years)
  
  setwd("..")}#set directory back to where I started

View(uniq.years.df)

