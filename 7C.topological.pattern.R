#reads a transmission pair tree or a set of trees and computes the topology class 
#(monophyletic-monophyletic [MM], paraphyletic-monophyletic [PM] or paraphyletic-polyphyletic [PP]) and phylogenetic metrics.

topological.pattern <-
    function(tr,
             states,
             outgroup,
             p.ancestral.logs,
             logs,
             pair.id) {
        #treefile: tree file
        #states:
        #load libraries:
        require(caper)
        require(ape)
        require(phytools) #reroot
        require(adephylo) #listDD
        require(dplyr) #pipes
        require(tibble) #column_to_rownames
        require(phangorn)
        
        names(states) <- c("accession", "state")
        
        #useful mini functions (operator inverse of %in%)
        '%!in%' <- function(x, y)
            ! ('%in%'(x, y))
        
        #create states variables for individuals
        if (!exists('states.source'))
            states.source = 0
        if (!exists('states.recipient'))
            states.recipient = 1
        if (missing(logs))
            logs = NA
        pair.id=pair.id
        #vectors to save data
        vector.lineages.source <- integer()
        vector.lineages.recipient <- integer()
        vector.topology <- character()
        vector.p.ancestral.pars <- numeric()
        vector.p.ancestral.ml.ER <- numeric()
        vector.p.ancestral.bayes <- numeric()
        vector.pars.direction <- character()
        vector.ml.ER.direction <- character()
        vector.bayes.direction <- character()
        vector.PD.TBL.source <- numeric()
        vector.PD.SBL.source <- numeric()
        vector.PD.UEH.source <- numeric()
        vector.PD.TIP.source <- numeric()
        vector.PD.MST.source <- numeric()
        vector.PD.TBL.recipient <- numeric()
        vector.PD.SBL.recipient <- numeric()
        vector.PD.UEH.recipient <- numeric()
        vector.PD.TIP.recipient <- numeric()
        vector.PD.MST.recipient <- numeric()
        vector.tree.distinctness.source <- numeric()
        vector.tree.distinctness.recipient <- numeric()
        vector.tree.MinRootToTip.source <- numeric()
        vector.tree.MaxRootToTip.source <- numeric()
        vector.tree.MinRootToTip.recipient <- numeric()
        vector.tree.MaxRootToTip.recipient <- numeric()
        vector.tree.closerToRoot <- character()
        vector.identities <- character()
        LANLdb_cluster_name <- character()
        if (outgroup) {
            states <- states[2:nrow(states),] #remove HXB2
            rownames(states) <- NULL
            states$state <- factor(states$state)
        }
        
        for (i in 1:length(tr)) {
            if (outgroup) {
                troot <- reroot(tr[[i]], 1) #phytools #root to outgroup
                troot <-
                    drop.tip(troot, 1) #remove HXB2
                dt <- troot
                troot$edge.length[troot$edge.length==0]<-max(nodeHeights(dt))*1e-6
                
            } else
                troot <- tr[[i]]
            if (i %in% c(1, seq(length(tr) * 0.1, length(tr), length(tr) * 0.1)))
                print(paste("reading tree number ", i, sep = ""))
            #caculate evolutionary distinctness and rooToTip
            pair.ed.calc <- ed.calc(clade.matrix(troot))
            distRoot.calc <-
                as.data.frame(distRoot(troot))
            distRoot.calc$species <-
                row.names(distRoot.calc)
            names(distRoot.calc) <-
                c('values', 'species')
            
            #get closer to root identity
            closerToRoot.calc <-
                as.data.frame(distRoot(troot,method="nNodes"))
            closerToRoot.calc$species <-
                row.names(closerToRoot.calc)
            names(closerToRoot.calc) <-
                c('values', 'species')
            
            identities <- closerToRoot.calc$species[which(closerToRoot.calc$values== min(closerToRoot.calc$values))]
            
            identities<-unique(gsub("\\..*","",identities))
            if(length(identities)>1) identities <- "both"
            
            #read tree
            taxa.source <-
                as.character(states$accession[which(states$state == states.source)]) #vector of source names
            taxa.rec <-
                as.character(states$accession[which(states$state == states.recipient)]) #vector of recipients names
            nodes <-
                1:troot$Nnode + Ntip(troot) #create vector of internal nodes
            nodes.eval <-
                nodes #create vector of internal nodes that will be iteratively updated
            lineages.recipient <-
                0 #create vector of lineage counts (in recipient) to be iteratively updated
            lineages.source <-
                0 #create ector of lineage counts (in source) to be iteratively updated
            tips <- 1:Ntip(troot) #create vector of tips
            DD <-
                listDD(troot, nameBy = c("label", "number")) #require(adephylo) #List direct descendants for all nodes of a tree
            for (node in nodes) {
                #Loop internal nodes
                if (node %in% nodes.eval) {
                    #Check if to be evaluated
                    Children <-
                        as.integer(unlist(DD[as.character(node)])) #retrieve children of node (direct descendants)
                    for (child in Children) {
                        #evaluate each child separately
                        numDescendants <-
                            getDescendants(troot, child) #list descendants
                        numnodes <-
                            numDescendants[which(numDescendants %in% nodes)] #keep internal nodes, remove tips
                        charDescendants <-
                            troot$tip.label[c(numDescendants)] #get names of descendants
                        chartips <-
                            charDescendants[!is.na(charDescendants)] #remove internal nodes, only keep tips
                        taxtips <-
                            as.character(states$state[which(states$accession %in% chartips)]) #get states (source/recipient) of tips
                        if (child %in% tips) {
                            #the child is a tip
                            if (troot$tip.label[child] %in% taxa.rec) {
                                lineages.recipient <-
                                    lineages.recipient + 1
                            } else {
                                lineages.source <- lineages.source + 1
                            }
                        }
                        else{
                            if (length(unique(taxtips)) == 1) {
                                #if all descedants are the same
                                nodes.eval <-
                                    nodes.eval[!nodes.eval %in% c(child,
                                                                  numnodes)] #remove nodes (they are monophyletic) from nodes.eval
                                if ((unique(
                                    taxtips
                                )) == states.recipient) {
                                    lineages.recipient <-
                                        lineages.recipient + 1
                                } else {
                                    lineages.source <- lineages.source + 1
                                }
                            }
                        }
                    }
                    rm(
                        Children,
                        child,
                        numDescendants,
                        numnodes,
                        charDescendants,
                        chartips,
                        taxtips
                    )
                }
            }
            if (lineages.recipient == 1 &
                lineages.source == 1)
                topology = "MM"
            else {
                #MM, PM, PP?
                if (lineages.recipient >= 2 &
                    lineages.source == 1 |
                    lineages.recipient == 1 &
                    lineages.source >= 2)
                    topology = "PM"
                else
                    topology = "PP"
            }
            #direction
            states.pars <-
                states  %>% column_to_rownames("accession") %>% as.matrix %>% as.phyDat(type =
                                                                                            "USER", levels = c("0", "1")) %>% ancestral.pars(troot, ., type = "MPR")
            p.ancestral.pars <-
                states.pars[[Ntip(troot) + 1]][1]
            states.ml.ER <-
                setNames(factor(states$state),
                         factor(states$accession)) %>% ace(
                             .,
                             troot,
                             method = "ML",
                             model = "ER",
                             type = "discrete"
                         )
            p.ancestral.ml.ER <- states.ml.ER$lik.anc[1]
            if (p.ancestral.pars >= 0.9) {
                pars.direction <- "consistent"
            } else{
                if (p.ancestral.pars <= 0.1) {
                    pars.direction <- "inconsistent"
                } else{
                    pars.direction <- "equivocal"
                }
            }
            if (p.ancestral.ml.ER >= 0.9) {
                ml.ER.direction <- "consistent"
            } else{
                if (p.ancestral.ml.ER <= 0.1) {
                    ml.ER.direction <- "inconsistent"
                } else{
                    ml.ER.direction <- "equivocal"
                }
            }
            #append data to vectors
            vector.lineages.source <-
                c(vector.lineages.source, lineages.source)
            vector.lineages.recipient <-
                c(vector.lineages.recipient,
                  lineages.recipient)
            vector.topology <-
                c(vector.topology, topology)
            vector.p.ancestral.pars <-
                c(vector.p.ancestral.pars, p.ancestral.pars)
            vector.p.ancestral.ml.ER <-
                c(vector.p.ancestral.ml.ER, p.ancestral.ml.ER)
            vector.pars.direction <-
                c(vector.pars.direction, pars.direction)
            vector.ml.ER.direction <-
                c(vector.ml.ER.direction, ml.ER.direction)
            
            #p.BI.signal #this required a corresponding loaded log file
            if (p.ancestral.logs) {
                p.ancestral.bayes <- logs[i, 19]
                if (p.ancestral.bayes >= 0.9) {
                    bayes.direction <- "consistent"
                } else{
                    if (p.ancestral.bayes <= 0.1) {
                        bayes.direction <- "inconsistent"
                    } else{
                        bayes.direction <- "equivocal"
                    }
                }
                vector.p.ancestral.bayes <-
                    c(vector.p.ancestral.bayes,
                      p.ancestral.bayes)
                vector.bayes.direction <-
                    c(vector.bayes.direction,
                      bayes.direction)
            }
            vector.PD.TBL.source <-
                c(vector.PD.TBL.source,
                  pd.calc(clade.matrix(troot), taxa.source, "TBL"))
            vector.PD.SBL.source <-
                c(vector.PD.SBL.source,
                  pd.calc(clade.matrix(troot), taxa.source, "SBL"))
            vector.PD.UEH.source <-
                c(vector.PD.UEH.source,
                  pd.calc(clade.matrix(troot), taxa.source, "UEH"))
            vector.PD.TIP.source <-
                c(vector.PD.TIP.source,
                  pd.calc(clade.matrix(troot), taxa.source, "TIP"))
            vector.PD.MST.source <-
                c(vector.PD.MST.source,
                  pd.calc(clade.matrix(troot), taxa.source, "MST"))
            vector.PD.TBL.recipient <-
                c(vector.PD.TBL.recipient,
                  pd.calc(clade.matrix(troot), taxa.rec, "TBL"))
            vector.PD.SBL.recipient <-
                c(vector.PD.SBL.recipient,
                  pd.calc(clade.matrix(troot), taxa.rec, "SBL"))
            vector.PD.UEH.recipient <-
                c(vector.PD.UEH.recipient,
                  pd.calc(clade.matrix(troot), taxa.rec, "UEH"))
            vector.PD.TIP.recipient <-
                c(vector.PD.TIP.recipient,
                  pd.calc(clade.matrix(troot), taxa.rec, "TIP"))
            vector.PD.MST.recipient <-
                c(vector.PD.MST.recipient,
                  pd.calc(clade.matrix(troot), taxa.rec, "MST"))
            vector.tree.distinctness.source <-
                c(vector.tree.distinctness.source,
                  sum(pair.ed.calc$spp$ED[pair.ed.calc$spp$species %in% taxa.source]))
            vector.tree.distinctness.recipient <-
                c(vector.tree.distinctness.recipient,
                  sum(pair.ed.calc$spp$ED[pair.ed.calc$spp$species %in% taxa.rec]))
            vector.tree.MinRootToTip.source <-
                c(vector.tree.MinRootToTip.source,
                  min(distRoot.calc$values[distRoot.calc$species %in% taxa.source]))
            vector.tree.MaxRootToTip.source <-
                c(vector.tree.MaxRootToTip.source,
                  max(distRoot.calc$values[distRoot.calc$species %in% taxa.source]))
            vector.tree.MinRootToTip.recipient <-
                c(vector.tree.MinRootToTip.recipient,
                  min(distRoot.calc$values[distRoot.calc$species %in% taxa.rec]))
            vector.tree.MaxRootToTip.recipie t <-
                c(vector.tree.MaxRootToTip.recipient,
                  max(distRoot.calc$values[distRoot.calc$species %in% taxa.rec]))
            vector.identities <- c(vector.identities, identities)
            LANLdb_cluster_name <- c(LANLdb_cluster_name, pair.id)
        }
        if (p.ancestral.logs) {
            logs[, "p(0).pars"] <- vector.p.ancestral.pars
            logs[, "p(0).ml.ER"] <-
                vector.p.ancestral.ml.ER
        }
        
        return(
            list(
                'vector.lineages.source' = vector.lineages.source,
                'vector.lineages.recipient' = vector.lineages.recipient,
                'vector.topology' = vector.topology,
                'vector.p.ancestral.pars' = vector.p.ancestral.pars,
                'vector.p.ancestral.ml.ER' = vector.p.ancestral.ml.ER,
                'vector.p.ancestral.bayes' = vector.p.ancestral.bayes,
                'vector.pars.direction' = vector.pars.direction,
                'vector.ml.ER.direction' = vector.ml.ER.direction,
                'vector.bayes.direction' = vector.bayes.direction,
                'vector.PD.TBL.source' = vector.PD.TBL.source,
                'vector.PD.SBL.source' = vector.PD.SBL.source,
                'vector.PD.UEH.source' = vector.PD.UEH.source,
                'vector.PD.TIP.source' = vector.PD.TIP.source,
                'vector.PD.MST.source' = vector.PD.MST.source,
                'vector.PD.TBL.recipient' = vector.PD.TBL.recipient,
                'vector.PD.SBL.recipient' = vector.PD.SBL.recipient,
                'vector.PD.UEH.recipient' = vector.PD.UEH.recipient,
                'vector.PD.TIP.recipient' = vector.PD.TIP.recipient,
                'vector.PD.MST.recipient' = vector.PD.MST.recipient,
                'vector.tree.distinctness.source' = vector.tree.distinctness.source,
                'vector.tree.distinctness.recipient' = vector.tree.distinctness.recipient,
                'vector.tree.MinRootToTip.source' =
                    vector.tree.MinRootToTip.source,
                'vector.tree.MaxRootToTip.source' =
                    vector.tree.MaxRootToTip.source,
                'vector.tree.MinRootToTip.recipient' =
                    vector.tree.MinRootToTip.recipient,
                'vector.tree.MaxRootToTip.recipient' =
                    vector.tree.MaxRootToTip.recipient,
                'vector.identities' = vector.identities,
                'LANLdb_cluster_name' = LANLdb_cluster_name,
                'logs' = logs
            )
        )
        
    }



#run function
topological.pattern(tr=trees, states=states, outgroup=TRUE, p.ancestral.logs=TRUE, logs=logs, pair.id=ids)

ids.top <- list.dirs(path = '.', full.names = FALSE, recursive = FALSE)
ids.top1 <- ids.top[5]

states <- read.csv("11113-11213.790.1692.states", sep ="")
