# Ashley Mae Conard
# get_tf_relations.r
# Last Mod. 2/19/2020
# Determines network and relationship between TFs across clusters

# Checking input arguments
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Type: /usr/local/bin/Rscript get_tf_relations.r 
            1)/PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/) 
            2) ORGANISM_NCBI_ID (# NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
            3) APP_DIR
            4) WINDOW (how many time steps)", call.=FALSE)
} else if (length(args) == 4) {
    cat("Passed in:", args,"\n")
} else{
    stop("Pass in 
            1) /PATH/TO/INPUT_OUTPUT_DIR/ (e.g. /timeor/results/primary/) 
            2) ORGANISM_NCBI_ID (# NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
            3) APP_DIR
            4) WINDOW (how many time steps)")
}

# Inputting Arguments
IN_OUTPUT <- args[1]
NCBI_TAXO <- strtoi(args[2]) # NCBI Taxonomy: 7227 (Drosophila melanogaster), 9606 (Homo sapiens), 10090 (Mus musculus)
APP_DIR <- as.character(args[3])
WINDOW <- strtoi(args[4])

# Installing libraries
library("STRINGdb")
library(dplyr)
library(stringr)
library(data.table)
library(tidyr)

# Create tuples of timepoints of size window 
get_windows <- function(timep, WIN){
    lst_time_tup = list() # full time + tuple window list
    for(i in 0:length(timep)){
        tuple_wind = list()
        for(j in 1:WIN){
            s = j+i
            tuple_wind = c(tuple_wind, timep[s])
        }
        if(!anyNA(tuple_wind)){
            lst_time_tup = c(lst_time_tup, tuple_wind)
        }
    }
    return(lst_time_tup)
}

# Get list of TFs in current timepoint regulation change
get_tfs_cur_flip <- function(sdb_info_tbl, c_flips, tf_n_exp){
    
    # Initialize subset of tfs and interactions dataframes
    subset_tfs <- data.frame(matrix(ncol=length(names(sdb_info_tbl)), nrow=0))
    colnames(subset_tfs) <- names(sdb_info_tbl)
    subset_interactions <- data.frame(matrix(ncol=length(names(tf_n_exp)), nrow=0))
    colnames(subset_interactions) <- names(tf_n_exp)
    
    # Create subset df of TFs (observed and predicted) cluster numbers
    for(i in rownames(sdb_info_tbl)){
        tf_clust_num <- floor(as.numeric(sdb_info_tbl[i,"cluster_number"]))
        if(toString(tf_clust_num) %in% c_flips[["cluster"]]){
            
            # Adding to subset of TFs to characterize in this time window
            subset_tfs[nrow(subset_tfs) + 1, ] <- sdb_info_tbl[i,]
        }        
    }    
    
    # Extracting just the list of TFs from subset tf df
    list_tfs <- subset_tfs[["gene"]]

    # Adding to subset of TF interactions to characterize in this time window
    subset_interactions <- tf_n_exp[tf_n_exp$A_gene_name %in% list_tfs & tf_n_exp$B_gene_name %in% list_tfs,]
    return(list(tfs=subset_tfs, inters=subset_interactions))
}

# Filling in the temporal relations table
fill_df_tg_temp_rel <- function(stringdb_in_gene_table, tf_net_exp, df_cl_flips, t_windows, WIN){
    
    # Initialize df_tf_temp_rel
    colu_names <- c("cluster", "regulator", "regulated", "regulation_type", "edge_type")
    df_tf_temp_rel <- data.frame(matrix(ncol = length(colu_names), nrow = 0))
    colnames(df_tf_temp_rel) <- colu_names
    
    # For each time window
    for(wi in seq(1,length(t_windows),WIN)){ # iterate in window chunks
        tup = t_windows[wi:(wi+(WIN-1))]
        print("Window: ")
        print(paste(tup,sep=" "))
        print(df_cl_flips$regTime)
        tmp_cur_fli <- df_cl_flips[df_cl_flips$regTime %in% tup,]
        
        # If regulation change from at least one timepoint to another 
        if(dim(tmp_cur_fli)[1]>1){ # no change between timepoints or only 1 TF (i.e. not interaction possible)

            # Get flips in window
            curr_flips <- tmp_cur_fli[order(tmp_cur_fli$regTime),]
            rownames(curr_flips) <- NULL

            # Get TFs and interactions within window
            window = get_tfs_cur_flip(stringdb_in_gene_table, curr_flips, tf_net_exp)
            list_tfs <- window$tfs[["gene"]]
            
            # Classify each interaction between TFs (4 types), and check if time changed is the same (same cluster), TF_A changed before TF_B, or TF_B changed TF_A
            ctf_inters <- t(combn(list_tfs, 2)) # generate all pairs of TFs
            for(row in 1:nrow(ctf_inters)){
                TF_A <- (ctf_inters[row,])[1]
                TF_B <- (ctf_inters[row,])[2]

                # Get TF_A and TF_B cluster, time, and type
                clust_TF_A <- as.numeric(window$tfs[window$tfs$gene==TF_A,"cluster_number"])
                clust_TF_B <- as.numeric(window$tfs[window$tfs$gene==TF_B,"cluster_number"])

                time_TF_A <- strtoi(gsub("[^0-9.-]", "", curr_flips[curr_flips$cluster==floor(clust_TF_A), "regTime"]))
                time_TF_B <- strtoi(gsub("[^0-9.-]", "", curr_flips[curr_flips$cluster==floor(clust_TF_B), "regTime"]))

                type_TF_A <- curr_flips[curr_flips$cluster==floor(clust_TF_A), "regType"]
                type_TF_B <- curr_flips[curr_flips$cluster==floor(clust_TF_B), "regType"]

                print("TF_A, cluster, time and reg. type")
                print(TF_A)
                print(clust_TF_A)
                print(time_TF_A)
                print(type_TF_A)
                
                print("TF_B, cluster, time and reg. type")
                print(TF_B)
                print(clust_TF_B)
                print(time_TF_B)
                print(type_TF_B)
                
                
                #### TF_A changes at same time TF_B (==)
                if(time_TF_A == time_TF_B){
                
                    # Get any interaction between TF_A and TF_B
                    inter_tf_A_tf_B <- window$inters[(window$inters$A_gene_name==TF_A & window$inters$B_gene_name==TF_B) | (window$inters$A_gene_name==TF_B & window$inters$B_gene_name==TF_A),]
                    
                    # If both TFs observed
                    if(clust_TF_A%%1==0 & clust_TF_B%%1==0){
                        
                        # T1: interact: obs_A => obs_B, obs_B => obs_A
                        if(nrow(inter_tf_A_tf_B)!=0){
                            print("equal, obs_obs_known_int") 
                            print(inter_tf_A_tf_B)
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "obs_obs_known_int")
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "obs_obs_known_int")
                            
                        # T2: no interact: obs_A - -> obs_B, obs_B - -> obs_A
                        } else{
                            print("equal, no interaction")
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "obs_obs_pred_int")
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "obs_obs_pred_int")
                        }
                    
                    # If one TF predicted
                    } else if((clust_TF_A%%1!=0 & clust_TF_B%%1==0) | (clust_TF_A%%1==0 & clust_TF_B%%1!=0)){
                        # T3: interact: pred => obs
                        if(nrow(inter_tf_A_tf_B)!=0){
                            print("equal, pred, obs interaction") 
                            print(inter_tf_A_tf_B)
                            
                            # If TF_A is pred
                            if(clust_TF_A%%1!=0 & clust_TF_B%%1==0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "pred_obs_known_int")
                            
                            # If TF_B is pred
                            } else if(clust_TF_A%%1==0 & clust_TF_B%%1!=0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "pred_obs_known_int")
                            }
                                
                        # T4: no interact: pred - -> obs
                        } else{
                            print("equal, pred, obs no interaction")  
                            
                            # If TF_A is pred
                            if(clust_TF_A%%1!=0 & clust_TF_B%%1==0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "pred_obs_pred_int")
                            
                            # If TF_B is pred
                            } else if(clust_TF_A%%1==0 & clust_TF_B%%1!=0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "pred_obs_pred_int")
                            }
                        }
                    }
                }
                
                
                #### TF_B changes before TF_A (>)
                if(time_TF_A > time_TF_B){
                    
                    # If both observed
                    if(clust_TF_A%%1==0 & clust_TF_B%%1==0){
                        
                        # Get any interaction between TF_A and TF_B
                        inter_tf_A_tf_B <- window$inters[(window$inters$A_gene_name==TF_A & window$inters$B_gene_name==TF_B) | (window$inters$A_gene_name==TF_B & window$inters$B_gene_name==TF_A),]
                        
                        # T1: interact: obs => obs
                        if(nrow(inter_tf_A_tf_B)!=0){
                            print("TF_B -> TF_A interaction")
                            print(inter_tf_A_tf_B)
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "obs_obs_known_int")
                            
                        # T2: no interact: obs - -> obs
                        } else{
                            print("TF_B -> TF_A no interaction")   
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "obs_obs_pred_int")
                        }
                    
                    # If one TF predicted
                    } else if((clust_TF_A%%1!=0 & clust_TF_B%%1==0) | (clust_TF_A%%1==0 & clust_TF_B%%1!=0)){
                        # T3: interact: pred => obs
                        if(nrow(inter_tf_A_tf_B)!=0){
                            
                            # If TF_B is pred
                            if(clust_TF_A%%1==0 & clust_TF_B%%1!=0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "pred_obs_known_int")
                            }
                            
                        # T4: no interact: pred - -> obs
                        } else{
                            print("TF_B -> TF_A, pred, obs no interaction")      
                            
                            # If TF_B is pred
                            if(clust_TF_A%%1==0 & clust_TF_B%%1!=0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_A), TF_B, TF_A, type_TF_A, "pred_obs_pred_int")
                            }  
                        }
                    }
                }

                #### TF_A changes before TF_B (<)
                if(time_TF_A < time_TF_B){
                    
                    # If both observed
                    if(clust_TF_A%%1==0 & clust_TF_B%%1==0){
                        
                        # Get any interaction between TF_A and TF_B
                        inter_tf_A_tf_B <- window$inters[(window$inters$A_gene_name==TF_A & window$inters$B_gene_name==TF_B) | (window$inters$A_gene_name==TF_B & window$inters$B_gene_name==TF_A),]
                        
                        # T1: interact: obs => obs
                        if(nrow(inter_tf_A_tf_B)!=0){
                            print("TF_A -> TF_B interaction")
                            print(inter_tf_A_tf_B)
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "obs_obs_known_int")
                            
                        # T2: no interact: obs - -> obs
                        } else{
                            print("TF_A -> TF_B no interaction")
                            df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "obs_obs_pred_int")
                        }
                    
                    # If one TF predicted
                    } else if((clust_TF_A%%1!=0 & clust_TF_B%%1==0) | (clust_TF_A%%1==0 & clust_TF_B%%1!=0)){
                        # T3: interact: pred => obs
                        if(nrow(inter_tf_A_tf_B)!=0){
                            print("TF_A -> TF_B, pred, obs interaction") 
                            print(inter_tf_A_tf_B)
                            
                            # If TF_A is pred
                            if(clust_TF_A%%1!=0 & clust_TF_B%%1==0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "pred_obs_known_int")
                            }
                            
                        # T4: no interact: pred - -> obs
                        } else{
                            print("TF_A -> TF_B, pred, obs no interaction")     
                            
                            # If TF_A is pred
                            if(clust_TF_A%%1!=0 & clust_TF_B%%1==0){
                                df_tf_temp_rel[nrow(df_tf_temp_rel) + 1,] = c(floor(clust_TF_B), TF_A, TF_B, type_TF_B, "pred_obs_pred_int")
                            }
                        }
                    }
                }
            }         
        }
    }
    print(paste(tup,sep=" "))
    #print("df_tf_temp_rel")
    #print(df_tf_temp_rel)
    #stop("blah")
    return(distinct(df_tf_temp_rel))
}

# Finding all flips between up and down regulation across clusters
find_flips <- function(df_c_lev_traj){
    
    # Create transition dataframe of regulation change: i.e. cluster by tuples list (action_time, regulation_type)
    df_transitions = data.frame(matrix(ncol=2, nrow=nrow(df_c_lev_traj)))
    colnames(df_transitions) <- c("cluster", "actT_regType")
    
    # For each cluster
    for(r in 1:nrow(df_c_lev_traj)){
        df_transitions$cluster[r] <- df_c_lev_traj$cluster[r]
        skip_gene_name = 1 # skip first column with gene name
        tp_1 = 1 # first timepoint
        
        for(c in names(df_c_lev_traj)){
            if(skip_gene_name){
                skip_gene_name = 0
                next
            } else{
                val <- df_c_lev_traj[[c]][r] # get mean expression value for that cluster and timepoint
                if(tp_1){# if first timepoint
                    
                    if(val<0){ # if negative mean expression value for that cluster and timepoint
                        tmp = -1
                    } else if(val>0){ # if positive mean expression value for that cluster and timepoint
                        tmp = 1
                    } else{ # if zero mean expression value for that cluster and timepoint
                        tmp = 0
                    }
                    tp_1 = 0 # no longer first timepoint
                } else { # if not first timepoint
                    # If first regulation change
                    if(is.na(df_transitions$actT_regType[r])){
                        if(tmp <= 0 && val > 0){ # detect upregulation
                            tmp = 1 # current state change to upregulated
                            df_transitions$actT_regType[r] <-  paste("a", c, sep=";")
                        } else if(tmp >= 0 && val < 0){# detect downregulation
                            tmp = -1 # current state change to downregulated
                            df_transitions$actT_regType[r] <- paste("r", c, sep=";")
                        } else{
                            next
                        }
                    # Append second regulation change
                    } else{ 
                        if(tmp <= 0 && val > 0){ # detect upregulation
                            tmp = 1 # current state change to upregulated
                            next_reg_change = paste("a", c, sep=";")
                            df_transitions$actT_regType[r] <-  list(append(df_transitions$actT_regType[r], next_reg_change))
                        } else if(tmp >= 0 && val < 0){# detect downregulation
                            tmp = -1 # current state change to downregulated
                            next_reg_change = paste("r", c, sep=";")
                            df_transitions$actT_regType[r] <- list(append(df_transitions$actT_regType[r], next_reg_change))
                        }else{
                            next
                        }
                    }
                }
            }
        }  
    }
    
    # Initializing list of first timepoint regulation types. If one is different (e.g. all are downregulated except one, propose that 'one' is initiating the cascade
    # Get frequency of up and downregulation for each cluster at 1st timepoint
    df_frq_reg_type_tp1 <- as.data.frame(table(sign(df_c_lev_traj[,2])))
    
    # If one cluster is different, propose that cluster as the cascade initator
    if(1 %in% df_frq_reg_type_tp1$Freq){
        
        # Get regulation sign and set regulation type
        reg_sign <- strtoi(df_frq_reg_type_tp1[df_frq_reg_type_tp1$Freq==1,]$Var1)
        if(reg_sign > 1){
            reg_type_tp1 <- "a"
        } else{
            reg_type_tp1 <- "r"
        }
        
        for(t in 1:length(df_frq_reg_type_tp1)){
            if(sign(df_c_lev_traj[t,2]) == reg_sign){
                clust_tp1 <- df_c_lev_traj[t,1]
                df_tp1 <-data.frame("cluster"=clust_tp1, "actT_regType"=paste(reg_type_tp1, names(df_c_lev_traj)[2],sep=";"))
            }
        }
    }

    df_transitions <- rbind(df_transitions, df_tp1)
    df_trans <- df_transitions %>% separate(actT_regType, c("regType", "regTime"))
  
    return(df_trans)
}

# Finding cluster level mean expression to determine overall expression changes of each cluster
get_cluster_level_mean_traj <- function(filens, all_tp_names){
    
    # Create dataframe of cluster by time points where each value is the mean expression for that cluster
    df_cluster_level_traj = data.frame(matrix(ncol=length(all_tp_names), nrow=length(filens)))
    colnames(df_cluster_level_traj) <- all_tp_names
    names(df_cluster_level_traj)[names(df_cluster_level_traj) == "X"] <- "cluster"

    # Fill df_cluster_level_traj
    for(dir_file in filens){
        geneNumClust=as.integer((sub(".*_([^.]+)\\.csv.*", "\\1", dir_file)))
        clust_table <- read.table(dir_file, sep=",", stringsAsFactors=FALSE, header=TRUE)

        skip_gene_name = 1 # skip first column with gene name
        for(col_na in names(clust_table)){
            if(skip_gene_name){
                skip_gene_name = 0
                next
            } else{
                col_vector <- clust_table[[col_na]]
                mode_col <- mean(col_vector)   
                df_cluster_level_traj$cluster[geneNumClust] <- geneNumClust
                df_cluster_level_traj[[col_na]][geneNumClust] <- mode_col
            }
        }
    }    
    return(df_cluster_level_traj)       
}

# Determining the timepoint(s) when more than 50% of genes (i.e. mode) when genes in cluster flip 
convert_to_stringdb_input <- function(obs_n_put){
    col_names <- c("gene",  "cluster_number")
    stringdb_input_table <- data.frame(matrix(ncol = length(col_names), nrow = 0))
    colnames(stringdb_input_table) <- col_names 

    # Adding observed TFs to information table
    for(i in 1:nrow(obs_n_put)){
        if (obs_n_put$observed_TFs[i] != ''){
            obs_col_genes <- unlist(strsplit(obs_n_put$observed_TFs[i], "\\|"))
            if (length(obs_col_genes)>1){
                for(j in 1:length(obs_col_genes)){
                    stringdb_input_table[nrow(stringdb_input_table) + 1,] = c(obs_col_genes[j],as.integer(obs_n_put$cluster[i]))
                }
            } else{
                stringdb_input_table[nrow(stringdb_input_table) + 1,] = c(obs_col_genes,as.integer(obs_n_put$cluster[i]))
            }
        }
        
        # If there are predicted genes in any of the clusters
        if (!is.na(obs_n_put$top_TF_1[i]) & obs_n_put$top_TF_1[i]!=''){
            
            # Adding putative gene to stringdb_input_table for that cluster
            stringdb_input_table[nrow(stringdb_input_table) + 1,] = c(obs_n_put$top_TF_1[i],as.double(obs_n_put$cluster[i])+0.5)
        }
    }
    return(stringdb_input_table)
}

main <- function(){
    
    # Window of 1 time step (window) is considering the change between 2 timepoints and so on
    WIN = WINDOW+1 
    
    # Getting results folder name for user's experiment 
    exp <- basename(IN_OUTPUT)
    print(paste("Results folder name: ",exp, sep=" "))
    
    # Iterating through all experiments in list_exp
    if(!grepl("_results", exp, fixed=T)){
        stop("Please pass in a named results folder (e.g. test_results)")
    }else{
    
#     # Getting all the experiment types 
#     list_exp = (list.dirs(path = IN_OUTPUT, recursive = FALSE, full.names=FALSE))   
    
#     # For all experiments in list_exp
#     for(exp in list_exp){    
        # Getting tf_reg_network folder
        dir_tf_rn <- paste(IN_OUTPUT, "temporal_relations", sep="/")
        dir.create(file.path(IN_OUTPUT, "temporal_relations"), showWarnings = FALSE)
        
        # Finding the cluster directories for each experiment (in list_exp)
        dir_clusts <- paste(IN_OUTPUT, "clusters", sep="/")

        # Finding only filenames that are *geneList*.csv or *geneName*.csv to get list of genes to calc. enrichment
        filenames <- Sys.glob(file.path(dir_clusts, "/*/", "*geneList*.csv"))
        print(dir_clusts)
        print(filenames)
        
        # Getting the number of timepoints
        all_timepoint_names <- names(read.table(filenames[1], sep=",", stringsAsFactors=FALSE, header=TRUE))

        # Getting cluster level trajectories
        df_clust_lev_traj <- get_cluster_level_mean_traj(filenames, all_timepoint_names)
        
        # Finding flips across each cluster (note 'flip' is a cluster-wide mean log2 fold change)
        df_clust_flips <- find_flips(df_clust_lev_traj)

        # Getting observed_predicted_tfs_n_encode
        obs_n_put <- fread(paste(dirname(dir_tf_rn),"factor_binding","observed_predicted_tfs_n_encode.csv",sep="/"), sep=",", select=c("cluster","observed_TFs","top_TF_1"))
        
        # Creating tuples of timepoints of size window 
        list_tup_windows <- get_windows(all_timepoint_names[2:length(all_timepoint_names)], WIN) # from 2 because first column of observed_putative_tfs_n_encode.csv is gene name
        
        # Converting this to stringdb_input list of TFs
        stringdb_input_gene_table <- convert_to_stringdb_input(obs_n_put)
        fwrite(stringdb_input_gene_table, file = paste(dir_tf_rn, "string_db_input.csv", sep="/"))
        print("stringdb input gene table")
        print(stringdb_input_gene_table)

        # Calling stringdb output 
        system(paste("Rscript", paste(APP_DIR, "/scripts/stringdb_top_tfs.r", sep="/"), IN_OUTPUT, NCBI_TAXO, "1", APP_DIR))
        
        # Filtering stringdb output to keep only experiments and experiments_transferred
        tf_network <- fread(paste(dirname(dir_tf_rn),"temporal_relations","stringdb_info_table.tsv",sep="/"), sep="\t", select=c("A_gene_name","B_gene_name", "experiments", "experiments_transferred"))
        tf_network_exp <- tf_network[tf_network$experiments>0 | tf_network$experiments_transferred>0 ]
        print("tf network")
        print(tf_network_exp)

        # Creating TF temporal relationship dataframe
        df_tf_temp_relations <- fill_df_tg_temp_rel(stringdb_input_gene_table, tf_network_exp, df_clust_flips, list_tup_windows, WIN)
        print("df tf temp relations")
        print(df_tf_temp_relations)
        fwrite(df_tf_temp_relations, file = paste(dir_tf_rn, "TF_temp_rel.csv", sep="/"))

        print("all_timepoint_names")
        print(all_timepoint_names)

        print("df_clust_lev_traj")
        print(df_clust_lev_traj)

        print("obs_n_put")
        print(obs_n_put)

        print("list_tup_windows")
        print(list_tup_windows)

        print("stringdb input gene table")
        print(stringdb_input_gene_table)

        print("tf network")
        print(tf_network_exp)

        print("df tf temp relations")
        print(df_tf_temp_relations)
    }
}

if(!interactive()) {
    main()
}