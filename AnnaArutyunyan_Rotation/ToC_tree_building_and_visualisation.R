# BUILDING AND VISUALISING TREES 
# FROM COSINE SIMILARITY DISTANT MATRICES

library(ape)
library(rpart)
library(maps)
library(phytools) # requires magick
library(quadprog)
library(phangorn)
library(qgraph) 
#packageVersion("phytools")

# reading all the cosine similarity distance matrices and creating NJ trees based on all DE genes

read_matrix_create_nj_tree <- function(filename_matrix){
  
  # reading the distance matrix
  curr_table = read.table(file=filename_matrix,header = TRUE,fill = TRUE,sep = '\t')
  curr_table = curr_table[c(2:ncol(curr_table))]
  curr_table_dist = as.dist(curr_table)
  
  # creating the neighbour joining tree
  nj.tree <- bionj(curr_table_dist)
  curr_tree_all_DE_init <- nj.tree
  
  return(curr_tree_all_DE_init)
}

# PLANARIA was comma separated so rewriting it to be tab separated as other datasets 
planaria_table = read.table(file='/home/jovyan/notebooks/Tree_of_Cells/data/cosine_matrix_DE_planaria.csv',header = TRUE,fill = TRUE,sep = ',')
# writing with tab separation to have all matrices tab separated 
write.table(as.matrix(planaria_table), '/home/jovyan/notebooks/Tree_of_Cells/data/cosine_matrix_DE_planaria_tab.csv', append = FALSE, sep = "\t",
            row.names = TRUE, col.names = TRUE)

# reading all the matrices and building NJ trees (all DE genes, cosine similarity)
mouse_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/data/cosine_matrix_DE_mouse.csv')
planaria_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/data/cosine_matrix_DE_planaria_tab.csv')
Sim_1112_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/Sim_1112_cosine_dist_matrix_DE.csv')
Sim_3814_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/Sim_3814_cosine_dist_matrix_DE.csv')
Sim_528_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/Sim_528_cosine_dist_matrix_DE.csv')
Sim_71640_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/Sim_71640_cosine_dist_matrix_DE.csv')
Sim_9901_tree_all_DE_init = read_matrix_create_nj_tree('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/Sim_9901_cosine_dist_matrix_DE.csv')

# Robustness test

read_trees_into_list <- function(tag_word, N_trees){
  tables = list()
  trees <- list()
  for (i in c(1:N_trees)) { 
    print(i)
    curr_address = paste('/home/jovyan/notebooks/Tree_of_Cells/Simulations_Tallulah/',tag_word,
                         i, '_cosine_dist_matrix_DE_subset0.8.csv',
                         sep='')
    curr_table <- read.table(file=curr_address, header = TRUE, fill = TRUE, sep = '\t')
    curr_table = curr_table[c(2:ncol(curr_table))]
    curr_table_dist = as.dist(curr_table)
    #print(curr_table)
    tables[[i]] = list(curr_table)
    trees[[i]] =  bionj((curr_table_dist))
    
  }
  return(trees)
}

generate_consensus_tree <- function(tree_list){
  #start.time <- Sys.time()
  consensus <- averageTree(tree_list)
  #end.time <- Sys.time()
  #time_cons <- end.time - start.time
  return(consensus)
  #print('runtime',time_cons)
}

# reading trees for all datasets
mouse_trees <- read_trees_into_list('mouse',100)
planaria_trees <- read_trees_into_list('planaria',33) # only 33 trees for planaria - big dataset
Sim_3814_trees <- read_trees_into_list('Sim_3814',100)
Sim_1112_trees <- read_trees_into_list('Sim_1112',100)
Sim_528_trees <- read_trees_into_list('Sim_528',100)
Sim_71640_trees <- read_trees_into_list('Sim_71640',100)
Sim_9901_trees <- read_trees_into_list('Sim_9901',100)


# creating consensus trees, overall takes about 2-3 hours
mouse_tree_consensus <- generate_consensus_tree(mouse_trees) 
planaria_tree_consensus <- generate_consensus_tree(planaria_trees)
Sim_3814_consensus <- generate_consensus_tree(Sim_3814_trees)
Sim_1112_consensus <- generate_consensus_tree(Sim_1112_trees)
Sim_528_consensus <- generate_consensus_tree(Sim_528_trees)
Sim_71640_consensus <- generate_consensus_tree(Sim_71640_trees)
Sim_9901_consensus <- generate_consensus_tree(Sim_9901_trees)


# RF distances between initial trees (all DE genes) and consensus tree across 100 80%-trees (80%-DE genes)
mouse_delta <- RF.dist(tree1 = mouse_tree_all_DE_init, tree2 = mouse_tree_consensus, normalize = TRUE, check.labels = TRUE,
                       rooted = TRUE)

planaria_delta <- RF.dist(tree1 = planaria_tree_all_DE_init, tree2 = planaria_tree_consensus, normalize = TRUE, check.labels = TRUE,
                          rooted = TRUE)

Sim_1112_delta <- RF.dist(tree1 = Sim_1112_tree_all_DE_init, tree2 = Sim_1112_consensus, normalize = TRUE, check.labels = TRUE,
                          rooted = TRUE)

Sim_3814_delta <- RF.dist(tree1 = Sim_3814_tree_all_DE_init, tree2 = Sim_3814_consensus, normalize = TRUE, check.labels = TRUE,
                          rooted = TRUE)

Sim_528_delta <- RF.dist(tree1 = Sim_528_tree_all_DE_init, tree2 = Sim_528_consensus, normalize = TRUE, check.labels = TRUE,
                         rooted = TRUE)

Sim_71640_delta <- RF.dist(tree1 = Sim_71640_tree_all_DE_init, tree2 = Sim_71640_consensus, normalize = TRUE, check.labels = TRUE,
                           rooted = TRUE)

Sim_9901_delta <- RF.dist(tree1 = Sim_9901_tree_all_DE_init, tree2 = Sim_9901_consensus, normalize = TRUE, check.labels = TRUE,
                          rooted = TRUE)



# generating 100 more RF distances per sample (between consensus and 80%-trees - there are 100 of them)
# to estimate confidence interval


generate_RF_list <- function(tree_list,consensus_tree){
  #start.time <- Sys.time()
  RF_list = list()
  for (tree in tree_list) { 
    newelem <- RF.dist(tree1 = tree, tree2 = consensus_tree, normalize = TRUE, check.labels = TRUE,
                       rooted = TRUE) 
    RF_list <- c(RF_list, newelem)
  }
  #end.time <- Sys.time()
  #time_cons <- end.time - start.time
  return(RF_list)
  #print('runtime',time_cons)
}

RF_list_mouse <- generate_RF_list(tree_list = mouse_trees, consensus_tree = mouse_tree_consensus)
RF_list_planaria <- generate_RF_list(tree_list = planaria_trees, consensus_tree = planaria_tree_consensus)
RF_list_3814 <- generate_RF_list(tree_list = Sim_3814_trees, consensus_tree = Sim_3814_consensus)
RF_list_528 <- generate_RF_list(tree_list = Sim_528_trees, consensus_tree = Sim_528_consensus)
RF_list_1112 <- generate_RF_list(tree_list = Sim_1112_trees, consensus_tree = Sim_1112_consensus)
RF_list_71640 <- generate_RF_list(tree_list = Sim_71640_trees, consensus_tree = Sim_71640_consensus)
RF_list_9901 <- generate_RF_list(tree_list = Sim_9901_trees, consensus_tree = Sim_9901_consensus)



# Calculating some statistics

tree_stats_mean_95CI <- function(dist_list){
  dist_list = as.numeric(dist_list)
  
  # Looking at distribution of RF distances between 100 80%-DE trees and consensus tree, optional
  #hist(dist_list)
  
  # Calculating means
  curr_mean = mean(dist_list)
  
  # Calculating 95% confidence interval (CI) for RF distances per sample, no assumptions about the distributions
  curr_CI = quantile(dist_list, probs=c(0.025, 0.975))
  return(curr_mean, curr_CI)
}

# Calculating means and 95% CI for all datasets
tree_stats_mean_95CI(RF_list_mouse)
tree_stats_mean_95CI(RF_list_planaria)
tree_stats_mean_95CI(RF_list_1112)
tree_stats_mean_95CI(RF_list_3814)
tree_stats_mean_95CI(RF_list_528)
tree_stats_mean_95CI(RF_list_71640)
tree_stats_mean_95CI(RF_list_9901)






