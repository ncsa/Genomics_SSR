## 1. Installing Packages
install_packages <- function(packages = c()){
  for (p in packages){
    install.packages(p)
    library(p, character.only = T)
  }
}

## 2. Reading all input files
read_input_files <- function(input_files = c()){
  Test_Seeds <- read.csv(input_file[1],header=FALSE, sep="\t", dec=".",stringsAsFactors = FALSE)
  Test_Seeds <- Test_Seeds$V1
  
  Parameters_file <- read.csv(input_file[2],header=TRUE, sep="\t", dec=".",stringsAsFactors = FALSE)
  
  r <- as.numeric(Parameters_File[1,2])
  delta <- as.numeric(Parameters_File[2,2])
  tau <- as.numeric(unlist(strsplit(Parameters_File[3,2],",")))
  k <- as.numeric(Parameters_File[4,2])
  lambda <- as.numeric(Parameters_File[5,2])
  eta <- as.numeric(Parameters_File[6,2])
  
  parameters <- list(Test_Seeds, r, delta, tau, k, lambda, eta)
  names(parameters) <- c("testseed", "r", "delta", "tau", "k", "lambda", "eta")
  return(parameters)
}

## 3. Reading Network files and building graphs
read_layers <- function(networks = c()){
  Layers <- vector("list", 3)
  PPI <- read.table(networks[1], sep=" ")
  PPI_Network <- graph.data.frame(PPI,directed=FALSE)
  PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
  
  Layers[[1]] <- PPI_Network
  names(Layers)[1] <- "PPI_NETWORK"
  
  Pathway <- read.table(networks[2], sep=" ")
  Pathway_Network <- graph.data.frame(Pathway,directed=FALSE)
  Pathway_Network <- igraph::simplify(Pathway_Network, remove.multiple = TRUE, remove.loops = TRUE)
  
  Layers[[2]] <- Pathway_Network
  names(Layers)[2] <- "PATHWAY_NETWORK"
  
  Coex <- read.table(networks[3], sep=" ")
  Coex_Network <- graph.data.frame(Coex,directed=FALSE)
  Coex_Network <- igraph::simplify(Coex_Network, remove.multiple = TRUE, remove.loops = TRUE)
  
  Layers[[3]] <- Coex_Network
  names(Layers)[3] <- "COEXPRESION_NETWORK"
  
  return(Layers)
}

## 4. All unique Nodes together
Total_nodes <- function(layers){
  Total_layers <- length(layers)
  
  All_nodes <- character()
  for (i in 1:Total_layers){
    Node_Names_Layer <- V(Layers[[i]])$name
    All_nodes <-c(All_nodes, Node_Names_Layer)
  }
  
  All_nodes <- unique(All_nodes)
  
  return(All_nodes)
}

## 5. Generating Adjacency matrix for multiplex layers
generating_multiplex <- function(layers, delta, N){
  eye_matrix <- Diagonal(N, x = 1)
  L <- length(layers)
  multi_adj_mat <- Matrix(0, ncol = N*L, nrow = N*L, sparse = TRUE)  # N*L = 17559 * 3 = 52677 => 52667 by 52667 matrix is created
  
  Col_NodeNames <- character()
  Row_NodeNames <- character()
  
  for (i in 1:L)
  {
    Adjacency_layer <- as_adjacency_matrix(layers[[i]], sparse = TRUE)
    
    #order function arranges all the rownames and colnames in increasingly manner
    
    Adjacency_layer <- Adjacency_layer[order(rownames(Adjacency_layer)),order(colnames(Adjacency_layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_layer),i,sep="_")
    
    Col_NodeNames <- c(Col_NodeNames,Layer_Col_Names)
    Row_NodeNames <- c(Row_Node_Names,Layer_Row_Names)
    
    ## We fill the diagonal blocks with the adjacencies matrix of each layer.
    ini_position_row <- 1 + (i-1)*N    # becoz we have three layer and we need to set for each layers diagonal and hence this
    end_position_row <- N + (i-1)*N
    multi_adj_mat[(ini_position_row:end_position_row),(ini_position_row:end_position_row)] <- (1-delta)*(Adjacency_layer)
    
    for (j in 1:L){
      ini_position_col <- 1 + (j-1)*N
      end_position_col <- N + (j-1)*N
      if (j != i){
        multi_adj_mat[(ini_position_row:end_position_row),(ini_position_col:ini_position_col)] <- (delta/(L-1))*eye_matrix
      }
    }
  }
  
  rownames(multi_adj_mat) <- Row_NodeNames
  colnames(multi_adj_mat) <- Col_NodeNames
  
  return(multi_adj_mat)
}

## 6. Generation of Bipartite Graph
bipartite_network <- function(all_nodes_sort, disease_names, Gene_Phenotype_relation,N,M){
  Bipartite_matrix <- Matrix(data=0, nrow=N, ncol=M)
  rownames(Bipartite_matrix) <- all_nodes_sort
  colnames(Bipartite_matrix) <- disease_names
  report <- character()
  
  for (i in 1:N){
    current_node <- all_nodes_sort[i]
    current_mim <- Gene_Phenotype_relation$mim_morbid[which(Gene_Phenotype_relation$hgnc_symbol == current_node)]
    
    if (length(current_mim) > 0){
      for (j in 1:length(current_mim)){
        if (!is.na(current_mim[j])){
          # We need to identify the phenotypes position on the matrix.
          index_disease <- which(colnames(Bipartite_matrix) %in%  current_mim[j])
          # We have to check if that index is present in the matrix.
          if (length(index_disease) == 1){ 
            Bipartite_matrix[i,index_disease] <- 1
          } else {
            error_message <- paste("MIM_CODE", current_mim[j], length(index_disease), "No phenotype found",sep=";", collapse = NULL)
            report <- c(report, error_message)
          }           
        }
      }  
    } else {
      error_message <- paste("HGNC_Symbol", current_node, "No HGNC found",sep=";", collapse = NULL)
      report <- c(report,error_message)
    }
  }
  Bipartite_and_error <- list(Bipartite_matrix,report)
  return(Bipartite_and_error)
}

## 7. Adding Bipartite graph to multiplex system
expanding_bipartite <- function(N, L, M, Bipartite_matrix){
  Super_Matrix <- Matrix(0,nrow=N*L, ncol=M, sparse = TRUE)
  Row_NodeNames <- character()
  
  for (i in 1:L){
    Layer_RowNames <- paste(rownames(Bipartite_matrix),i,sep="_")
    Row_NodeNames <- c(Row_NodeNames,Layer_RowNames)
    ini_position_row <- 1 + (i-1)*N
    end_position_row <- N + (i-1)*N
    Super_Matrix[(ini_position_row:end_position_row),] <- Bipartite_matrix
  }  
  
  rownames(Super_Matrix) <- Row_NodeNames
  colnames(Super_Matrix) <- colnames(Bipartite_matrix)
}



packages <- c('igraph', 'Matrix')
install_packages(packages)

input_files <- c('Input_Files//Test_Seeds.txt', 'Input_Files//Parameters_Example.txt')
parameters <- read_input_files(input_files)

networks <- c("Networks/PPI_2016-11-23.gr", "Networks/AllPathways_2016-11-24.gr", "Networks/Co-Expression_2016-11-23.gr")
list_layers <- read_layers(networks)
L <- length(list_layers)

All_nodes <- Total_nodes(list_layers)
All_nodes_sort <- sort(All_nodes)
N <- length(All_nodes)

multi_adj_mat <- generating_multiplex(list_layers, parameters$delta, N)

Disease <- read.table("Networks//DiseaseSimilarity_2016-12-06.gr",sep=" ")
Disease_Network <- graph.data.frame(Disease,directed=FALSE)
Disease_Network <- igraph::simplify(Disease_Network, remove.multiple = TRUE, remove.loops = TRUE)
Adjacency_mat_Diseases <- as_adjacency_matrix(Disease_Network,sparse = TRUE)

## number of diseases in network
M <- nrow(Adjacency_matrix_Diseases)

Gene_Phenotype_relation<-read.table("Input_Files//Gene_Phenotype_relation.txt", sep="\t", header=TRUE,stringsAsFactors = FALSE)
Gene_Phenotype_relation <- Gene_Phenotype_relation[which(Gene_Phenotype_relation$hgnc_symbol %in% Node_Names_all), ]

Bipartite_matrix_and_error <- bipartite_network(All_nodes_sort, colnames(Adjacency_mat_Diseases), Gene_Phenotype_relation,N,M)
Bipartite_matrix <- Bipartite_matrix_and_error[[1]]
Error_report <- Bipartite_matrix_and_error[[2]]

Super_Matrix <- expanding_bipartite(N,L,M,Bipartite_matrix)