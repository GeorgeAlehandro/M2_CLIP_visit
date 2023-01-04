library(dplyr)
library(tictoc)
library(monocle3)

#' \code{evaluate_speed}
#' @param data numeric matrix. Flow, CyTOF or scRNASeq data in matrix of shape n cells x m markers.
#' @param labels character vector. Contains the annotations for each cell of data.
#' @param dimred numeric matrix. Low-dimensional embeddings of data. 
#' @param origin character string (default NULL). Label of the cell population to be used as origin to compute trajectories. 
#' @param ori_id numeric (default 1). Numeric index of the cell to be used as origin to compute trajectories. Specify either origin or ori_id but not both. 
#' @param method character vector (default c('monocle')). Names of the trajectory inference methods to evaluate. Currently available: 'monocle'
#' @param sizes vector (default c(1000, 'all')). Size of the dataset(s) to feed to the method. Can accept multiple values to test on multiple datasets. Datasets for a value smaller then the original size of data are prepared via uniform random sampling. Origin cell is then added to the dataset. If 'all', takes the entire dataset without any sampling. 
#' @param rep numeric (default 1). Amount of times to repeat the entire process. Set to 5-10 for a more accurate measure of speed (will take longer).  
#' @param is.artificial logical (default FALSE). Flag to pass if the dataset has been artificially generated, if TRUE values will be fixed to avoid some issues with common pre-processing functions.
#' @details Function used to test one or multiple trajectory inference methods, returning running times, post-analysis objects, and plotting the result.
#' @return returns a named list (times, obj). Each element is a list and has a key corresponding to a method that was used. times contains the average running time for each dataset tested, and obj contains the object returned by the test function associated with the tested method. 
#' @example 
#' res <- evaluate_speed(data, labels, dimred_vae, sizes = c(1000, 5000, 'all'), method = 'monocle', origin='M4')
#' @export


evaluate_speed <- function(data, labels, dimred,
                           origin = NULL, ori_id = 1,
                           method=c('monocle'), sizes=c(1000, 'all'), rep = 1,
                           is.artificial = F){
    objects <- list()
    n <- nrow(data)
    times <- list()[1:length(sizes)] %>% setNames(sizes)
    if(is.artificial) data <- cell.data.like(data)
    rownames(data) <- 1:nrow(data)
    for(s in sizes){
        # sample datasets and find and append origin
        if(is.null(ori_id)){
            ori_id = which(labels==origin)[1]
        }
        if(s=='all' || s>n){
            ss <- 1:n
        }else{
            ss <- sample(1:n, s, replace=F)
        }
        
        if(!(ori_id %in% ss)){
            ss <- c(ss, ori_id)
            ori_id_new <- length(ss)
        }else{
            ori_id_new <- which(ss==ori_id)
        }
        sdata <- data[ss,]
        sdimred <- dimred[ss,]
        slabels <- labels[ss]
        
        ###
        # test methods
        ###
        
        st = 0
        for(i in 1:rep){
            print(paste("testing ",as.character(length(ss))," events, iter ",i, sep=""))
            
            if("monocle" %in% method){
                start.time = Sys.time()
                print("starting monocle")
                
                cds <- test_monocle(sdata, sdimred)
                
                end.time <- Sys.time()
                tm <- round(end.time-start.time, 2)
                # Beware: execution time will be printed without the unit (seconds, minutes or hours)
                print(paste("monocle exec time for ", 
                            as.character(length(ss)), 
                            " events: ",
                            tm,
                            sep=""))
                times[[as.character(s)]] <- tm
                # Plot doesn't always start automatically
                # If no plot is displaying, do the following outside of the function:
                # monocle3::plot_cells(res$obj$monocle[['insert size here']])
                plot_cells(cds)
                objects[[method]][[as.character(s)]] <- cds
                
            }
        }
    }
    return(list(times=times,obj=objects))
}


test_monocle <- function(data, dimred, method='UMAP'){
    cds <- new_cell_data_set(t(as.matrix(data)))
    ## Step 1: Normalize and pre-process the data
    cds <- preprocess_cds(cds, num_dim = 100)
    ## Step 2: Remove batch effects with cell alignment
    # cds <- align_cds(cds, alignment_group = "batch")
    
    ## Step 3: Reduce the dimensions using UMAP
    # cds <- reduce_dimension(cds)
    
    # manually add dimred
    # even if the method used isn't UMAP, name should be left as 'UMAP' for monocle
    row.names(dimred) <- colnames(cds)
    reducedDims(cds)[[method]] <- dimred
    ## Clear out any old graphs:
    cds@principal_graph_aux[[method]] <- NULL
    cds@principal_graph[[method]] <- NULL
    cds@clusters[[method]] <- NULL
    
    ## Step 4: Cluster the cells
    cds <- cluster_cells(cds)
    
    ## Step 5: Learn a graph
    cds <- learn_graph(cds)
    
    ## Step 6: Order cells
    # cds <- order_cells(cds)
    
    return(cds)    
}

cell.data.like <- function(data){
    epsi = 1e-7
    data <- apply(data, 2, function(col){
        mini = min(col)
        col <- col - mini + 1e-7
        return(col)
    })
    return(data)
}
