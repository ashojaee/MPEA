
MPEA= function(expressions,pathway,y,num_simulation=9999,ncore=6,seed.number=NA){
require(energy)
require(plyr)
require(foreach)
require(doParallel)
  
individual_pathway_length <- sapply(pathway_collapsed,length)
gene_number <- length(gene_names)
pathway_number <- length(individual_pathway_length)

dcor_results <- list()
univariate_results <- list()
univariate_results_adjusted <- list()

for (ii in 1:length(pathway)){
  pathway_expressions <- expressions[,colnames(expressions)%in% pathway[[ii]]]
  
  dcor_results[ii] <- dcor((pathway_expressions),(y))
  
}
names(dcor_results)  <- names( pathway)
if (!is.na(seed.number)){
set.seed= seed.number
}
# surrogate_dcor <- list()
cl<-makeCluster(ncore)
registerDoParallel(cl)

surrogate_dcor <- foreach (ii = 1:pathway_number,  .export = c("num_simulation","gene_number","individual_pathway_length","y"), .packages = "energy") %dopar%{
  if(ii%%100==0) print(ii)
  aux_dcor <- NULL
  for ( jj in 1: num_simulation)
  {
    aux_genes <- sample(1:gene_number,individual_pathway_length[ii])
    aux_dcor[jj] <- dcor(expressions[,aux_genes],y)
    
  }  
  # surrogate_dcor[[ii]] <- aux_dcor
  aux_dcor
}
stopCluster(cl)
nominal_dcor_pvalues <- sapply(1:pathway_number, function(x) (1+(sum(surrogate_dcor[[x]]>dcor_results[x])))/(1+num_simulation))
names(nominal_dcor_pvalues) <- names(pathway)
return(nominal_dcor_pvalues)
}



MPEA_subsampling= function(expressions,y,
                             L=dim(expressions)[2],
                             pathway,
                             subsampling_size=10,
                             numsim1 = 1000,
                             numsim2 = 10000-1,
                             parallel_cores= detectCores()-2,
                           dcor_test_repetition=1000,
                             sim_mode = c("mean")){
  
  
  # expressions == expression matrix
  # y == vector of experiment design
  # pathway == vector of pathway gene indeces ( columns of expression matrix)
  # subsampling_size == Subsampling window size
  #numsim1 == pathway subsampling simulation runs
  #numsim2 == expression array subsampling simulation runs
  # dcor_test_repetition == R in dcor.test
  #  sim_mode== "mean" or"max". Either maximum or average dcor chosen for analysis.
  
  
  
  library(energy)
  library(doParallel)
  library(ggplot2)
  packageVersion("energy")

  # print( paste("Actual Pathway Dcor = ", dcor(expressions[,pathway],y)))

  pathway_dcor=dcor(expressions[,pathway],y)
  pathway_dcor_test=dcor.test(expressions[,pathway],y,R = dcor_test_repetition)

  pathway_dcor_sim <- NULL
  for (ii in 1:numsim1){
    pathway_dcor_sim[ii] <-  dcor(expressions[,sample(pathway, size = subsampling_size)],y)
  }

  pathway_dcor_sim_output <- ifelse(sim_mode=="max",max(pathway_dcor_sim),mean(pathway_dcor_sim))
  pathway_dcor_sim_output
  
  
  
  cl<-makeCluster(parallel_cores)
  registerDoParallel(cl)
  background_dcor_sim <- foreach (ii = 1:numsim2, .combine = c, .export = c("L", "pathway", "expressions","subsampling_size","y","sim_mode"), .packages = "energy" ) %dopar% {
    # print(ii)
    aux_pathway <- sample(1:L, length(pathway)) # generate an aux pathway the same size as the actual pathway
    aux_dcor <- NULL
    for (jj in 1: numsim1){
      aux_dcor[jj] <-  dcor(expressions[,sample(aux_pathway, size = subsampling_size)],y)
    }
    ifelse(sim_mode=="max", max(aux_dcor), mean(aux_dcor))
    
  }
  stopCluster(cl)
 

  p_value <- (sum(background_dcor_sim>pathway_dcor_sim_output)+1)/(numsim2+1)
  return(list(pathway_dcor=pathway_dcor,pathway_dcor_test=pathway_dcor_test, simulated_p_value=p_value))
}


