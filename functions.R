#function - make table to export spillovers
spillovers_to_dataframe <- function(spillovers){
  main <- cbind(as.data.frame(spillovers[[1]])*100, from(spillovers)[[1]])
  last <- c(to(spillovers)[[1]], overall(spillovers)[[1]])
  spillover_table <- rbind(main, last)
  names(spillover_table)[length(names(spillover_table))] <- 'FROM others'
  rownames(spillover_table)[length(rownames(spillover_table))] <- 'TO others'
  return(spillover_table)
}


big_var_est <- function(data) {
  Model1 = constructModel(as.matrix(data), p = 1, struct = "Basic", gran = c(50, 50), 
                          VARX = list(), verbose = T, ownlambdas = FALSE)
  Model1Results = cv.BigVAR(Model1)
}

#gran: Two arguments that characterize the grid of penalty parameters. The first denotes 
#the depth of the grid and the second the number of gridpoints.


big_var_est_fixed_lambda <- function(data) {
  Model1 = constructModel(as.matrix(data), p = 1, struct = "Basic", gran = c(optimal_lambda), 
                          VARX = list(), verbose = T, ownlambdas = TRUE)
  Model1Results = cv.BigVAR(Model1)
}
