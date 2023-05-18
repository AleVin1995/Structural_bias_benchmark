# Calculate weighted mean sensitivity for Local Drop Out correction (LDO) [single sample and single chromosome assumed]
LDO_weighted_mean_subunit <- function(data,
                              params = list(
                                omega = 1e6,
                                decay_func = (function(x)exp(x)),
                                in_subunit_weight = 0,
                                side = "both"
                              )
){
  
  require(dplyr)
  require(magrittr)
  
  missing_colnames <- setdiff(c("DEPENDENCY_SCORE", "POSITION", "GENE_NAME"), colnames(data))
  if( length(missing_colnames) > 0) stop( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n" ) )
  if( !is.numeric(   data$DEPENDENCY_SCORE ) ) stop( "The class of data$DEPENDENCY_SCORE is not \'numeric\'. \n")
  if( !is.numeric(   data$POSITION ) )         stop( "The class of data$POSITION is not \'numeric\'. \n")
  if( !is.character( data$GENE_NAME ) )        stop( "The class of data$GENE_NAME is not \'character\'. \n")
  if( "ESSENTIAL" %in% colnames(data) ){
    if(!is.logical( data$ESSENTIAL)) stop( "The class of data$ESSENTIAL is not \'logical\'. \n")
  } else {
    warning("No ESSENTIAL column was found in the input data frame. All guides will be used for the weighted mean calculation. \n")
    data$ESSENTIAL <- FALSE
  }
  
  if( length(setdiff(c("omega","decay_func","in_subunit_weight","side"),names(params))) != 0){
    stop( paste0( "Missing parameter variables: ", paste0(setdiff(c("omega","decay_func","in_subunit_weight","side"),names(params)), collapse = ", "), ". \n") )
  }
  
  
  weight.matrix <- params$decay_func(matrix(rep(data$POSITION, length(data$POSITION)), length(data$POSITION)) %>% add(-t(.)) %>% abs %>% divide_by(-params$omega))
  if( params$side == "left"){
    weight.matrix[upper.tri(weight.matrix, diag = FALSE)] <- 0 
  } else if( params$side == "right"){
    weight.matrix[lower.tri(weight.matrix, diag = FALSE)] <- 0
  }
  
  groups.mat <- matrix(rep(data$GENE_NAME, length(data$GENE_NAME)), length(data$GENE_NAME))
  weight.matrix[groups.mat == t(groups.mat)] <- params$in_subunit_weight * weight.matrix[groups.mat == t(groups.mat)]
  
  weight.matrix[,data$ESSENTIAL] <- 0
  
  weighted.mean.dependencyScore <- (weight.matrix %*% data$DEPENDENCY_SCORE) %>% divide_by( rowSums(weight.matrix))
  
  weighted.mean.dependencyScore[is.nan(weighted.mean.dependencyScore)] <- 0
  
  return(data.frame(data, LDO_weighted_mean = weighted.mean.dependencyScore))
  
}


# Calculate weighted mean sensitivity for Local Drop Out correction (LDO)
LDO_weighted_mean <- function(data,
                              params = list(
                                omega = 1e6,
                                decay_func = (function(x)exp(x)),
                                in_subunit_weight = 0,
                                side = "both"
                              ),
                              verbose = TRUE
){
  
  missing_colnames <- setdiff( c("SAMPLE_NAME", "CHROMOSOME"), colnames(data))
  if( length(missing_colnames) == 2){
    warning( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n",
                    "All data points will be assumed to be coming from a single sample and a single chromosome.") )
    out <- LDO_weighted_mean_subunit(data, params)
    return(out)
  }
  if( length(missing_colnames) == 1){
    warning( paste0("The following column is missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n") )
    if( missing_colnames == "CHROMOSOME")  warning("All data points will be assumed to be coming from a single chromosome.")
    if( missing_colnames == "SAMPLE_NAME") warning("All data points will be assumed to be coming from a single sample.")
  }
  
  data_split <- data %>% split(., data[, c("SAMPLE_NAME", "CHROMOSOME")])
  
  res <- list()
  idx <- 1:length(data_split)
  
  for (i in idx){
    res[[i]] <- LDO_weighted_mean_subunit(data_split[[i]], params)
  }
  
  return(res)

}


# Calculate LDO regression tree (subunit, i.e. sample, chromsome)
LDO_regression_subunit <- function( data,
                            params = list(minbucket = 2, 
                                          minsubunits = 3, 
                                          cp_init = -3, 
                                          cp_iteration = 0.1)){
  
  require(rpart)
  require(tidyverse)
  
  params.default = list(minbucket = 2, 
                        minsubunits = 3, 
                        cp_init = -3, 
                        cp_iteration = 0.1)
  
  missing_colnames <- setdiff(c("DEPENDENCY_SCORE", "POSITION", "GENE_NAME"), colnames(data))
  if( length(missing_colnames) > 0) stop( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n" ) )
  if( !is.numeric(   data$DEPENDENCY_SCORE ) ) stop( "The class of data$DEPENDENCY_SCORE is not \'numeric\'. \n")
  if( !is.numeric(   data$POSITION ) )         stop( "The class of data$POSITION is not \'numeric\'. \n")
  if( !is.character( data$GENE_NAME ) )        stop( "The class of data$GENE_NAME is not \'character\'. \n")
  if( "ESSENTIAL" %in% colnames(data) ){
    if(!is.logical( data$ESSENTIAL)) stop( "The class of data$ESSENTIAL is not \'logical\'. \n")
    # remove_bin_tmp <- !data$ESSENTIAL
  } else {
    warning("No ESSENTIAL column was found in the input data frame. All guides will be used for the regression tree calculation. \n")
    # remove_bin_tmp <- rep( TRUE, nrow(data))
    data$ESSENTIAL <- FALSE
  }
  
  remove_bin_tmp <- !data$ESSENTIAL
  k <- params$cp_init
  iteration <- params$cp_iteration
  mean.guides.per.gene <- data$GENE_NAME %>% table %>% mean(na.rm = T)
  
  check.while <- TRUE
  while(check.while){
    
    # estimate regression tree on "non-essential" guides
    tree_tmp <- rpart(DEPENDENCY_SCORE ~ POSITION, data = data, subset = remove_bin_tmp, control = rpart.control(cp = 10^k, minbucket = params$minbucket * mean.guides.per.gene))
    LDO_fix <- predict(tree_tmp, data)
    
    # calculate number of genes in the smallest <leaf>
    tmp <- data[,c("POSITION", "GENE_NAME")] %>% mutate(LDO_fix = LDO_fix) %>% as.data.frame
    tmp.rle <- rle(tmp$LDO_fix)
    tmp.rle$values <- 1:length(tmp.rle$values) %>% set_names(names(tmp.rle$values))
    tmp$LDO_nEvents <- inverse.rle(tmp.rle)
    tmp <- tmp %>% group_by(LDO_nEvents) %>% 
      mutate(SUBUNIT = GENE_NAME) %>% 
      mutate( LDO_nGenes = (table(SUBUNIT) > 0.75*mean.guides.per.gene) %>% sum) %>%
      as.data.frame
    
    check.while <- min(tmp$LDO_nGenes, na.rm = T) <= params$minsubunits
    
    k <- k + iteration
  }
  
  data$LDO_fix <- LDO_fix
  data$final_cp_param <- k - iteration
  
  return(data)
}


# Calculate LDO regression tree
LDO_regression <- function( data,
                            params = list(minbucket = 2, 
                                          minsubunits = 3, 
                                          cp_init = -3, 
                                          cp_iteration = 0.1),
                            verbose = TRUE){

  require(rpart)
  require(tidyverse)

  missing_colnames <- setdiff( c("SAMPLE_NAME", "CHROMOSOME"), colnames(data))
  if( length(missing_colnames) == 2){
    warning( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n",
                    "All data points will be assumed to be coming from a single sample and a single chromosome.") )
    out <- LDO_regression_subunit(data, params)
    return(out)
  }
  if( length(missing_colnames) == 1){
    warning( paste0("The following column is missing from the input data frame: ", paste0(missing_colnames, collapse = ", "), ". \n") )
    if( missing_colnames == "CHROMOSOME")  warning("All data points will be assumed to be coming from a single chromosome.")
    if( missing_colnames == "SAMPLE_NAME") warning("All data points will be assumed to be coming from a single sample.")
  }
  
  group_byvar <- c("SAMPLE_NAME", "CHROMOSOME") 
  data <- data %>% group_by_at(group_byvar)
  
  data <- LDO_regression_subunit(data, params) %>%
    as.data.frame()
    
  return(data)
}


#' Calculate LDO copy number correction
#'
#'@param data A data frame in which to interpret the variables. Column names must include: DEPENDENCY_SCORE, POSITION, GENE_NAME, SAMPLE_NAME (optional), CHROMOSOME (optional).
#'@param params_weigthed_mean a list of parameters for the LDO weighted mean regression.
#'omega: width of the window used for the decay_func (based on POSITION values).
#'decay_func: the monotonuously decreasing weight function as a function of the POSITION values.
#'in_subunit_weight: a priori weight set to the guides targeting the same gene (i.e. GENE_NAME).
#'side: decay_function should be applied on the left, right or both side spanning from each guide.
#'@param params_regression a list of parameters for the LDO regression tree.
#'cp_init: Initial complexity parameter. Any split in the regression tree that does not decrease the overall lack of fit by a factor of cp is not attempted.
#'cp_iteration: complexity parameter iteration step applied to the initial cp_init.
#'minbucket: the minimum number of observations in any terminal <leaf> node.
#'minsubunits: the minimum number of subunits (e.g. genes) in each terminal <leaf> node.
#'@param verbose the gene level identifier, e.g. gene symbols
#'
#'@examples
#'data2 <- LDO(data)
LDO <- function( data,
                 params_weigthed_mean = list(
                   omega = 1e6,
                   decay_func = (function(x)exp(x)),
                   in_subunit_weight = 0,
                   side = "both"
                 ),
                 params_regression = list(minbucket = 2, 
                               minsubunits = 3, 
                               cp_init = -3, 
                               cp_iteration = 0.1),
                 verbose = FALSE){
  ##########
  colnames.data <- colnames(data)
  
  # Expand the list of potential pan-lethal using a weighted average method
  cat("Expand the list of potential pan-lethal using a weighted average method \n")
  data <- LDO_weighted_mean(data = data,
                            params = params_weigthed_mean,
                            verbose)
  
  data <- bind_rows(data)
  data <- data %>%
    mutate( DEPENDENCY_SCORE_tmp = DEPENDENCY_SCORE - LDO_weighted_mean ) %>%
    rename(ESSENTIAL1 = ESSENTIAL) %>%
    mutate( ESSENTIAL = ESSENTIAL1 | (abs(DEPENDENCY_SCORE_tmp) > quantile(abs(DEPENDENCY_SCORE_tmp),0.85)) ) %>%
    arrange(SAMPLE_NAME, CHROMOSOME, POSITION)
  cat("\n")
  
  # LDO regression tree
  cat("Calculate regression tree \n")
  data <- LDO_regression(data,
                          params = params_regression,
                         verbose) %>%
    rename(ESSENTIAL_LDO = ESSENTIAL) %>%
    mutate(DEPENDENCY_SCORE_LDO = DEPENDENCY_SCORE - LDO_fix)
  cat("\n")
  
  # remove temporary columns
  if( "ESSENTIAL" %in% colnames.data){
    data <- data %>% rename(ESSENTIAL = ESSENTIAL1) %>% 
      .[,c(colnames.data, "ESSENTIAL_LDO", "LDO_fix", "DEPENDENCY_SCORE_LDO")]
  } else {
    data <- data %>%
      .[,c(colnames.data, "ESSENTIAL_LDO", "LDO_fix", "DEPENDENCY_SCORE_LDO")]
  }
  
  data <- data %>% as.data.frame
  
  return(data)
  
}

######
#' Calculate GAM copy number correction
#'
#'@param data A data frame in which to interpret the variables, the columns used in formula and subunit need to be included
#'@param formula a formula indicating the variables to be used in the gam regression. By default: DEPENDENCY_SCORE ~ CNA
#'@param subunit the subunits in which the gam regression should be performed, by default SAMPLE_NAME is used
#'
#'@examples
#'data2 <- LDO(data, "SCORE", "POSITION", "GENESYMBOL", "ESSENTIAL", params = list( cp = 10^-3, minbucket = 10, minsubunits = 3))
GAM <- function(data, formula = "DEPENDENCY_SCORE ~ CNA + EXP", subunit = "SAMPLE_NAME", verbose = FALSE){
  ##########
  require(tidyverse)
  require(mgcv)
  
  formula <- formula %>% as.formula
  formula.terms <- attr( terms(formula), "factors") %>% rownames
  formula.terms.right <- attr( terms(formula), "factors") %>% colnames
  missing_colnames1 <- setdiff( formula.terms, colnames(data))
  if( length(missing_colnames1) > 0) stop( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames1, collapse = ", "), ". \n",
                                                  "or the format of the formula is faulty.") )
  
  missing_colnames2 <- setdiff( subunit, colnames(data))
  if( length(missing_colnames2) > 0){
    warning( paste0("The following columns are missing from the input data frame: ", paste0(missing_colnames2, collapse = ", "), ". \n",
                    "All data points will be assumed to be coming from the same model.") )
  } else {
    data <- data %>% group_by(.dots = subunit)
  }
  
  gam.formula <- paste0( formula.terms[1], " ~ ", paste0( "s(", formula.terms.right, ")", collapse = " + ")) %>% as.formula
  
  for( k in formula.terms.right ){
    data[ is.na(data[,k]), k] <- median( t( data[,k]), na.rm = T)
  }
  
  depout <- data %>% as.data.frame %>% gam( gam.formula, data = .)
  
  tmp <- depout %>% coef %>% names
  right.remove.idx <- formula.terms.right[-1] %>% lapply( function(x){ tmp %>% grep(x,.)}) %>% unlist %>% union(1)
  
  GAM_fix <- predict(depout, type = "lpmatrix")[, -right.remove.idx] %*% coef(depout)[-right.remove.idx] %>% as.vector
  data <- data.frame(data, GAM_fix)
    
  data <- data %>% mutate(DEPENDENCY_SCORE_GAM = DEPENDENCY_SCORE - GAM_fix) %>% as.data.frame
  
  return(data)
}