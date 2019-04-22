


2^10


k <- 6
replicas <- 1
params <- paste("X", 1:k, sep = "")
setting <- computations(16, 6)


model <- A1





  




















new_estimator <- function(sample.size, type, k, params, replicas, model) {
  col_names <- c("A", paste("X", 1:k, sep = ""))
  # Define the settings
  setting <- computations(x = sample.size, k)
  # Create the matrices
  A <- scrambled_replicas(setting$Nb, k, replicas)
  names(A) <- 1:replicas
  # Separate the matrices
  X <- lapply(A, function(x) separate_matrices(x))
  # Name the matrices
  for(i in names(X)) {
    names(X[[i]]) <- col_names
  }
  # Name the columns
  for(i in names(X)) {
    for(j in names(X[[i]])) {
      colnames(X[[i]][[j]]) <- params
    }
  }
  # Compute model output
  Y <- lapply(X, function(x) lapply(x, function(y) model(y)))
  if(type == "old") { # RUN THE TRADITIONAL APPROACH --------------
    ###############################################################
    # Compute STi following the traditional approach
    out <- lapply(Y, function(x) sobol_Ti(do.call(c, x), params))
    # Finalize
    final <- rbindlist(out, idcol = "replica") %>%
      .[, sample.size := setting$Nc] %>%
      .[, algorithm:= "old"]
  } 
  if(type == "new") {# RUN THE NEW ALGORITHM ---------------------
    ###############################################################
    
    # EXTRACT THE WARM-UP SAMPLES ----------------------------------
    
    # Retrieve model output for the warm-up samples
    Y.warmup <- lapply(Y, function(x) lapply(x, function(y) y[1:setting$Nts]))
    # Compute Sobol' Ti for the warm-up samples
    out <- lapply(Y.warmup, function(x) sobol_Ti(do.call(c, x), params))
    # Retrieve the fixed parameters
    STi.fixed <- lapply(out, function(x) 
      x[V1 < quantile(V1, probs = 1 - 75 / 100)])
    # Retrieve 1/4 of the non-fixed parameters
    STi.non.fixed <- lapply(out, function(x) 
      x[V1 > quantile(V1, probs = 1 - 75 / 100)]) %>%
      lapply(., function(x) c("A", x[, parameters]))
    # Retrieve model runs of the non.fixed parameters I
    Y.saved <- list()
    for(i in names(Y)) {
      Y.saved[[i]] <-  Y[[i]][STi.non.fixed[[i]]]
    }
    # RUN THE EXTRA RUNS ------------------------------------------
    
    # Retrieve the sample matrices with the non-fixed parameters that 
    # will receive the extra runs
    X.extra.runs <- list()
    for(i in names(X)) {
      X.extra.runs[[i]] <- X[[i]][-1]
      X.extra.runs[[i]] <- X.extra.runs[[i]][STi.non.fixed[[i]][-1]]
    }
    # Extract the rows where the extra model runs will be computed
    AX <- lapply(X.extra.runs, function(x) 
      lapply(x, function(y) y[c(5:((5+setting$Nextra) - 1)), ]))
    # Create the new points in each column of the non-fixed parameters
    for(i in seq_along(AX)) {
      for(j in names(AX[[i]])) {
        AX[[i]][[j]][, j] <- runif(setting$Nextra)
      }
    }
    # Run the extra runs
    Y.extra.runs <- lapply(AX, function(x) lapply(x, function(y) model(y)))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # Run the extra runs
    Y.extra.runs <- lapply(AX, function(x) lapply(x, function(y) model(y)))
    
    
    
    
    
    da <- lapply(Y.extra.runs, function(x) do.call(cbind, x))[[1]] %>%
      data.table() %>%
      melt(., measure.vars = 1:ncol(.), 
           variable.name = "parameters",
           value.name = "Y_AX")
  
    
    # unite the vectors of the same simulation
    Y <- lapply(Y.saved, function(x) c(do.call(cbind, x)))[[1]]
    
    # Retrieve vector with the name of the parameters
    # with the extra runs
    params <- lapply(STi.non.fixed, function(x) x[-1])[[1]]
    
   
    
     # RE-ARRANGE COMPUTATION OF STI
    k <- length(new[[1]])
    # Calculate the length of the A matrix
    p <- length(1:(length(Y) / (k + 1)))
    # Extract the model output of the A matrix
    Y_A <- Y[1:p]
    # Extract the model output of the AB matrix
    Y_AB <- Y[(p+1):length(Y)]
    # Create vector with parameters
    parameters <- rep(params, each = length(Y_A))
    # merge vector with data table
    vec <- cbind(Y_A, Y_AB)
    out <- data.table::data.table(vec, parameters)
    
    
    
    merge(out, da, by = "parameters", all = TRUE)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    # MERGE THE WARM UP, THE SAVED AND THE EXTRA RUNS -----------------
    
    # Concatenate the warm-up, the saved and the extra model runs
    # for the non-fixed parameters
    all <- list()
    for(i in names(Y.saved)) {
      for(j in names(Y.saved[[i]])) {
        all[[i]][[j]] <- c(Y.saved[[i]][[j]],
                           Y.extra.runs[[i]][[j]])
      }
    }
    # Add NA values to the A matrix until it reaches length AB
    pop <- list()
    for(i in names(all)) {
      pop[[i]] <- lapply(all[[1]], function(x) { 
        length(x) <- max(lengths(all[[1]])); x })
    }
    # unite the vectors of the same simulation
    Y <- lapply(pop, function(x) c(do.call(cbind, x)))
    
    # COMPUTE STI FOR THE NON-FIXED PARAMETERS -----------------------
    
    # Retrieve vector with the name of the parameters
    # with the extra runs
    new <- lapply(STi.non.fixed, function(x) x[-1])
    # Compute Sobol' Ti for the extra runs
    out <- list()
    for(i in names(Y)) {
      out[[i]] <- sobol_Ti(Y[[i]], params = new[[i]])
    }
    
    # MERGE STI FOR THE FIXED AND NON-FIXED PARAMETERS ---------------
    
    # Cbind the fixed parameters
    all.together <- list()
    for(i in names(out)) {
      all.together[[i]] <- rbind(out[[i]], STi.fixed[[i]]) %>%
        .[order(parameters)]
    }
    # Finalize
    final <- rbindlist(all.together, idcol = "replica") %>%
      .[, sample.size := setting$Nc] %>%
      .[, algorithm:= "new"]
  }
  return(final)
}






























computations(16, 6)

# DEFINE GENERAL SETTINGS -----------------------------------------------------

# Create vector with the name of the test functions
test_functions <- c("A1", "A2", "B1", "B2", "B3", "C1", "C2")

# Vector power of two
x <- seq(4, 13, 1)

# Get the initial sample sizes
Nb <- sapply(x, function(x) 2 ^ x)

# Set number of factors
k <- 6
params <- paste("X", 1:k, sep = "")

# Set number of sample matrix replicas
replicas <- 50
    
# RUN THE MODEL ---------------------------------------------------------------

out <- list()
run_model <- as.list(c(test_functions))
names(run_model) <- test_functions
estimators <- c("new", "old")
for(i in names(run_model)) {
  if(i == "A1") {
    test_F <- A1
  } else if(i == "A2") {
    test_F <- A2
  } else if(i == "B1") {
    test_F <- B1
  } else if(i == "B2") {
    test_F <- B1
  } else if(i == "B3") {
    test_F <- B3
  } else if(i == "C1") {
    test_F <- C1
  } else if(i == "C2") {
    test_F <- C2
  }
  out[[i]] <- lapply(estimators, function(x)
    lapply(Nb, function(Nb) new_estimator(sample.size = Nb, 
                                          type = x, 
                                          k = k, 
                                          params = params, 
                                          replicas = replicas, 
                                          model = test_F)))
}   
    
# COMPUTE MAE -----------------------------------------------------------------

# Arrange data
temp <- lapply(out, function(x) 
  lapply(x, function(y) rbindlist(y))) %>%
  lapply(., function(x) rbindlist(x))

# Merge indices with analytical values
for(i in names(temp)) {
  temp[[i]] <- cbind(temp[[i]], analytical[[i]])
}

# Compute the MAE
MAE <- rbindlist(temp, idcol = "Function") %>%
  setnames(., c("V1", "value"), c("estimated", "analytical")) %>%
  .[, .(MAE = mean(abs(estimated - analytical))), 
    .(Function, sample.size, algorithm)] %>%
  .[, sample.size:= as.numeric(sample.size)]   

# PLOT MAE --------------------------------------------------------------------

ggplot(MAE, aes(sample.size, MAE, color = algorithm)) +
  geom_point() +
  geom_line() +
  scale_color_discrete(name = expression(S[Ti]), 
                       labels = c("New algorithm", "Jansen 1999")) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Total cost", 
       y = "MAE") +
  facet_wrap(~Function, 
             ncol = 4) +
  theme_bw() +
  theme(aspect.ratio = 1, 
        legend.position = "top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))

























  
  
  
  
  # EXTRACT THE WARM-UP SAMPLES ----------------------------------
  
  # Retrieve model output for the warm-up samples
  Y.warmup <- lapply(Y, function(x) lapply(x, function(y) y[1:setting$Nts]))
  # Compute Sobol' Ti for the warm-up samples
  out <- lapply(Y.warmup, function(x) sobol_Ti(do.call(c, x), params))
  # Retrieve the fixed parameters
  STi.fixed <- lapply(out, function(x) 
    x[V1 < quantile(V1, probs = 1 - 75 / 100)])
  # Retrieve 1/4 of the non-fixed parameters
  STi.non.fixed <- lapply(out, function(x) 
    x[V1 > quantile(V1, probs = 1 - 75 / 100)]) %>%
    lapply(., function(x) c("A", x[, parameters]))
  # Retrieve model runs of the non.fixed parameters I
  Y.saved <- list()
  for(i in names(Y)) {
    Y.saved[[i]] <-  Y[[i]][STi.non.fixed[[i]]]
  }
  # RUN THE EXTRA RUNS ------------------------------------------
  
  # Retrieve the sample matrices with the non-fixed parameters that 
  # will receive the extra runs
  X.extra.runs <- list()
  for(i in names(X)) {
    X.extra.runs[[i]] <- X[[i]][-1]
    X.extra.runs[[i]] <- X.extra.runs[[i]][STi.non.fixed[[i]][-1]]
  }
  # Extract the rows where the extra model runs will be computed
  AX <- lapply(X.extra.runs, function(x) 
    lapply(x, function(y) y[c(5:((5+setting$Nextra) - 1)), ]))
  # Create the new points in each column of the non-fixed parameters
  for(i in seq_along(AX)) {
    for(j in names(AX[[i]])) {
      AX[[i]][[j]][, j] <- runif(setting$Nextra)
    }
  }
  # Run the extra runs
  Y.extra.runs <- lapply(AX, function(x) lapply(x, function(y) A1(y)))
  
  # MERGE THE WARM UP, THE SAVED AND THE EXTRA RUNS -----------------
  
  # Concatenate the warm-up, the saved and the extra model runs
  # for the non-fixed parameters
  all <- list()
  for(i in names(Y.saved)) {
    for(j in names(Y.saved[[i]])) {
      all[[i]][[j]] <- c(Y.saved[[i]][[j]],
                         Y.extra.runs[[i]][[j]])
    }
  }
  # Add NA values to the A matrix until it reaches length AB
  pop <- list()
  for(i in names(all)) {
    pop[[i]] <- lapply(all[[1]], function(x) { 
      length(x) <- max(lengths(all[[1]])); x })
  }
  # unite the vectors of the same simulation
  Y <- lapply(pop, function(x) c(do.call(cbind, x)))
  
  # COMPUTE STI FOR THE NON-FIXED PARAMETERS -----------------------
  
  # Retrieve vector with the name of the parameters
  # with the extra runs
  new <- lapply(STi.non.fixed, function(x) x[-1])
  # Compute Sobol' Ti for the extra runs
  out <- list()
  for(i in names(Y)) {
    out[[i]] <- sobol_Ti(Y[[i]], params = new[[i]])
  }
  
  # MERGE STI FOR THE FIXED AND NON-FIXED PARAMETERS ---------------
  
  # Cbind the fixed parameters
  all.together <- list()
  for(i in names(out)) {
    all.together[[i]] <- rbind(out[[i]], STi.fixed[[i]]) %>%
      .[order(parameters)]
  }
  # Finalize
  final <- rbindlist(all.together, idcol = "replica") %>%
    .[, sample.size := setting$Nc] %>%
    .[, algorithm:= "new"]
}













# Create vector for the A and AB matrices
  col_names <- c("A", paste("X", 1:k, sep = ""))
  # Define the settings
  setting <- computations(x = sample.size, k)
  # Create the matrices
  A <- scrambled_replicas(setting$Nb, k, replicas)
  names(A) <- 1:replicas
  # Separate the matrices
  X <- lapply(A, function(x) separate_matrices(x))
  # Name the matrices
  for(i in names(X)) {
    names(X[[i]]) <- col_names
  }
  # Name the columns
  for(i in names(X)) {
    for(j in names(X[[i]])) {
      colnames(X[[i]][[j]]) <- params
    }
  }
  # Compute model output
  Y <- lapply(X, function(x) lapply(x, function(y) A1(y)))
  
  # EXTRACT THE WARM-UP SAMPLES
  ###############################################################
  
  # Retrieve model output for the warm-up samples
  Y.warmup <- lapply(Y, function(x) lapply(x, function(y) y[1:setting$Nts]))
  # Compute Sobol' Ti for the warm-up samples
  out <- lapply(Y.warmup, function(x) sobol_Ti(do.call(c, x), params))
  # Retrieve the fixed parameters
  STi.fixed <- lapply(out, function(x) 
    x[V1 < quantile(V1, probs = 1 - 75 / 100)])
  # Retrieve 1/4 of the non-fixed parameters
  STi.non.fixed <- lapply(out, function(x) 
    x[V1 > quantile(V1, probs = 1 - 75 / 100)]) %>%
    lapply(., function(x) c("A", x[, parameters]))
  # Retrieve model runs of the non.fixed parameters I
  Y.saved <- list()
  for(i in names(Y)) {
    Y.saved[[i]] <-  Y[[i]][STi.non.fixed[[i]]]
  }
  # RUN THE EXTRA RUNS
  ###############################################################
  
  # Retrieve the sample matrices with the non-fixed parameters that 
  # will receive the extra runs
  X.extra.runs <- list()
  for(i in names(X)) {
    X.extra.runs[[i]] <- X[[i]][-1]
    X.extra.runs[[i]] <- X.extra.runs[[i]][STi.non.fixed[[i]][-1]]
  }
  # Extract the rows where the extra model runs will be computed
  AX <- lapply(X.extra.runs, function(x) 
    lapply(x, function(y) y[c(5:((5+setting$Nextra) - 1)), ]))
  # Create the new points in each column of the non-fixed parameters
  for(i in seq_along(AX)) {
    for(j in names(AX[[i]])) {
      AX[[i]][[j]][, j] <- runif(setting$Nextra)
    }
  }
  # Run the extra runs
  Y.extra.runs <- lapply(AX, function(x) lapply(x, function(y) A1(y)))
  
  # MERGE THE WARM UP, THE SAVED AND THE EXTRA RUNS 
  ###################################################################
  
  # Concatenate the warm-up, the saved and the extra model runs
  # for the non-fixed parameters
  all <- list()
  for(i in names(Y.saved)) {
    for(j in names(Y.saved[[i]])) {
      all[[i]][[j]] <- c(Y.saved[[i]][[j]],
                         Y.extra.runs[[i]][[j]])
    }
  }
  # Add NA values to the A matrix until it reaches length AB
  pop <- list()
  for(i in names(all)) {
    pop[[i]] <- lapply(all[[1]], function(x) { 
      length(x) <- max(lengths(all[[1]])); x })
  }
  # unite the vectors of the same simulation
  Y <- lapply(pop, function(x) c(do.call(cbind, x)))
  
  # COMPUTE STI FOR THE NON-FIXED PARAMETERS
  #####################################################################
  
  # Retrieve vector with the name of the parameters
  # with the extra runs
  new <- lapply(STi.non.fixed, function(x) x[-1])
  # Compute Sobol' Ti for the extra runs
  out <- list()
  for(i in names(Y)) {
    out[[i]] <- sobol_Ti(Y[[i]], params = new[[i]])
  }
  
  # MERGE STI FOR THE FIXED AND NON-FIXED PARAMETERS
  #####################################################################
  
  # Cbind the fixed parameters
  all.together <- list()
  for(i in names(out)) {
    all.together[[i]] <- rbind(out[[i]], STi.fixed[[i]]) %>%
      .[order(parameters)]
  }
  # Finalize
  final <- rbindlist(all.together, idcol = "replica") %>%
    .[, sample.size := setting$Nc] %>%
    .[, algorithm:= "new"]
  
  
  


  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
# Create vector for the A and AB matrices
col_names <- c("A", paste("X", 1:k, sep = ""))

# Define the settings
setting <- computations(x = sample.size, k)

# Create the matrices
A <- scrambled_replicas(setting$Nb, k, replicas)

names(A) <- 1:replicas

# Separate the matrices
X <- lapply(A, function(x) separate_matrices(x))

# Name the matrices
for(i in names(X)) {
  names(X[[i]]) <- col_names
}

# Name the columns
for(i in names(X)) {
  for(j in names(X[[i]])) {
    colnames(X[[i]][[j]]) <- params
  }
}

# Compute model output
Y <- lapply(X, function(x) lapply(x, function(y) A1(y)))

# EXTRACT THE WARM-UP SAMPLES
###############################################################

# Retrieve model output for the warm-up samples
Y.warmup <- lapply(Y, function(x) lapply(x, function(y) y[1:setting$Nts]))

# Compute Sobol' Ti for the warm-up samples
out <- lapply(Y.warmup, function(x) sobol_Ti(do.call(c, x), params))

# Retrieve the fixed parameters
STi.fixed <- lapply(out, function(x) 
  x[V1 < quantile(V1, probs = 1 - 75 / 100)])

# Retrieve 1/4 of the non-fixed parameters
STi.non.fixed <- lapply(out, function(x) 
  x[V1 > quantile(V1, probs = 1 - 75 / 100)]) %>%
  lapply(., function(x) c("A", x[, parameters]))

# Retrieve model runs of the non.fixed parameters I
Y.saved <- list()
for(i in names(Y)) {
  Y.saved[[i]] <-  Y[[i]][STi.non.fixed[[i]]]
}

# RUN THE EXTRA RUNS
###############################################################

# Retrieve the sample matrices with the non-fixed parameters that 
# will receive the extra runs
X.extra.runs <- list()
for(i in names(X)) {
  X.extra.runs[[i]] <- X[[i]][-1]
  X.extra.runs[[i]] <- X.extra.runs[[i]][STi.non.fixed[[i]][-1]]
}

# Extract the rows where the extra model runs will be computed
AX <- lapply(X.extra.runs, function(x) 
  lapply(x, function(y) y[c(5:((5+setting$Nextra) - 1)), ]))

# Create the new points in each column of the non-fixed parameters
for(i in seq_along(AX)) {
  for(j in names(AX[[i]])) {
    AX[[i]][[j]][, j] <- runif(setting$Nextra)
  }
}

# Run the extra runs
Y.extra.runs <- lapply(AX, function(x) lapply(x, function(y) A1(y)))

# MERGE THE WARM UP, THE SAVED AND THE EXTRA RUNS 
###################################################################

# Concatenate the warm-up, the saved and the extra model runs
# for the non-fixed parameters
all <- list()
for(i in names(Y.saved)) {
  for(j in names(Y.saved[[i]])) {
    all[[i]][[j]] <- c(Y.saved[[i]][[j]],
                        Y.extra.runs[[i]][[j]])
  }
}

# Add NA values to the A matrix until it reaches length AB
pop <- list()
for(i in names(all)) {
  pop[[i]] <- lapply(all[[1]], function(x) { 
    length(x) <- max(lengths(all[[1]])); x })
}

# unite the vectors of the same simulation
Y <- lapply(pop, function(x) c(do.call(cbind, x)))

# COMPUTE STI FOR THE NON-FIXED PARAMETERS
#####################################################################

# Retrieve vector with the name of the parameters
# with the extra runs
new <- lapply(STi.non.fixed, function(x) x[-1])

# Compute Sobol' Ti for the extra runs
out <- list()
for(i in names(Y)) {
  out[[i]] <- sobol_Ti(Y[[i]], params = new[[i]])
}

# Cbind the fixed parameters
all.together <- list()
for(i in names(out)) {
  all.together[[i]] <- rbind(out[[i]], STi.fixed[[i]]) %>%
    .[order(parameters)]
}

# Finalize
final <- rbindlist(all.together, idcol = "replica") %>%
  .[, sample.size := setting$Nc] %>%
  .[, algorithm:= "new"]











sobol_compute_Ti <- function(Y_A, Y_AB) {
  n <- length(Y_A[!is.na(Y_A)])
  f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
  VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
  STi <- ((1 / (2 * n)) * sum((Y_A - Y_AB) ^ 2, na.rm = TRUE)) / VY
  return(STi) 
}


sobol_Mapply_Ti <- function(d) {
  return(mapply(sobol_compute_Ti,
                d[, "Y_A"],
                d[, "Y_AB"]))
}

sobol_Ti <- function(Y, params) {
  # Calculate the number of parameters
  k <- length(params)
  # Calculate the length of the A matrix
  p <- length(1:(length(Y) / (k + 1)))
  # Extract the model output of the A matrix
  Y_A <- Y[1:p]
  # Extract the model output of the AB matrix
  Y_AB <- Y[(p+1):length(Y)]
  # Create vector with parameters
  parameters <- rep(params, each = length(Y_A))
  # merge vector with data table
  vec <- cbind(Y_A, Y_AB)
  out <- data.table::data.table(vec, parameters)
  out.1 <- out %>%
    # remove rows with NA
    stats::na.omit()
  # Compute Sobol'indices
  output <- out.1[, sobol_Mapply_Ti(.SD), by = parameters]
  return(output)
}


N <- 1000
k <- 6
params <- paste("X", 1:k, sep = "")
A <- scrambled_replicas(N, k, 1)[[1]]
Y <- B1(A)
sobol_Ti(Y, params)




sample.size <- 2^4

settings <- computations(sample.size, 6)
k <- 6
params <- paste("X", 1:k, sep = "")
model <- A1
replicas <- 10

A <- scrambled_replicas(settings$Nb, k = k, 2)





new_algorithm <- function(sample.size, params, type, replicas, model) {
  settings <- computations(sample.size, length(params))
  A <- scrambled_replicas(settings$Nb, k = length(params), replicas)
  out <- list()
  for(i in seq_along(A)) {
    out[[i]] <- new_estimator(A = A[[i]], 
                              type = type, 
                              params = params, 
                              model = model)
  }
return(out)
}



sample.size <- 2^4
params <- paste("X", 1:k, sep = "")
model <- A1
replicas <- 10

new_algorithm(sample.size = sample.size, 
              params = params, 
              type = "new", 
              replicas = replicas, 
              model = model)






sample.size <- 2^4
params <- paste("X", 1:k, sep = "")
model <- A1
replicas <- 10




new_estimator <- function(A, sample.size, type, params, model) {
  settings <- computations(sample.size, length(params))
  # Compute model output on the initial sample size
  Y <- model(A) 
  if(type == "old") { # RUN THE TRADITIONAL APPROACH --------------
    ###############################################################
    output <- sobol_Ti(Y = Y, params = params)
  }
  if(type == "new") { # RUN THE NEW ALGORITHM ---------------------
    ###############################################################
    # Add a column with the parameters and the model output
    A <- data.table(A)[, Y:= cbind(Y)] %>%
      .[, parameters:= rep(c("A", params), each = settings$Nb)] %>%
      setnames(., paste("V", 1:k, sep = ""), paste("X", 1:k, sep = ""))
    
    # STi ON THE WARM UP SAMPLE ---------------------------------------------------
    STi.warmup <- sobol_Ti(A[, .SD[1:settings$Nts], parameters][, Y], params)
    # Retrieve a vector with the non-fixed parameters
    STi.non.fixed <- STi.warmup[V1 > quantile(V1, probs = 1 - 75 / 100)][, parameters]
    # Vector with the A and the non-fixed parameters
    A.STi.non.fixed <- c("A", STi.non.fixed)
    # Retrieve the STi indices of the fixed parameters
    STi.fixed <- STi.warmup[V1 < quantile(V1, probs = 1 - 75 / 100)]
    
    # RETRIEVE ALL RUNS FOR THE MOST IMPORTANT PARAMETERS---------------------------
    YA <- A[parameters %in% A.STi.non.fixed]
    Y <- YA[, Y]
    
    # RUN EXTRA RUNS ---------------------------------------------------------------
    # Retrieve 1/4 of the rows of the AB matrices of the most
    # important parameters to compute the extra model runs
    AX <- A[parameters %in% STi.non.fixed][
      , .SD[1:((1+settings$Nextra) - 1)], parameters][, !"Y"] 
    # Substitute the column of the j parameter for a random numer
    set.seed(666)
    for(j in colnames(AX[, 2:6])) {
      AX[parameters == j, j] <- runif(settings$Nextra)
    }
    # Compute model output
    Y.extra <- model(as.matrix(AX[, -1]))
    # Add model output
    AX <- data.table(AX)[, Y_AX:= cbind(Y.extra)]
    # Create a vector with the new points located at the proper place
    AY.extra <- AX[, c(Y_AX, rep(NA, settings$Nts+settings$Nextra)), parameters]
    # ADD NA to cover the A matrix
    temp <- rbind(data.table(parameters = "A", V1 = rep(NA, settings$Nb)), 
                  AY.extra) %>%
      .[, .(V1)] %>%
      setnames(., "V1", "Y_AX")
    # Bind
    Y_AX <- cbind(YA, temp)[!parameters == "A"][, Y_AX]
    Yfinal <- c(Y, Y_AX)
    
    # COMPUTE STI -------------------------------------------------------------
    # Calculate the number of parameters
    k <- length(STi.non.fixed)
    # Extract the model output of the A matrix
    p <- length(1:(length(Yfinal) / ((2 * k) + 1)))
    # Extract the model output of the A matrix
    Y_A <- Yfinal[1:p]
    # Extract the model output of the AB matrix
    Y_AB <- Yfinal[(p+1): (p * (k + 1))]
    # Extract the model output of the AX matrix
    Y_AX <- Yfinal[(p * (k + 1) + 1): length(Yfinal)]
    # Create vector with parameters
    parameters <- rep(STi.non.fixed, each = length(Y_A))
    # merge vector with data table
    vec <- cbind(Y_A, Y_AB, Y_AX)
    out <- data.table::data.table(vec, parameters)
    
    # Define parameters
    f0 <- (1 / length(Y_A[!is.na(Y_A)])) * sum(Y_A, na.rm = TRUE)
    VY <- 1 / length(Y_A[!is.na(Y_A)]) * sum((Y_A - f0) ^ 2, na.rm = TRUE)
    
    # Compute
    first <- out[, .(parameters, Y_A, Y_AB)] %>%
      setnames(., c("Y_A", "Y_AB"), c("one", "two"))
    second <- out[, .(parameters, Y_A, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_A", "Y_AX"), c("one", "two"))
    third <- out[, .(parameters, Y_AB, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_AB", "Y_AX"), c("one", "two"))
    
    temp <- rbind(first, second, third) %>%
      .[order(parameters)] %>%
      .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters]
    # Finalize
    output <- rbind(STi.fixed, temp)[order(parameters)]
  }
  return(output)
}




A <- scrambled_replicas(settings$Nb, k = k, 1)[[1]]
sample.size <- 16
k <- 6
params <- paste("X", 1:k, sep = "")
model <- A1

new_estimator(A, sample.size = sample.size, type = "old",  params = params, model = B1)












sample.size <- 2^10

settings <- computations(sample.size, 6)
k <- 6
params <- paste("X", 1:k, sep = "")
model <- A1
replicas <- 10

A <- scrambled_replicas(settings$Nb, k = k, 1)[[1]]



new_estimator(A, sample.size = 1000, type = "new",  params = params, model = B1)



  





























############################## THIS WORKS


sample.size <- 2^8
k <- 6
params <- paste("X", 1:k, sep = "")
model <- A1
replicas <- 30

A <- scrambled_replicas(settings$Nb, k = k, 1)[[1]]

new_estimator(A, sample.size = sample.size, params = params, model = model, type = "new")
































new_estimator <- function(A, sample.size, params, model, type) {
  settings <- computations(sample.size, length(params))
  # Compute model output on the initial sample size
  Y.initial <- model(A) 
  if(type == "old") { # RUN THE TRADITIONAL APPROACH --------------
    ###############################################################
    output <- sobol_Ti(Y.initial, params) %>%
      .[, model.runs := settings$Nc] %>%
      .[, algorithm:= "old"]
  }
  if(type == "new") { # RUN THE NEW ALGORITHM ---------------------
    ###############################################################
    # Add a column with the parameters and the model output
    A <- data.table(A)[, Y:= cbind(Y.initial)] %>%
      .[, parameters:= rep(c("A", params), each = settings$Nb)] %>%
      setnames(., paste("V", 1:k, sep = ""), paste("X", 1:k, sep = ""))
    
    # STi ON THE WARM UP SAMPLE ---------------------------------------------------
    STi.warmup <- sobol_Ti(A[, .SD[1:settings$Nts], parameters][, Y], params)
    # Retrieve a vector with the non-fixed parameters
    STi.non.fixed <- STi.warmup[V1 > quantile(V1, probs = 1 - 75 / 100)][, parameters]
    # Vector with the A and the non-fixed parameters
    A.STi.non.fixed <- c("A", STi.non.fixed)
    # Retrieve the STi indices of the fixed parameters
    STi.fixed <- STi.warmup[V1 < quantile(V1, probs = 1 - 75 / 100)]
    
    # RETRIEVE ALL RUNS FOR THE MOST IMPORTANT PARAMETERS---------------------------
    YA <- A[parameters %in% A.STi.non.fixed]
    Y <- YA[, Y]
    
    # RUN EXTRA RUNS ---------------------------------------------------------------
    # Retrieve 1/4 of the rows of the AB matrices of the most
    # important parameters to compute the extra model runs
    AX <- A[parameters %in% STi.non.fixed][
      , .SD[1:((1+settings$Nextra) - 1)], parameters][, !"Y"] 
    # Substitute the column of the j parameter for a random numer
    set.seed(666)
    for(j in colnames(AX[, 2:6])) {
      AX[parameters == j, j] <- runif(settings$Nextra)
    }
    # Compute model output
    Y.extra <- model(as.matrix(AX[, -1]))
    # Add model output
    AX <- data.table(AX)[, Y_AX:= cbind(Y.extra)]
    # Create a vector with the new points located at the proper place
    AY.extra <- AX[, c(Y_AX, rep(NA, settings$Nts+settings$Nextra)), parameters]
    # ADD NA to cover the A matrix
    temp <- rbind(data.table(parameters = "A", V1 = rep(NA, settings$Nb)), 
                  AY.extra) %>%
      .[, .(V1)] %>%
      setnames(., "V1", "Y_AX")
    # Bind
    Y_AX <- cbind(YA, temp)[!parameters == "A"][, Y_AX]
    Yprove <- c(Y, Y_AX)
    
    # Calculate the number of parameters
    k <- length(STi.non.fixed)
    # Extract the model output of the A matrix
    p <- length(1:(length(Yprove) / ((2 * k) + 1)))
    # Extract the model output of the A matrix
    Y_A <- Y[1:p]
    # Extract the model output of the AB matrix
    Y_AB <- Y[(p+1): (p * (k + 1))]
    # Extract the model output of the AX matrix
    Y_AX <- Yprove[(p * (k + 1) + 1): length(Yprove)]
    # Create vector with parameters
    parameters <- rep(STi.non.fixed, each = length(Y_A))
    # merge vector with data table
    vec <- cbind(Y_A, Y_AB, Y_AX)
    out <- data.table::data.table(vec, parameters)
    
    # Define parameters
    n <- length(Y_A[!is.na(Y_A)])
    f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
    VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
    
    # Compute
    first <- out[, .(parameters, Y_A, Y_AB)] %>%
      setnames(., c("Y_A", "Y_AB"), c("one", "two"))
    second <- out[, .(parameters, Y_A, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_A", "Y_AX"), c("one", "two"))
    third <- out[, .(parameters, Y_AB, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_AB", "Y_AX"), c("one", "two"))
    # Calculate output
    output <- rbind(first, second, third) %>%
      .[order(parameters)] %>%
      .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters] %>%
      rbind(., STi.fixed) %>%
      .[order(parameters)] %>%
      .[, model.runs := settings$Nc] %>%
      .[, algorithm:= "new"]
  }
  return(output)
}






new_algorithm <- function(sample.size, params, model, type, replicas) {
  settings <- computations(sample.size, length(params))
  A <- scrambled_replicas(settings$Nb, k = length(params), replicas)
  out <- lapply(A, function(x) new_estimator(x, sample.size, params, model, type))
  return(out)
}



new_algorithm(16, params = params, model = A1, type = "old", replicas)





# Vector power of two
x <- seq(4, 7, 1)

# Get the initial sample sizes
Nb <- sapply(x, function(x) 2 ^ x)


lapply(Nb, function(Nb) new_algorithm(Nb,  params = params, model = A1, type = "new", replicas))







# Create vector with the name of the test functions
test_functions <- c("A1", "A2", "B1", "B2", "B3", "C1", "C2")

# Vector power of two
x <- seq(4, 13, 1)

# Get the initial sample sizes
Nb <- sapply(x, function(x) 2 ^ x)

# Set number of factors
k <- 6
params <- paste("X", 1:k, sep = "")

# Set number of sample matrix replicas
replicas <- 50




# RUN THE MODEL ---------------------------------------------------------------

out <- list()
run_model <- as.list(c(test_functions))
names(run_model) <- test_functions
estimators <- c("new", "old")
for(i in names(run_model)) {
  if(i == "A1") {
    test_F <- A1
  } else if(i == "A2") {
    test_F <- A2
  } else if(i == "B1") {
    test_F <- B1
  } else if(i == "B2") {
    test_F <- B1
  } else if(i == "B3") {
    test_F <- B3
  } else if(i == "C1") {
    test_F <- C1
  } else if(i == "C2") {
    test_F <- C2
  }
  out[[i]] <- mclapply(estimators, function(x) 
    lapply(Nb, function(Nb) new_algorithm(Nb, 
                                          params = params, 
                                          model = test_F, 
                                          type = x, 
                                          replicas = replicas)), 
    mc.cores = detectCores() - 1)
}

# COMPUTE MAE -----------------------------------------------------------------


temp <- out

for(i in names(out)) {
  for(j in seq_along(out[[i]])){
    names(temp[[i]][[j]]) <- Nb
    for(k in seq_along(out[[i]][[j]])) {
      names(temp[[i]][[j]][[k]]) <- 1:replicas
    }
  }
}


temp <- lapply(temp, function(x) lapply(x, function(y) 
  lapply(y, function(z) rbindlist(z, idcol = "replicas")))) %>%
  lapply(., function(x) lapply(x, function(y) rbindlist(y, idcol = "N"))) %>%
  lapply(., function(x) rbindlist(x)) 

# Merge indices with analytical values
for(i in names(temp)) {
  temp[[i]] <- cbind(temp[[i]], analytical[[i]])
}

# Compute the MAE
MAE <- rbindlist(temp, idcol = "Function") %>%
  setnames(., c("V1", "value"), c("estimated", "analytical")) %>%
  .[, .(MAE = mean(abs(estimated - analytical))), 
    .(Function, sample.size, algorithm)] %>%
  .[, sample.size:= as.numeric(sample.size)]


# PLOT MAE --------------------------------------------------------------------

ggplot(MAE, aes(sample.size, MAE, color = algorithm)) +
  geom_point() +
  geom_line() +
  scale_color_discrete(name = expression(S[Ti]), 
                       labels = c("New algorithm", "Jansen 1999")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  labs(x = "Total cost", 
       y = "MAE") +
  facet_wrap(~Function, 
             ncol = 4) +
  theme_bw() +
  theme(aspect.ratio = 1, 
        legend.position = "top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(fill = "transparent",
                                         color = NA),
        legend.key = element_rect(fill = "transparent",
                                  color = NA))





































# Compute model output on the initial sample size
Y.initial <- model(A) 
# Add a column with the parameters and the model output
A <- data.table(A)[, Y:= cbind(Y.initial)] %>%
  .[, parameters:= rep(c("A", params), each = settings$Nb)] %>%
  setnames(., paste("V", 1:k, sep = ""), paste("X", 1:k, sep = ""))

# STi ON THE WARM UP SAMPLE ---------------------------------------------------
STi.warmup <- sobol_Ti(A[, .SD[1:settings$Nts], parameters][, Y], params)
# Retrieve a vector with the non-fixed parameters
STi.non.fixed <- STi.warmup[V1 > quantile(V1, probs = 1 - 75 / 100)][, parameters]
# Vector with the A and the non-fixed parameters
A.STi.non.fixed <- c("A", STi.non.fixed)
# Retrieve the STi indices of the fixed parameters
STi.fixed <- STi.warmup[V1 < quantile(V1, probs = 1 - 75 / 100)]

# RETRIEVE ALL RUNS FOR THE MOST IMPORTANT PARAMETERS---------------------------
YA <- A[parameters %in% A.STi.non.fixed]
Y <- YA[, Y]

# RUN EXTRA RUNS ---------------------------------------------------------------
# Retrieve 1/4 of the rows of the AB matrices of the most
# important parameters to compute the extra model runs
AX <- A[parameters %in% STi.non.fixed][
  , .SD[1:((1+settings$Nextra) - 1)], parameters][, !"Y"] 
# Substitute the column of the j parameter for a random numer
set.seed(666)
for(j in colnames(AX[, 2:6])) {
  AX[parameters == j, j] <- runif(settings$Nextra)
}
# Compute model output
Y.extra <- model(as.matrix(AX[, -1]))
# Add model output
AX <- data.table(AX)[, Y_AX:= cbind(Y.extra)]
# Create a vector with the new points located at the proper place
AY.extra <- AX[, c(Y_AX, rep(NA, settings$Nts+settings$Nextra)), parameters]
# ADD NA to cover the A matrix
temp <- rbind(data.table(parameters = "A", V1 = rep(NA, settings$Nb)), 
      AY.extra) %>%
  .[, .(V1)] %>%
  setnames(., "V1", "Y_AX")
# Bind
Y_AX <- cbind(YA, temp)[!parameters == "A"][, Y_AX]
Yprove <- c(Y, Y_AX)

# Calculate the number of parameters
k <- length(STi.non.fixed)
# Extract the model output of the A matrix
p <- length(1:(length(Yprove) / ((2 * k) + 1)))
# Extract the model output of the A matrix
Y_A <- Y[1:p]
# Extract the model output of the AB matrix
Y_AB <- Y[(p+1): (p * (k + 1))]
# Extract the model output of the AX matrix
Y_AX <- Yprove[(p * (k + 1) + 1): length(Yprove)]
# Create vector with parameters
parameters <- rep(STi.non.fixed, each = length(Y_A))
# merge vector with data table
vec <- cbind(Y_A, Y_AB, Y_AX)
out <- data.table::data.table(vec, parameters)

# Define parameters
n <- length(Y_A[!is.na(Y_A)])
f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)

# Compute
first <- out[, .(parameters, Y_A, Y_AB)] %>%
  setnames(., c("Y_A", "Y_AB"), c("one", "two"))
second <- out[, .(parameters, Y_A, Y_AX)] %>%
  na.omit() %>%
  setnames(., c("Y_A", "Y_AX"), c("one", "two"))
third <- out[, .(parameters, Y_AB, Y_AX)] %>%
  na.omit() %>%
  setnames(., c("Y_AB", "Y_AX"), c("one", "two"))

output <- rbind(first, second, third) %>%
  .[order(parameters)] %>%
  .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters] %>%
  rbind(., STi.fixed) %>%
  .[order(parameters)]
















out[, sobol_Mapply_Ti(.SD), by = parameters]































sobol_Ti <- function(Y, params, type) {
  # Calculate the number of parameters
  k <- length(params)
  if(type == "normal") {
    # Extract the model output of the A matrix
    Y_A <- Y[1:p]
    # Extract the model output of the AB matrix
    Y_AB <- Y[(p+1):length(Y)]
    # Create vector with parameters
    parameters <- rep(params, each = length(Y_A))
    # merge vector with data table
    vec <- cbind(Y_A, Y_AB)
    out <- data.table::data.table(vec, parameters)
    # Compute Sobol'indices
    output <- out[, sobol_Mapply_Ti(.SD), by = parameters]
  }
  else if(type == "new") {
    # Extract the model output of the A matrix
    p <- Y:A[1:length(1:(length(Y) / ((2 * k) + 1)))]
    # Extract the model output of the A matrix
    p <- Y:A[1:length(1:(length(Y) / ((2 * k) + 1)))]
    # Extract the model output of the AB matrix
    Y_AB <- Y[(p+1): (p * (k + 1))]
    # Extract the model output of the AX matrix
    Y_AX <- Y[(p * (k + 1) + 1): length(Y)]
    # Create vector with parameters
    parameters <- rep(params, each = length(Y_A))
    # merge vector with data table
    vec <- cbind(Y_A, Y_AB, Y_AX)
    out <- data.table::data.table(vec, parameters)
    # Define parameters
    n <- length(Y_A[!is.na(Y_A)])
    f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
    VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
    # Compute
    first <- out[, .(parameters, Y_A, Y_AB)] %>%
      setnames(., c("Y_A", "Y_AB"), c("one", "two"))
    second <- out[, .(parameters, Y_A, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_A", "Y_AX"), c("one", "two"))
    third <- out[, .(parameters, Y_AB, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_AB", "Y_AX"), c("one", "two"))
    output <- rbind(first, second, third) %>%
      .[order(parameters)] %>%
      .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters]
  }
  return(output)
}





sobol_Ti <- function(Y, params, type) {
  # Calculate the number of parameters
  k <- length(params)
  # Extract the model output of the A matrix
  p <- Y:A[1:length(1:(length(Y) / ((2 * k) + 1)))]
  # Extract the model output of the AB matrix
  Y_AB <- Y[(p+1): (p * (k + 1))]
  # Extract the model output of the AX matrix
  Y_AX <- Y[(p * (k + 1) + 1): length(Y)]
  # Create vector with parameters
  parameters <- rep(params, each = length(Y_A))
  # merge vector with data table
  vec <- cbind(Y_A, Y_AB, Y_AX)
  out <- data.table::data.table(vec, parameters)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # Extract the model output of the A matrix
  Y_A <- Y[1:p]
  # Extract the model output of the AB matrix
  Y_AB <- Y[(p+1): (p * (k + 1))]
  # Extract the model output of the AX matrix
  Y_AX <- Y[(p * (k + 1) + 1): length(Y)]
  # Create vector with parameters
  parameters <- rep(params, each = length(Y_A))
  # merge vector with data table
  vec <- cbind(Y_A, Y_AB, Y_AX)
  out <- data.table::data.table(vec, parameters)
  # Compute Sobol'indices
  if(type == "normal") {
    output <- out[, sobol_Mapply_Ti(.SD), by = parameters]
  } else if(type == "new") {
    # Define parameters
    n <- length(Y_A[!is.na(Y_A)])
    f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
    VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
    
    # Compute
    first <- out[, .(parameters, Y_A, Y_AB)] %>%
      setnames(., c("Y_A", "Y_AB"), c("one", "two"))
    second <- out[, .(parameters, Y_A, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_A", "Y_AX"), c("one", "two"))
    third <- out[, .(parameters, Y_AB, Y_AX)] %>%
      na.omit() %>%
      setnames(., c("Y_AB", "Y_AX"), c("one", "two"))
    
    output <- rbind(first, second, third) %>%
      .[order(parameters)] %>%
      .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters]
  }
  return(output)
}


sobol_Ti(Y = Y, params = STi.non.fixed, type = "old")






































# COMPUTE STI -----------------------
# Calculate the length of the A matrix
p <- length(1:(length(Y) / (length(STi.non.fixed) + 1)))
# Extract the model output of the A matrix
Y_A <- Y[1:p]
# Extract the model output of the AB matrix
Y_AB <- Y[(p+1):length(Y)]
# Create vector with parameters
parameters <- rep(STi.non.fixed, each = length(Y_A))
# merge vector with data table
vec <- cbind(Y_A, Y_AB, Y_AX)
out <- data.table::data.table(vec, parameters)



n <- length(Y_A[!is.na(Y_A)])
f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
 

out[, sobol_Mapply_Ti(.SD), by = parameters]


first <- out[, .(parameters, Y_A, Y_AB)] %>%
  setnames(., c("Y_A", "Y_AB"), c("one", "two"))
second <- out[, .(parameters, Y_A, Y_AX)] %>%
  na.omit() %>%
  setnames(., c("Y_A", "Y_AX"), c("one", "two"))
third <- out[, .(parameters, Y_AB, Y_AX)] %>%
  na.omit() %>%
  setnames(., c("Y_AB", "Y_AX"), c("one", "two"))

rbind(first, second, third) %>%
  .[order(parameters)] %>%
  .[, ((1 / (2 * .N)) * sum((one - two) ^ 2, na.rm = TRUE)) / VY, parameters]



analytical




first <- out[, first:= (Y_A - Y_AB) ^2][, .(first, parameters)]
second <- out %>%
  na.omit() %>% 
  .[, first:= (Y_A - Y_AX) ^2]
second <- second[, .(first, parameters)]

third <- out %>%
  na.omit() %>% 
  .[, first:= (Y_AB - Y_AX) ^2]
third <- third[, .(first, parameters)]




rbind(first, second, third) %>%
  .[, ((1 / (2 * nrow(.SD))) * sum(first, na.rm = TRUE)) / VY, parameters]



















out %>%
  .[, effect1 := (Y_A - Y_AB) ^2 ] %>%
  .[, effect2:= (Y_A - Y_AX) ^2] %>%
  .[, effect3:= (Y_AB - Y_AX) ^2]
  
  
out[, ((1 / (2 * n)) * sum(effect1)) / VY, parameters]

out %>%
  na.omit() %>%
  .[, ((1 / (2 * .N)) * sum(effect2)) / VY, parameters]

out %>%
  na.omit() %>%
  .[, ((1 / (2 * .N)) * sum(effect3)) / VY, parameters]









out[, ((1 / (2 * n)) * sum((Y_A - Y_AB) ^ 2, na.rm = TRUE)), parameters]



out[, effect1:= Y_A - Y_AB]
out[, effect2:= Y_A - Y_AX]
out[, effect3:= Y_AB - Y_AX]

da <- out %>%
  na.omit() %>%
  .[, .(one = (1 / (2 * .N)) * sum(effect1) ^ 2, 
        two = (1 / (2 * .N)) * sum(effect2) ^ 2, 
        three = (1 / (2 * .N)) * sum(effect3) ^ 2), parameters] %>%
  .[, mean:= (one + two + three) / VY]
da
da

A[!parameters %in% "A"]








sobol_compute_Ti <- function(Y_A, Y_AB) {
  n <- length(Y_A[!is.na(Y_A)])
  f0 <- (1 / n) * sum(Y_A, na.rm = TRUE)
  VY <- 1 / n * sum((Y_A - f0) ^ 2, na.rm = TRUE)
  STi <- ((1 / (2 * n)) * sum((Y_A - Y_AB) ^ 2, na.rm = TRUE)) / VY
  return(STi) 
}

sobol_Mapply_Ti <- function(d) {
  return(mapply(sobol_compute_Ti,
                d[, "Y_A"],
                d[, "Y_AB"]))
}





















