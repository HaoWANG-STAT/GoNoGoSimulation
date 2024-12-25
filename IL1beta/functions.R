#' Summarise the Operating Characteristics into a 3-Outcome (Go/NoGo/Consider) Table
#' 
#' @param d (`tibble`)  The input data frame
#' @param NoGoThreshold (`integer`) The vector of No-Go threshold
#' @param GoThreshold (`integer`) The vector of Go threshold


createOCSummaryTable_3ODM <- function(
    d, 
    NoGoThreshold = 20,
    GoThreshold = 50,
    DeltaControl = 10, # effect of placebo used in simulations
    nSim = 2000 # number of simulations for each scenario
) {
  
  #header <-  paste0("N=", 2 * (d %>% distinct(NPerGroup) %>% pull(NPerGroup))) # table header
  col_labels  <-  c(
    "Effect of Selnoflast (%)", 
    rep(c("Go   ", "NoGo  ", "Consider"), d %>% distinct(NPerGroup) %>% nrow())
  )
  
  # Create the table
  d %>% 
    filter(Threshold %in% c(NoGoThreshold, GoThreshold) ) %>% # parameter Go and NoGo threshold
    mutate(NTotal = as.factor(2 * NPerGroup),
           Outcome = ifelse(Threshold==GoThreshold, "Go", "NoGo"),
           Prob = ifelse(Threshold==GoThreshold, GoProb, 1-GoProb)
    ) %>% 
    arrange(Threshold) %>% 
    select(-NPerGroup, -Threshold, -GoProb) %>% 
    mutate(Prob = Prob * 100) %>%
    pivot_wider(
      values_from = Prob,
      names_from = c(Outcome)
    ) %>% 
    mutate(Consider=100-Go-NoGo) %>% 
    pivot_wider(
      values_from = c(Go, NoGo, Consider),
      names_from = c(NTotal)
    ) %>% 
    select(Delta,
           Go_30, NoGo_30, Consider_30,
           Go_60, NoGo_60, Consider_60,
           Go_90, NoGo_90, Consider_90) %>% 
    kable(
      digits = 1,
      col.names = col_labels,
      full.width = FALSE
    )  %>% 
    add_header_above(
      c(
        " ", 
        "N=30" = 3,
        "N=60" = 3,
        "N=90" = 3
      )
    ) %>% 
    add_header_above(
      c("",
        setNames(9, paste0("Go if >=", GoThreshold, "%, " ,"No Go if <", NoGoThreshold, "%"))
    )
    )
}


#' Summarise the Operating Characteristics into a 3-Outcome (Go/NoGo/Consider) Figure
#' 
#' @param d (`tibble`)  The input data frame
#' @param NoGoThreshold (`integer`) The vector of No-Go threshold
#' @param GoThreshold (`integer`) The vector of Go threshold

createOCSummaryFigure_3ODM <- function(
    d, 
    NoGoThreshold = 20,
    GoThreshold = 50,
    DeltaControl = 10, # effect of placebo used in simulations
    nSim = 2000 # number of simulations for each scenario
) {
  
  # Create the figure 
  f <- d %>% 
    filter(Threshold %in% c(NoGoThreshold, GoThreshold) ) %>% # parameter Go and NoGo threshold
    mutate(NTotal = as.factor(2 * NPerGroup),
           Outcome = ifelse(Threshold==GoThreshold, "Go", "NoGo"),
           Prob = ifelse(Threshold==GoThreshold, GoProb, 1-GoProb)
    ) %>% 
    arrange(Threshold) %>% 
    select(-NPerGroup, -Threshold, -GoProb) %>% 
    mutate(Prob = Prob * 100) %>%
    pivot_wider(
      values_from = Prob,
      names_from = c(Outcome)
    ) %>% 
    mutate(Consider=100-Go-NoGo) %>%
    pivot_longer(cols = c(Go, NoGo, Consider),
                 values_to = "Prob",
                 names_to= "Decision"
    ) %>% 
    ggplot() +
    geom_line(aes(x = Delta, y = Prob, colour = Decision, linetype = NTotal)) +
    scale_color_manual(values = c(Go = "darkgreen", NoGo = "red", Consider = "darkblue")) +
    #facet_wrap(vars(NTotal), labeller = as_labeller(function(z) paste0("N=",z), label_both)) +
    labs(
      x = "True response rate for selnoflast",
      colour = "Decision",
      y = "Probability",
      title = paste0("Based on ", nSim, " simulations per scenario.  Pbo response rate: ", DeltaControl, "%")
    ) +
    theme_minimal()
  ggplotly(f, tooltip = c("y", "x")) %>% layout(hovermode = 'x')
  
}

#' Summarise the Operating Characteristics
#' 
#' @param d (`tibble`)  The input data frame
#' @param col_labels (`character`) The vector of column labels
#' @param group_label_css (`character`) The CSS used to format the group rows
createOCSummaryTable <- function(
    d, 
    col_labels = c(
      "Effect of Selnoflast (%)", 
      as.character(2 * (d %>% distinct(NPerGroup) %>% pull(NPerGroup)))
    ),
    group_label_css = "background-color: #666; color: #fff;"
) {
  # Derive the group labels and their positions
  groups <- d %>% 
    group_by(Threshold) %>%  
    summarise(N = n() / length(unique(NPerGroup))) 
  idx <- groups %>% pull(N)
  names(idx) <- paste0("Threshold: ", groups %>% pull(Threshold), "%")

  # Create the table
  d %>% 
    mutate(NTotal = as.factor(2 * NPerGroup)) %>% 
    select(-NPerGroup) %>% 
    arrange(Threshold) %>% 
    mutate(GoProb = GoProb * 100) %>% 
    pivot_wider(
      values_from = GoProb,
      names_from = NTotal
    ) %>% 
    select(-Threshold) %>% 
    kable(
      digits = 1,
      col.names = col_labels,
      full.width = FALSE
    ) %>% 
    pack_rows(
      index = idx, 
      label_row_css = group_label_css
    ) %>% 
    add_header_above(
      c(
        " " = 1, 
        "Trial size" = d %>% distinct(NPerGroup) %>% nrow()
      )
    )
}


#' Simulate the Operating Characteristics (Must be used along with function generateBiomarkerData2)
#' Modification of getOperatingCharacteristics1: not replicating data generation/simulation 
#' for each threshold as this parameter should be analysed based on same group of datasets 
#' generated under every combination of delta and sample size
#' 
#' #' @param nSim (`positive_integer`) the number of trials to simulate for each
#' combination of `nPerGroup`, `deltaTest` and `thesholdTest`
#' @param nPerGroup (`positive_integer`) a vector containing the sample sizes
#' @param deltaTest (`positive_integer`) a vector containing the deltas of selnoflast
#' @param deltaControl (`positive_integer`) a value defining the placebo effect
#' @param threshold (`positive_integer`) a vector containing the thresholds

#' @param generateData (`function`) the function used to generate participant-level data
#' @param fileName (`character`) the name of the data file from which the simulations,
#' if they exist, should be read
#' @param seed (`positive_integer`) the seed passed to `generateData` to initialise the
#' PRNG
#' @param dataFile (`character`) the name of the data file to read or create
#' @param ... passed to `generateData`
#' @return a tible containg the following columns:
#'   GoProb - The probability of getting a Go decision for this combination of
#'     parameters
#'   NPerGroup - The number of participants per group
#'   Delta - The true treatment effect, expressed as a percentage
#'   Threshold - The observed treatment effect, expressed as a percentage, required for
#'     a Go decision
#'     
#' @section Usage Notes:
#' The `generateData` function should take parameters `nPerGroup`, `seed` and `...` as
#' defined above, together with `delta` as described in `generateBiomarkerData2`.  
#' For a more detailed description of the data generation method, see
#' `generateBiomarkerData2`.  
#' If `dataFile` exists, its contents are read and returned.  if not, the 
#' simulations are performed and the result returned.
#' 
#' Parallelisation is used where possible.
getOperatingCharacteristics4 <- function(
    nSim = 1000,
    nPerGroup = c(30, 60, 90) / 2,     # varying parameters
    deltaTest = seq(0, 95, 5),         # varying parameters
    deltaControl = 0,
    threshold = seq(20, 50, 10),   # varying parameters
    generateData = generateBiomarkerData2,
    seed = floor(.Machine$integer.max * runif(1)),
    dataFile,
    ...
) {
  # Get a list containing all permutations of values in deltaTest and threshold
  # There must be a better way of doing this!
  perms <- tibble() %>% 
    expand(deltaT = deltaTest, nPerGroup = nPerGroup) %>% 
    rowwise() %>% 
    group_map(
      function(.x, .y) list(deltaT = .x$deltaT, nPerGroup = .x$nPerGroup)
    )
  
  set.seed(seed)
  
  if (file.exists(dataFile)) {
    d <- readRDS(dataFile)
    message(paste0("Data read from ", dataFile, "."))
  } else {
    startTime = Sys.time()
    coresAvailable <- detectCores()
    coresNeeded <- min(coresAvailable - 1, length(perms))
    
    d <- mclapply(
      perms,
      function(perm) {
        deltaT <- perm$deltaT
        nPerGroup <- perm$nPerGroup
        #lapply(
        #nPerGroup,
        #function(n) {
        lapply(
          seq_len(nSim),
          function(dummy) {
            x <- generateData(
              nPerGroup = nPerGroup,
              delta = list(
                "Test" = deltaT,           # varying parameters
                "Control" = deltaControl
              ),
              seed = floor(.Machine$integer.max * runif(1)),
              ...
            ) 
            x %>% 
              pivot_wider(
                values_from = Value,
                names_from = TimePoint
              ) 
          }
        ) %>% 
          bind_rows(.id = "Trial") %>%
          group_by(Trial) %>%
          group_map(
            function(.x, .y) {
              #.x %>%
              lapply(
                threshold, 
                function(threshold_each) analyseTrial(d = .x, threshold = threshold_each) %>% 
                  add_column(Threshold=threshold_each) 
              ) %>% 
                bind_rows() %>%  # convert listing to a tibble/dataframe, this step should not be removed
                bind_cols(.y)  # add the group variable .y (ie, Trial) back to the tibble
            }
          ) %>%
          bind_rows() %>%  # convert listing to a tibble/dataframe, this step should not be removed
          group_by(Threshold) %>% 
          summarise(GoProb = mean(Outcome == "Go")) %>% 
          add_column(Delta = deltaT) %>% 
          add_column(NPerGroup = nPerGroup)
      },
      mc.cores = coresNeeded
    ) 
    
    d <- d %>% bind_rows() # convert listing to a tibble/dataframe, this step should not be removed
    
    saveRDS(d, dataFile)
    endTime = Sys.time()
    t
    cat(
      paste0(
        "Data generated in ", 
        str_extract(
          capture.output(print(endTime - startTime, digits = 5)), 
          "\\d+\\.\\d+ .+"
        ), 
        " and results saved to ", 
        dataFile,
        "."
      )
    )
  }
  d
}

#' Simulate the Operating Characteristics (Must be used along with function generateBiomarkerDat or generateBiomarkerData1)
#' Modification of getOperatingCharacteristics1: not replicating data generation/simulation 
#' for each threshold as this parameter should be analysed based on same group of datasets 
#' generated under every combination of delta and sample size
#' 
#' @param nSim (`positive_integer`) the number of trials to simulate for each
#' combination of `nPerGroup` and `deltaTest`
#' @param nPerGroup (`positive_integer`) a vector containing the sample sizes to be
#' simulated
#' @param deltaTest (`positive_integer`) a vector containing the deltas to be
#' simulated
#' @param thresholdTest (`positive_integer`) a vector containing the thresholds to be
#' analysed
#' @param generateData (`function`) the function used to generate participant-level data
#' @param fileName (`character`) the name of the data file from which the simulations,
#' if they exist, should be read
#' @param seed (`positive_integer`) the seed passed to `generateData` to initialise the
#' PRNG
#' @param dataFile (`character`) the name of the data file to read or create
#' @param ... passed to `generateData`
#' @return a tible containg the following columns:
#'   GoProb - The probability of getting a Go decision for this combination of
#'     parameters
#'   NPerGroup - The number of participants per group
#'   Delta - The true treatment effect, expressed as a percentage
#'   Threshold - The observed treatment effect, expressed as a percentage, required for
#'     a Go decision
#'     
#' @section Usage Notes:
#' The `generateData` function should take parameters `nPerGroup`, `seed` and `...` as
#' defined above, together with `delta` as described in `generateBiomarkerData1`.  
#' For a more detailed description of the data generation method, see
#' `generateBiomarkerData1`.  
#' If `dataFile` exists, its contents are read and returned.  if not, the 
#' simulations are performed and the result returned.
#' 
#' Parallelisation is used where possible.
getOperatingCharacteristics3 <- function(
    nSim = 1000,
    nPerGroup = c(30, 60, 90) / 2,
    deltaTest = seq(0, 95, 5),
    threshold = seq(20, 50, 10),
    generateData = generateBiomarkerData1,
    seed = floor(.Machine$integer.max * runif(1)),
    dataFile,
    ...
) {
  # Get a list containing all permutations of values in deltaTest and thresholdTest
  # There must be a better way of doing this!
  perms <- tibble() %>% 
    expand(delta = deltaTest, nPerGroup = nPerGroup) %>% 
    rowwise() %>% 
    group_map(
      function(.x, .y) list(delta = .x$delta, nPerGroup = .x$nPerGroup)
    )
  
  set.seed(seed)
  
  if (file.exists(dataFile)) {
    d <- readRDS(dataFile)
    message(paste0("Data read from ", dataFile, "."))
  } else {
    startTime = Sys.time()
    coresAvailable <- detectCores()
    coresNeeded <- min(coresAvailable - 1, length(perms))
    
    d <- mclapply(
      perms,
      function(perm) {
        delta <- perm$delta
        nPerGroup <- perm$nPerGroup
        #lapply(
        #nPerGroup,
        #function(n) {
        lapply(
          seq_len(nSim),
          function(dummy) {
            x <- generateData(
              nPerGroup = nPerGroup,
              deltaTest = delta,
              seed = floor(.Machine$integer.max * runif(1)),
              ...
            ) 
            x %>% 
              pivot_wider(
                values_from = Value,
                names_from = TimePoint
              ) 
          }
        ) %>% 
          bind_rows(.id = "Trial") %>%
          group_by(Trial) %>%
          group_map(
            function(.x, .y) {
              #.x %>%
              lapply(
                threshold, 
                function(threshold_each) analyseTrial(d = .x, threshold = threshold_each) %>% 
                  add_column(Threshold=threshold_each) 
              ) %>% 
                bind_rows() %>%  # convert listing to a tibble/dataframe, this step should not be removed
                bind_cols(.y)  # add the group variable .y (ie, Trial) back to the tibble
            }
          ) %>%
          bind_rows() %>%  # convert listing to a tibble/dataframe, this step should not be removed
          group_by(Threshold) %>% 
          summarise(GoProb = mean(Outcome == "Go")) %>% 
          add_column(Delta = delta) %>% 
          add_column(NPerGroup = nPerGroup)
      },
      mc.cores = coresNeeded
    ) 
    
    d <- d %>% bind_rows() # convert listing to a tibble/dataframe, this step should not be removed
    
    saveRDS(d, dataFile)
    endTime = Sys.time()
    t
    cat(
      paste0(
        "Data generated in ", 
        str_extract(
          capture.output(print(endTime - startTime, digits = 5)), 
          "\\d+\\.\\d+ .+"
        ), 
        " and results saved to ", 
        dataFile,
        "."
      )
    )
  }
  d
}

#getOperatingCharacteristics3(  nSim=1,
#                               generateData = generateBiomarkerData1,
#                               dataFile = here("IL1beta", "test.Rds"),)



#' Simulate the Operating Characteristics(Must be used along with function generateBiomarkerData2)
#' 
#' @param nSim (`positive_integer`) the number of trials to simulate for each
#' combination of `nPerGroup`, `deltaTest` and `thesholdTest`
#' @param nPerGroup (`positive_integer`) a vector containing the sample sizes
#' @param deltaTest (`positive_integer`) a vector containing the deltas of selnoflast
#' @param deltaControl (`positive_integer`) a value defining the placebo effect
#' @param threshold (`positive_integer`) a vector containing the thresholds

#' @param generateData (`function`) the function used to generate participant-level data
#' @param fileName (`character`) the name of the data file from which the simulations,
#' if they exist, should be read
#' @param seed (`positive_integer`) the seed passed to `generateData` to initialise the
#' PRNG
#' @param dataFile (`character`) the name of the data file to read or create
#' @param ... passed to `generateData`
#' @return a tible containg the following columns:
#'   GoProb - The probability of getting a Go decision for this combination of
#'     parameters
#'   NPerGroup - The number of participants per group
#'   Delta - The true treatment effect, expressed as a percentage
#'   Threshold - The observed treatment effect, expressed as a percentage, required for
#'     a Go decision
#'     
#' @section Usage Notes:
#' The `generateData` function should take parameters `nPerGroup`, `seed` and `...` as
#' defined above, together with `delta` as described in `generateBiomarkerData2`.  
#' For a more detailed description of the data generation method, see
#' `generateBiomarkerData2`.  
#' If `dataFile` exists, its contents are read and returned.  if not, the 
#' simulations are performed and the result returned.
#' 
#' Parallelisation is used where possible.
getOperatingCharacteristics2 <- function(
    nSim = 1000,
    nPerGroup = c(30, 60, 90) / 2,     # varying parameters
    deltaTest = seq(0, 95, 5),         # varying parameters
    deltaControl = 0,
    threshold = seq(20, 50, 10),   # varying parameters
    generateData = generateBiomarkerData2,
    seed = floor(.Machine$integer.max * runif(1)),
    dataFile,
    ...
) {
  # Get a list containing all permutations of values in deltaTest and threshold
  # There must be a better way of doing this!
  perms <- tibble() %>% 
    expand(deltaTest = deltaTest, threshold = threshold) %>% 
    rowwise() %>% 
    group_map(
      function(.x, .y) list(deltaTest = .x$deltaTest, threshold = .x$threshold)
    )
  
  set.seed(seed)
  
  if (file.exists(dataFile)) {
    d <- readRDS(dataFile)
    message(paste0("Data read from ", dataFile, "."))
  } else {
    startTime = Sys.time()
    coresAvailable <- detectCores()
    coresNeeded <- min(coresAvailable - 1, length(perms))
    
    d <- mclapply(
      perms,
      function(perm) {
        deltaTest <- perm$deltaTest
        threshold <- perm$threshold
        lapply(
          nPerGroup,
          function(n) {
            lapply(
              seq_len(nSim),
              function(dummy) {
                x <- generateData(
                  nPerGroup = n,
                  delta = list(
                    "Test" = deltaTest,           # varying parameters
                    "Control" = deltaControl
                  ),
                  seed = floor(.Machine$integer.max * runif(1)),
                  ...
                ) 
                x %>% 
                  pivot_wider(
                    values_from = Value,
                    names_from = TimePoint
                  ) 
              }
            ) %>% 
              bind_rows(.id = "Trial") %>%
              group_by(Trial) %>%
              group_map(
                function(.x, .y) {
                  .x %>%
                    analyseTrial(threshold = threshold) %>%
                    bind_cols(.y)
                }
              ) %>%
              bind_rows(.id = "Trial") %>%
              summarise(GoProb = mean(Outcome == "Go")) %>%
              add_column(NPerGroup = n)
          }
        ) %>% 
          bind_rows() %>% 
          add_column(Delta = deltaTest) %>% 
          add_column(Threshold = threshold)
      },
      mc.cores = coresNeeded
    ) 
    
    d <- d %>% bind_rows()
    
    saveRDS(d, dataFile)
    endTime = Sys.time()
    t
    cat(
      paste0(
        "Data generated in ", 
        str_extract(
          capture.output(print(endTime - startTime, digits = 5)), 
          "\\d+\\.\\d+ .+"
        ), 
        " and results saved to ", 
        dataFile,
        "."
      )
    )
  }
  d
}


#' Simulate the Operating Characteristics of a Design (Must be used along with function generateBiomarkerDat or generateBiomarkerData1)
#' Modification of getOperatingCharacteristics to permit change the response rate of control arm
#' 
#' @param nSim (`positive_integer`) the number of trials to simulate for each
#' combination of `nPerGroup`, `deltaTest` and `theshold`
#' @param nPerGroup (`positive_integer`) a vector containing the sample sizes to be
#' simulated
#' @param deltaTest (`positive_integer`) a vector containing the deltas to be
#' simulated
#' @param threshold (`positive_integer`) a vector containing the thresholds to be
#' simulated
#' @param generateData (`function`) the function used to generate participant-level data
#' @param fileName (`character`) the name of the data file from which the simulations,
#' if they exist, should be read
#' @param seed (`positive_integer`) the seed passed to `generateData` to initialise the
#' PRNG
#' @param dataFile (`character`) the name of the data file to read or create
#' @param ... passed to `generateData`
#' @return a tible containg the following columns:
#'   GoProb - The probability of getting a Go decision for this combination of
#'     parameters
#'   NPerGroup - The number of participants per group
#'   Delta - The true treatment effect, expressed as a percentage
#'   Threshold - The observed treatment effect, expressed as a percentage, required for
#'     a Go decision
#'     
#' @section Usage Notes:
#' The `generateData` function should take parameters `nPerGroup`, `seed` and `...` as
#' defined above, together with `delta` as described in `generateBiomarkerData1`.  
#' For a more detailed description of the data generation method, see
#' `generateBiomarkerData`.  
#' If `dataFile` exists, its contents are read and returned.  if not, the 
#' simulations are performed and the result returned.
#' 
#' Parallelisation is used where possible.
getOperatingCharacteristics1 <- function(
    nSim = 1000,
    nPerGroup = c(30, 60, 90) / 2,
    deltaTest = seq(0, 95, 5),
    threshold = seq(20, 50, 10),
    generateData = generateBiomarkerData1,
    seed = floor(.Machine$integer.max * runif(1)),
    dataFile,
    ...
) {
  # Get a list containing all permutations of values in deltaTest and thresholdTest
  # There must be a better way of doing this!
  perms <- tibble() %>% 
    expand(delta = deltaTest, threshold = threshold) %>% 
    rowwise() %>% 
    group_map(
      function(.x, .y) list(delta = .x$delta, threshold = .x$threshold)
    )
  
  set.seed(seed)
  
  if (file.exists(dataFile)) {
    d <- readRDS(dataFile)
    message(paste0("Data read from ", dataFile, "."))
  } else {
    startTime = Sys.time()
    coresAvailable <- detectCores()
    coresNeeded <- min(coresAvailable - 1, length(perms))
    
    d <- mclapply(
      perms,
      function(perm) {
        delta <- perm$delta
        threshold <- perm$threshold
        lapply(
          nPerGroup,
          function(n) {
            lapply(
              seq_len(nSim),
              function(dummy) {
                x <- generateData(
                  nPerGroup = n,
                  deltaTest = delta,
                  seed = floor(.Machine$integer.max * runif(1)),
                  ...
                ) 
                x %>% 
                  pivot_wider(
                    values_from = Value,
                    names_from = TimePoint
                  ) 
              }
            ) %>% 
              bind_rows(.id = "Trial") %>%
              group_by(Trial) %>%
              group_map(
                function(.x, .y) {
                  .x %>%
                    analyseTrial(threshold = threshold) %>%
                    bind_cols(.y)
                }
              ) %>%
              bind_rows(.id = "Trial") %>%
              summarise(GoProb = mean(Outcome == "Go")) %>%
              add_column(NPerGroup = n)
          }
        ) %>% 
          bind_rows() %>% 
          add_column(Delta = delta) %>% 
          add_column(Threshold = threshold)
      },
      mc.cores = coresNeeded
    ) 
    
    #%>% bind_rows()
    d <- d %>% bind_rows()
    
    saveRDS(d, dataFile)
    endTime = Sys.time()
    t
    cat(
      paste0(
        "Data generated in ", 
        str_extract(
          capture.output(print(endTime - startTime, digits = 5)), 
          "\\d+\\.\\d+ .+"
        ), 
        " and results saved to ", 
        dataFile,
        "."
      )
    )
  }
  d
}


#' Simulate the Operating Characteristics of a Design
#' 
#' @param nSim (`positive_integer`) the number of trials to simulate for each
#' combination of `nPerGroup`, `deltaTest` and `thesholdTest`
#' @param nPerGroup (`positive_integer`) a vector containing the sample sizes to be
#' simulated
#' @param deltaTest (`positive_integer`) a vector containing the deltas to be
#' simulated
#' @param threshold (`positive_integer`) a vector containing the thresholds to be
#' simulated
#' @param generateData (`function`) the function used to generate participant-level data
#' @param fileName (`character`) the name of the data file from which the simulations,
#' if they exist, should be read
#' @param seed (`positive_integer`) the seed passed to `generateData` to initialise the
#' PRNG
#' @param dataFile (`character`) the name of the data file to read or create
#' @param ... passed to `generateData`
#' @return a tible containg the following columns:
#'   GoProb - The probability of getting a Go decision for this combination of
#'     parameters
#'   NPerGroup - The number of participants per group
#'   Delta - The true treatment effect, expressed as a percentage
#'   Threshold - The observed treatment effect, expressed as a percentage, required for
#'     a Go decision
#'     
#' @section Usage Notes:
#' The `generateData` function should take parameters `nPerGroup`, `seed` and `...` as
#' defined above, together with `delta` as described in `generateBiomarkerData`.  
#' For a more detailed description of the data generation method, see
#' `generateBiomarkerData`.  
#' If `dataFile` exists, its contents are read and returned.  if not, the 
#' simulations are performed and the result returned.
#' 
#' Parallelisation is used where possible.
getOperatingCharacteristics <- function(
    nSim = 1000,
    nPerGroup = c(30, 60, 90) / 2,
    deltaTest = seq(0, 95, 5),
    threshold = seq(20, 50, 10),
    generateData = generateBiomarkerData,
    seed = floor(.Machine$integer.max * runif(1)),
    dataFile,
    ...
) {
  # Get a list containing all permutations of values in deltaTest and threshold
  # There must be a better way of doing this!
  perms <- tibble() %>% 
    expand(deltaTest = deltaTest, threshold = threshold) %>% 
    rowwise() %>% 
    group_map(
      function(.x, .y) list(deltaTest = .x$deltaTest, threshold = .x$threshold)
    )
  
  set.seed(seed)
  
  if (file.exists(dataFile)) {
    d <- readRDS(dataFile)
    message(paste0("Data read from ", dataFile, "."))
  } else {
    startTime = Sys.time()
    coresAvailable <- detectCores()
    coresNeeded <- min(coresAvailable - 1, length(perms))
    
    d <- mclapply(
      perms,
      function(perm) {
        deltaTest <- perm$deltaTest
        threshold <- perm$threshold
        lapply(
          nPerGroup,
          function(n) {
            lapply(
              seq_len(nSim),
              function(dummy) {
                x <- generateData(
                  nPerGroup = n,
                  deltaTest = deltaTest,
                  seed = floor(.Machine$integer.max * runif(1)),
                  ...
                ) 
                x %>% 
                  pivot_wider(
                    values_from = Value,
                    names_from = TimePoint
                  ) 
              }
            ) %>% 
              bind_rows(.id = "Trial") %>%
              group_by(Trial) %>%
              group_map(
                function(.x, .y) {
                  .x %>%
                    analyseTrial(threshold = threshold) %>%
                    bind_cols(.y)
                }
              ) %>%
              bind_rows(.id = "Trial") %>%
              summarise(GoProb = mean(Outcome == "Go")) %>%
              add_column(NPerGroup = n)
          }
        ) %>% 
          bind_rows() %>% 
          add_column(Delta = deltaTest) %>% 
          add_column(Threshold = threshold)
      },
      mc.cores = coresNeeded
    ) 
    
    d <- d %>% bind_rows()
    
    saveRDS(d, dataFile)
    endTime = Sys.time()
    t
    cat(
      paste0(
        "Data generated in ", 
        str_extract(
          capture.output(print(endTime - startTime, digits = 5)), 
          "\\d+\\.\\d+ .+"
        ), 
        " and results saved to ", 
        dataFile,
        "."
      )
    )
  }
  d
}

analyseTrial <- function(d, threshold) {
d %>% 
    mutate(logBL = log(Baseline),
           logEP = log(Endpoint)) %>% 
  group_by(Treatment) %>% 
  summarise(
    across(
      c(Baseline, Endpoint, Delta, logBL, logEP),
      list(Mean = \(x) mean(x, na.rm = TRUE),
           SD = \(x) sd(x, na.rm = TRUE)  # output SD of raw and log value for further evaluation of standard deviation of endpoint and change from baseline
      )
    ),
    .groups = "keep"
  ) %>% 
  mutate(PctMeanDelta = Delta_Mean / Baseline_Mean * 100,
         MeanChangeLog = logBL_Mean - logEP_Mean) %>% 
  select(Treatment, PctMeanDelta, Baseline_SD, Endpoint_SD, Delta_SD, logBL_SD, logEP_SD, MeanChangeLog) %>% 
  pivot_wider(
    names_from = Treatment,
    values_from = c(PctMeanDelta, Baseline_SD, Endpoint_SD, Delta_SD, logBL_SD, logEP_SD, MeanChangeLog)
  ) %>% 
  mutate(Outcome = ifelse(PctMeanDelta_Test < (PctMeanDelta_Control - threshold), "Go", "NoGo"))
}

#' Generate Biomarker Data with Intra-Patient Correlation
#' 
#' This function generates IL1-beta (or any other biomarker) data with intra-patient correlation.
#' 
#' @param nPerGroup (`integer`) the number of participants per group
#' @param mu (`numeric`) the mean of the biomarker before transformation
#' @param sigma2 (`numeric`) the variance of the biomarker
#' @param rho (`numeric`) the correlation between visits
#' @param deltaTest (`numeric`) the selnoflast effect.
#' @param deltaControl (`numeric`) the placebo effect.
#' @param effect_limit (`numeric`)\cr The endpoint value that defines the maximum
#' possible (100%) treatment effect.  See Usage Notes below.
#' @param seed (`integer`)  the seed used to initialise the PRNG
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @section Usage Notes:
#' The default values of `delta` and `effect_limit` lower the endpoint values in the 
#' Test group by a proportion that has a Beta distribution with parameters 
#' alpha = 7 and beta = 3 (for a mean of 70%) and those in the Control group by 
#' a proportion with a Beta distribution with parameters alpha = 1 and beta = 9
#' {for a mean of 10%)}.
#' 
#' The formula for calculating the treatment-adjusted endpoint value is
#' x_adj = x_unadj - p * (x_unadj - limit)
#' where x_unadj is the unadjusted endpoint value, limit is the value of
#' `effect_limit` above and p is a random number between 0 and 1 generated 
#' according to the rules in `delta`.  For example, if x_unadj is 100, limit is 50
#' and p is 0.2, then x_adj is 100 - 0.2 * (100 - 50) = 90.
#' 
#' The effect limit should be specified on the raw scale, even though endpoints 
#' are generated on the log scale and then back transformed to the raw scale.
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @section Usage Notes:
#' 
generateBiomarkerData <- function(
    nPerGroup = 10,
    mu = c("il1beta" = 5.6),
    sigma2 = 0.54 * 0.54,
    rho = 0.5,
    deltaTest = 70,
    deltaControl = 0,
    effect_limit = 0,
    seed = 172910
) {
  # Prepare
  set.seed(seed)
  # Execute
  v <- matrix(c(sigma2, rho * sigma2, rho * sigma2, sigma2), ncol = 2)
  m <- rep(mu, 2)
  names(m) <- c("Baseline", "Endpoint")
  nTotal <- 2 * nPerGroup
  delta = list(
    "Test" = function() {
      if (deltaTest==0) {1} else {1 - rbeta(1, deltaTest/10, (100 - deltaTest)/10) } # for zero placebo effect
    },
    "Control" = function() {
      if (deltaControl==0) {1} else { 1 - rbeta(1, deltaControl/10, (100 - deltaControl)/10) }
    } # for zero placebo effect
  )
  d <- as_tibble(rmvnorm(nTotal, m, v))%>% 
    add_column(PID = 1:nTotal, .before = 1) %>% 
    mutate(Treatment = rep(c("Test", "Control"), each = nPerGroup)) %>% 
    rowwise() %>% 
    mutate(
      Baseline = exp(Baseline),
      UnadjustedEndpoint = exp(Endpoint),
      Factor = delta[[Treatment]](),
      Endpoint = ifelse(
        UnadjustedEndpoint < effect_limit,
        UnadjustedEndpoint,
        effect_limit + Factor * (exp(Endpoint) - effect_limit)
      ),
      Delta = Endpoint - Baseline
    ) %>%
    pivot_longer(
      -c(PID, Treatment, UnadjustedEndpoint, Factor),
      names_to = "TimePoint",
      values_to = "Value"
    )  %>% 
    select(-UnadjustedEndpoint, -Factor)
  d
}
#' Generate Biomarker Data with Intra-Patient Correlation
#' apply a different definition of "partial" benefit from generateBiomarkerData
#' 
#' This function generates IL1-beta (or any other biomarker) data with intra-patient correlation.
#' 
#' @param nPerGroup (`integer`) the number of participants per group
#' @param mu (`numeric`) the mean of the biomarker before transformation
#' @param sigma2 (`numeric`) the variance of the biomarker
#' @param rho (`numeric`) the correlation between visits
#' @param deltaTest (`numeric`) the selnoflast effect.
#' @param deltaControl (`numeric`) the placebo effect.
#' @param effect_limit (`numeric`)\cr The threshold below which the active treatment
#' has no effect.  See Usage Notes below.
#' @param seed (`integer`)  the seed used to initialise the PRNG
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @param seed (`integer`)  the seed used to initialise the PRNG
#' @section Usage Notes:
#' The default values of `delta` and `effect_limit` lower the endpoint values in the 
#' Test group by a proportion that has a Beta distribution with parameters 
#' alpha = 7 and beta = 3 (for a mean of 70%) and those in the Control group by 
#' a proportion with a Beta distribution with parameters alpha = 1 and beta = 9
#' {for a mean of 10%)}.
#' 
#' The effect of active treatment depends on the size of the unadjusted endpoint 
#' value relative to the threshold.  If the unadjusted endpoint is above the 
#' threshold, the active treatment effect applies.  Otherwise, the control 
#' treatment effect applies.
#' 
#' The effect limit should be specified on the raw scale, even though endpoints 
#' are generated on the log scale and then back transformed to the raw scale.
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @section Usage Notes:
#' 
generateBiomarkerData1 <- function(
    nPerGroup = 10,
    mu = c("il1beta" = 5.6),
    sigma2 = 0.54 * 0.54,
    rho = 0.5,
    deltaTest = 70,
    deltaControl = 0, # to specify placebo effect
    effect_limit = 0,
    seed = 172910
) {
  # Prepare
  set.seed(seed)
  # Execute
  v <- matrix(c(sigma2, rho * sigma2, rho * sigma2, sigma2), ncol = 2)
  m <- rep(mu, 2)
  names(m) <- c("Baseline", "Endpoint")
  nTotal <- 2 * nPerGroup
  delta = list(
    "Test" = function() {
      if (deltaTest==0) {1} else {1 - rbeta(1, deltaTest/10, (100 - deltaTest)/10) } # for zero placebo effect
    },
    "Control" = function() {
      if (deltaControl==0) {1} else { 1 - rbeta(1, deltaControl/10, (100 - deltaControl)/10) }
    } # for zero placebo effect
  )
  d <- as_tibble(rmvnorm(nTotal, m, v))%>% 
    add_column(PID = 1:nTotal, .before = 1) %>% 
    mutate(Treatment = rep(c("Test", "Control"), each = nPerGroup)) %>% 
    rowwise() %>% 
    mutate(
      Baseline = exp(Baseline),
      UnadjustedEndpoint = exp(Endpoint),
      Factor = ifelse(
        UnadjustedEndpoint > effect_limit,
        delta[[Treatment]](),
        delta[["Control"]]()
      ),
      Endpoint = Factor * exp(Endpoint),
      Delta = Endpoint - Baseline
    ) %>%
    pivot_longer(
      -c(PID, Treatment, UnadjustedEndpoint, Factor),
      names_to = "TimePoint",
      values_to = "Value"
    )  %>% 
    select(-UnadjustedEndpoint, -Factor)
  d
}


#' Generate Biomarker Data with Intra-Patient Correlation
#' use a different data generation strategy from generateBiomarkerData and generateBiomarkerData1
#' apply similar definition of "partial" benefit as generateBiomarkerData1
#' 
#' @param nPerGroup (`integer`) the number of participants per group
#' @param mu (`numeric`) the mean of the biomarker in log-scale at baseline
#' the mean of endpoint in each arm can be derived mu2= mu1+1/2*sigma1^2 -1/2*sigma2^2 + log(1-delta)
#' refer to https://www.statlect.com/probability-distributions/log-normal-distribution
#' @param sigma {`list`} a list of values that define standard deviation in log-scale
#' The list should contain three elements, "Baseline", "Test", and "Control". 
#' @param rho (`numeric`) the correlation between visits
#' @param delta {`list`} a list of values that define the % of mean change from baseline
#' in raw scale [(mean(Baseline)-mean(Endpoint)/mean(Baseline)*100%].  
#' The list should contain two elements, one labelled "Test", the other "Control". 
#' The default value lowers the endpoint values in the 
#' Test group by 70% and those in the Control group by 0%.
#' @param effect_limit (`numeric`)\cr The threshold below which the active treatment
#' has no effect.  See Usage Notes below.
#' @param seed (`integer`)  the seed used to initialise the PRNG
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @section Usage Notes:
#' The effect of active treatment depends on the size of the baseline 
#' value relative to the threshold.  If the baseline value is above the 
#' threshold, the active treatment effect applies.  Otherwise, the "zero" 
#' treatment effect applies.
#' #' The effect limit should be specified on the raw scale, even though endpoints 
#' are generated on the log scale and then back transformed to the raw scale.
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#'
generateBiomarkerData2 <- function(
    nPerGroup = 10,
    rho = 0.5, # baseline-endpoint correlation coefficient
    delta = list(
      "Test" = 70,
      "Control" = 0
    ),
    sigma = list(
      "Baseline" = 1, # standard deviation of log-baseline for test/control arm
      "Control" = 1, # standard deviation of log-endpoint for control arm
      "Test" = 1 # standard deviation of log-endpoint for test arm
    ),
    mu = c("Baseline" = 4),
    effect_limit = 0,
    seed = 172910
) {
  # Prepare
  set.seed(seed)
  # Execute
  
  # set variance-covariance matrix
  v <- function(Treatment) {
    matrix(c(sigma[["Baseline"]]^2, rho*sigma[["Baseline"]]*sigma[[Treatment]], 
             rho*sigma[["Baseline"]]*sigma[[Treatment]], sigma[[Treatment]]^2), 
           ncol = 2)
  }
  
  # set mean vector: mu2= mu1 + 1/2*sigma1^2 -1/2*sigma2^2 + log(1-delta)
  m <- function(Treatment) {
    c(mu, mu + 1/2*sigma[["Baseline"]]^2 - 1/2*sigma[[Treatment]]^2 + log(1-delta[[Treatment]]/100) ) 
  } 
  
  
  nTotal <- 2 * nPerGroup
  d <- tibble(PID = 1:nTotal,
              Treatment = rep(c("Test", "Control"), each = nPerGroup)) %>% 
    rowwise() %>%  
    mutate(logValue=rmvnorm(1, m(Treatment), v(Treatment)),
           Baseline = exp(logValue[,1]),
           Endpoint = ifelse(Treatment =="Test" & Baseline<=effect_limit, # situation of partial benefit
                             exp(logValue[,2]- (1/2*sigma[["Baseline"]]^2 - 1/2*sigma[[Treatment]]^2 + log(1-delta[[Treatment]]/100))),
                             exp(logValue[,2])
           ),
           Delta = Endpoint-Baseline # patient level change from baseline
    ) %>% 
    select(-c(logValue)) %>% 
    pivot_longer(
      -c(PID, Treatment),
      names_to = "TimePoint",
      values_to = "Value"
    )
  d
}


# Do not use this function without revision to match generateDataIL2beta above

#' Generate IL1beta, Neutrophil and CRP Data based on UBIOPRED data
#' 
#' This function generates visit level biomarker data based on the UBIOPRED dataset
#' in UBIOPREDjohn.csv supplied by Richard Welford.  It honours the observed 
#' correlations between biomarkers within visit, but observations across visits
#' are uncorrelated.
#' @param nPerGroup (`integer`) the number of participants per group
#' @param mu (`numeric`) a named vector of means, one for each biomarker
#' @param sigma2 (`matrix`) the covariance matrix for the biomarkers
#' @param delta {`list`} a list of functions that define how to modify the endpoint
#' values of each biomarker.  The list should contain two elements, one labelled
#' "Test", the other "Control". The default value lowers the endpoint values in the 
#' Test group by 70% and those in the Control group by 10%.
#' @param seed (`integer`)  the seed used to initialise the PRNG
#' @return a tibble in long format, with columns `PID`, `TimePoint`, `Variable`
#' and `Value`. The values of `TimePoint` are `"Baseline"`, `"Endpoint"` and 
#' `"Delta"`, the latter containing the Endpoint - Baseline differences.  The values
#' of `Variable` are names of `mu`.
#' @section Usage Notes:
#' This function has not been tested with anything other than three biomarkers.
generateDataBiomarkers <- function(
    nPerGroup = 10,
    mu = c(il1beta = 256.89, neutrophils = 48.64, hsCRP = 4.66),
    sigma2 = matrix(
      c(
        1.0000 * 107.74 * 107.74, 0.8682 * 107.74 * 25.75, 0.3789 * 107.74 * 11.81,
        0.8682 * 107.74 * 25.75,  1.0000 * 25.75 * 25.75, -0.1570 * 25.75 * 11.81,
        0.3789 * 107.74 * 11.81, -0.1570 * 25.75 * 11.81,  1.0000 * 11.81 * 11.81
      ),
      nrow = 3
    ),
    delta = list(
      "Test" = function(x) x * rnorm(1, mean = 0.3, sd = 0.1),
      "Control" = function(x) x * rnorm(1, mean = 0.9, sd = 0.1)
    ),
    seed = 561158
) {
  # Prepare
  set.seed(seed)
  nTotal <- 2 * nPerGroup
  muBase <- mu
  names(muBase) <- paste0("b_", names(mu))
  muEnd <- mu
  names(muEnd) <- paste0("e_", names(mu))
  # Execute
  bind_cols(
    rmvnorm(nTotal, muBase, sigma2), 
    rmvnorm(nTotal, muEnd, sigma2)
  ) %>% 
    add_column(PID = 1:nTotal, .before = 1) %>% 
    mutate(Treatment = rep(c("Test", "Control"), each = nPerGroup)) %>% 
    rowwise() %>% 
    # This step is incorrect: It calculates the revised endpoint value, 
    # not the change from baseline
    mutate(
      across(
        starts_with("e_"), 
        \(.) delta[[Treatment]](.), 
        .names = "d{.col}"
      )
    ) %>% 
    pivot_longer(
      -c(PID, Treatment),
      names_to = c("TimePoint", "Variable"),
      values_to = "Value",
      names_sep = "_"
    ) %>% 
    mutate(
      TimePoint = case_when(
        TimePoint == "b" ~ "Baseline",
        TimePoint == "e" ~ "Endpoint",
        TimePoint == "de" ~ "Delta",
      )
    )
}
