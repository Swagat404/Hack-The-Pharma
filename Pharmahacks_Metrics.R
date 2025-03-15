# Could implement both relative MSE and MAE for K_d and S_inhibitor

# MSE: Relative MSE < 1 => Good
get_rel_mse <- function(prediction_vector, true_vector) {
  num_predictions <- length(prediction_vector)
  mse <- sum((prediction_vector - true_vector)^2) / num_predictions
  rel_mse <- mse / var(true_vector)
  paste("Relative MSE:", rel_mse)
}

# MAE: Relative MAE < 1 => Model is better than the naive model
get_mae <- function(prediction_vector, true_vector) {
  num_predictions <- length(prediction_vector)
  mae <- sum(abs(prediction_vector - true_vector)) / num_predictions
  
  # Naive model assumes that the previous value is the prediction for the 
  # current time step
  naive_forecast <- abs(diff(true_vector))
  mae_naive <- mean(naive_forecast)
  
  rel_mae <- mae_model / mae_naive
  paste("Relative MAE:", rel_mae)
}


# Compute selectivity of inhibitor
num_kinases <- 442
# if threshold_K_d is set at 3000 nM => low affinity 
#                           => kinases likely to NOT be inhibited by that inhibitor
# if threshold_K_d is set at 300 nM => high affinity 
#                           => kinases likely to be inhibited by that inhibitor

get_S_inhibitor_vector <- function(pred_K_d_matrix, threshold_K_d) {
  n_matrix <- ifelse(pred_K_d_matrix < threshold_K_d, 1, 0) 
  # n_ij = 1 if the ith kinase is bound to the jth inhibitor, = 0 is not bound
  n_kinase_vector <- apply(n_matrix, 2, sum)
  pred_S_inhibitor_vector <- n_kinase / num_kinases
  return(pred_S_inhibitor_vector) # vector for all the inhibitors
}

# for now, set selectivity threshold to classify them as highly, medium or lowly selective
# Highly selective => picky => inhibitor binds to less kinases => n_kinase decreases
#                    => selectivity decreases
# selectivity scores range from 0 to 1 by definition as proportion of target kinases
# the inhibitor is bound to, so set at 0.5
S_threshold = 0.5
get_df_classify_inhibitors <- function() {
  S_inhib_bound <- get_S_inhibitor(pred_K_d_matrix, 300) # at two thresholds
  S_inhib_loose <- get_S_inhibitor(pred_K_d_matrix, 3000) # loose => not bound
  df_S_inhibitors <- dataframe(S_inhib_(300) = S_inhib_bound
                      , S_inhib_(3000) = S_inhib_loose) %>% 
                      mutate(Selectivity_Level = ifelse(S_inhib_(3000) < S_threshold, "Low",
                                                ifelse(S_inhib_(300) < S_threshold, "Medium",
                                                       "High")))
}
