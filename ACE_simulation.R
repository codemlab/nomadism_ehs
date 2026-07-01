# simulate true probability of events
true_probs <- runif(1000, 0,1)

# simulate true outcomes as binary
true_outcomes <- rbinom(1000, 1, true_probs)

# simulate human forecasts as probabilities that differ from true probs with some normal error
humans <- (true_probs + rnorm(1000, 0, 0.1)) |> pmax(0) |> pmin(1) # constrain between 0 and 1

# Simulate larger number of LLMs
N_llms <- 100
llm_error <- runif(N_llms,0.1,1) # each llm differs in accuracy

# initialise list of vectors of forecasts for each llm
llms <- list()
for (i in 1:N_llms){
  llms[[i]] <- (true_probs + rnorm(1000, 0, llm_error[i])) |> pmax(0) |> pmin(1) #each LLM has different forecast depending on their accuracy, reflected in the std dev of the error
}

# Calculate accuracy and correlation according to text of article
# "Forecasting accuracy was measured as (1 - Brier Score), where the Brier Score (BS) represents the mean squared error between probabilistic forecasts and binary outcomes
brier_scores <- unlist(lapply(llms, function(x) mean((x - true_outcomes)^2)))
llm_accuracy <- 1 - brier_scores

# "As for Human-AI correlation, we calculated the Pearson correlation coefficient between each LLM's forecast set and the aggregated human forecasts"
hai_correlation <- unlist(lapply(llms,function(x) cor(x, humans)))


plot(llm_accuracy, hai_correlation)


