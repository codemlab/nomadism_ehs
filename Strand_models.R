# From Alami et al 2026 nomadism paper, attempt to use STRAND model.


### try strand:
library(devtools)
install_github('ctross/STRAND@phosphorescent_desert_buttons')
library(STRAND)


strandsupport <- create_network("strand",data = alt, census = dvil, directed = T, multiple = F,
                                in_question = c("Natural Hazards IN", "Loans IN"), 
                                out_question = c("Natural Hazards OUT", "Loans OUT"),
                                outcome_mode = "bernoulli",
                                link_mode = "logit", check_standardization = F)


fitstrand = fit_block_plus_social_relations_model(data=strandsupport, data = alt, census = dvil,
  block_regression = ~ bone + sex + hid + former_nomad,
  focal_regression = ~ sex + age + edu + ln_income_total + ln_wealth + externals + former_nomad,
  target_regression = ~ sex + age + edu + ln_income_total + ln_wealth + externals + former_nomad,
  dyad_regression = ~ relatedness,
  mode="mcmc",
  stan_mcmc_parameters = list(chains = 2, parallel_chains = 2, refresh = 1,
                                iter_warmup = 500, iter_sampling = 500,
                                max_treedepth = NULL, adapt_delta = .98))
sumstrand <- summarize_strand_results(fitstrand)
vis_1 = strand_caterpillar_plot(sumstrand, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE, only_slopes=TRUE)

save.image("/home/edseab/github/codem/nomadism_models_01Sep2025.RData")
vis_1

# Now the latent model
latentsupport <- create_network("latent", data = alt, census = dvil, directed = T, multiple = F,
in_question = c("Natural Hazards IN", "Loans IN"), 
out_question = c("Natural Hazards OUT", "Loans OUT"),
outcome_mode = "bernoulli",
link_mode = "logit", check_standardization = F)

fitlatent = fit_latent_network_model(data=latentsupport,
  block_regression = ~ bone + sex + hid + former_nomad,
  focal_regression = ~ sex + age + edu + ln_income_total + ln_wealth + externals+ former_nomad,
  target_regression = ~ sex + age + edu + ln_income_total + ln_wealth + externals + former_nomad,
  dyad_regression = ~ relatedness,
  fpr_regression = ~ age + ln_wealth,
  rtt_regression = ~ age + ln_wealth,
  theta_regression = ~ 1,
  mode="mcmc",
  stan_mcmc_parameters = list(chains = 2, parallel_chains = 2, refresh = 1,
                                iter_warmup = 500, iter_sampling = 500,
                                max_treedepth = NULL, adapt_delta = .98))
sumlatent <- summarize_strand_results(fitlatent)
vis_2 = strand_caterpillar_plot(sumlatent, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), normalized=TRUE, only_slopes=TRUE)

save.image("/home/edseab/github/codem/nomadism_models_02Sep2025.RData")

