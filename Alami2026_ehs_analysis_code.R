######################################################################
# Data analysis EHBEA 2024 - Sarah
# coder: Ed Seabright, Sarah Alami
# commit: 
######################################################################

# Set Seed
set.seed(10101)

# Load packages
library(lme4)
library(broom.helpers)
library(ggdist)
library(ggpattern)
library(DiagrammeR)
library(sna)
library(dplyr)
library(tibble)
library(statnet)
library(ergm)
library(network)
library(glmmTMB)
library(kinship2)
library(ggplot2)
library(gt)
library(gtsummary)
library(broom.mixed)
library(brms)
library(tidyr)
library(janitor)
library(boot)
library(kableExtra)
library(scales)

dir.create("model_output/", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/", showWarnings = FALSE, recursive = TRUE)

# Load data
d <- read.csv("manuscript_data/main_dataset.csv")
alt <- read.csv("manuscript_data/social_networks.csv")
rel_mat <- read.csv("manuscript_data/relatedness_matrix.csv",row.names = 1)
nom <- read.csv("manuscript_data/nomadism_qual.csv")


# Collapse nomad sub-types 
d$group2 <- d$group
d$group2[which(d$group2=="part-time nomad")] <- "nomad"
d$group2[which(d$group2=="full-time nomad")] <- "nomad"
d$nomadic <- as.integer(d$group2 == "nomad")

nom$group2 <- nom$group
nom$group2[which(nom$group2=="part-time nomad")] <- "nomad"
nom$group2[which(nom$group2=="full-time nomad")] <- "nomad"
nom$nomadic <- as.integer(nom$group2 == "nomad")


# Separate village interview from nomad interview
dvil <- d[which(d$interviewed_econ==1),]
dnom <- d[which(d$interviewed_nomadism==1),]

# Income model
dinc <- dvil[!is.na(dvil$net_income & dvil$nomadic==0),]
dinc$income_std <- dinc$net_income/sd(dinc$net_income)
dinc$ln_income <- (dinc$net_income + 1)
dinc$age_std <- scale(dinc$age)
dinc$time_com_std <- scale(dinc$years_in_com)
dinc$yrs_edu_std <- scale(dinc$yrs_edu)
dinc$former_nomad <- as.numeric(dinc$group2 == "former nomad")

prior_list <- c(
  prior(normal(0, 1), class = "sd"),  # for group-level SDs
  prior(normal(0, 5), class = "Intercept")
)

inc_model_brm <- brm(
  bf(
    income_std ~ sex + age_std + former_nomad + (1 | householdID),
    hu ~ 1 + age_std + sex + former_nomad + (1 | householdID)
  ),
  family = hurdle_lognormal(),
  data = dinc,
  prior = prior_list,
  backend = "cmdstanr",
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 1,
  threads = threading(1),
  init = "0",
  save_model = T,
)
draws_inc <- as_draws_df(inc_model_brm)

dincfull <- dinc[!is.na(dinc$yrs_edu_std)& !is.na(dinc$time_com_std),]
inc_model_brm_full <- brm(
  bf(
    income_std ~ sex + age_std + former_nomad + time_com_std + darija + yrs_edu_std + (1 | householdID),
    hu ~ 1 + age_std + sex + former_nomad + time_com_std + darija + yrs_edu_std + (1 | householdID)
  ),
  family = hurdle_lognormal(),
  data = dincfull,
  prior = prior_list,
  backend = "cmdstanr",
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 2,
  threads = threading(1),
  init = "0",
  save_model = T,
)
draws_inc_full <- as_draws_df(inc_model_brm_full)

# wealth models
dw <- dinc[!is.na(dinc$wealth),]
dw$wealth_std <- (dw$wealth+1)/sd(dw$wealth)
dw$ln_wealth <- log(dw$wealth + 1)

wealth_model_brm <- brm(
  bf(
    wealth_std ~ sex + age_std + former_nomad + (1 | householdID)
  ),
  family = lognormal(),
  data = dw,
  prior = prior_list,
  backend = "cmdstanr",
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 2,
  threads = threading(1),
  init = "0",
  save_model = T,
)


draws_wealth <- as_draws_df(wealth_model_brm)

dwfull <- dw[!is.na(dw$yrs_edu_std)& !is.na(dw$time_com_std),]

wealth_model_brm_full <- brm(
  bf(
    wealth_std ~ sex + age_std + former_nomad + time_com_std + darija + yrs_edu_std + income_std +(1 | householdID)
  ),
  family = lognormal(),
  data = dwfull,
  prior = prior_list,
  backend = "cmdstanr",
  chains = 4,
  warmup = 1000,
  iter = 2000,
  cores = 2,
  threads = threading(1),
  init = "0",
  save_model = T,
)


draws_wealth_full <- as_draws_df(wealth_model_brm_full)
##==================== DIAGNOSTICS (brms models) ====================##

samples <- nuts_params(inc_model_brm)
table(samples$Value[samples$Parameter == "divergent__"])
sum(rhat(inc_model_brm)>=1.01)
sum(neff_ratio(inc_model_brm)< 0.5)
sum(neff_ratio(inc_model_brm)< 0.1)

samples <- nuts_params(inc_model_brm_full)
table(samples$Value[samples$Parameter == "divergent__"])
sum(rhat(inc_model_brm_full)>=1.01)
sum(neff_ratio(inc_model_brm_full)< 0.5)
sum(neff_ratio(inc_model_brm_full)< 0.1)

samples <- nuts_params(wealth_model_brm)
table(samples$Value[samples$Parameter == "divergent__"])
sum(rhat(wealth_model_brm)>=1.01)
sum(neff_ratio(wealth_model_brm)< 0.5)
sum(neff_ratio(wealth_model_brm)< 0.1)

samples <- nuts_params(wealth_model_brm_full)
table(samples$Value[samples$Parameter == "divergent__"])
sum(rhat(wealth_model_brm_full)>=1.01)
sum(neff_ratio(wealth_model_brm_full)< 0.5)
sum(neff_ratio(wealth_model_brm_full)< 0.1)

## network creation
source('code/network_creation_function.R')
### socialising:
socnet <- create_network("network", data = alt,
                                    census = dvil[dvil$nomadic==0,],
                                    in_question = c('Festivities IN','Social MEN',"Social WOMEN",'Couscous'), 
                                    out_question = 'Festivities OUT', 
                                    directed = T)

mod_soc<-  ergm(socnet$net~edges + mutual + nodefactor('Bone') +nodematch('Bone') + nodecov('Age') +
                nodefactor('Sex') + nodeifactor('Sex') + nodematch('Sex')+ 
                nodematch('HouseholdID')  + nodecov('Externals') + nodeicov('Externals') +
                nodecov('Yrs schooling') + nodeicov('Yrs schooling') +
                nodecov('Income (ln)') + nodeicov('Income (ln)') + absdiff('Income (ln)') + 
                nodecov('Wealth (ln)') + nodeicov('Wealth (ln)') + absdiff('Wealth (ln)') + 
                nodefactor("Former nomad") + nodeifactor("Former nomad") + nodematch("Former nomad") + 
                edgecov(socnet$relmat))


socgraph <- create_network("graph", data = alt, census= dvil[dvil$nomadic==0,],
in_question = c('Festivities IN', 'Couscous','Social MEN',"Social WOMEN"), out_question = 'Festivities OUT')

render_graph(socgraph)

### work:
worknet <- create_network("network", data = alt, census = dvil[dvil$nomadic==0,], in_question = c('Work IN', 'Daily Labour',"Childcare"), out_question = 'Work OUT')

mod_work<-   ergm(worknet$net~edges + mutual + nodefactor('Bone') +nodematch('Bone') + nodecov('Age') +
                nodefactor('Sex') + nodeifactor('Sex') + nodematch('Sex')+ 
                nodematch('HouseholdID')  + nodecov('Externals') + nodeicov('Externals') +
                nodecov('Yrs schooling') + nodeicov('Yrs schooling') +
                nodecov('Income (ln)') + nodeicov('Income (ln)') + absdiff('Income (ln)') + 
                nodecov('Wealth (ln)') + nodeicov('Wealth (ln)') + absdiff('Wealth (ln)') + 
                nodefactor("Former nomad") + nodeifactor("Former nomad") + nodematch("Former nomad") + 
                edgecov(worknet$relmat))
workgraph <- create_network("graph", data = alt, census = dvil[dvil$nomadic==0,],in_question = c('Work IN', 'Daily Labour',"Childcare"), out_question = c("Work OUT"))

render_graph(workgraph)

### support:
supportnet <- create_network("network", data = alt, census = dvil[dvil$nomadic==0,], directed = T, multiple = F,in_question = c("Natural Hazards IN", "Loans IN"), out_question = c("Natural Hazards OUT", "Loans OUT"))

mod_support<-  ergm(supportnet$net~edges + mutual + nodefactor('Bone') +nodematch('Bone') + nodecov('Age') +
                nodefactor('Sex') + nodeifactor('Sex') + nodematch('Sex')+ 
                nodematch('HouseholdID')  + nodecov('Externals') + nodeicov('Externals') +
                nodecov('Yrs schooling') + nodeicov('Yrs schooling') +
                nodecov('Income (ln)') + nodeicov('Income (ln)') + absdiff('Income (ln)') + 
                nodecov('Wealth (ln)') + nodeicov('Wealth (ln)') + absdiff('Wealth (ln)') + 
                nodefactor("Former nomad") + nodeifactor("Former nomad") + nodematch("Former nomad") + 
                edgecov(supportnet$relmat))


supportgraph <- create_network("graph",data = alt, census = dvil[dvil$nomadic==0,],in_question = c("Natural Hazards IN", "Loans IN"), out_question = c("Natural Hazards OUT", "Loans OUT"))

render_graph(supportgraph)



# Age at settlement categories (needed by inline scalars)
nom$age_cats <- cut(
  nom$age_settled,
  breaks = c(-Inf, 17, 25, 40, Inf),
  labels = c("Child/Teenager", "Young adult", "Adult", "Middle-aged/Older adult"),
  right = TRUE
)

# Aspiration to settle — freq_why_want referenced in qmd
search_terms <- c("tired", "livestock", "protection", "school", "infrastructure", "engaged", "bored")
freq_why_want <- sapply(search_terms, function(x) sum(grepl(x, nom$want_settle_why_cat)))

##==================== DEMOGRAPHY (manuscript inline scalars) ====================##



# d$group2 already created above
n <- nrow(d)
sex_table <- prop.table(table(d$group2, d$sex), margin = 1) |> data.frame()
sex_table <- sex_table[sex_table$Var2 == "female", ]

never_nom   <- d[d$group2 == "never nomad", ]
nn_prop     <- nrow(never_nom) / n
nn_sex      <- sex_table$Freq[sex_table$Var1 == "never nomad"]

former_nom  <- d[d$group2 == "former nomad", ]
fn_prop     <- nrow(former_nom) / n
fn_sex      <- sex_table$Freq[sex_table$Var1 == "former nomad"]

current_nom <- d[d$group2 == "nomad", ]
cn_prop     <- nrow(current_nom) / n
cn_sex      <- sex_table$Freq[sex_table$Var1 == "nomad"]

d |> group_by(group2) |>
  summarise(
    mean_age = round(mean(age), 1),
    sd_age   = round(sd(age), 1),
    mean_edu = round(mean(yrs_edu, na.rm=T), 1),
    sd_edu   = round(sd(yrs_edu, na.rm=T), 1),
    dar_pct  = round(sum(darija == "fluent") * 100 / n(), 1),
    .groups  = "drop"
  ) -> des_s

mean_ages_range <- range(des_s$mean_age) |> round(1)
sex_props_range <- range(sex_table$Freq * 100) |> round(1)

# Water/food insecurity from d (for bootstrap CIs)
likert_levels_des <- c("never", "rarely", "sometimes", "often", "always")

ws_des <- d[, c("nomadic", "water_need", "water_interruption", "water_change_plans")]
ws_des$water_need         <- as.numeric(factor(ws_des$water_need,         levels = likert_levels_des)) - 1
ws_des$water_interruption <- as.numeric(factor(ws_des$water_interruption, levels = likert_levels_des)) - 1
ws_des$water_change_plans <- as.numeric(factor(ws_des$water_change_plans, levels = likert_levels_des)) - 1
ws_des$water_insecurity   <- ws_des$water_need + ws_des$water_interruption + ws_des$water_change_plans
ws_des <- ws_des[!is.na(ws_des$water_insecurity), ]

fs_des <- d[, c("nomadic", "fewer_meals", "no_food", "hungry_alldaynight")]
fs_des$fewer_meals        <- as.numeric(factor(fs_des$fewer_meals,        levels = likert_levels_des)) - 1
fs_des$no_food            <- as.numeric(factor(fs_des$no_food,            levels = likert_levels_des)) - 1
fs_des$hungry_alldaynight <- as.numeric(factor(fs_des$hungry_alldaynight, levels = likert_levels_des)) - 1
fs_des$food_insecurity    <- fs_des$fewer_meals + fs_des$no_food + fs_des$hungry_alldaynight
fs_des <- fs_des[!is.na(fs_des$food_insecurity), ]

mean_stat <- function(data, indices) mean(data[indices])

ws_des |> group_by(nomadic) |>
  summarise(mean = round(mean(water_insecurity)), sd = round(sd(water_insecurity)),
            ci_lower = NA_real_, ci_upper = NA_real_, .groups = "drop") |>
  as.data.frame() -> ws_s

for (i in 0:1) {
  ci <- boot(ws_des$water_insecurity[ws_des$nomadic == i], mean_stat, R = 5000) |>
    boot.ci(type = "perc")
  ws_s[i + 1, c("ci_lower", "ci_upper")] <- round(ci$percent[4:5], 1)
}

fs_des |> group_by(nomadic) |>
  summarise(mean = round(mean(food_insecurity), 1), sd = round(sd(food_insecurity), 1),
            ci_lower = NA_real_, ci_upper = NA_real_, .groups = "drop") |>
  as.data.frame() -> fs_s

for (i in 0:1) {
  ci <- boot(fs_des$food_insecurity[fs_des$nomadic == i], mean_stat, R = 5000) |>
    boot.ci(type = "perc")
  fs_s[i + 1, c("ci_lower", "ci_upper")] <- round(ci$percent[4:5], 1)
}

##==================== SEDENTARISATION (manuscript inline scalars) ====================##
# Sedentarisation process: age at settlement
n_nom_sub         <- sum(nom$group == "former nomad")
n_age_adult       <- sum(nom$age_cats == "Adult",                   na.rm = T)
pct_age_adult     <- round(n_age_adult  / n_nom_sub * 100, 1)
n_age_young       <- sum(nom$age_cats == "Young adult",             na.rm = T)
pct_age_young     <- round(n_age_young  / n_nom_sub * 100, 1)
n_age_older       <- sum(nom$age_cats == "Middle-aged/Older adult", na.rm = T)
pct_age_older     <- round(n_age_older  / n_nom_sub * 100, 1)
n_age_child       <- sum(nom$age_cats == "Child/Teenager",          na.rm = T)
pct_age_child     <- round(n_age_child  / n_nom_sub * 100, 1)

nomads_raw    <- nom[nom$group2 == "nomad", ]
ft_nomads_raw <- nom[nom$group == "full-time nomad", ]
fnom_raw      <- nom[nom$group == "former nomad", ]

n_former_nomad  <- sum(d$group == "former nomad")
n_part_nomad    <- sum(d$group == "part-time nomad")
n_all_non_nomad <- n_former_nomad + n_part_nomad

search_terms    <- c("tired", "school", "livestock", "marriage", "conflict", "family", "defend", "overcrowding")
search_term     <- search_terms
freq_why_settle <- sapply(search_terms, function(x) sum(grepl(x, nom$why_settle_cats), na.rm = TRUE))
names(freq_why_settle) <- search_terms

pct_why_settle_tired <- round(freq_why_settle["tired"] / n_nom_sub * 100, 1)

# Sedentarisation process: who decided
own_pat     <- "^Me$|^Me,|^I am|^Me and"
n_decided   <- sum(!is.na(fnom_raw$who_decided))
n_not_own   <- sum(!grepl(own_pat, fnom_raw$who_decided) & !is.na(fnom_raw$who_decided))
pct_not_own <- round(n_not_own / n_decided * 100, 1)

# Sedentarisation process: why this village
swwy                 <- nom$settle_where_why_cat[nom$group == "former nomad"]
n_settle_immediate   <- sum(swwy == "support_from_immediate_family", na.rm = T)
pct_settle_immediate <- round(n_settle_immediate / n_nom_sub * 100, 1)
n_settle_inlaws      <- sum(swwy == "support_from_in_laws",          na.rm = T)
pct_settle_inlaws    <- round(n_settle_inlaws     / n_nom_sub * 100, 1)
n_settle_unrelated   <- sum(swwy == "support_from_unrelated",        na.rm = T)
pct_settle_unrelated <- round(n_settle_unrelated  / n_nom_sub * 100, 1)
n_settle_purchase    <- sum(swwy == "purchase_house_land",           na.rm = T)
pct_settle_purchase  <- round(n_settle_purchase   / n_nom_sub * 100, 1)
n_settle_inherited   <- sum(swwy == "inherited_land_house",          na.rm = T)
pct_settle_inherited <- round(n_settle_inherited  / n_nom_sub * 100, 1)
n_settle_marriage    <- sum(swwy == "marriage_to_villager",          na.rm = T)
pct_settle_marriage  <- round(n_settle_marriage   / n_nom_sub * 100, 1)

# Aspirations: want to settle (full-time nomads only, from nom_raw)
n_ft_nomads  <- sum(nom$group == "full-time nomad")
n_want_yes   <- sum(ft_nomads_raw$want_settle == "yes", na.rm = T)
pct_want_yes <- round(n_want_yes / n_ft_nomads * 100, 1)
n_want_no    <- sum(ft_nomads_raw$want_settle == "no",  na.rm = T)
# freq_why_want computed earlier from nom$want_settle_why_cat

# Social support: current nomads (qual.csv)
n_nom_support      <- sum(!is.na(nom$nomad_support[nom$group2 == "nomad"]))
pct_support_never  <- round(mean(nom$nomad_support[nom$group2=="nomad"] == "never",  na.rm=T)*100, 1)
pct_support_rarely <- round(mean(nom$nomad_support[nom$group2=="nomad"] == "rarely", na.rm=T)*100, 1)
pct_support_often  <- round(mean(nom$nomad_support[nom$group2=="nomad"] == "often",  na.rm=T)*100, 1)
pct_support_always <- round(mean(nom$nomad_support[nom$group2=="nomad"] == "always", na.rm=T)*100, 1)

# Social support: change since settling (former + part-time nomads, qual.csv)
cs_fnom            <- nom$change_support[nom$group != "full-time nomad"]
n_support_change   <- sum(!is.na(cs_fnom))
n_support_easier   <- sum(grepl("easier", cs_fnom), na.rm = T)
pct_support_easier <- round(n_support_easier / n_support_change * 100, 1)
n_support_same     <- sum(cs_fnom == "same", na.rm = T)
pct_support_same   <- round(n_support_same   / n_support_change * 100, 1)
n_support_harder   <- sum(grepl("harder", cs_fnom), na.rm = T)
pct_support_harder <- round(n_support_harder / n_support_change * 100, 1)

# Respect: current nomads (nom_raw)
n_nom_resp           <- nrow(nomads_raw)
n_resp_oth_always    <- sum(nomads_raw$nomad_respect_others  == "always", na.rm=T)
pct_resp_oth_always  <- round(n_resp_oth_always  / n_nom_resp * 100, 1)
n_resp_oth_often     <- sum(nomads_raw$nomad_respect_others  == "often",  na.rm=T)
pct_resp_oth_often   <- round(n_resp_oth_often   / n_nom_resp * 100, 1)
n_resp_oth_rarely    <- sum(nomads_raw$nomad_respect_others  == "rarely", na.rm=T)
pct_resp_oth_rarely  <- round(n_resp_oth_rarely  / n_nom_resp * 100, 1)
n_resp_oth_never     <- sum(nomads_raw$nomad_respect_others  == "never",  na.rm=T)
pct_resp_oth_never   <- round(n_resp_oth_never   / n_nom_resp * 100, 1)
n_resp_set_always    <- sum(nomads_raw$nomad_respect_settled == "always", na.rm=T)
pct_resp_set_always  <- round(n_resp_set_always  / n_nom_resp * 100, 1)
n_resp_set_often     <- sum(nomads_raw$nomad_respect_settled == "often",  na.rm=T)
pct_resp_set_often   <- round(n_resp_set_often   / n_nom_resp * 100, 1)
n_resp_set_rarely    <- sum(nomads_raw$nomad_respect_settled == "rarely", na.rm=T)
pct_resp_set_rarely  <- round(n_resp_set_rarely  / n_nom_resp * 100, 1)
n_resp_set_never     <- sum(nomads_raw$nomad_respect_settled == "never",  na.rm=T)
pct_resp_set_never   <- round(n_resp_set_never   / n_nom_resp * 100, 1)

# Respect: former nomads past (nom_raw) and present
n_fnom_resp              <- nrow(fnom_raw)
n_fnom_resp_oth_always   <- sum(fnom_raw$nomad_respect_others  == "always", na.rm=T)
pct_fnom_resp_oth_always <- round(n_fnom_resp_oth_always / n_fnom_resp * 100, 1)
n_fnom_resp_set_always   <- sum(fnom_raw$nomad_respect_settled == "always", na.rm=T)
pct_fnom_resp_set_always <- round(n_fnom_resp_set_always / n_fnom_resp * 100, 1)
n_fnom_resp_oth_low      <- sum(fnom_raw$nomad_respect_others  %in% c("rarely","never"), na.rm=T)
pct_fnom_resp_oth_low    <- round(n_fnom_resp_oth_low / n_fnom_resp * 100, 1)
n_fnom_resp_set_low      <- sum(fnom_raw$nomad_respect_settled %in% c("rarely","never"), na.rm=T)
pct_fnom_resp_set_low    <- round(n_fnom_resp_set_low / n_fnom_resp * 100, 1)

# Conflicts: current nomads (qual.csv)
n_conf_nom_never  <- sum(nom$nomad_conflicts[nom$group2=="nomad"] == "never",  na.rm=T)
pct_conf_nom_never <- round(n_conf_nom_never / n_nom_support * 100, 1)
n_conf_nom_rarely <- sum(nom$nomad_conflicts[nom$group2=="nomad"] == "rarely", na.rm=T)
pct_conf_nom_rarely <- round(n_conf_nom_rarely / n_nom_support * 100, 1)
n_conf_nom_often  <- sum(nom$nomad_conflicts[nom$group2=="nomad"] == "often",  na.rm=T)
pct_conf_nom_often <- round(n_conf_nom_often / n_nom_support * 100, 1)

# Conflicts: former nomads past (qual.csv)
n_conf_fnom_never   <- sum(nom$nomad_conflicts[nom$group=="former nomad"] == "never",  na.rm=T)
pct_conf_fnom_never <- round(n_conf_fnom_never  / n_nom_sub * 100, 1)
n_conf_fnom_rarely  <- sum(nom$nomad_conflicts[nom$group=="former nomad"] == "rarely", na.rm=T)
pct_conf_fnom_rarely <- round(n_conf_fnom_rarely / n_nom_sub * 100, 1)
n_conf_fnom_often   <- sum(nom$nomad_conflicts[nom$group=="former nomad"] == "often",  na.rm=T)
pct_conf_fnom_often <- round(n_conf_fnom_often  / n_nom_sub * 100, 1)

# Conflicts: former nomads settled (qual.csv)
n_sconf_fnom_never     <- sum(nom$settled_conflicts[nom$group=="former nomad"] == "never",     na.rm=T)
pct_sconf_fnom_never   <- round(n_sconf_fnom_never   / n_nom_sub * 100, 1)
n_sconf_fnom_rarely    <- sum(nom$settled_conflicts[nom$group=="former nomad"] == "rarely",    na.rm=T)
pct_sconf_fnom_rarely  <- round(n_sconf_fnom_rarely  / n_nom_sub * 100, 1)
n_sconf_fnom_sometimes <- sum(nom$settled_conflicts[nom$group=="former nomad"] == "sometimes", na.rm=T)
pct_sconf_fnom_sometimes <- round(n_sconf_fnom_sometimes / n_nom_sub * 100, 1)

# Freedom/values: current nomads and former nomads (nom_raw)
pct_nom_free_often   <- round(mean(nomads_raw$nomad_freedom == "often",  na.rm=T)*100, 1)
pct_nom_free_always  <- round(mean(nomads_raw$nomad_freedom == "always", na.rm=T)*100, 1)
pct_fnom_free_often  <- round(mean(fnom_raw$nomad_freedom   == "often",  na.rm=T)*100, 1)
pct_fnom_free_always <- round(mean(fnom_raw$nomad_freedom   == "always", na.rm=T)*100, 1)
pct_fnom_sfree_often  <- round(mean(fnom_raw$settled_freedom == "often",  na.rm=T)*100, 1)
pct_fnom_sfree_always <- round(mean(fnom_raw$settled_freedom == "always", na.rm=T)*100, 1)
pct_fnom_sfree_rarely <- round(mean(fnom_raw$settled_freedom == "rarely", na.rm=T)*100, 1)
n_nom_free          <- sum(!is.na(nomads_raw$nomad_freedom))
n_fnom_free         <- sum(!is.na(fnom_raw$nomad_freedom) & !is.na(fnom_raw$settled_freedom))
pct_nom_free_high   <- round(mean(nomads_raw$nomad_freedom %in% c("always","often"), na.rm=T)*100, 1)
pct_fnom_free_high  <- round(mean(fnom_raw$nomad_freedom   %in% c("always","often"), na.rm=T)*100, 1)
pct_fnom_sfree_high <- round(mean(fnom_raw$settled_freedom %in% c("always","often"), na.rm=T)*100, 1)

# Regrets: would you go back? (nom_raw, former nomads)
n_go_back_total    <- sum(!is.na(fnom_raw$go_back))
n_no_go_back       <- sum(fnom_raw$go_back == "no",     na.rm=T)
pct_no_go_back     <- round(n_no_go_back     / n_go_back_total * 100, 1)
n_yes_go_back      <- sum(fnom_raw$go_back == "yes",    na.rm=T)
pct_yes_go_back    <- round(n_yes_go_back    / n_go_back_total * 100, 1)
n_unsure_go_back   <- sum(fnom_raw$go_back == "unsure", na.rm=T)
pct_unsure_go_back <- round(n_unsure_go_back / n_go_back_total * 100, 1)

##==================== DISCUSSION INLINE SCALARS ====================##

# Education and Darija convenience vars (from des_s)
edu_nom   <- des_s$mean_edu[des_s$group2 == "nomad"]
edu_fnom  <- des_s$mean_edu[des_s$group2 == "former nomad"]
edu_never <- des_s$mean_edu[des_s$group2 == "never nomad"]
dar_nom   <- des_s$dar_pct[des_s$group2 == "nomad"]
dar_fnom  <- des_s$dar_pct[des_s$group2 == "former nomad"]
dar_never <- des_s$dar_pct[des_s$group2 == "never nomad"]

# ERGM: labour network outgoing ties OR + 95% Wald CI
se_work        <- sqrt(diag(vcov(mod_work)))
or_work_out    <- round(exp(coef(mod_work)["nodefactor.Former nomad.Yes"]), 2)
ci_work_out_lo <- round(exp(coef(mod_work)["nodefactor.Former nomad.Yes"] - 1.96 * se_work["nodefactor.Former nomad.Yes"]), 2)
ci_work_out_hi <- round(exp(coef(mod_work)["nodefactor.Former nomad.Yes"] + 1.96 * se_work["nodefactor.Former nomad.Yes"]), 2)
pct_work_out   <- round((or_work_out - 1) * 100)

# ERGM: homophily (% greater odds) across three networks
or_work_hom  <- round(exp(coef(mod_work)["nodematch.Former nomad"]), 2)
or_soc_hom   <- round(exp(coef(mod_soc)["nodematch.Former nomad"]), 2)
or_sup_hom   <- round(exp(coef(mod_support)["nodematch.Former nomad"]), 2)
pct_work_hom <- round((or_work_hom - 1) * 100)
pct_soc_hom  <- round((or_soc_hom  - 1) * 100)
pct_sup_hom  <- round((or_sup_hom  - 1) * 100)

##==================== WEALTH/INCOME INLINE SCALARS ====================##

exchange <- 0.09838  # exchange rate June 2023
d$income     <- round(d$income     * exchange)
d$net_income <- round(d$net_income * exchange)
d$wealth     <- round(d$wealth     * exchange)

##==================== SAVE ====================##


##================================ Relative Need ===============================##
# COLORS
#Nomad#D55E00
#Former nomad#FFC800
#Never nomad#0072B2
#Combined villager#009E73

labels_need <- c("less", "same","more")
nom$nomad_relative_need <- factor(nom$nomad_relative_need, levels = labels_need)
nom$settled_relative_need <- factor(nom$settled_relative_need, levels = labels_need)

freq_need <- data.frame(
  "current nomads: relative to nomads"                  = as.numeric(table(nom$nomad_relative_need[nom$nomadic==1])),
  "former nomads(past): relative to nomads"             = as.numeric(table(nom$nomad_relative_need[nom$nomadic==0])),
  "former nomads(present): relative to other villagers" = as.numeric(table(nom$settled_relative_need[nom$nomadic==0])),
  row.names = labels_need, check.names = FALSE)

df_need <- freq_need |>
  rownames_to_column("category") |>
  pivot_longer(-category, names_to = "group", values_to = "freq") |>
  group_by(group) |>
  mutate(prct = freq / sum(freq) * 100) |>
  ungroup()

      des_perceived_need <- ggplot(df_need, aes(x = factor(category, levels = c("less", "same", "more")), y = prct, fill = group)) +
        geom_bar_pattern(stat = "identity", position = "dodge", aes(pattern = group), pattern_density = 0.05, pattern_spacing = 0.01, colour = "black") +
        labs( x = "Response categories",
             y = "Percent participant in each group") +
        theme_bw() +
        scale_fill_manual(values = c("current nomads: relative to nomads" = "#D55E00", "former nomads(past): relative to nomads" = "#FFC800", "former nomads(present): relative to other villagers" = "#FFC800"))+
        scale_pattern_manual(name = "group", values = c("current nomads: relative to nomads" = "none", "former nomads(past): relative to nomads" = "none", "former nomads(present): relative to other villagers" = "stripe"))+
          theme(
            axis.title.x = element_text(size = 18),
            axis.title.y = element_text(size = 18),
            axis.text.x = element_text(size = 16),
            axis.text.y = element_text(size = 16),
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 16)) 
      
            pdf("perceived_need.pdf",width=12,height=9)
            des_perceived_need
            dev.off()
ggsave("figures/perceived_need.png", des_perceived_need, width = 12, height = 9)

fnom_set_relneed <- round(df_need$prct[9],1)
fnom_nom_relneed <- round(df_need$prct[8],1)
nom_nom_relneed <- round(df_need$prct[7],1)


##================= Perceived access to social support ==================##
# nomad and former nomads access to support 
nom$nomad_support <- factor(
  nom$nomad_support,
  levels = c("never", "rarely", "often", "always"),
  ordered = TRUE
)
nom %>%
  filter(!is.na(nomad_support)) %>%
  count(group, nomad_support) %>%
  group_by(group) %>%
  mutate(
    percent = round(100 * n / sum(n), 1)
  ) %>%
  ungroup()

# change in support since settling 
nom$change_support_simple <- case_when(grepl("harder",nom$change_support) ~ "harder",
                           grepl("easier", nom$change_support) ~ "easier")

prop.table(table(nom$change_support_simple)) * 100

##================= Respect ==================##
table(nom$group, nom$nomad_respect_others)

#nomads respect by other nomads (now)
# ALWAYS : 59/75 (78.67%)
# OFTEN : 8/75 (10.67%) 
# RARELY : 1/75 (1.33%)
# NEVER : 7/75 (9.33%)

# former nomads respect by other nomads (past)
# ALWAYS : 17/35 (45.71%)
# OFTEN : 9/35 (25.71%) 
# RARELY : 3/35 (8.57%)
# NEVER : 6/35 (17.14%)

table(nom$group, nom$nomad_respect_settled)

#nomads respect by settled neighbours (now)
# ALWAYS : 55/75 (73.33%)
# OFTEN : 12/75 (16.00%)
# RARELY : 3/75 (4.00%)
# NEVER : 5/75 (6.67%)

# former nomads respect by settled neighbours (now)
# ALWAYS : 10/35 (28.57%)
# OFTEN : 10/35 (28.57%)
# RARELY : 8/35 (22.86%)
# NEVER : 7/35 (20.00%)

##================= Conflicts ==================##

table(nom$group,nom$nomad_conflicts)

# current nomad
# OFTEN: 11/75(14.7%)
# SOMETIMES: 0%
# RARELY: 3/75 (4%)
# NEVER: 61/75 (81.3%)

# former nomad in past nomadic life
# OFTEN: 5/35 (14.3%)
# SOMETIMES: 0%
# RARELY: 9/35 (25.7%)
# NEVER: 21/35 (60.0%)

table(nom$group,nom$settled_conflicts)

# Former nomads in present village
# OFTEN: 0%
# SOMETIMES: 2/35 (5.7%)
# RARELY: 5/35 (14.3%)
# NEVER: 27/35 (77.1%)

### Advantages of nomadism
labels_advnom <- c("Environment: space, quiet, nature", "Freedom, self-sufficiency","Mental wellbeing, ability to avoid social problems", "Caring for livestock","Good life contigent on rainfall", "No advantages of nomadic life","Food and cultural traditions","Identity and pride in being a nomad")
search_terms <- c("atmosphere", "freedom", "peace", "livestock", "conditional", "nothing", "food", "identity")
freq_advnom <- sapply(search_terms,function(x)c(sum(grepl(x,nom$advnom_cat[nom$group2=="nomad"])),sum(grepl(x,nom$advnom_cat[nom$group2=="former nomad"]))))

n_nom_adv_total  <- sum(nom$group2 == "nomad")
pct_nom_peace    <- round(freq_advnom[1, "peace"] / n_nom_adv_total * 100, 1)
pct_fnom_peace   <- round(freq_advnom[2, "peace"] / n_nom_sub * 100, 1)
pct_all_peace    <- round(sum(freq_advnom[, "peace"]) / (n_nom_adv_total + n_nom_sub) * 100, 1)

df_advnom <- data.frame(perceived_advantages = rep(labels_advnom,2),
                        count = c(freq_advnom[1,], freq_advnom[2,]),
                        total_group = rep(c(sum(nom$group2=="nomad"),sum(nom$group2=="former nomad")), each = length(labels_advnom)),
                        group = rep(c("nomad","former nomad"), each = length(labels_advnom)))


df_advnom$perc <- (df_advnom$count) * 100 /df_advnom$total_group

# PLOTS
# COLORS
col_nomad           <- "#D55E00"
col_former_nomad    <- "#FFC800"
col_never_nomad     <- "#0072B2"
col_combined_villager <- "#009E73"

pal_group <- c(
  "nomad"             = col_nomad,
  "former nomad"      = col_former_nomad,
  "never nomad"       = col_never_nomad,
  "combined villager" = col_combined_villager
)

order_nomads <- df_advnom %>%
  filter(group == "nomad") %>%
  select(perceived_advantages, perc)

df_advnom_ord <- df_advnom %>%
  left_join(order_nomads, by = "perceived_advantages",
            suffix = c("", "_nomads")) %>%
  mutate(
    perceived_advantages = reorder(perceived_advantages, perc_nomads)
  )


p_advnom <- ggplot(df_advnom_ord,
    aes(x = perceived_advantages,
        y = perc,
        fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "nomad" = col_nomad,
      "former nomad" = col_former_nomad
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = "Perceived advantages of nomadic life",
    y = "Percentage within group (%)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 13)

# Advantages of settling

labels_advset <- c("Comfort/lower physical effort", "Schooling", "Healthcare", "Infrastructure", "Leisure", "Religion", "Safety and protection from elements", "Sanitation and cleanliness", "Services and markets", "Social life and support")
search_terms <- c("comfort", "edu", "health", "infra", "leisure", "reli", "safety", "sanitation", "market", "support")
freq_advset <- sapply(search_terms, function(x) c(sum(grepl(x, nom$advset_cat[nom$group2=="nomad"])), sum(grepl(x, nom$advset_cat[nom$group2=="former nomad"]))))

df_advset <- data.frame(
  perceived_advantages = rep(labels_advset, 2),
  count = c(freq_advset[1,], freq_advset[2,]),
  total_group = rep(c(sum(nom$group2=="nomad"), sum(nom$group2=="former nomad")), each = length(labels_advset)),
  group = rep(c("nomad", "former nomad"), each = length(labels_advset))
)

df_advset$perc <- df_advset$count * 100 / df_advset$total_group

order_fnom <- df_advset %>%
  filter(group == "former nomad") %>%
  select(perceived_advantages, perc)

df_advset_ord <- df_advset %>%
  left_join(order_fnom, by = "perceived_advantages",
            suffix = c("", "_fnom")) %>%
  mutate(
    perceived_advantages = reorder(perceived_advantages, perc_fnom)
  )

p_advset <- ggplot(df_advset_ord,
  aes(x = perceived_advantages,
      y = perc,
      fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "nomad" = col_nomad,
      "former nomad" = col_former_nomad
    ),
    guide = guide_legend(reverse = TRUE)
  ) +
  labs(
    x = "Perceived advantages of settled life",
    y = "Percentage within group (%)",
    fill = "Group"
  ) +
  theme_minimal(base_size = 13)

library(patchwork)
(p_advnom / p_advset) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom")

ggsave("figures/perceived_advantages.png", width = 12, height = 14)


##==================== POSTERIOR PLOTS ====================##

p_income <- bayesplot::mcmc_areas(
  as.matrix(inc_model_brm),
  regex_pars = "^b_",
  prob = 0.9
)
ggsave("figures/income_model.png", p_income, width = 8, height = 6)

p_wealth <- bayesplot::mcmc_areas(
  as.matrix(wealth_model_brm),
  regex_pars = "^b_",
  prob = 0.9
)
ggsave("figures/wealth_model.png", p_wealth, width = 8, height = 4)

##==================== MODELS OVER TIME ====================##

wealth_model <- glmmTMB(ln_wealth ~ sex + age + former_nomad + (1|householdID), data = dw)
inc_nomads <- lm(ln_income ~ years_in_com + sex + age, data = dinc[which(dinc$group2=='former nomad'),])
wealth_nomads <- lm(ln_wealth ~ years_in_com*former_nomad + sex + age, data = dw)

dw$wealth_resid <- residuals(lmer(ln_wealth ~ sex + age +(1|householdID), data = dw))
dw$predicted_wealth <- predict(wealth_nomads, newdata = dw)

ggplot(dw[which(dw$wealth>0 & !is.na(dw$years_in_com)),], aes(x = years_in_com, y = wealth_resid, color = factor(former_nomad), fill = factor(former_nomad))) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", aes(group = factor(former_nomad)), se = T) +
  labs(
    title = "Interaction between Years in Community and Former Nomad Status",
    x = "Years in Community",
    y = "Adjusted Log Wealth (Residuals from Age & Sex)",
    color = "Family formerly nomadic", fill = "Family formerly nomadic"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("firebrick", "forestgreen")) +
  scale_fill_manual(values = c("firebrick", "forestgreen"))
ggsave("figures/wealth_interaction.png", width = 12, height = 6)

# income nomad comparison
dall_inc <- d[!is.na(d$income),]
ggplot(dall_inc, aes(x = group2, y = income, fill = group2)) +
  geom_boxplot(width = 0.3, na.rm = TRUE) +
  stat_slab(side = "right", scale = 0.4, na.rm = TRUE) +
  theme_minimal(base_size = 25) +
  labs(y = "Annual income (USD)") +
  scale_fill_manual(values = c("former nomad" = "slateblue", "never nomad" = "forestgreen", "nomad" = "coral"))
ggsave("figures/income_comparison.png", width = 12, height = 6)


# wealth nomad comparison
dall_w <- d[!is.na(d$wealth),]
ggplot(dall_w, aes(x = group2, y = wealth, fill = group2)) +
  geom_boxplot(na.rm = TRUE) +
  labs(y = "Wealth (USD)") +
  theme_minimal()

##==================== SUPPLEMENTARY DESCRIPTIVES ====================##

models <- list(
  "Labour (agriculture + childcare)" = mod_work,
  "Socialisation and festivities"     = mod_soc,
  "Material support (including loans)" = mod_support
)

source("code/dag_plot.R")

des_education_by_group <- ggplot(d, aes(x = factor(group2, levels = c("nomad", "former nomad", "never nomad")), y = yrs_edu, fill = group2)) +
  geom_boxplot(alpha = 0.8, outlier.shape = 21, na.rm = TRUE) +
  labs(x = NULL, y = "Years of education") +
  theme_bw(base_size = 16) +
  scale_fill_manual(values = c("nomad" = "#D55E00", "former nomad" = "#FFC800", "never nomad" = "#0072B2")) +
  theme(legend.position = "none")
ggsave("figures/des_education_by_group.png", des_education_by_group, width = 8, height = 6)

# education table
d %>%
  group_by(group2) %>%
  summarise("Years of education" = round(mean(yrs_edu, na.rm=T), 1),
            "% Arabic speakers"  = round(mean(darija, na.rm=T), 2)*100) %>%
  arrange(match(group2, c("nomad", "former nomad", "never nomad"))) %>%
  kbl(booktabs = T) %>%
  kable_styling(latex_options = "striped")


## ======================= DEMOGRAPHICS SUMMARY TABLE =======================##

n_age    <- sum(!is.na(d$age))
n_edu    <- sum(!is.na(d$yrs_edu))
n_female <- sum(!is.na(d$sex))

numeric_summary <- d %>%
  summarise(
    Age     = paste0(round(mean(age, na.rm=TRUE), 1), " (", median(age, na.rm=TRUE),
                     ", ", min(age, na.rm=TRUE), "-", max(age, na.rm=TRUE), "; N=", n_age, ")"),
    Yrs_Edu = paste0(round(mean(yrs_edu, na.rm=TRUE), 1), " (", median(yrs_edu, na.rm=TRUE),
                     ", ", min(yrs_edu, na.rm=TRUE), "-", max(yrs_edu, na.rm=TRUE), "; N=", n_edu, ")"),
    Female  = paste0(round(mean(sex == "female", na.rm=TRUE)*100, 1),
                     "% (", sum(sex == "female", na.rm=TRUE), "/", n_female, ")")
  ) %>%
  pivot_longer(cols = everything(), names_to = "Variable", values_to = "Summary") %>%
  mutate(Variable = case_when(
    Variable == "Age"     ~ "Age (mean, median, range)",
    Variable == "Yrs_Edu" ~ "Years of Education (mean, median, range)",
    Variable == "Female"  ~ "Female (%)"
  ))

darija_summary <- d %>%
  count(darija) %>%
  mutate(Summary  = paste0(round(n / sum(n) * 100, 1), "% (", n, ")"),
         Variable = paste0("Darija: ", darija, " (%)")) %>%
  select(Variable, Summary)

group_summary <- d %>%
  count(group) %>%
  mutate(Summary  = paste0(round(n / sum(n) * 100, 1), "% (", n, ")"),
         Variable = paste0("Group: ", group, " (%)")) %>%
  select(Variable, Summary)

full_summary <- bind_rows(numeric_summary, darija_summary, group_summary) %>%
  slice(match(c("Age (mean, median, range)", "Years of Education (mean, median, range)", "Female (%)",
                "Darija: none (%)", "Darija: little (%)", "Darija: fluent (%)",
                "Group: never nomad (%)", "Group: former nomad (%)",
                "Group: part-time nomad (%)", "Group: full-time nomad (%)"), Variable))

## ======================= Descriptives by group =======================##

aggregate(d$yrs_edu, list(d$group2), mean)
aggregate(d$yrs_edu, list(group2 = d$group2),
          function(x) c(n = length(x), mean = mean(x), sd = sd(x), min = min(x), max = max(x)))

p.edu <- ggplot(d, aes(x = group2, y = yrs_edu)) +
  geom_boxplot(outlier.shape = NA, aes(fill = group2), alpha = 0.4, na.rm = TRUE) +
  geom_jitter(aes(color = group2), width = 0.15, size = 2, alpha = 0.6) +
  scale_color_manual(values = c("never nomad" = "#0072B2", "former nomad" = "#FFC800", "nomad" = "#D55E00")) +
  scale_fill_manual(values  = c("never nomad" = "#0072B2", "former nomad" = "#FFC800", "nomad" = "#D55E00")) +
  labs(x = "Group", y = "Years of Education") +
  theme_minimal(base_size = 16) +
  theme(legend.position = "none",
        axis.text.x  = element_text(size = 16), axis.text.y  = element_text(size = 16),
        axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18))

summary_stats <- d %>%
  group_by(group2) %>%
  summarise(
    `N (% of sample)`                    = sprintf("%d (%.1f%%)", n(), n()/nrow(d)*100),
    `Age (mean, SD, range)`              = sprintf("%.1f (%.1f), %d–%d",
                                                    mean(age, na.rm=TRUE), sd(age, na.rm=TRUE),
                                                    min(age, na.rm=TRUE),  max(age, na.rm=TRUE)),
    `% female`                           = sprintf("%.1f%%", mean(sex == "female", na.rm=TRUE)*100),
    `Years of education (mean, SD, range)` = sprintf("%.1f (%.1f), %d–%d",
                                                    mean(yrs_edu, na.rm=TRUE), sd(yrs_edu, na.rm=TRUE),
                                                    min(yrs_edu, na.rm=TRUE),  max(yrs_edu, na.rm=TRUE)),
    .groups = "drop"
  ) %>% rename(Group = group2)

darija_long <- d %>%
  group_by(group2, darija) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group2) %>%
  mutate(`%` = sprintf("%.1f%%", round(n / sum(n) * 100, 1))) %>%
  select(-n) %>%
  mutate(Variable = paste("Fluency in Darija: %", darija)) %>%
  select(group2, Variable, `%`) %>%
  pivot_wider(names_from = group2, values_from = `%`) %>%
  mutate(across(everything(), as.character))

summary_by_group <- bind_rows(
  summary_stats %>% pivot_longer(-Group, names_to = "Variable", values_to = "Value") %>%
    pivot_wider(names_from = Group, values_from = Value),
  darija_long
)

##==================== WATER AND FOOD INSECURITY PLOTS ====================##

likert_levels <- c("never", "rarely", "sometimes", "often", "always")

water_security_cols <- c("indivID", "water_need", "water_interruption", "water_change_plans", "group", "group2")
ws <- d[, water_security_cols]
ws$group2[ws$group2 != "nomad"] <- "villager\n(including former nomad)"
ws$water_need         <- as.numeric(factor(ws$water_need,         levels = likert_levels)) - 1
ws$water_interruption <- as.numeric(factor(ws$water_interruption, levels = likert_levels)) - 1
ws$water_change_plans <- as.numeric(factor(ws$water_change_plans, levels = likert_levels)) - 1
ws$water_insecurity   <- ws$water_need + ws$water_interruption + ws$water_change_plans
ws <- ws[complete.cases(ws[c("water_insecurity", "group2")]), ]

des_water_insecurity <- ggplot(ws, aes(x = water_insecurity, fill = group2)) +
  stat_bin(aes(y = after_stat(density)), position = 'dodge', binwidth = 1, na.rm = TRUE) +
  labs(x = "Water insecurity score", y = "Proportion", fill = "group") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),  axis.text.y  = element_text(size = 16),
        legend.text  = element_text(size = 16),   legend.title = element_text(size = 16)) +
  scale_fill_manual(values = c("nomad" = "#D55E00", "villager\n(including former nomad)" = "#009E73"))
ggsave("figures/des_water_insecurity.png", des_water_insecurity, width = 8, height = 6)

ws %>% group_by(group2) %>%
  summarise(n = n(), mean_score = mean(water_insecurity), sd_score = sd(water_insecurity),
            min_score = min(water_insecurity), max_score = max(water_insecurity))

food_security_cols <- c("indivID", "fewer_meals", "no_food", "hungry_alldaynight", "group", "group2")
fs <- d[, food_security_cols]
fs$group2[fs$group2 != "nomad"] <- "villager\n(including former nomad)"
fs$fewer_meals        <- as.numeric(factor(fs$fewer_meals,        levels = likert_levels)) - 1
fs$no_food            <- as.numeric(factor(fs$no_food,            levels = likert_levels)) - 1
fs$hungry_alldaynight <- as.numeric(factor(fs$hungry_alldaynight, levels = likert_levels)) - 1
fs$food_insecurity    <- fs$fewer_meals + fs$no_food + fs$hungry_alldaynight
fs <- fs[complete.cases(fs[c("food_insecurity", "group2")]), ]
fs$food_insecurity <- as.numeric(fs$food_insecurity)

des_food_insecurity <- ggplot(fs, aes(x = food_insecurity, fill = group2)) +
  stat_bin(aes(y = after_stat(density)), position = 'dodge', binwidth = 1, na.rm = TRUE) +
  labs(x = "Food insecurity score", y = "Proportion", fill = "group") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18),
        axis.text.x  = element_text(size = 16),  axis.text.y  = element_text(size = 16),
        legend.text  = element_text(size = 16),   legend.title = element_text(size = 16)) +
  scale_fill_manual(values = c("nomad" = "#D55E00", "villager\n(including former nomad)" = "#009E73")) +
  scale_x_continuous(breaks = min(fs$food_insecurity):max(fs$food_insecurity))
ggsave("figures/des_food_insecurity.png", des_food_insecurity, width = 8, height = 6)

fs %>% group_by(group2) %>%
  summarise(n = n(), mean_score = mean(food_insecurity), sd_score = sd(food_insecurity),
            min_score = min(food_insecurity), max_score = max(food_insecurity))

##==================== SEDENTARISATION DESCRIPTIVE PLOTS ====================##

labels_why_settle <- c("Physically tired or ill", "Children's schooling", "Livestock loss",
                       "Marriage to villager", "Conflict", "Joined family",
                       "Defense of material wealth", "Overcrowding in tent")
search_terms_settle <- c("tired", "school", "livestock", "marriage", "conflict", "family", "defend", "overcrowding")
freq_why_settle_plot <- sapply(search_terms_settle, function(x) sum(grepl(x, nom$why_settle_cats)))

df_why_settle <- data.frame(labels_why_settle, freq_why_settle = freq_why_settle_plot)
df_why_settle$labels_why_settle <- factor(
  df_why_settle$labels_why_settle,
  levels = df_why_settle$labels_why_settle[order(df_why_settle$freq_why_settle, decreasing = FALSE)]
)

des_why_settle <- ggplot(df_why_settle, aes(x = labels_why_settle, y = freq_why_settle)) +
  geom_bar(stat = "identity", show.legend = FALSE, fill = "#FFC800") +
  coord_flip() +
  labs(x = "Reason for settling", y = "Frequency of mentions") +
  theme_bw(base_size = 16) +
  theme(axis.text.y = element_text(size = 16), axis.text.x = element_text(size = 16)) +
  ylim(0, max(df_why_settle$freq_why_settle) + 3)
ggsave("figures/des_why_settle.png", des_why_settle, width = 10, height = 6)

des_age_settled <- ggplot(nom, aes(x = age_settled)) +
  geom_histogram(binwidth = 2, fill = "#FFC800", color = "white") +
  coord_cartesian(xlim = c(8, 60)) +
  labs(x = "Age when settled", y = "Frequency") +
  theme_minimal(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"))
ggsave("figures/des_age_settled.png", des_age_settled, width = 8, height = 6)

df_why_village <- as.data.frame(table(nom$settle_where_why_cat))
des_why_tz <- ggplot(df_why_village, aes(x = reorder(Var1, Freq), y = Freq)) +
  geom_col(fill = "#FFC800") +
  labs(x = "Why choose the study village to settle?", y = "Frequency of mentions") +
  theme_minimal(base_size = 16) +
  coord_flip()
ggsave("figures/des_why_tz.png", des_why_tz, width = 10, height = 6)

labels_why_want <- c("nomadism is physically tiring",
                     "declining livestock/drought/increasing expenses",
                     "protection from the elements/danger to small children",
                     "children's schooling",
                     "access to infrastructure (water, electricity, markets)",
                     "engaged to villager",
                     "boredom and solitude of nomadic life")
df_reasons_want_settling <- data.frame(labels_why_want, freq_why_want)

des_reasons_want_settling <- ggplot(df_reasons_want_settling, aes(x = reorder(labels_why_want, freq_why_want), y = freq_why_want)) +
  geom_bar(stat = "identity", color = "#D55E00", fill = "#D55E00") +
  labs(x = "Reasons for wanting to settle in the future", y = "Frequency of mentions") +
  theme_bw() +
  coord_flip() +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x  = element_text(size = 16, angle = 45, hjust = 1),
        axis.text.y  = element_text(size = 16),
        legend.text  = element_text(size = 16), legend.title = element_text(size = 16))
ggsave("figures/des_reasons_want_settling.png", des_reasons_want_settling, width = 10, height = 6)

tbl_not_want_settle <- subset(nom, nom$want_settle == "no")
table(tbl_not_want_settle$want_settle_why_cat)


##=========  MAP PLOT ==========
source("code/map_plot.R")

### Save
save(list = ls(), file = "model_output/alami_et_al_2026.Rdata")
