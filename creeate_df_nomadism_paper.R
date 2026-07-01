library(kinship2)
library(readxl)
library(tidyverse)

d <- read.csv('clean_data/censuses/census_tz_2023.csv')
nom <- read.csv("clean_data/nomadism_20251022.csv")
econ <- read.csv('clean_data/econ_tz_20251022.csv')
reg <- read.csv('clean_data/registers/register.csv')
dem <- read_excel("raw_data/demography_v20230416_-_all_versions_-_False_-_2023-08-24-15-20-34.xlsx", sheet=1)
alt <- read.csv('clean_data/sn_alters_2024_03_08.csv')


# Demography file pid
dem$pid <- d$interviewed_pid[match(dem$'_id', d$iid)]

# AGE
## Calculate age as of June 1st, 2024
ref_date <- as.Date("2023-06-01")
d$age_june2023 <- as.numeric(floor((ref_date - as.Date(d$dob)) / 365.25))
d$age_june2023[is.na(d$age_june2023)] <- d$age[is.na(d$age_june2023)]
econ$age <- d$age_june2023[match(econ$pid,d$pid)]
d$age <- d$age_june2023

# Income
d$income <- econ$total_net_income_mad[match(d$pid,econ$pid)] + econ$expenses_mad[match(d$pid,econ$pid)]
d$wealth <- econ$wealth[match(d$pid,econ$pid)]
d$wealth_w_savings <- econ$wealth_w_savings[match(d$pid,econ$pid)]
d$loans <- econ$loans[match(d$pid,econ$pid)]
d$net_income <- econ$total_net_income_mad[match(d$pid,econ$pid)]

# Darija fluency
d$darija <- as.numeric(d$darija == "yes")

### Different ways of classifying origin
d$former_nomad <- as.numeric(d$origin=="Nomad") |> factor(labels = c("No", "Yes"))
d$origin <- factor(d$origin, levels = c('Study village', 'Nomad', 'Neighboring village'))
d$group <- ifelse(d$origin == "Nomad", "former nomad", "never nomad")
d$group[d$pid=="WVQ"] <- "part-time nomad"

# Water and food insecurity
d$water_need <- econ$water_security1[match(d$pid,econ$pid)]
d$water_interruption <- econ$water_security2[match(d$pid,econ$pid)]
d$water_change_plans <- econ$water_security3[match(d$pid,econ$pid)]

d$fewer_meals <- econ$fewer_meals[match(d$pid,econ$pid)]
d$no_food <- econ$no_food[match(d$pid,econ$pid)]
d$hungry_alldaynight <- econ$hungry_alldaynight[match(d$pid,econ$pid)]

# Bone
nom$bone <- nom$bone_tribe

# MOROCCAN ARABIC
nom$darija <- ifelse(nom$nomad_darija =="none",0,1)

# SUBSET NOM
nom <- nom[which(nom$lives_where!="other"),] # Remove non-nomad, non TZ folks.
nom$group <- case_when(nom$lives_where=="nomad" ~ "full-time nomad",
                       grepl("and", nom$lives_where) ~ "part-time nomad",
                       nom$lives_where=="tz" ~ "former nomad"
                    )

nom$nomadic <- ifelse(nom$group=="former nomad", 0, 1)

# Nom income and wealth
nom$income <- nom$nomad_total_income_mad
nom$wealth <- nom$nomad_total_wealth_mad
nom$years_in_com <- NA
nom$net_income <- NA

nom$fewer_meals <- nom$nomad_food_insecurity1
nom$no_food <- nom$nomad_food_insecurity2
nom$hungry_alldaynight <- nom$nomad_food_insecurity3

nom$loans <- NA
# IDs
nom$householdID <- nom$hid
d$householdID <- d$hid <- d$HID

nom$indivID <- nom$pid
d$indivID <- d$pid

# Create data frame 
descriptive_cols <- c("indivID","householdID","age","sex","group","bone","yrs_edu","darija","income","net_income","wealth","loans","years_in_com",
                      "water_need", "water_interruption", "water_change_plans",
                     "fewer_meals", "no_food", "hungry_alldaynight", "pid", "hid")
desc_df <- d[d$pid %in% econ$pid,descriptive_cols] |> rbind(nom[!nom$pid %in% d$pid,descriptive_cols])

desc_df<- desc_df[desc_df$pid != "KNZ",] #KNZ has no data across dbs, dropped

# Change ids
desc_df$indivID <- paste0("ID",sample(1:nrow(desc_df)))
unique_hids <- unique(desc_df$householdID)
desc_df$householdID <- paste0("ID",sample(1:length(unique_hids)))[match(desc_df$householdID,unique_hids)]

# interviewed columns
desc_df$interviewed_econ <- as.numeric(desc_df$pid %in% econ$pid)
desc_df$interviewed_networks <- as.numeric(desc_df$pid %in% econ$pid)
desc_df$interviewed_nomadism <- as.numeric(desc_df$pid %in% nom$pid)

write.csv(desc_df, "./manuscripts/nomadism_ehs/manuscript_data/nomadism_paper_main_df.csv", row.names=F)
desc_df_anon <- desc_df[,!colnames(desc_df) %in% c("pid","hid")]


write.csv(desc_df_anon, "./manuscripts/nomadism_ehs/manuscript_data/main_dataset.csv", row.names=F)

# Relationship matrix
ids <- desc_df$pid[desc_df$pid %in% reg$pid]

rmat <- kinship(reg$pid, dadid= reg$father_pid, momid = reg$mother_pid, sex = ifelse(reg$sex == 'male',1,ifelse(reg$sex=='female'),2,NA))
relmat <- as.matrix(rmat[ids,ids]*2)
colnames(relmat) <- rownames(relmat) <- desc_df$indivID[match(colnames(relmat), desc_df$pid)]

write.csv(relmat, "./manuscripts/nomadism_ehs/manuscript_data/relatedness_matrix.csv")

# Alters
others <- unique(c(alt$ego_pid, alt$pid))
others <- others[!others %in% desc_df$pid]
others_new <- paste0("OID", 1:length(others))
all_pids <- c(desc_df$pid, others)
all_newids <- c(desc_df$indivID, others_new)

alt <- alt |> mutate(egoID = all_newids[match(ego_pid,all_pids)], 
                     alterID = all_newids[match(pid,all_pids)]) |>
              select(question, egoID, alterID, in_out) 

write.csv(alt, "./manuscripts/nomadism_ehs/manuscript_data/social_networks.csv", row.names=F)

# Qualitative data
qual <- nom[,c("age_settled", "group", "why_settle_cats","settle_where_why_cat","want_settle","want_settle_why_cat","nomad_relative_need","settled_relative_need","nomad_support","change_support","nomad_respect_settled","nomad_respect_others","nomad_conflicts","settled_conflicts","advnom_cat", "advset_cat", "go_back", "nomad_freedom", "settled_freedom")]
write.csv(qual, "./manuscripts/nomadism_ehs/manuscript_data/nomadism_qual.csv", row.names=F)
