library(brms)
library(tidyverse)
library(posterior)
library(ggdist)
library(ggridges)
library(bayestestR)

library(ggplot2)
library(ggdist)
library(dplyr)
library(tibble)

# Function to extract and label parameters
extract_draws <- function(fit, model_name, outcome = "Income") {
  draws <- as_draws_df(fit)
  priors <- fit$prior

  # Define the parameters of interest in order
  parameters <- paste0("b_", c("former_nomadYes", "hu_former_nomadYes",
                   "darija", "hu_darija",
                   "yrs_edu_std", "hu_yrs_edu_std"))
  if(outcome == "Wealth")  parameters <- c(parameters, paste0("b_", c("time_com_std", "income_std")))

  draws_long <- draws |>
    select(all_of(parameters[parameters %in% names(draws)])) |>
    pivot_longer(everything(), names_to = "parameter", values_to = "value") |>
    mutate(
      model = model_name,
      effect = case_when(
        outcome == "Wealth" ~ "Log Wealth",
        str_detect(parameter, "^b_hu_") ~ "Log-odds of no income",
        TRUE ~ "Log Income (if >0)"
      ),
      clean_param = parameter |>
        str_remove("^b_hu_") |>
        str_remove("^b_") |>
        recode(
          "former_nomadYes" = "Former nomad (1/0)",
          "yrs_edu_std" = "Years of education",
          "darija" = "Darija speaker (1/0)",
          "timeincom" = "Time in community",
          "income" = "Income (Ln)",
          "income_std" = "Income (Ln)",
          "time_com_std" = "Time in community"
        )
    )
}

model_names <- c("Total effect", "Direct effect + mediators")
# Get draws from both models
draws1 <- extract_draws(inc_model_brm, model_names[1])
draws2 <- extract_draws(inc_model_brm_full, model_names[2])

# Combine and assign y-axis row numbers for plotting
all_draws <- bind_rows(draws1, draws2) %>%
  mutate(
    y_group = factor(clean_param, levels = c(
      "Darija speaker (1/0)",
      "Years of education",
      "Former nomad (1/0)"
    )),
    row = case_when(
      parameter == "b_former_nomadYes" ~ 1,
      parameter == "b_hu_former_nomadYes" ~ 2,
      parameter == "b_yrs_edu_std" ~ 4,
      parameter == "b_hu_yrs_edu_std" ~ 5,
      parameter == "b_darija" ~ 7,
      parameter == "b_hu_darija" ~ 8
    ),
    row = factor(row, levels = 8:1),
    model = factor(model, levels = model_names)
  ) 

# Plot  
label_df <- data.frame(
  x = -5,
  y = c(7.5, 4.5, 1.5),
  label = c("Former nomad (1/0)", "Years of education", "Darija speaker (1/0)"),
  model = factor(model_names[1], levels = model_names)
)

ggplot(all_draws, aes(x = value, y = row, fill = effect)) +
  stat_halfeye(
    .width = 0.95,
    point_interval = median_qi,
    slab_size = 0.7,
    slab_alpha = 0.8,
    normalize = "groups"
  ) +
    facet_grid(cols = vars(model))+
  scale_y_discrete(drop = F,
    labels = ""
  ) +
  scale_fill_manual(
    name = NULL,
    values = c("Log Income (if >0)" = "#2166ac", "Log-odds of no income" = "#b2182b")
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.spacing.x = unit(1.5, "cm") ,
    strip.text.y.left = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 110) ,
    panel.grid.major.y = element_blank(),
    legend.background = element_rect(fill = "white", color = "gray80")
  ) + coord_cartesian(xlim = c(-2.3, 2.3), clip = "off") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.8) + 

     geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 4,
      group = NA,
      hjust = 0
    ) +
  labs(x = "Posterior distribution", y = NULL)

ggsave("/home/edseab/github/codem/manuscripts/nomadism_ehs/figures/income_model.png", width = 10, height = 6, bg = 'white', dpi = 300)

#wealth plot
# Get draws from both models
draws3 <- extract_draws(wealth_model_brm, model_names[1], "Wealth")
draws4 <- extract_draws(wealth_model_brm_full, model_names[2], "Wealth")

# Combine and assign y-axis row numbers for plotting
all_draws <- bind_rows(draws3, draws4) %>%
  mutate(
    y_group = factor(clean_param, levels = c(
      "Darija speaker (1/0)",
      "Years of education",
      "Former nomad (1/0)",
      "Time in community",
      "Income (Ln)"
    )),
    row = case_when(
      parameter == "b_former_nomadYes" ~ 1,
      parameter == "b_yrs_edu_std" ~ 2,
      parameter == "b_darija" ~ 3,
      parameter == "b_time_com_std" ~ 4,
      parameter == "b_income_std" ~ 5
    ),
    row = factor(row, levels = 5:1),
    model = factor(model, levels = model_names)
  ) 

# Plot  
label_df <- data.frame(
  x = -1,
  y = c(5, 4, 3, 2, 1),
  label = c("Former nomad (1/0)", "Years of education", "Darija speaker (1/0)","Time in community", "Income (Ln)"),
  model = factor(model_names[1], levels = model_names)
)

ggplot(all_draws, aes(x = value, y = row, fill = effect)) +
  stat_halfeye(
    .width = 0.95,
    point_interval = median_qi,
    slab_size = 0.7,
    slab_alpha = 0.8,
    normalize = "groups"
  ) +
    facet_grid(cols = vars(model))+
  scale_y_discrete(drop = F,
    labels = ""
  ) +
    scale_fill_manual(
      name = NULL,
      values = "#2166ac"
      ) +
  theme_minimal(base_size = 16) +
  theme(
    strip.text.y.left = element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 10, r = 5, b = 10, l = 130) ,
    panel.grid.major.y = element_blank(),
    legend.background = element_rect(fill = "white", color = "gray80")
  ) + coord_cartesian(xlim = c(-0.5, 0.3), clip = "off") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.6, alpha = 0.8) + 

     geom_text(
      data = label_df,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      size = 4,
      group = NA,
      hjust = 0
    ) +
  labs(x = "Posterior distribution", y = NULL)

ggsave("/home/edseab/github/codem/manuscripts/nomadism_ehs/figures/wealth_model.png", width = 10, height = 6, bg = 'white', dpi = 300)