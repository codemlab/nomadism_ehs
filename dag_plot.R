library(ggplot2)
library(grid)  # for arrow()

# ---- Box dimensions ----

make_var <- function(name, x, y, box_width = 1.8, box_height = 0.5) {
  data.frame(
    name, x, y,
    xmin = x - box_width/2,
    xmax = x + box_width/2,
    ymin = y - box_height/2,
    ymax = y + box_height/2
  )
}

# ---- Variables (with Mediators stacked vertically) ----
vars <- rbind(
  make_var("Former nomad (1/0)", -4, 0),
  make_var("Age", -1, 2),
  make_var("Sex",  1, 2),
  make_var("Time in community",      0, -1.5),
  make_var("Darija speaker (1/0)",   0, -2.3),
  make_var("Education",              0, -3.1),
  make_var("Unobserved causal pathway", 0, 0, box_width = 2.2),
  make_var("Income", 4,  1),
  make_var("Wealth", 4, -1)
)

# ---- Group boxes (Mediators label goes below) ----
groups <- rbind(
  make_var("Confounds", 0, 2, box_width = 5, box_height = 1),
  make_var("Mediators", 0, -2.3, box_width = 2, box_height = 2.2),
  make_var("Outcomes", 4, 0, box_width = 2, box_height = 3)
)

# ---- Arrow helper: connect midpoints of inner-facing edges ----
arrow_between_edge <- function(from_df, to_df, from_name, to_name,
  from_edge = "inner", to_edge = "inner") {
f <- from_df[from_df$name == from_name, ]
t <- to_df[to_df$name == to_name, ]
stopifnot(nrow(f) == 1, nrow(t) == 1)

edge_coords <- function(box, edge, other_box) {
if (edge == "inner") {
if (other_box$x > box$x) return(c(box$xmax, (box$ymin + box$ymax) / 2)) # right
if (other_box$x < box$x) return(c(box$xmin, (box$ymin + box$ymax) / 2)) # left
if (other_box$y > box$y) return(c((box$xmin + box$xmax) / 2, box$ymax)) # top
return(c((box$xmin + box$xmax) / 2, box$ymin))                          # bottom
}
switch(edge,
left   = c(box$xmin, (box$ymin + box$ymax) / 2),
right  = c(box$xmax, (box$ymin + box$ymax) / 2),
top    = c((box$xmin + box$xmax) / 2, box$ymax),
bottom = c((box$xmin + box$xmax) / 2, box$ymin),
stop("Invalid edge name"))
}

from_xy <- edge_coords(f, from_edge, t)
to_xy   <- edge_coords(t, to_edge, f)

data.frame(x = from_xy[1], y = from_xy[2],
xend = to_xy[1], yend = to_xy[2])
}


# ---- Build arrow data ----
arrows_straight <- rbind(
  arrow_between_edge(vars, groups, "Former nomad (1/0)", "Mediators"),
  arrow_between_edge(vars, vars,   "Former nomad (1/0)", "Unobserved causal pathway"),
  arrow_between_edge(groups, vars, "Confounds", "Former nomad (1/0)", "left","top"),
  arrow_between_edge(groups, groups,"Confounds","Outcomes", "right", "top"),
  arrow_between_edge(groups, vars, "Confounds", "Unobserved causal pathway"),
  arrow_between_edge(groups, groups,"Mediators","Outcomes", "right", "bottom"),
  arrow_between_edge(vars, groups, "Unobserved causal pathway", "Outcomes"),
  arrow_between_edge(vars, vars,   "Income", "Wealth")
)

# The one curved arrow: Confounds -> Mediators (route around Unobserved)
arrows_curved <- arrow_between_edge(groups, groups, "Confounds", "Mediators")

# ---- Plot ----
ggplot() +
  # Grouping boxes
  geom_rect(data = groups,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", linewidth = 1) +
  # Group labels (Confounds, Outcomes above; Mediators below)
  geom_text(data = subset(groups, name != "Mediators"),
            aes(x = x, y = ymax + 0.35, label = name),
            fontface = "bold", size = 6, vjust = 0) +
  geom_text(data = subset(groups, name == "Mediators"),
            aes(x = x, y = ymin - 0.35, label = name),
            fontface = "bold", size = 6, vjust = 1) +

  # Variable boxes + labels
  geom_rect(data = vars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black") +
  geom_text(data = vars, aes(x = x, y = y, label = name), size = 4) +

  # Straight arrows
  geom_segment(data = arrows_straight,
               aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.6, "cm"), type = "closed"),
               linewidth = 0.7) +

  # Curved arrow (single curvature value for the whole layer)
  geom_curve(data = arrows_curved,
             aes(x = x, y = y, xend = xend, yend = yend),
             curvature = -0.9,  # bend leftwards around the Unobserved box
             arrow = arrow(length = unit(0.6, "cm"), type = "closed"),
             linewidth = 0.7) +

  coord_fixed() +
  theme_void()

ggsave("manuscripts/nomadism_ehs/figures/dag.png", width = 12, height = 6)
