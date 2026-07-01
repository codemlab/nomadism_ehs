create_network <- function(output = c("network","graph", "strand","latent"), data, census, rmat= rel_mat,in_question = NA, out_question = NA, directed = F, minage = 18, multiple = F, impute_missing = T,...){
output = match.arg(output)


    data |> filter(question %in% in_question) |>
      mutate(from=alterID,to=egoID)|>
      select(from,to) -> df_in
    data |> filter(question %in% out_question) |>
      mutate(from=egoID,to=alterID) |>
      select(from,to) -> df_out
  
  d_full <- rbind(df_in,df_out)

cen <- census

nodes <- unique(c(d_full$from,d_full$to,cen$indivID[which(cen$age>=minage & !is.na(cen$householdID))]))
nodes <- nodes[which(cen$interviewed_networks[match(nodes,cen$indivID)]==1)]

d_in <- df_in[df_in$from %in% nodes & df_in$to %in% nodes,]
d_out <- df_out[df_out$from %in% nodes & df_out$to %in% nodes,]
d <- d_full[d_full$from %in% nodes & d_full$to %in% nodes,]

bone <- cen$bone[match(nodes,cen$indivID)]
bone[is.na(bone)] <- 'other'
age <- cen$age[match(nodes, cen$indivID)]
sex <- cen$sex[match(nodes, cen$indivID)]
hid <- cen$householdID[match(nodes, cen$indivID)]
former_nomad <- as.character(ifelse(cen$group[match(nodes, cen$indivID)]=="former nomad", "Yes","No"))
edu <- cen$yrs_edu[match(nodes, cen$indivID)]
edu[is.na(edu)] <- mean(edu, na.rm=T)

income <- cen$income[match(nodes, cen$indivID)]  
income_total <- aggregate(cen$income,list(cen$householdID),mean,na.rm=T)
inc_total <- income_total$x[match(hid,income_total$Group.1)]
if(impute_missing == T) inc_total[is.na(inc_total)]<- mean(inc_total,na.rm=T)
ln_income_total <- log(inc_total+1)/10 
ln_income <-log(income+1)/10

wealth_total <- aggregate(cen$wealth, list(cen$householdID), mean, na.rm=T)
wealth_total <- wealth_total$x[match(hid, wealth_total$Group.1)]
if(impute_missing == T) wealth_total[is.na(wealth_total)]<- mean(wealth_total,na.rm=T)
ln_wealth <-log(wealth_total+1)/10 

loans <- cen$loans[match(nodes, cen$indivID)]
if(impute_missing == T) loans[is.na(loans)] <- "no"
  
externals <- sapply(nodes,function(x) length(c(which(!d_full$from %in% cen$indivID & d_full$to==x),which(!d_full$to %in% cen$indivID & d_full$from==x))))
age[is.na(age)]<- mean(age,na.rm=T)
  
if (!multiple){
  df <- d[!duplicated(paste(d$from,d$to)),]
  }
    
relmat <- as.matrix(rmat[nodes,nodes])  


if (output %in% c("strand","latent")){
    
  
  if (output == "strand"){
    outcome <-  matrix(0, nrow = length(nodes), length(nodes))
    rownames(outcome) <- colnames(outcome) <- nodes
    for (i in 1:nrow(d)) {
      ego_id <- as.character(d$from[i])
      alter_id <- as.character(d$to[i])
      outcome[ego_id, alter_id] <- 1
    }
    outcome <- list(outcome = outcome)
  }

  if (output == "latent"){
    TransferIn <- TransferOut <- matrix(0, nrow = length(nodes), length(nodes))
    rownames(TransferIn) <- colnames(TransferIn) <- rownames(TransferOut) <- colnames(TransferOut) <- nodes
    for (i in 1:nrow(d_in)) {
      ego_in <- as.character(d_in$to[i])
      alter_in <- as.character(d_in$from[i])
      TransferIn[ego_in, alter_in] <- 1
    }
    for (i in 1:nrow(d_out)) {
      ego_out <- as.character(d_out$to[i])
      alter_out <- as.character(d_out$from[i])
      TransferOut[ego_out, alter_out] <- 1
    }
    outcome <- list(TransferIn = TransferIn , TransferOut = TransferOut)
    }
  
  dyad <- list(relatedness=relmat)
  groups <- data.frame(sex = factor(sex),
                       bone = factor(bone),
                       hid = factor(hid),
                       loans = factor(loans),
                      former_nomad = factor(former_nomad))
  indiv <- data.frame(
        age = STRAND::standardize(age),
        sex = sex,
        former_nomad = factor(former_nomad),
        externals = STRAND::standardize(externals),
        income_total = STRAND::standardize(inc_total),
        ln_income_total = STRAND::standardize(ln_income_total),
        wealth = STRAND::standardize(wealth_total),
        ln_wealth = STRAND::standardize(ln_wealth),
        edu = STRAND::standardize(edu))
  
rownames(indiv) <- rownames(groups) <- nodes
dat = make_strand_data(outcome = outcome,
  block_covariates = groups, 
  individual_covariates = indiv, 
  dyadic_covariates = dyad, ...)

return(dat)
}  
  
  
if(output %in% c( "network", "graph")) {

node_df <- create_node_df(n=length(nodes),
bone=bone,
age=scale(age),
sex=sex,
hid=hid,
former_nomad=former_nomad,
externals=scale(externals),
income_total=inc_total, 
ln_income_total=scale(log(inc_total+1)),
wealth=wealth_total,
ln_wealth=scale(log(wealth_total+1)),
loans=loans,
edu=edu,
shape = c('diamond','circle'))
  


edge_df <- create_edge_df(from=(1:length(nodes))[match(df$from,nodes)],to=(1:length(nodes))[match(df$to,nodes)],use_labels=F, directed = directed)
net <- as.network(df,multiple=multiple)
net_names <- net %v% "vertex.names"
add.vertices(net, sum(!nodes %in% net_names))

set.vertex.attribute(net, 'vertex.names', c(net_names,nodes[!nodes %in% net_names]))
set.vertex.attribute(net, 'Bone', bone[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Former nomad', former_nomad[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Age', age[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Sex', sex[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Externals', externals[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'HouseholdID', hid[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Income (ln)', ln_income_total[match(net %v% "vertex.names",nodes)]) 
set.vertex.attribute(net, 'Wealth (ln)', ln_wealth[match(net %v% "vertex.names",nodes)]) 
set.vertex.attribute(net, 'Loans', loans[match(net %v% "vertex.names",nodes)])
set.vertex.attribute(net, 'Yrs schooling', edu[match(net %v% "vertex.names",nodes)])
set.edge.value(net, 'Relatedness', relmat)

## graph 
  if (output == "graph"){
colour_func <- scales::gradient_n_pal(colours=c('white', 'LightBlue', 'DarkBlue'),values=c(0,6,11.5))

create_graph(nodes_df=node_df,edges_df=edge_df,directed=directed) %>%
   join_node_attrs(df = get_eigen_centrality(.)) %>%
    set_node_attrs(
    node_attr = fillcolor,
    values = colour_func(log(inc_total+1)),
  ) %>%
  set_node_attrs(
    node_attr = label,
    values = "",
  ) %>%
  set_node_attrs(
    node_attr = shape,
    values = c('square','circle')[match(former_nomad,unique(former_nomad))],
      )%>%
   rescale_node_attrs(
    node_attr_from = eigen_centrality,
    to_lower_bound = 0.1,
    to_upper_bound = 1.5,
    node_attr_to = height
  ) %>%
   rescale_node_attrs(
    node_attr_from = eigen_centrality,
    to_lower_bound = 0.1,
    to_upper_bound = 1.5,
    node_attr_to = width
  )-> graph
}

if(output == "graph") return  (graph) 
if(output == "network") return (list(net = net, relmat = relmat))
}
}

