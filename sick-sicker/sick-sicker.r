####################
# Setup
####################
library(tidyverse)
library(DiagrammeR)
library(here)
library(glue)
library(Matrix)
library(expm)
library(ggrepel)    # For plotting
library(ellipse)    # For plotting
library(scales)     # For dollar signs and commas
library(randtoolbox)
library(igraph)
library(tidygraph)
library(ggraph)
library(progress)
options("scipen" = 100, "digits" = 5)

source("sick-sicker/darth_functions.r")
source("sick-sicker/darth_functions2.r")
sub_in_expression <- function(x) map2_chr(names(dist_lut),dist_lut, ~ifelse(grepl(glue('^{.x}'),x),gsub(glue('^{.x}'),.y,x),NA)) %>% na.omit() %>% as.vector()
map2_progress <- function(.x,.y ,.f, ..., .id = NULL) {
  # Source: https://www.jamesatkins.net/post/2020/progress-bar-in-purrrmap_df/
  .f <- purrr::as_mapper(.f, ...)
  pb <- progress::progress_bar$new(total = length(.x), force = TRUE)

  f <- function(...) {
    pb$tick()
    .f(...)
  }
  purrr::map2(.x, .y,f, ..., .id = .id)
}
map_progress <- function(.x, .f, ..., .id = NULL) {
  # Source: https://www.jamesatkins.net/post/2020/progress-bar-in-purrrmap_df/
  .f <- purrr::as_mapper(.f, ...)
  pb <- progress::progress_bar$new(total = length(.x), force = TRUE)

  f <- function(...) {
    pb$tick()
    .f(...)
  }
  purrr::map(.x, f, ..., .id = .id)
}

#############################################
# Read in model structure and parameter info
#############################################

input_file <- normalizePath(here("sick-sicker/sick-sicker.xlsx")); input_file
params_raw <- readxl::read_xlsx(input_file,sheet="parameters"); params_raw
rates_raw <- readxl::read_xlsx(input_file,sheet = "transition_rates");  rates_raw
probs_raw <- readxl::read_xlsx(input_file,sheet="transition_probabilities"); probs_raw
payoffs_raw <- readxl::read_xlsx(input_file, sheet = "payoffs"); payoffs_raw

###########################
## Metaparameters
###########################
v_names_states <- unique(probs_raw$from); v_names_states
v_n_states <- length(v_names_states)
v_names_str <- unique(rates_raw$strategy) %>% na.omit() %>% as.vector(); v_names_str
n_strategies <- length(v_names_str)
n_cycles <- params_raw %>% filter(param=="n_cycles") %>% pull(baseline) ; n_cycles
n_sim <- params_raw %>% filter(param == "n_sim") %>% pull(baseline); n_sim

##############################
## Step 1: Process Parameters
##############################

# Parameters with scalar values
dist_lut <- list("gamma" = "qgamma",
                 "lognormal" = "qlnorm",
                 "beta" = "qbeta",
                 "normal" = "qnorm")
params_sc <-
  params_raw  %>%
  filter(is.na(formula)) %>%
  select(param,baseline) %>%
  na.omit() %>%
  deframe() %>%
  as.list(); params_sc

params_psa_ <-
  params_raw  %>%
  filter(!is.na(distribution)) %>%
  select(param,distribution) %>%
  na.omit() %>%
  deframe() %>%
  as.list(); params_psa_

psa_ <- halton(n=params_sc$n_sim, dim = length(params_psa_))
colnames(psa_) <- unlist(params_psa_)
colnames(psa_) <- lapply(colnames(psa_),function(x) sub_in_expression(x))
params_psa <-
  colnames(psa_) %>%
  map(~(glue("{gsub(')$','',.x)}, p={psa_[,.x]})")))  %>%
  set_names(names(params_psa_)) %>%
  bind_cols() %>% #glimpse()
  rowwise() %>%
  mutate_all(~eval(str2expression(.))) %>%
  ungroup() %>%
  as.list()

params_sc_psa <-
  params_sc[setdiff(names(params_sc),names(params_psa))] %>%
  map(~(rep(.x,params_sc$n_sim)))

params_psa <-
  append(params_psa, params_sc_psa)

# Parameters that are functions of other parameters
params_fmla <-
  params_raw %>%
  filter(!is.na(formula)) %>%
  select(param,formula) %>%
  na.omit() %>%
  deframe() %>%
  as.list(); params_fmla

# Execute the formulas and create the master params object
params <-
  append(params_sc,map(params_fmla,~(eval(str2expression(.x),envir=params_sc))))

params_psa <-
  append(params_psa,map(params_fmla,~(eval(str2expression(.x),envir=params_psa))))

  ### Discount weight for costs and effects ----
  v_dwc  <- 1 / ((1 + (params$d_e * params$cycle_length)) ^ (0:params$n_cycles))
  v_dwe  <- 1 / ((1 + (params$d_c * params$cycle_length)) ^ (0:params$n_cycles))
  ## Within-cycle correction (WCC) using Simpson's 1/3 rule ----
  v_wcc <- gen_wcc(n_cycles = params$n_cycles,  # Function included in "R/Functions.R". The latest version can be found in `darthtools` package
                   method = "Simpson1/3") # vector of wcc

################################
## Step 2: Transition Matrices
################################

  .x = "SoC"

  get_transition_matrices <- function(strategies, params, jumpover = NULL) {
    m_P <-
      strategies %>% map(~({
        n <- params %>% map_dbl(~(length(.x))) %>% max()
        m_Q <- P_ <- list()

        m_Q_ <-
          rates_raw %>%
          filter(strategy == .x) %>%
          filter(from %in% all_of(v_names_states))   %>%
          select_at(vars(from,v_names_states)) %>%
          column_to_rownames(var = "from") %>%
          as.matrix(); m_Q_

        if (!is.null(jumpover)) {
          for (jj in names(jumpover)) {
            ff <- jumpover[[jj]][[1]]; ff
            tt <- jumpover[[jj]][[2]]; tt
            tmp_ <- rbind(cbind(m_Q_,rep("0",nrow(m_Q_))),rep("0",ncol(m_Q_)+1)); tmp_
            colnames(tmp_) <- c(colnames(m_Q_),jj)
            rownames(tmp_) <- c(rownames(m_Q_),jj)
            tmp_[ff,jj] <- m_Q_[ff,tt]; tmp_
            m_Q_ <- tmp_
          }
        }

        for (s in 1:n) {
          params_ <- params %>% transpose() %>% pluck(s); params_
          m_Q[[s]] <- apply(m_Q_,c(1,2),function(x) eval(str2expression(x), envir=params_)); m_Q[[s]]
          P_[[s]] <-   expm(m_Q[[s]]) %>%
            as.matrix(); P_[[s]]


        }

        P_

      })) %>%
      set_names(strategies)
    m_P %>% transpose()
  }

#!!!!!!!!
#!!!!!!!! Removing jumpover version for now until we get it figured out
#!!!!!!!!

  # m_P <- get_transition_matrices(strategies = v_names_str, params = params, jumpover = list("acc_Sick" = list("Healthy","Sick"))) %>% pluck(1)
  # m_P_psa <-  get_transition_matrices(strategies = v_names_str, params = params_psa,jumpover = list("acc_Sick" = list("Healthy","Sick")))

  m_P <- get_transition_matrices(strategies = v_names_str, params = params, jumpover = NULL) %>% pluck(1)
  m_P_psa <-  get_transition_matrices(strategies = v_names_str, params = params_psa,jumpover = NULL)

# Visualize Model Structure, including implied jumpover states

  #brew install graphviz

  m_P_raw <- get_transition_matrices(strategies = v_names_str, params = params) %>% pluck(1)
  g <- from_adj_matrix(m_P_raw[[1]]!=0,mode="directed")
  #render_graph(g, layout = "nicely")

  g_ <- probs_raw %>%
    column_to_rownames("from")  %>%
    as.matrix()
  g_[g_!="0"] <- 1; g_
  g_ <- from_adj_matrix(g_, mode="directed")
  #render_graph(g_, layout = "nicely")

  jumpover <- get_edge_df(g) %>% anti_join(get_edge_df(g_),c("from","to")) %>%
    mutate(jumpover=TRUE) %>%
    select(from,to,jumpover)

  graph_from_adjacency_matrix(m_P[[1]]>0, mode= "directed") %>%
    tidygraph::as_tbl_graph() %>%
    activate(edges) %>%
    left_join(jumpover,c("from","to")) %>%
    mutate(jumpover = ifelse(!is.na(jumpover),1,0)) %>%
    ggraph( 'matrix') +
    #geom_edge_link(aes(colour = factor(jumpover)),edge_width=1.5) +
    #geom_node_point(size=10) +
    geom_node_label(aes(label = name),size=4)+
    theme_graph() +
    scale_edge_color_manual(name = "Transition Type",values = c("black","red"),labels = c("Direct","Jumpover")) +
    geom_edge_arc(aes(colour = factor(jumpover),edge_linetype = factor(jumpover)), arrow = arrow(length = unit(4, 'mm')),
                  start_cap = circle(8, 'mm'),
                  end_cap = circle(8, 'mm'), edge_width = .75) +
    theme(legend.position = "bottom") +
    scale_edge_linetype_manual(name = "Transition Type",values = c(1,5),labels = c("Direct","Jumpover"))

################################
## Step 3: Prepare Payoffs
################################

get_payoffs <- function(params, jumpover = NULL) {
  n <- params %>% map_dbl(~(length(.x))) %>% max()

  out <-
    map2_progress(rep(list(payoffs_raw),n),transpose(params),~({
      tmp_ <- .y

      payoff_ <-
      .x %>%
        mutate(costs = map_dbl(costs,~(eval(str2expression(.x),envir=tmp_))),
               qalys = map_dbl(qalys,~(eval(str2expression(.x),envir=tmp_))))

      if (!is.null(jumpover)) {
        for (jj in names(jumpover)) {
          ff <- jumpover[[jj]][[1]]; ff
          tt <- jumpover[[jj]][[2]]; tt

          payoff_j <- payoff_ %>% filter(label == tt) %>%
            mutate(label = jj,
                   description = jj)

          payoff_ <- payoff_ %>% bind_rows(payoff_j)

          v_names_states_ <- c(v_names_states,jj)
        }
      } else {
        v_names_states_ <- v_names_states
      }
      payoff_ %>%
        group_by(strategy) %>%
        nest() %>%
        deframe() %>%
        map(~(
          .x %>%
            select(label,costs,qalys)  %>%
            gather(type,value,-label) %>%
            group_by(type) %>%
            nest() %>%
            deframe() %>%
            map(~(.x %>%
                    spread(label,value) %>%
                    select_at(vars(v_names_states_)) %>%
                    as.matrix()))
        ))
    }))
  if (n==1) out <- out %>% pluck(1)
  return(out)
}

  #!!!!!!!!
  #!!!!!!!! Removing jumpover version for now until we get it figured out
  #!!!!!!!!

# payoffs <- get_payoffs(params,jumpover = list("acc_Sick" = list("Healthy","Sick")))
# payoffs_psa <- get_payoffs(params_psa,jumpover = list("acc_Sick" = list("Healthy","Sick")))

payoffs <- get_payoffs(params,jumpover = NULL)
payoffs_psa <- get_payoffs(params_psa,jumpover = NULL)

################################
## Step 4: Markov Traces
################################

# Construct shell markov trace
l_m_M <-
  m_P %>%
  map(~({
    tmp <- .x[1,] %>% as.matrix() %>% t()
    tmp <- matrix(0,ncol=(ncol(tmp)),nrow=(n_cycles+1),dimnames = list(paste0(0:n_cycles),colnames(.x)))
    tmp[1,1] <- 1
    tmp
  }))

l_m_M_psa <- rep(list(l_m_M),n_sim)

for (t in 1:n_cycles) {
  res <-
    map2(l_m_M,m_P,~({
      .x[paste0(t-1),] %*% .y
    }))
  l_m_M <-
    map2(l_m_M,res,~({
      .x[paste0(t),] <- .y
      .x
    }))
}

plot_trace(l_m_M[[2]])

# TK This could also be done much faster

for (k in 1:n_sim) {
  for (t in 1:n_cycles) {
    res <-
      map2(l_m_M_psa[[k]],m_P_psa[[k]],~({
        .x[paste0(t-1),] %*% .y
      }))
    l_m_M_psa[[k]] <-
      map2(l_m_M_psa[[k]],res,~({
        .x[paste0(t),] <- .y
        .x
      }))
  }
}

p <- plot_trace(l_m_M[[2]])
p_ <- ggplot_build(p)[["data"]][[1]] %>%
  mutate(iteration=0)

p_psa <- l_m_M_psa %>%
  map_progress(~({
    tmp <- plot_trace(.x[[2]])
    ggplot_build(tmp)[["data"]][[1]]
  }
  )) %>%
  bind_rows(.id = "iteration") %>%
  mutate(iteration = factor(iteration))

ggplot(data = p_, aes(x=x,y=y, group = factor(iteration))) + geom_line() +
  theme_light() +
  facet_wrap(~group) +
  geom_line(data = p_psa, alpha =0.01 )

################################
## Step 5: Calculate Costs and QALYs
################################

v_tot_qaly <-
  map2(l_m_M, payoffs,~({
    v_u_str <- .y$qalys %>% as.vector()
    t(.x %*% v_u_str) %*% (v_dwe * v_wcc)
})) %>%
  unlist()

v_tot_cost <-
  map2(l_m_M, payoffs,~({
    v_c_str <- .y$costs %>% as.vector()
    t(.x %*% v_c_str) %*% (v_dwe * v_wcc)
  })) %>%
  unlist()

df_cea <- calculate_icers(cost       = v_tot_cost,
                          effect     = v_tot_qaly,
                          strategies = v_names_str)
df_cea

v_tot_qaly_psa <-
  map2(l_m_M_psa, payoffs_psa, ~({
    map2(.x, .y,~({
      v_u_str <- .y$qalys %>% as.vector()
      t(.x %*% v_u_str) %*% (v_dwe * v_wcc)
    })) %>%
      unlist()
  }))

v_tot_cost_psa <-
  map2(l_m_M_psa, payoffs_psa, ~({
    map2(.x, .y,~({
      v_c_str <- .y$costs %>% as.vector()
      t(.x %*% v_c_str) %*% (v_dwe * v_wcc)
    })) %>%
      unlist()
  }))

df_cea_psa <-
  map2_progress(v_tot_cost_psa, v_tot_qaly_psa,~({
    calculate_icers(cost       = .x,
                    effect     = .y,
                    strategies = v_names_str)}))
#df_cea_psa

## CEA table in proper format ----
table_cea <- format_table_cea(df_cea) # Function included in "R/Functions.R"; depends on the `scales` package
table_cea

## CEA frontier -----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot(df_cea, label = "all", txtsize = 16) +
  expand_limits(x = max(table_cea$QALYs) + 0.1) +
  theme(legend.position = c(0.8, 0.2))

##################
# PSA
##################
df_c <- v_tot_cost_psa %>% bind_rows %>% data.frame()
df_e <- v_tot_qaly_psa %>% bind_rows() %>% data.frame()
df_psa_input <- params_psa %>% bind_cols() %>% data.frame()
l_psa <- make_psa_obj(cost          = df_c,
                      effectiveness = df_e,
                      parameters    = df_psa_input,
                      strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness)<- v_names_str
colnames(l_psa$cost)<- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 200000, by = 5000)

### Cost-Effectiveness Scatter plot ----
#* Function included in "R/Functions.R"; depends on `tidyr` and `ellipse` packages.
#* The latest version can be found in `dampack` package
plot.psa(l_psa) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")


### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
df_out_ce_psa <- summary.psa(l_psa)

#* Function included in "R/Functions.R"; depends on the `dplyr` package
#* The latest version can be found in `dampack` package
df_cea_psa <- calculate_icers(cost       = df_out_ce_psa$meanCost,
                              effect     = df_out_ce_psa$meanEffect,
                              strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
#* Function included in "R/Functions.R"; depends on the `ggplot2`  and `ggrepel` packages.
#* The latest version can be found in `dampack` package
plot.icers(df_cea_psa)

## Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
#* Functions included in "R/Functions.R". The latest versions can be found in `dampack` package
ceac_obj <- ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary.ceac(ceac_obj)
#* CEAC & CEAF plot
plot.ceac(ceac_obj) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.82, 0.5))

## Expected Loss Curves (ELCs) ----
#* Function included in "R/Functions.R".The latest version can be found in `dampack` package
elc_obj <- calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj
#* ELC plot
plot.exp_loss(elc_obj, log_y = FALSE,
              txtsize = 16, xlim = c(0, NA), n_x_ticks = 14,
              col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)",
                     breaks = number_ticks(10),
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7))

## Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
plot.evpi(evpi, effect_units = "QALY")


# https://www.youtube.com/watch?v=6L5LPkzvbtI
# http://www.kibohut.com/download/index.php
# https://clas.ucdenver.edu/marcelo-perraillon/teaching/cost-effectiveness-analysis-hsmp-6609
# https://www.youtube.com/watch?v=r_-PLJsN8ng


