library(signnet)
library("Rcpp")
library(ergm)
library(ggraph)
library(igraph)
library(tidygraph)
library(ggpubr)
library(coda)
library(data.table)
library(ggmcmc)
library(ergm.sign)
library(gridExtra)
library(grid)
library(Hmisc)
library(data.table)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())
source(file = "../other_functions.R")
dca_data = fread("DCAD-v1.0-dyadic.csv")
library(peacesciencer)
download_extdata()

tmp_years = 2000:2010
debugonce(create_dataset_dca)
tmp_data = create_dataset_dca(years = tmp_years,
                              only_biggest_cluster =T,
                              aggregate_years =5,mid_var = "any_mid",
                              alliance_var = "any_cow",dca_data = dca_data)
information = as.data.frame(round(tmp_data$actor_df[,.("Mean GDP" = mean(exp(gdp)),"Variance GDP" = var(exp(gdp)),
                                                       "Mean CINC" = mean(exp(cinc)),"Variance CINC" = var(exp(cinc)),
                                                       "Mean Population" = mean(exp(pop)), "Variance Population" = var(exp(pop)),
                                                       "Mean Polity Score" = mean(polity), "Variance Polity Score" = var(polity)),
                                                    by = .(year)], digits = 3))
# information$`Variance CINC` = "<0.001"
# info_1 = as.matrix(information[,c(1:5)])
# info_2 = as.matrix(information[,c(1, 6:9)])
#
# latex(info_2,label = "tbl:res_cov",first.hline.double= F,
#       cgroupTexCmd = " ",rgroupTexCmd = " ",
#       caption = "Decriptive information about the exogenous cpvariates.")
#
#
aggre_data = tmp_data

tmp_network = aggre_data$target_matrix

alpha = matrix(log(2),nrow = 1,ncol = 1)
alpha = lapply(seq_len(length(tmp_network)), function(X) alpha)
repetition = aggre_data$cov_repetition
repetition_pos =aggre_data$cov_repetition_pos
repetition_neg = aggre_data$cov_repetition_neg
repetition_neut = aggre_data$cov_repetition_neut
contiguity = aggre_data$cov_contiguity
polity = aggre_data$cov_polity
gdp = aggre_data$cov_gdp
distance = aggre_data$cov_distance
gdppc = aggre_data$cov_gdppc
cinc = aggre_data$cov_cinc
common_friends = aggre_data$common_friends
common_enemies = aggre_data$common_enemies


library(parallel)
no_cores <- 5
# Setup cluster
file.remove("outfile.txt")
clust <- makeForkCluster(no_cores, outfile= "outfile.txt")
clusterCall(clust, function() {library(ergm.sign);library(Rcpp);library(RcppArmadillo)})

model_1 = tsergm(formula = tmp_network ~ edges_pos+
                   isolates_pos +
                   isolates_neg +
                   cov_dyad_pos(data = repetition_pos) +
                   cov_dyad_neg(data = repetition_neg) +
                   cov_dyad_pos(data = polity) +
                   cov_dyad_neg(data = polity) +
                   cov_dyad_pos(data = cinc) +
                   cov_dyad_neg(data = cinc) +
                   cov_dyad(data = distance) +
                   cov_dyad_pos(data = gdp) +
                   cov_dyad_neg(data = gdp) +
                   # gwdasp(data = alpha) +
                   # gwdsp_pos(data = alpha)+
                   # gwdsp_neg(data = alpha)+
                   gwese_pos(data = alpha)+
                   gwese_neg(data = alpha) +
                   gwesf_pos(data = alpha) +
                   gwesf_neg(data = alpha)+
                   # # gwesp_neg(data = alpha) +
                   gwdegree_pos(data = alpha) +
                   gwdegree_neg(data = alpha) +
                   # degree_pos(type = c(3))
                   edges_neg,
                 control =control.sergm(method_est = c("Stepping"),tol = 0.2,
                                        method_var = c("MPLE"),max_it = 30,
                                        eval_likelihood = F, cluster = clust,
                                        sampler_var = sampler.sergm(number_networks = 3000,
                                                                    init_empty = T,
                                                                    n_proposals_burn_in = 100000,
                                                                    n_proposals = 10000,mh = F,seed = 13),
                                        sampler_path_sampling  = sampler.sergm(number_networks = 3000,
                                                                               init_empty = T,
                                                                               n_proposals_burn_in = 100000,
                                                                               n_proposals = 1000,mh = F,seed = 13),
                                        sampler_est = sampler.sergm(number_networks = 2000,
                                                                    init_empty = F,
                                                                    n_proposals_burn_in = 1000*10,
                                                                    n_proposals = 1000,mh = F,seed = 13),
                                        Stepping_number_grids = 1000))

summary(model_1)
save(model_1, file = "models/model_1.Rds")

model_1_assessment = t_model_assessment(model_1,cluster = clust,sampler.sergm(number_networks = 3000,
                                                                              init_empty = T,
                                                                              n_proposals_burn_in = 10000,
                                                                              n_proposals = 1000,mh = F,seed = 13))
save(model_1_assessment, file = "models/model_dependence_assessment.Rds")

load("models/model_dependence_assessment.Rds")

plots_1 = t_plot_model_assessment(mod_ass = model_1_assessment)
rm(model_1_assessment)
gc()
# Small function to have particular major and minor labels in a graph generated with ggplot2
minor_label = function(from = 0, to = 75, label = c(0,seq(10,70,10))) {
  res = character(length = length(from:to))
  res[label+1] = label
  return(res)
}

for(i in 1:10){
  cat(i,"\n")
  tmp_plot = plots_1[i,3][[1]]
  plots_1[i,3][[1]] = tmp_plot + scale_x_discrete(breaks = 0:75,label = minor_label(0,75,c(0,seq(10,70,10))))
  tmp_plot = plots_1[i,2][[1]]
  plots_1[i,2][[1]] = tmp_plot + scale_x_discrete(breaks = 0:75,label = minor_label(0,75,c(0,seq(5,35,5))))
  ggsave(plot = plots_1[i,3][[1]],filename = paste0("plots/positive_degree_gof_scale_",i,".png"),width = 6,height = 3.5)
  ggsave(plot = plots_1[i,1][[1]],filename = paste0("plots/ese_gof_scale_",i,".png"),width = 3.5,height = 3.5)
  ggsave(plot = plots_1[i,4][[1]],filename = paste0("plots/negative_degree_gof_scale_",i,".png"),width = 3.5,height = 3.5)
  ggsave(plot = plots_1[i,2][[1]],filename = paste0("plots/esf_gof_scale_",i,".png"),width = 6,height = 3.5)
}

clusterCall(clust, function() {gc(full = T)})

model_2 = tsergm(formula = tmp_network ~ edges_pos+
                   isolates_pos +
                   isolates_neg +
                   cov_dyad_pos(data = repetition_pos) +
                   cov_dyad_neg(data = repetition_neg) +
                   cov_dyad_pos(data = polity) +
                   cov_dyad_neg(data = polity) +
                   cov_dyad_pos(data = cinc) +
                   cov_dyad_neg(data = cinc) +
                   cov_dyad(data = distance) +
                   cov_dyad_pos(data = gdp) +
                   cov_dyad_neg(data = gdp) +
                   cov_dyad_pos(data = common_friends) +
                   cov_dyad_neg(data = common_friends) +
                   cov_dyad_pos(data = common_enemies) +
                   cov_dyad_neg(data = common_enemies) +
                   gwdegree_pos(data = alpha) +
                   gwdegree_neg(data = alpha) +
                   edges_neg,
                 control.sergm(method_est = c("Stepping"),tol = 0.2,
                               method_var = c("MPLE"),max_it = 30,
                               eval_likelihood = F, cluster = clust,
                               sampler_var = sampler.sergm(number_networks = 3000,
                                                           init_empty = T,
                                                           n_proposals_burn_in = 100000,
                                                           n_proposals = 10000,mh = F,seed = 13),
                               sampler_path_sampling  = sampler.sergm(number_networks = 3000,
                                                                      init_empty = T,
                                                                      n_proposals_burn_in = 100000,
                                                                      n_proposals = 1000,mh = F,seed = 13),
                               sampler_est = sampler.sergm(number_networks = 2000,
                                                           init_empty = F,
                                                           n_proposals_burn_in = 1000*10,
                                                           n_proposals = 1000,mh = F,seed = 13),
                               Stepping_number_grids = 1000))

clusterCall(clust, function() {gc(full = T)})

save(model_2, file = "models/model_2.Rds")

model_2_assessment = t_model_assessment(model_2,cluster = clust,sampler.sergm(number_networks = 3000,
                                                                              init_empty = T,
                                                                              n_proposals_burn_in = 10000,
                                                                              n_proposals = 1000,mh = F,seed = 13))

clusterCall(clust, function() {gc(full = T)})

plots_2 = t_plot_model_assessment(mod_ass = model_2_assessment)
rm(model_2_assessment)

for(i in 1:10){
  tmp_plot = plots_2[i,3][[1]]
  plots_2[i,3][[1]] = tmp_plot + scale_x_discrete(breaks = 0:75,label = minor_label(0,75,c(0,seq(10,70,10))))
  tmp_plot = plots_2[i,2][[1]]
  plots_2[i,2][[1]] = tmp_plot + scale_x_discrete(breaks = 0:75,label = minor_label(0,75,c(0,seq(5,35,5))))
  ggsave(plot = plots_2[i,3][[1]],filename = paste0("plots/alt_positive_degree_gof_scale_",i,".png"),width = 6,height = 3.5)
  ggsave(plot = plots_2[i,1][[1]],filename = paste0("plots/alt_ese_gof_scale_",i,".png"),width = 3.5,height = 3.5)
  ggsave(plot = plots_2[i,4][[1]],filename = paste0("plots/alt_negative_degree_gof_scale_",i,".png"),width = 3.5,height = 3.5)
  ggsave(plot = plots_2[i,2][[1]],filename = paste0("plots/alt_esf_gof_scale_",i,".png"),width = 6,height = 3.5)
}

# MCMC Plots ----
colnames(model_1$mcmc_chain) = c("Edges +", "Edges -", "Isolates +",
                                 "Isolates -", "Stability +", "Stability -",
                                 "Abs. Polity Diff. +", "Abs. Polity Diff. -",
                                 "CINC Ratio +", "CINC Ratio -",
                                 "Abs. GDP Diff. +", "Abs. GDP Diff. -","Abs. Distance ±",
                                 "GWESE +", "GWESE -",
                                 "GWESF +", "GWESF -",
                                 "GWD +","GWD -")


S <- ggs(model_1$mcmc_chain)
# Plot the densities of the samples
density_1 = ggs_density(S[S$Parameter %in% c("Edges +","Edges -","Isolates +","Isolates -"),]) + theme_pubr() + ylab("Density")
density_2 = ggs_density(S[S$Parameter %in% c("Stability +", "Stability -", "Abs. Polity Diff. +","Abs. Polity Diff. -"),]) + theme_pubr() + ylab("Density")
density_3 = ggs_density(S[S$Parameter %in% c("CINC Ratio +", "CINC Ratio -", "Abs. GDP Diff. +","Abs. GDP Diff. -"),]) + theme_pubr() + ylab("Density")
density_4 = ggs_density(S[S$Parameter %in% c("Abs. Distance ±", "GWD +","GWD -"),]) + theme_pubr() + ylab("Density")
density_5 = ggs_density(S[S$Parameter %in% c("GWESE +", "GWESE -", "GWESF +","GWESF -"),]) + theme_pubr() + ylab("Density")

# Plot the chains of the samples
trace_1 = ggs_traceplot(S[S$Parameter %in% c("Edges +","Edges -","Isolates +","Isolates -"),]) + theme_pubr()+ ylab("Trace plot")
trace_2 = ggs_traceplot(S[S$Parameter %in% c("Stability +", "Stability -", "Abs. Polity Diff. +","Abs. Polity Diff. -"),]) + theme_pubr()+ ylab("Trace plot")
trace_3 = ggs_traceplot(S[S$Parameter %in% c("CINC Ratio +", "CINC Ratio -", "Abs. GDP Diff. +","Abs. GDP Diff. -"),]) + theme_pubr() + ylab("Trace plot")
trace_4 = ggs_traceplot(S[S$Parameter %in% c("Abs. Distance ±", "GWD +","GWD -"),]) + theme_pubr() + ylab("Trace plot")
trace_5 = ggs_traceplot(S[S$Parameter %in% c("GWESE +", "GWESE -", "GWESF +","GWESF -"),]) + theme_pubr() + ylab("Trace plot")


res_1 = grid.arrange(trace_1,density_1, ncol = 2)
res_2 = grid.arrange(trace_2,density_2, ncol = 2)
res_3 = grid.arrange(trace_3,density_3, ncol = 2)
res_4 = grid.arrange(trace_4,density_4, ncol = 2)
res_5 = grid.arrange(trace_5,density_5, ncol = 2)

ggsave(res_1,filename = "plots/dca_mcmc_diag_1.pdf",width = 7,height = 6.66)
ggsave(res_2,filename = "plots/dca_mcmc_diag_2.pdf",width = 7,height = 6.66)
ggsave(res_3,filename = "plots/dca_mcmc_diag_3.pdf",width = 7,height = 6.66)
ggsave(res_4,filename = "plots/dca_mcmc_diag_4.pdf",width = 7,height = 6.66*0.75)
ggsave(res_5,filename = "plots/dca_mcmc_diag_5.pdf",width = 7,height = 6.66)
