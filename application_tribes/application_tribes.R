library(signnet)
library("Rcpp")
library(ergm)
library(ggraph)
library(igraph)
library(tidygraph)
library(ggpubr)
library(coda)
library(ergm.sign)
library(stringr)
library(data.table)
library(ggmcmc)
model_assessment

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source(file = "../other_functions.R")
rm(list=ls())

data("tribes")
save(tribes, file = paste0("tribes.RData"))
g = tribes
network_tmp = as.matrix(as_adjacency_matrix(g,sparse = F,attr = "sign"))

library(parallel)
no_cores <- 20
# Setup cluster
clust <- makeForkCluster(no_cores)
clusterCall(clust, function() {library(ergm.sign);library(Rcpp);library(RcppArmadillo)})

model_mple = sergm(network_tmp ~ edges_pos+edges_neg,
                   control = control.sergm(method_est = "MPLE",method_var = "MPLE",eval_likelihood =T,
                                           cluster = clust,
                                           sampler_var = sampler.sergm(number_networks = 1000),
                                           sampler_path_sampling = sampler.sergm(seed = 12, number_networks = 1000,
                                                                                 n_proposals_burn_in =  10000,
                                                                                 n_proposals = 2000,mh = F,init_empty = T)))

alpha = matrix(1.5,nrow = 1,ncol = 1)

model_mle = sergm(formula = network_tmp ~ edges_pos+edges_neg+
                    gwese_pos(data = alpha) +
                    gwese_neg(data = alpha) +
                    gwesf_pos(data = alpha) +
                    gwdegree_pos(data = alpha),
                  control = control.sergm(method_est = c("Stepping"),tol = 0.1,n_bridges = 16,
                                          mple_init = T,cluster = clust,eval_likelihood = T,
                                          max_it =20,Stepping_number_grids = 2000,
                                          sampler_est = sampler.sergm(seed = 123,number_networks = 1000,
                                                                      n_proposals_burn_in =  10000,
                                                                      n_proposals = 1000,mh = F,init_empty = T), 
                                          sampler_var = sampler.sergm(seed = 123,number_networks = 3000,
                                                                      n_proposals_burn_in =  10000,init_empty = T,
                                                                      n_proposals = 1000,mh = F),
                                          sampler_path_sampling = sampler.sergm(seed = 12, number_networks = 1000,
                                                                        n_proposals_burn_in =  10000,
                                                                        n_proposals = 2000,mh = F,init_empty = T)))


stopCluster(clust)
mod_ass_1 = model_assessment(model_mle,sampler = sampler.sergm(seed = 123,number_networks = 1000,
                                                               n_proposals_burn_in =  10000,
                                                               n_proposals = 1000,mh = F))

mod_ass_2 = model_assessment(model_mple,sampler = sampler.sergm(seed = 123,number_networks = 1000,
                                                         n_proposals_burn_in =  10000,
                                                         n_proposals = 1000,mh = F))

library(cowplot)
require(gridExtra)
# Figures 22 and 23 in 4.2 of the SM)
pdf("tribes_plots/tribes_model_assessment_1.pdf")
res = plot_model_assessment(mod_ass = mod_ass_1)
grid.arrange(res[[1]], res[[2]], res[[3]], res[[4]], ncol=2)
dev.off()
pdf("tribes_plots/tribes_model_assessment_2.pdf")
res = plot_model_assessment(mod_ass_2)
grid.arrange(res[[1]], res[[2]], res[[3]], res[[4]], ncol=2)
dev.off()


# Table 4 from 4.2 in the SM
res = t_generate_table(digits = 3, model_mle, model_mple)
library(Hmisc)
old_names = res[,1]
new_names = character(length(old_names))
new_names = c("Edges Positive", "Edges Negative", "GWESE", "GWESF^+", "GWESF^-", 
                           "GWD^+", "AIC")


# MCMC Plots (Figures 24 and 25 in 4.3 of the SM)----
colnames(model_mle$mcmc_chain) = c("Edges +", "Edges -","GWESE +", "GWESE -", "GWESF +", 
                                 "GWD +")

S <- ggs(model_mle$mcmc_chain)
# Plot the densities of the samples 
density_1 = ggs_density(S[S$Parameter %in% c("Edges +","Edges -","GWD +"),]) + theme_pubr() + ylab("Density")
density_2 = ggs_density(S[S$Parameter %in% c("GWESE +", "GWESE -", "GWESF +"),]) + theme_pubr() + ylab("Density")

# Plot the chains of the samples 
trace_1 = ggs_traceplot(S[S$Parameter %in%  c("Edges +","Edges -","GWD +"),]) + theme_pubr() + ylab("Trace plot")
trace_2 = ggs_traceplot(S[S$Parameter %in% c("GWESE +", "GWESE -", "GWESF +"),]) + theme_pubr()+ ylab("Trace plot")

res_1 = grid.arrange(trace_1,density_1, ncol = 2)
res_2 = grid.arrange(trace_2,density_2, ncol = 2)

ggsave(res_1,filename = "tribes_plots/tribes_mcmc_diag_1.pdf",width = 7,height = 6.66)
ggsave(res_2,filename = "tribes_plots/tribes_mcmc_diag_2.pdf",width = 7,height = 6.66)
