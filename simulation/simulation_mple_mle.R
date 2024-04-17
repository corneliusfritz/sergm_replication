library(ergm.sign)
library(data.table)
library(ggpubr)
n_actors = 100
n_sim = 100
n_obs_grid = seq(from = 5, to = 50, by = 5)
gt = c(-5,-3, 0.5, -0.5, 0.5)
seed = 123

library(parallel)
no_cores <- 10
# Setup cluster
clust <- makeForkCluster(no_cores, outfile = "output.txt")
clusterCall(clust, function() {library(ergm.sign);library(Rcpp);library(RcppArmadillo)})
n_sim = 100

simulation_mple_mle = function(i,n_actors, n_sim, n_obs_grid, gt, seed) {
  n_obs = n_obs_grid[i]
  network = matrix(data = 0,nrow = n_actors, ncol = n_actors)
  tmp_network = replicate(n_obs,network, simplify=FALSE)
  alpha =  replicate(n_obs, matrix(log(2)), simplify=FALSE)
  formula_tmp = tmp_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+ gwdegree_neg(data = alpha) + gwese_pos(data = alpha)
  environment(formula_tmp) = environment()
  simulated_nets = tsergm_simulate(formula_tmp, coef = gt, only_stats = F,
                                   sampler = sampler.sergm(number_networks = n_sim,seed = seed,mh = F))
  simulated_nets[[1]]$stats
  res_stepping = list()  
  res_mple = list()  
  
  for(k in 1:n_sim) {
    simulated_network = list()
    for(j in 1:n_obs){
      simulated_network[[j]] = map_to_mat(edges_pos = simulated_nets[[j]]$networks[[k]]$edges_pos, 
                                          edges_neg = simulated_nets[[j]]$networks[[k]]$edges_neg, n_actors = n_actors)
    }
    now = Sys.time()
    model_stepping =  tsergm(formula = simulated_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+ gwdegree_neg(data = alpha) + gwese_pos(data = alpha),
                             control = control.sergm(method_est = "Stepping",eval_likelihood = F,method_var = "MPLE",tol =0.01, 
                                                     sampler_est = sampler.sergm(number_networks = 150,init_empty = T)))
    # sampler_est = sampler.sergm(number_networks = 150,seed = seed*k,mh = F,init_empty = F)
    time_taken  = as.numeric(difftime(Sys.time() , now), units="hours")
    
    res_stepping[[k]] = list(model = model_stepping, time = time_taken)
  }
  save(res_stepping, file = paste0("simulation/results/simulation_1_res_stepping_n_obs_",n_obs_grid[i],".RData"))
  cat("Finished with n = ", n_obs_grid[i], "\n")
}


parLapply(X = 1:length(n_obs_grid),cl = clust,fun = simulation_mple_mle,
          n_actors = n_actors, n_sim = n_sim, n_obs_grid = n_obs_grid,
          gt = gt, seed = seed) 

stopCluster(clust)

file.remove("output.txt")
complete_res_stepping = list()


for(i in 1:10){
  load(paste0("simulation/results/simulation_1_res_stepping_n_obs_",n_obs_grid[i],".RData"))
  tmp_res_stepping = rbindlist(lapply(res_stepping, FUN = function(x){as.data.table(t(c(x$model$coefficients,x$time)))}))
  names(tmp_res_stepping) = c("Edges +", "Edges -", "GWDegree +", "GWDegree -", "GWESE","time")
  # names(tmp_res_stepping) = c("Edges +", "Edges -", "GWDegree +", "time")
  tmp_res_stepping$number_networks = n_obs_grid[i]
  complete_res_stepping[[i]] = tmp_res_stepping 
}

complete_res_stepping = rbindlist(complete_res_stepping)
complete_res_stepping$time = complete_res_stepping$time*60

melt_res_stepping = melt.data.table(complete_res_stepping,id.vars = "number_networks")

# New facet label names 
supp.labs <- c("Edges +", "Edges -", "GWD +", "GWD -", "GWESE +", "Time in Min")
names(supp.labs) <- levels(melt_res_stepping$variable)

pdf("simulation/plot_mle.pdf")
ggplot(data = melt_res_stepping, aes(x = factor(number_networks), y = value)) +
  geom_boxplot(fill = "grey80") +
  geom_hline(data = dummy2, aes(yintercept = Z), color = "gold", lty = 2) +
  facet_grid(rows = vars(variable),scales="free",
             labeller = labeller(variable = supp.labs)) +
  stat_summary(fun =mean, geom="point",  color = c("darkred"), shape = 17, size = 2) +
  xlab("Number of Networks") +
  ylab("Estimates") + 
  theme_pubr()
dev.off()
