library(ergm.sign)
library(dplyr)
library(stringr)
n_actors = 100
n_sim = 100
# gt = c("edges_pos" = -5,"edges_neg" = -3,"gwdegree_pos" = 0.5, "gwdegree_neg" = -0.5,"gwese_pos" = 0.7,"gwesf_pos" = 1)
gt = c("edges_pos" = -2,"edges_neg" = -4,"gwdegree_pos" = -2, "gwdegree_neg" = -0.3,
       "gwese_pos" = 0.6,"gwesf_pos" = 0.5,"gwese_neg" = 0.5,"gwesf_neg" = 0.3)
seed = 123
n_obs = 10

network = matrix(data = 0,nrow = n_actors, ncol = n_actors)
tmp_network = replicate(n_obs,network, simplify=FALSE)
alpha =  replicate(n_obs, matrix(log(2)), simplify=FALSE)
# formula_tmp = tmp_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+ gwdegree_neg(data = alpha) + gwese_pos(data = alpha)+ gwesf_pos(data = alpha)
formula_tmp = tmp_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+ gwdegree_neg(data = alpha) + gwese_pos(data = alpha)+ gwesf_pos(data = alpha) +
  gwese_neg(data = alpha) + gwesf_neg(data = alpha)

simulated_nets = tsergm_simulate(formula_tmp, 
                                 coef = gt, 
                                 only_stats = F,
                                 sampler = sampler.sergm(number_networks = n_sim,
                                                         seed = seed,mh = F,
                                                         n_proposals_burn_in = 2000000,
                                                         n_proposals = 20000))

plot(simulated_nets[[1]]$stats[,8])

acf(simulated_nets[[1]]$stats[,6])

networks = list()
for(k in 1:n_sim){
  simulated_network = list()
  for(j in 1:n_obs){
    simulated_network[[j]] = map_to_mat(edges_pos = simulated_nets[[j]]$networks[[k]]$edges_pos, 
                                        edges_neg = simulated_nets[[j]]$networks[[k]]$edges_neg, n_actors = n_actors)
  }
  networks[[k]] = simulated_network
}

library(parallel)
no_cores <- 20
# Setup cluster
clust <- makeForkCluster(no_cores, outfile = "output.txt")
clusterCall(clust, function() {library(ergm.sign);library(Rcpp);library(RcppArmadillo)})

function_simulate_cp = function(k,networks, n_actors, alpha) {
  res = list()
  now = Sys.time()
  simulated_network = networks[[k]]
  model_stepping =  tsergm(formula = simulated_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+
                             gwdegree_neg(data = alpha) + gwese_pos(data = alpha)+ gwesf_pos(data = alpha),
                           control = control.sergm(method_est = "Stepping",eval_likelihood = F,tol =0.001, method_var = "Fisher",
                                                   sampler_est = sampler.sergm(number_networks = 250,init_empty = T,seed = 122),
                                                   sampler_var = sampler.sergm(number_networks = 500,init_empty = T,seed = 123)))
  time_taken  = as.numeric(difftime(Sys.time() , now), units="hours")
  
  # res_stepping = list(model = model_stepping, model_coef = model_stepping$coefficients,model_var = model_stepping$var, time = time_taken)
  res_stepping = list(model_coef = model_stepping$coefficients,model_var = model_stepping$var, time = time_taken)
  rm(model_stepping)
  gc(full = T)
  now = Sys.time()
  mdel_mple =  tsergm(formula = simulated_network ~ edges_pos + edges_neg + gwdegree_pos(data = alpha)+
                        gwdegree_neg(data = alpha) + gwese_pos(data = alpha)+ gwesf_pos(data = alpha),
                      control = control.sergm(method_est = "MPLE",eval_likelihood = F,method_var = "MPLE"))
  time_taken  = as.numeric(difftime(Sys.time() , now), units="hours")
  # res_mple = list(model = mdel_mple, model_coef = mdel_mple$coefficients,model_var = mdel_mple$var, time = time_taken)
  res_mple = list(model_coef = mdel_mple$coefficients,model_var = mdel_mple$var, time = time_taken)
  rm(mdel_mple)
  
  gc(full = T)
  save(res_stepping, file = paste0("simulation/results/simulation_2_res_stepping_",k,".RData"))
  save(res_mple, file = paste0("simulation/results/simulation_2_res_mple_",k,".RData"))
  cat("Save!\n")
} 

# debugonce(function_simulate_cp)
# res = lapply(X = 1:n_sim,FUN = function_simulate_cp,networks = networks,
#            n_actors = n_actors, alpha = alpha)

parLapply(X = 1:n_sim,cl = clust,fun = function_simulate_cp,n_actors = n_actors,networks = networks, alpha = alpha) 
# save(res, file = paste0("../sergm_replication/simulation/results/res_cp.RData"))
stopCluster(clust)
file.remove("output.txt")
gc(full = T)

res = list()

for(i in 1:n_sim){
  load(file = paste0("simulation/results/simulation_2_res_stepping_",i,".RData"))
  load(file = paste0("simulation/results/simulation_2_res_mple_",i,".RData"))
  res[[i]] = list()
  res[[i]]$mle = res_stepping
  res[[i]]$mple = res_mple
}
errors_ML = list()
errors_MPL = list()
coefs_ML = list()
coefs_MPL = list()
is_within_ML = list()
is_within_MPL = list()

for(i in 1:length(res)){
  tmp_mle_model = res[[i]]$mle
  tmp_mple_model = res[[i]]$mple
  
  coefs_ML[[i]] = tmp_mle_model$model_coef
  errors_ML[[i]] = (tmp_mle_model$model_coef - gt)^2
  confint_tmp = matrix(c(tmp_mle_model$model_coef - sqrt(diag(tmp_mle_model$model_var))*qnorm(p = 1-0.025), 
           tmp_mle_model$model_coef + sqrt(diag(tmp_mle_model$model_var))*qnorm(p = 1-0.025)),ncol = 2)
  res_tmp = c()
  for(j in 1:length(tmp_mle_model$model_coef)){
    res_tmp[j] =  between(  gt[j],left = confint_tmp[j,1], right = confint_tmp[j,2])
  }
  is_within_ML[[i]] = res_tmp
  #
  #
  coefs_MPL[[i]] = tmp_mple_model$model_coef
  errors_MPL[[i]] = (tmp_mple_model$model_coef - gt)^2
  confint_tmp = matrix(c(tmp_mple_model$model_coef - sqrt(diag(tmp_mple_model$model_var))*qnorm(p = 1-0.025), 
                         tmp_mple_model$model_coef + sqrt(diag(tmp_mple_model$model_var))*qnorm(p = 1-0.025)),ncol = 2)
  res_tmp = c()
  for(j in 1:length(tmp_mle_model$model_coef)){
    res_tmp[j] =  between(  gt[j],left = confint_tmp[j,1], right = confint_tmp[j,2])
  }
  is_within_MPL[[i]] = res_tmp
}

coefs_ML = do.call(rbind,coefs_ML)
coefs_MPL = do.call(rbind,coefs_MPL)
errors_MPL = do.call(rbind,errors_MPL)
errors_ML = do.call(rbind,errors_ML)
cp_ML = do.call(rbind,is_within_ML)
cp_MPL = do.call(rbind,is_within_MPL)
colMeans(cp_ML)
colMeans(cp_MPL)
colMeans(sqrt(errors_ML))
colMeans(errors_MPL)
result = list(coefs_ML = coefs_ML,
              coefs_MPL = coefs_MPL, 
              errors_MPL = errors_MPL,
              errors_ML = errors_ML, 
              cp_ML = cp_ML)


get_info_simulation = function(result,gt, digits = 3) {
  coefs_ml = result$coefs_ML
  ave_ml = apply(coefs_ml,2, mean )
  coefs_std = t((t(coefs_ml) - gt)^2)
  rmse_ml = sqrt(1/nrow(coefs_std)*colSums(coefs_std))
  cp_ml = colMeans(cp_ML)
  # Get info for rem
  coefs_mpl = result$coefs_MPL
  ave_mpl = apply(coefs_mpl,2, mean )
  coefs_std = t((t(coefs_mpl) - gt)^2)
  rmse_mpl = sqrt(1/nrow(coefs_std)*colSums(coefs_std))
  cp_mpl = colMeans(cp_MPL)
  
  res = cbind(coefs = gt,ave_ml,rmse_ml,cp_ml,ave_mpl,rmse_mpl,cp_mpl)
  rownames(res) = str_to_title(names(gt))
  tmp = rownames(res)
  tmp = gsub("\\_",replacement = " ",tmp)
  rownames(res) = tmp
  res = round(res,digits = digits)
  class(res) = "result_simulation"
  return(res)
}
result_simulation = get_info_simulation(result = result,gt = gt,digits = 3)
library(Hmisc)


latex(result_simulation, label = "tbl:simulation_res",    file = "",
      colheads = c("Coef.", rep(c("AVE", "RMSE", "CP"),times = 2)),
      rowname = c("Edges +", "Edges -", "GWD +", "GWD -", "GWESE +", "GWESF +"),
      cgroup = c("", "ML", "MPL"), n.cgroup = c(1, 3, 3), cgroupTexCmd = " ",
      rgroupTexCmd = " ",colnamesTexCmd = " ",
              caption = "Result of the simulation study for the EcREM and REM \n      with the three data-generating processes (DG 1, DG2, DG3). \n      For each DG (data-generating process) and covariate, we note the AVE (average estimate), RMSE (Root-Mean-Squared Error), \n      and CP (Coverage Probability). Further, for each Data Generating process (DG) the Percentage of True Events (PTE) is given.")

# 
# errors_ml = sweep(coefs_ML,MARGIN = 2,STATS = gt,FUN = "-")
# errors_mpl = sweep(coefs_MPL,MARGIN = 2,STATS = gt,FUN = "-")
# apply(errors_mpl^2,MARGIN = 2,FUN = sum)
# apply(errors_ml^2,MARGIN = 2,FUN = sum)
# colMeans(coefs_ML)
# apply(coefs_ML,MARGIN = 2,FUN = var)
# apply(coefs_MPL,MARGIN = 2,FUN = var)
# 
# colMeans(coefs_MPL)
