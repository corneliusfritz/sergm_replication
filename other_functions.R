
NVL = function (...) {
  for (e in eval(substitute(alist(...)))) {
    x <- eval(e, parent.frame())
    if (!is.null(x)) 
      break
  }
  x
}

.catchToList = function (expr) 
{
  val <- NULL
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, w$message)
    invokeRestart("muffleWarning")
  }
  myError <- NULL
  eHandler <- function(e) {
    myError <<- e$message
    NULL
  }
  val <- tryCatch(withCallingHandlers(expr, warning = wHandler), 
                  error = eHandler)
  list(value = val, warnings = myWarnings, error = myError)
}

ERRVL = function (...) 
{
  x <- NULL
  for (e in eval(substitute(alist(...)))) {
    x <- eval(if (inherits(x, "try-error")) 
      do.call(substitute, list(e, list(. = x)))
      else e, parent.frame())
    if (!inherits(x, "try-error")) 
      return(x)
  }
  stop("No non-error expressions passed.")
}

NVL2 = function (test, notnull, null = NULL) 
{
  if (is.null(test)) 
    null
  else notnull
}


sergm.MCMCse <- function (theta0, theta,statsmatrices,H) { 
  av <- colMeans(statsmatrices)
  xsims <- sweep(statsmatrices, 2, av)
  # to change
  gsims <--xsims
  xobs <- -av
  xsim <- as.matrix(xsims)
  gsim <- as.matrix(gsims)
  prob <- rep.int(1/nrow(xsim), nrow(xsim))
  cov.zbar <- spectrum0.mvar(gsims) * sum(prob^2)
  imp.factor <- sum(prob^2) * length(prob)
  novar <- rep(TRUE, nrow(H))
  novar <- diag(H) < sqrt(.Machine$double.eps)
  cov.zbar.obs <- cov.zbar
  cov.zbar.obs[, ] <- 0
  H.obs <- H
  H.obs[, ] <- 0
  imp.factor.obs <- NULL
  H <- H[!novar, !novar, drop = FALSE]
  H.obs <- H.obs[!novar, !novar, drop = FALSE]
  cov.zbar <- cov.zbar[!novar, !novar, drop = FALSE]
  cov.zbar.obs <- cov.zbar.obs[!novar, !novar, drop = FALSE]
  mc.cov.offset <- matrix(0, ncol = length(theta), nrow = length(theta))
  H <- H.obs - H
  cov.zbar <- cov.zbar + cov.zbar.obs
  mc.cov <- matrix(NA, ncol = length(novar), nrow = length(novar))
  mc.cov0 <- solve(H, cov.zbar, tol = 1e-20)
  mc.cov0 <- solve(H, t(mc.cov0), tol = 1e-20)
  mc.cov[!novar, !novar] <- mc.cov0
  attr(mc.cov, "imp.factor") <- imp.factor
  return(mc.cov)
}


spectrum0.mvar = function (x, order.max = NULL, aic = is.null(order.max), tol = .Machine$double.eps^0.5, ...) 
{ breaks <- if (is.mcmc.list(x)) 
    c(0, cumsum(sapply(x, niter)))
  else NULL
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  v <- matrix(0, p, p)
  novar <- abs(apply(x, 2, stats::sd)) < tol
  x <- x[, !novar, drop = FALSE]
  if (ncol(x) == 0) 
    stop("All variables are constant.")
  first_local_min <- function(x) {
    d <- diff(c(Inf, x, Inf))
    min(which(d >= 0)) - 1
  }
  e <- eigen(cov(x), symmetric = TRUE)
  Q <- e$vectors[, sqrt(pmax(e$values, 0)/max(e$values)) > 
                   tol * 2, drop = FALSE]
  xr <- x %*% Q
  ind.var <- cov(xr)
  xr <- if (!is.null(breaks)) 
    do.call(mcmc.list, lapply(lapply(seq_along(breaks[-1]), 
                                     function(i) xr[(breaks[i] + 1):(breaks[i + 1]), , 
                                                    drop = FALSE]), mcmc))
  else as.mcmc.list(mcmc(xr))
  ord <- NVL(order.max, ceiling(10 * log10(niter(xr))))
  xr <- do.call(rbind, c(lapply(unclass(xr)[-nchain(xr)], function(z) rbind(cbind(z), 
                                                                            matrix(NA, ord, nvar(z)))), unclass(xr)[nchain(xr)]))
  arfit <- .catchToList(ar(xr, aic = is.null(order.max), order.max = ord, 
                           na.action = na.pass, ...))
  while ((!is.null(arfit$error) || ERRVL(try(any(eigen(arfit$value$var.pred, 
                                                       only.values = TRUE)$values < 0), silent = TRUE), TRUE)) && 
         ord > 0) {
    ord <- ord - 1
    if (ord <= 0) 
      stop("Unable to fit ar() even with order 1; this is likely to be due to insufficient sample size or a trend in the data.")
    arfit <- .catchToList(ar(xr, aic = is.null(order.max), 
                             order.max = ord, na.action = na.pass, ...))
  }
  arfit <- arfit$value
  if (aic && arfit$order > (ord <- first_local_min(arfit$aic) - 
                            1)) {
    arfit <- ar(xr, aic = ord == 0, order.max = max(ord, 
                                                    1), na.action = na.pass)
  }
  arvar <- arfit$var.pred
  arcoefs <- arfit$ar
  arcoefs <- NVL2(dim(arcoefs), apply(arcoefs, 2:3, base::sum), 
                  sum(arcoefs))
  adj <- diag(1, nrow = ncol(xr)) - arcoefs
  iadj <- solve(adj)
  v.var <- iadj %*% arvar %*% t(iadj)
  infl <- exp((determinant(v.var)$modulus - determinant(ind.var)$modulus)/ncol(ind.var))
  v.var <- Q %*% v.var %*% t(Q)
  v[!novar, !novar] <- v.var
  attr(v, "infl") <- infl
  v
}

create_dataset_dca = function(years = 1995:2005, atop = T, only_biggest_cluster = T, 
                              aggregate_years = 3, mid_var = "any_mid", alliance_var = "any_atop", dca_data) {
  # Year t
  years_t = seq(to = max(years), by = 1, length.out = aggregate_years)
  years_t_minus = seq(to = max(years)- aggregate_years, by = 1, length.out = aggregate_years)
  edgeList = create_dyadyears(system = "cow", mry = TRUE, directed = TRUE, subset_years = years)
  edgeList = add_gwcode_to_cow(edgeList)
  edgeList = add_atop_alliance(edgeList)
  edgeList$any_atop = rowSums(edgeList[,grep("atop", names(edgeList))])>0
  edgeList = add_cow_alliance(edgeList)
  edgeList$any_cow = rowSums(edgeList[,grep("cow", names(edgeList))])>0
  edgeList = add_cow_mids(edgeList)
  edgeList$any_mid = (edgeList$cowmidongoing + edgeList$cowmidonset) > 0
  edgeList = add_cow_trade(edgeList)
  edgeList$acc_flow = edgeList$smoothflow1 + edgeList$flow2
  edgeList = add_democracy(edgeList)
  edgeList = add_minimum_distance(edgeList)
  edgeList = add_capital_distance(edgeList)
  edgeList = add_nmc(edgeList)
  edgeList = add_contiguity(edgeList)
  edgeList = add_sdp_gdp(edgeList)
  edgeList$id = paste(edgeList$ccode1, edgeList$ccode2,sep = "_")
  edgeList$years_t = edgeList$year %in% years_t
  edgeList$polity21[is.infinite(edgeList$polity21)] = NA
  setDT(edgeList)
  edgeList$year_country_dyad_id = paste(edgeList$ccode1,edgeList$ccode2,edgeList$year)
  dca_data$year_country_dyad_id = paste(dca_data$ccode1,dca_data$ccode2,dca_data$year)
  # edgeList$dcaAnyV1 = dca_data$dcaAnyV1[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  # edgeList$dcaAnyV2 = dca_data$dcaAnyV2[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  
  edgeList$dcaAnyV1 = dca_data$dcaGeneralV1[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  edgeList$dcaAnyV2 = dca_data$dcaAnyV2[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  edgeList$dcaGeneralV1 = dca_data$dcaGeneralV1[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  edgeList$dcaGeneralV2 = dca_data$dcaGeneralV2[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  edgeList$dcaSectorV1 = dca_data$dcaSectorV1[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  edgeList$dcaSectorV2 = dca_data$dcaSectorV2[match(edgeList$year_country_dyad_id,dca_data$year_country_dyad_id)]
  
  # edgeList$cinc1 = cinc$cinc[match(edgeList$year_country_1_id,cinc$year_country_id)]
  # edgeList$cinc2 = cinc$cinc[match(edgeList$year_country_2_id,cinc$year_country_id)]
  edgeList$cinc_1_bigger = edgeList$cinc1>edgeList$cinc2
  edgeList$cinc1_tmp = 0
  edgeList$cinc1_tmp[edgeList$cinc_1_bigger] = edgeList$cinc1[edgeList$cinc_1_bigger]
  edgeList$cinc1_tmp[!edgeList$cinc_1_bigger] = edgeList$cinc2[!edgeList$cinc_1_bigger]
  edgeList$cinc2_tmp = 0
  edgeList$cinc2_tmp[edgeList$cinc_1_bigger] = edgeList$cinc2[edgeList$cinc_1_bigger]
  edgeList$cinc2_tmp[!edgeList$cinc_1_bigger] = edgeList$cinc1[!edgeList$cinc_1_bigger]
  edgeList$cinc_rate = log(edgeList$cinc1_tmp/edgeList$cinc2_tmp)
  
  
  
  n_actors = length(unique(edgeList$ccode1))
  unique_years = unique(edgeList$year)
  tmp_matrix = matrix(data = 0, nrow = n_actors, ncol = n_actors)
  colnames(tmp_matrix) = unique(edgeList$ccode1)
  rownames(tmp_matrix) = unique(edgeList$ccode1)
  edgeList$id_1 = match(edgeList$ccode1, unique(edgeList$ccode1))
  edgeList$id_2 = match(edgeList$ccode2, unique(edgeList$ccode1))
  
  target_matrices = list()
  cov_repetition_list =  list()
  cov_repetition_pos_list =  list()
  cov_repetition_neg_list =  list()
  cov_repetition_neut_list = list()
  cov_contiguity_list = list()
  cov_distance_list =  list()
  cov_polity_list =  list()
  cov_gdp_list =  list()
  cov_gdppc_list =  list()
  cov_cinc_list =  list()
  cov_population_list =  list()
  common_friends_list =  list()
  common_enemies_list = list()
  n = 1
  # Include all countries that were active at least once
  active_countries = unique(c(edgeList[any_mid ==1 | dcaAnyV2 ==1]$ccode1, 
           edgeList[any_mid ==1 | dcaAnyV2 ==1]$ccode2))
  additional_exlude_missing = unique(c(unique(edgeList[is.na(edgeList$wbgdp2011est1),ccode1]),
           unique(edgeList[is.na(edgeList$wbgdppc2011est1),ccode1]),
           unique(edgeList[is.na(edgeList$polity21),ccode1])))
  
  tmp = edgeList[,.(include = length(unique(year)) == length(years)), by = .(ccode1)]
  
  active_countries = active_countries[!(active_countries %in%additional_exlude_missing) & !(active_countries %in%tmp$ccode1[tmp$include == F])]
  
  active_ind = rownames(tmp_matrix) %in% active_countries
  
  # edgeList = edgeList[(ccode1 %in% active_countries) & (ccode2 %in% active_countries)]
 
  
  
  # table(c(edgeList[any_mid ==1 | dcaAnyV2 ==1]$id_1, 
  #   edgeList[any_mid ==1 | dcaAnyV2 ==1]$id_2))
  # 
  # 
  # 
  # c(edgeList[any_mid ==1 | dcaAnyV2 ==1]$id_1, 
  #   edgeList[any_mid ==1 | dcaAnyV2 ==1]$id_2)
  # 
  for(t in unique_years[-1]){
    target_matrix = tmp_matrix
    cov_repetition = tmp_matrix
    edgeList_tmp = edgeList[year == t]
    edgeList_old = edgeList[year == (t-1)]
    
    target_matrix[as.matrix(edgeList_tmp[any_mid ==1,.(id_1,id_2)])] = -1
    target_matrix[as.matrix(edgeList_tmp[dcaAnyV2 ==1,.(id_1,id_2)])] = 1
    target_matrices[[n]] = target_matrix[active_ind,active_ind]
    
    cov_repetition[as.matrix(edgeList_old[any_mid ==1,.(id_1,id_2)])] = -1
    cov_repetition[as.matrix(edgeList_old[dcaAnyV2 ==1,.(id_1,id_2)])] = 1
    cov_repetition_list[[n]] = cov_repetition[active_ind,active_ind]
    
    cov_repetition_pos_list[[n]] = (cov_repetition==1)[active_ind,active_ind]
    cov_repetition_neg_list[[n]] =(cov_repetition==-1)[active_ind,active_ind]
    cov_repetition_neut_list[[n]] = (cov_repetition==0)[active_ind,active_ind]
    
    common_friends_list[[n]] = ((cov_repetition==1)%*%(cov_repetition==1))[active_ind,active_ind]
    common_enemies_list[[n]] = ((cov_repetition==-1)%*%(cov_repetition==-1))[active_ind,active_ind]
    
    
    cov_contiguity = tmp_matrix
    cov_contiguity[as.matrix(edgeList_old[conttype != 0, .(id_1, id_2)])] = 1
    cov_contiguity_list[[n]] = cov_contiguity[active_ind,active_ind]
    cov_distance = tmp_matrix
    cov_distance[as.matrix(edgeList_old[,.(id_1, id_2)])] = log1p(edgeList_old$capdist)
    cov_distance_list[[n]] = cov_distance[active_ind,active_ind]
    
     # ind_missing = rowSums(is.na(cov_distance)) > 10
    
    cov_cinc = tmp_matrix
    cov_cinc[as.matrix(edgeList_old[, .(id_1, id_2)])] = edgeList_old$cinc_rate
    cov_cinc_list[[n]] = cov_cinc[active_ind,active_ind]
    
    
    actor_df = data.table(id = match(active_countries,rownames(tmp_matrix)), 
                          cowid = active_countries)
    
    actor_df$gdppc =  log(edgeList_old$wbgdppc2011est1[match( actor_df$id, edgeList_old$id_1)])
    actor_df$gdp =  log(edgeList_old$wbgdp2011est1[match( actor_df$id, edgeList_old$id_1)])
    actor_df$pop =  log(edgeList_old$wbpopest1[match( actor_df$id, edgeList_old$id_1)])
    actor_df$polity =  edgeList_old$polity21[match( actor_df$id, edgeList_old$id_1)]
    cov_polity = outer(actor_df$polity,actor_df$polity, FUN = function(x,y) abs(x-y))
    cov_gdp = outer(actor_df$gdp,actor_df$gdp, FUN = function(x,y) abs(x-y))
    cov_gdppc =  outer(actor_df$gdppc,actor_df$gdppc, FUN = function(x,y) abs(x-y))
    cov_population = outer(actor_df$pop,actor_df$pop, FUN = function(x,y) abs(x-y))
    cov_polity_list[[n]] = cov_polity
    cov_gdp_list[[n]] = cov_gdp
    cov_gdppc_list[[n]] = cov_gdppc
    cov_population_list[[n]] = cov_population
    n = n +1
  }
  
  # if(only_biggest_cluster) {
  #   i = 1
  #   edges_pos = data.frame(which(target_matrices[[i]]== 1, arr.ind = T))
  #   edges_pos$sign = 1
  #   edges_neg = data.frame(which(target_matrices[[i]]== -1, arr.ind = T))
  #   edges_neg$sign = -1
  #   
  # # edges_neg$val = factor(edges_neg$val)
  #   n_actors = nrow(target_matrices[[i]])
  #   graph = tbl_graph(edges = rbind(edges_pos, edges_neg),nodes = data.frame(name = 1:n_actors),directed = F)
  #   Isolated = which(degree(graph)==0)
  #   graph = delete.vertices(graph, Isolated)
  #   ggsigned_new(graph)
  #   
  #   
  #   
  #   components <- igraph::clusters(graph, mode="weak")
  #   biggest_cluster_id <- which.max(components$csize)
  #   vert_ids <- V(graph)[components$membership != biggest_cluster_id]
  #   include[vert_ids] = F
  # }
  return(list(target_matrix = target_matrices,
              cov_repetition = cov_repetition_list,
              cov_repetition_pos = cov_repetition_pos_list,
              cov_repetition_neg = cov_repetition_neg_list,
              cov_repetition_neut = cov_repetition_neut_list,
              cov_contiguity = cov_contiguity_list, 
              cov_distance = cov_distance_list, 
              cov_polity = cov_polity_list, 
              cov_gdp = cov_gdp_list, 
              cov_gdppc = cov_gdppc_list, 
              cov_cinc = cov_cinc_list,
              cov_population = cov_population_list, 
              common_friends = common_friends_list, 
              common_enemies = common_enemies_list))
}


generate_table = function(digits = 3, ...){
  mod_objects <- list(...)
  # mod_objects <- c(as.list(environment()), list(...))
  coef_names = unique(unlist(lapply(mod_objects, FUN = function(x){names(x$coefficients)})))
  # Augment the coefficients of all  
  mod_objects_new = lapply(mod_objects, FUN = function(x, coef_names){
    x$coefficients_new = numeric(length = length(coef_names))
    names(x$coefficients_new) = coef_names
    x$coefficients_new[match(names(x$coefficients), names(x$coefficients_new))] = 
      round(x$coefficients, digits = digits)
    x$sd_new = numeric(length = length(coef_names))
    names(x$sd_new) = coef_names
    x$sd_new[match(names(x$coefficients), names(x$coefficients_new))] = 
      round(sqrt(diag(x$var)), digits = digits)
    x$z_new = round(x$coefficients_new/x$sd_new, digits = digits)
    x$lower_ci = round(x$coefficients_new - qnorm(p = 1-0.025)*x$sd_new, digits = digits)
    x$upper_ci = round(x$coefficients_new + qnorm(p = 1-0.025)*x$sd_new, digits = digits)
    x$ci = paste0("[",x$lower_ci,",",x$upper_ci,"]")
    # x$coefficients_new = paste0(x$coefficients_new, x$sd_new)
    x$ci[x$ci == "[0,0]"] = ""
    x$coefficients_new[x$coefficients_new == 0] = "-"
    x$z_new[is.nan(x$z_new)] = ""
    res = matrix(data = "", nrow = length(coef_names)*2 +1,ncol = 2)
    even_rows = seq(2,length(coef_names)*2,by = 2)
    uneven_rows = seq(1,length(coef_names)*2,by = 2)
    res[uneven_rows,1] = x$coefficients_new
    res[even_rows,1] =    x$ci
    res[uneven_rows,2] =    x$z_new
    res[length(coef_names)*2 +1,1] = round(AIC(x), digits = digits)
    return(res)
  }, coef_names = coef_names)
  
  df = (do.call("cbind", mod_objects_new))
  colnames(df) = rep(c("Coef./CI", "Z Val."), times = length(mod_objects_new))
  uneven_rows = seq(1,length(coef_names)*2,by = 2)
  coef_names_aug = character(length(coef_names)*2 +1)
  coef_names_aug[uneven_rows] = coef_names
  coef_names_aug[length(coef_names)*2 +1] = "AIC"
  return(cbind(coef_names_aug,df))
}

t_generate_table = function(digits = 3, ...){
  mod_objects <- list(...)
  # mod_objects <- c(as.list(environment()), list(...))
  coef_names = unique(unlist(lapply(mod_objects, FUN = function(x){names(x$coefficients)})))
  # Augment the coefficients of all
  mod_objects_new = lapply(mod_objects, FUN = function(x, coef_names){
    x$coefficients_new = numeric(length = length(coef_names))
    names(x$coefficients_new) = coef_names
    x$coefficients_new[match(names(x$coefficients), names(x$coefficients_new))] =
      round(x$coefficients, digits = digits)
    x$sd_new = numeric(length = length(coef_names))
    names(x$sd_new) = coef_names
    x$sd_new[match(names(x$coefficients), names(x$coefficients_new))] =
      round(sqrt(diag(x$var)), digits = digits)
    x$z_new = round(x$coefficients_new/x$sd_new, digits = digits)
    x$lower_ci = round(x$coefficients_new - qnorm(p = 1-0.025)*x$sd_new, digits = digits)
    x$upper_ci = round(x$coefficients_new + qnorm(p = 1-0.025)*x$sd_new, digits = digits)
    x$ci = paste0("[",x$lower_ci,",",x$upper_ci,"]")
    # x$coefficients_new = paste0(x$coefficients_new, x$sd_new)
    x$ci[x$ci == "[0,0]"] = ""
    x$coefficients_new[x$coefficients_new == 0] = "-"
    x$z_new[is.nan(x$z_new)] = ""
    res = matrix(data = "", nrow = length(coef_names) +1,ncol = 2)
    # even_rows = seq(2,length(coef_names)*2,by = 2)
    # uneven_rows = seq(1,length(coef_names)*2,by = 2)
    res[1:length(coef_names),1] = x$coefficients_new
    res[1:length(coef_names),2] =    x$ci
    # res[uneven_rows,2] =    x$z_new
    res[length(coef_names) +1,1] = round(AIC(x), digits = digits)
    return(res)
  }, coef_names = coef_names)

  df = (do.call("cbind", mod_objects_new))
  colnames(df) = rep(c("Coef.", "CI"), times = length(mod_objects_new))
  # uneven_rows = seq(1,length(coef_names)*2,by = 2)
  coef_names_aug = character(length(coef_names) +1)
  coef_names_aug[1:length(coef_names)] = coef_names
  coef_names_aug[length(coef_names) +1] = "AIC"
  return(cbind(coef_names_aug,df))
}

t_plot_model_assessment = function(mod_ass){
  plots_list = list()
  for(i in 1:length(mod_ass$df_ese_list)){
    plots= list()
    # 1. Edgewise-shared Enemies plots 
    setDT(mod_ass$df_ese_list[[i]])
    plot_data = melt.data.table(mod_ass$df_ese_list[[i]])
    
    line_data = data.table(x = as.numeric(names(mod_ass$count_ese[[i]])), y = as.numeric(mod_ass$count_ese[[i]]))
    line_data$x = line_data$x +1 
    if(max(line_data$x)<(max(as.numeric(as.character(plot_data$variable)))+1)){
      line_data = rbind(line_data,data.table(x = as.numeric(x = (max(line_data$x)+1):(max(as.numeric(as.character(plot_data$variable)))+1)), y = 0))
    }
    
    plots[[1]] = ggplot(data = plot_data) +
      geom_boxplot(aes(x = (variable), y = (value))) +
      geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
      theme_pubr() + 
      xlab("Edgewise-Shared Enemies") +
      ylab("Number of observations") +
      scale_y_continuous(trans='log1p', breaks = c(0,2,10,50,250, 1000))
    # 2. Edgewise-shared Partners plots 
    setDT(mod_ass$df_esp_list[[i]])
    plot_data = melt.data.table(mod_ass$df_esp_list[[i]])
    line_data = data.table(x = as.numeric(names(mod_ass$count_esp[[i]])), y = as.numeric(mod_ass$count_esp[[i]]))
    line_data$x = line_data$x +1 
    if(max(line_data$x)<(max(as.numeric(as.character(plot_data$variable)))+1)){
      line_data = rbind(line_data,data.table(x = as.numeric(x = (max(line_data$x)+1):(max(as.numeric(as.character(plot_data$variable)))+1)), y = 0))
    }
    plots[[2]] = ggplot(data = plot_data) +
      geom_boxplot(aes(x = (variable), y = value)) +
      geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
      theme_pubr() + 
      xlab("Edgewise-Shared Friends") +
      ylab("Number of observations") +
      scale_y_continuous(trans='log1p', breaks = c(0,2,10,50,250, 600))
    
    # 3. Positive degree distribution
    setDT(mod_ass$df_dp_list[[i]])
    plot_data = melt.data.table(mod_ass$df_dp_list[[i]])
    line_data = data.table(x = as.numeric(names(mod_ass$count_dp[[i]])), y = as.numeric(mod_ass$count_dp[[i]]))
    line_data$x = line_data$x +1 
    if(max(line_data$x)<(max(as.numeric(as.character(plot_data$variable)))+1)){
      line_data = rbind(line_data,data.table(x = as.numeric(x = (max(line_data$x)+1):(max(as.numeric(as.character(plot_data$variable)))+1)), y = 0))
    }
    plots[[3]] = ggplot(data = plot_data) +
      geom_boxplot(aes(x = (variable), y = value)) +
      geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
      theme_pubr() + 
      xlab("Positive Degree") +
      ylab("Number of observations")
    
    # 4. Negative degree distribution
    setDT(mod_ass$df_dn_list[[i]])
    plot_data = melt.data.table(mod_ass$df_dn_list[[i]])
    line_data = data.table(x = as.numeric(names(mod_ass$count_dn[[i]])), y = as.numeric(mod_ass$count_dn[[i]]))
    line_data$x = line_data$x +1 
    if(max(line_data$x)<(max(as.numeric(as.character(plot_data$variable)))+1)){
      line_data = rbind(line_data,data.table(x = as.numeric(x = (max(line_data$x)+1):(max(as.numeric(as.character(plot_data$variable)))+1)), y = 0))
    }
    plots[[4]] = ggplot(data = plot_data) +
      geom_boxplot(aes(x = (variable), y = value)) +
      geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
      theme_pubr() + 
      xlab("Negative Degree") +
      ylab("Number of observations")+
      scale_y_continuous(trans='log1p', breaks = c(0,2,10,50,150))
    
    
    # plot_data = melt.data.table(data.table(mod_ass$df_igraph[[i]][,-1]))
    # line_data = data.table(x = 1:length(mod_ass$igraph_data_obs[[i]][-1]), y = as.numeric(mod_ass$igraph_data_obs[[i]][-1]))
    # # line_data$x = line_data$x +1 
    # plots[[5]] =ggplot(data = plot_data) +
    #   geom_violin(aes(x = (variable), y = value)) +
    #   geom_point(data = line_data, aes(x = x, y = y), color = "red") + 
    #   theme_pubr() + 
    #   xlab("Signed triad census") +
    #   ylab("Number of triagles")
    plots_list[[i]] = plots
  }
  plots_list = do.call(rbind, plots_list)
  
  return(plots_list)
}


plot_model_assessment = function(mod_ass){
  plots= list()
  # 1. Edgewise-shared Enemies plots 
  setDT(mod_ass$df_ese)
  plot_data = melt.data.table(mod_ass$df_ese)
  line_data = data.table(x = as.numeric(names(mod_ass$count_ese)), y = as.numeric(mod_ass$count_ese))
  line_data$x = line_data$x +1
  if(max(line_data$x-1)<max(as.numeric(as.character(plot_data$variable)))){
    x = (max(line_data$x):(max(as.numeric(as.character(plot_data$variable)))))+1
    line_data = rbind(line_data,data.table(x = x, y = 0))
  }
  
  
  plots[[1]] = ggplot(data = plot_data) +
    geom_boxplot(aes(x = (variable), y = value)) +
    geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
    theme_pubr() + 
    xlab("Edgewise-Shared Enemies") +
    ylab("Number of observations")
  # 2. Edgewise-shared Partners plots 
  setDT(mod_ass$df_esp)
  plot_data = melt.data.table(mod_ass$df_esp)
  line_data = data.table(x = as.numeric(names(mod_ass$count_esp)), y = as.numeric(mod_ass$count_esp))
  line_data$x = line_data$x +1 
  if(max(line_data$x-1)<max(as.numeric(as.character(plot_data$variable)))){
    x = (max(line_data$x):(max(as.numeric(as.character(plot_data$variable)))))+1
    line_data = rbind(line_data,data.table(x = x, y = 0))
  }
  
  plots[[2]] = ggplot(data = plot_data) +
    geom_boxplot(aes(x = (variable), y = value)) +
    geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
    theme_pubr() + 
    xlab("Edgewise-Shared Friends") +
    ylab("Number of observations")
  # 3. Positive degree distribution
  setDT(mod_ass$df_dp)
  plot_data = melt.data.table(mod_ass$df_dp)
  line_data = data.table(x = as.numeric(names(mod_ass$count_dp)), y = as.numeric(mod_ass$count_dp))
  line_data$x = line_data$x +1 
  if(max(line_data$x-1)<max(as.numeric(as.character(plot_data$variable)))){
    x = (max(line_data$x):(max(as.numeric(as.character(plot_data$variable)))))+1
    line_data = rbind(line_data,data.table(x = x, y = 0))
  }
  
  plots[[3]] = ggplot(data = plot_data) +
    geom_boxplot(aes(x = (variable), y = value)) +
    geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
    theme_pubr() + 
    xlab("Positive Degree") +
    ylab("Number of observations")
  # 4. Negative degree distribution
  setDT(mod_ass$df_dn)
  plot_data = melt.data.table(mod_ass$df_dn)
  line_data = data.table(x = as.numeric(names(mod_ass$count_dn)), y = as.numeric(mod_ass$count_dn))
  line_data$x = line_data$x +1 
  if(max(line_data$x-1)<max(as.numeric(as.character(plot_data$variable)))){
    x = (max(line_data$x):(max(as.numeric(as.character(plot_data$variable)))))+1
    line_data = rbind(line_data,data.table(x = x, y = 0))
  }
  
  plots[[4]] = ggplot(data = plot_data) +
    geom_boxplot(aes(x = (variable), y = value)) +
    geom_line(data = line_data, aes(x = x, y = y), color = "red") + 
    theme_pubr() + 
    xlab("Negative Degree") +
    ylab("Number of observations")
  
  plot_data = melt.data.table(as.data.table(mod_ass$df_igraph[,-1]))
  line_data = data.table(x = 1:length(mod_ass$igraph_data_obs[-1]), y = as.numeric(mod_ass$igraph_data_obs[-1]))
  # line_data$x = line_data$x +1 
  # plots[[5]] =ggplot(data = plot_data) +
  #   geom_violin(aes(x = (variable), y = value)) +
  #   geom_point(data = line_data, aes(x = x, y = y), color = "red") + 
  #   theme_pubr() + 
  #   xlab("Signed triad census") +
  #   ylab("Number of triagles")
  
  return(plots)
}
