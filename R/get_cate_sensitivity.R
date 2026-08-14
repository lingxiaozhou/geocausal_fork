#' Sensitivity analysis for a Hajek CATE estimator
#'
#' @description Evaluates the sensitivity of a Hajek conditional average
#' treatment effect (CATE) estimator to bounded departures from the fitted
#' treatment-assignment mechanism.
#' @param obs observed density
#' @param cf1 counterfactual density 1
#' @param cf2 counterfactual density 2
#' @param treat column of a hyperframe that summarizes treatment data. In the form of `hyperframe$column`.
#' @param pixel_count_out column of a hyperframe that summarizes the number of outcome events in each pixel
#' @param lag integer that specifies lags to calculate causal estimates.
#' @param trunc_level the level of truncation for the weights (0-1).
#' @param time_after whether to include one unit time difference between treatment and outcome. By default = TRUE
#' @param entire_window owin object (the entire region of interest)
#' @param E_mat optional covariance matrix (excluding the intercept) for the effect modifier. If provided, then the regression model will be based on this matrix. If `intercept = TRUE`, then a column of 1 will be add to `E_mat`.
#' @param em treat column of a hyperframe that summarizes the effect modifier data. In the form of `hyperframe$column`. It can be NULL if E_mat is provided.
#' @param nbase number of bases for splines
#' @param spline_type type of splines. Either `"ns"` or `"bs"`. 
#' @param intercept whether to include intercept in the regression model. Default is TRUE. 
#' @param eval_values a vector of values of the effect modifier for which CATE will be evaluated. Default is a `seq(a,b,length.out=20)` where `a` and `b` are minimum and maximum values of the effect modifier.
#' @param eval_mat evaluated spline basis (excluding the intercept) matrix at `eval_values`.  If `intercept = TRUE`, then a column of 1 will be add to `eval_mat`.
#' @param save_weights retained for compatibility; currently has no effect.
#' @param gamma_vals numeric vector of sensitivity parameter values greater than
#' or equal to one. The default is `seq(1.05, 1.2, 0.05)`.
#' @param ... arguments passed onto the function 
#' 
#' @returns A list containing pointwise feasible lambda intervals, the lambda
#' values and slack values closest to attaining a zero effect, and indicators of
#' zero attainability. Corresponding objects with the `_beta` suffix contain
#' coefficient-level results. The list also contains the first attainable gamma
#' values, summaries by evaluation point and coefficient, and `eval_values`.
#'
#' @details `E_mat` should be a matrix or array of dimensions \eqn{n} by \eqn{m} where \eqn{n} is the product of image dimensions and number of time period,
#' and \eqn{m} is `nbase`-`intercept`. If you want to construct your own covariate matrix `E_mat`, you should use `get_em_vec()` to convert
#' the effect modifier (usually a column of a hyperframe) to a vector, and then construct the spline basis based on the vector. The covariate matrix `E_mat` should not
#' include the intercept column. For each value in `gamma_vals`, the function uses linear programming to determine whether a zero effect is attainable under the corresponding bounds.
#'
#' @importFrom Rglpk Rglpk_solve_LP
#' @export

get_cate_sensitivity <- function(obs, cf1, cf2, treat, pixel_count_out,lag, trunc_level=0.95, time_after=TRUE,entire_window = NULL,
                                 em = NULL,E_mat = NULL,
                                 nbase = 6, spline_type = "ns",intercept = TRUE,
                                 eval_values = NULL, eval_mat = NULL,save_weights = TRUE,gamma_vals = seq(1.05,1.2,0.05),...) {
  # pixel_ratio <- 8451/(128*128)
  chisq_stat <- NULL
  p.value <- NULL
  total_effect <- NULL
  mean_effect <- NULL
  E_mat_provided <- !is.null(E_mat)
  weights <- NULL
  
  # --------------------------check the format of the arguments-------------------
  is_positive_integer_within_range <- function(x, min_val, max_val) {
    all(x %% 1 == 0 & x > 0 & x >= min_val & x <= max_val)
  }
  
  
  if(is.null(em) & is.null(E_mat)){
    stop("Both em and E_mat are null.")
  }
  
  if(is.null(eval_mat) & !is.null(E_mat)){
    stop("eval_mat is missing when E_mat is provided")
  }
  
  if(is.null(eval_values) & !is.null(eval_mat)){
    stop("eval_values is missing when eval_mat is provided")
  }
  
  if(is.null(nbase)){
    nbase <- ncol(E_mat)+intercept
  }
  
  
  
  #-------------------------get the weighted surfaces------------------------------------
  message("Get weighted surfaces... \n")
  # CF1
  estimates_1 <- get_weighted_surf(obs_dens = obs,
                                   cf_dens = cf1,
                                   treatment_data = treat,
                                   smoothed_outcome = pixel_count_out,
                                   mediation = FALSE, cate = TRUE,
                                   obs_med_log_sum_dens = NA,
                                   cf_med_log_sum_dens = NA,
                                   lag = lag, entire_window = entire_window,
                                   time_after,
                                   truncation_level = trunc_level)
  
  ## CF2
  estimates_2 <- get_weighted_surf(obs_dens = obs,
                                   cf_dens = cf2,
                                   treatment_data = treat,
                                   smoothed_outcome = pixel_count_out,
                                   mediation = FALSE, cate = TRUE,
                                   obs_med_log_sum_dens = NA,
                                   cf_med_log_sum_dens = NA,
                                   lag = lag, entire_window = entire_window,
                                   time_after,
                                   truncation_level = trunc_level)
  
  
  dimyx <- dim(estimates_1$weighted_surface_arr_haj)[c(1,2)]
  message('In the analysis, pixel grid dimension is ', dimyx[1],"x",dimyx[2],"\n")
  time_points <- dim(estimates_1$weighted_surface_arr_haj)[3]
  # dimyx <- dim(estimates$av_surface[[1]])[c(1,2)]
  # time_points <- dim(estimates$av_surface[[1]])[3]
  
  
  
  #------------------------construct the basis matrix if not provided----------------------
  if(is.null(E_mat)){
    message("Generate spline basis...\n")
    if("im"%in%class(em[[1]])){
      em <- lapply(1:length(em), function(x) as.matrix(as.im(em[[x]],dimyx = dimyx)))
    }else if(!"array"%in%class(em[[1]])){
      stop("em is not a list of im or 2D arrays")
    }
    
    if(length(em)==1){
      em <- lapply(1:time_points, function(x) em[[1]])
    }else{
      em <- em[1:time_points]
    }
    
    if(spline_type=="ns"){
      E_mat <- splines::ns(unlist(em[1:time_points]),df = nbase-intercept, ...)
    }
    
    if(spline_type=="bs"){
      E_mat <- splines::bs(unlist(em[1:time_points]),df = nbase-intercept,...)
    }
    
  }
  knots <- attr(E_mat, "knots")
  
  if(is.null(eval_values)){
    eval_values <- seq(min(unlist(em),na.rm = TRUE),max(unlist(em),na.rm = TRUE),length.out = 20)
  }
  if(is.null(eval_mat)){
    eval_mat <- predict(E_mat,newx = eval_values)
  }
  
  
  if(intercept){
    E_mat <- cbind(1,E_mat)
    eval_mat <- cbind(1,eval_mat)
  }
  
  
  #---------------------------------fit the model-----------------------------------
  message("Fit the model...\n")
  # form a dataframe for the estimates and covariates
  df <- cbind.data.frame(E = E_mat,estimates1 = c(estimates_1$weighted_surface_arr_haj),estimates2 = c(estimates_2$weighted_surface_arr_haj))
  weights <- rbind(estimates_1$weights,estimates_2$weights)
  # df <- cbind.data.frame(E = E_mat,estimates1 = c(estimates$av_surface_haj[[1]]),estimates2 = c(estimates$av_surface_haj[[6]]))
  # weights <- rbind(estimates$weights[1,],estimates$weights[6,])
  
  
  df$est <- df$estimates2-df$estimates1
  df$time <- sort(rep(1:time_points,prod(dimyx)))
  df <- na.omit(df) # remove all the na values
  mean_effect <- mean(df$est,na.rm = TRUE)
  total_effect <- mean_effect*sum(df$time==1)
  
  
  valid_time_points <- 0
  p <- ncol(E_mat) # number of covariates
  
  beta_arr <- array(NA,c(time_points,2*p+2))
  
  
  
  for (tt in 1:time_points) {
    df_tt <- df[df$time==tt,-ncol(df)]
    
    # Fit the spline
    tryCatch({
      Z <- as.matrix(df_tt[,-c(ncol(df_tt)-0:2)])
      Q <- solve(t(Z)%*%Z)
      
      beta_arr[tt,] <- c(Q%*%t(Z)%*%df_tt$estimates1,Q%*%t(Z)%*%df_tt$estimates2,weights[,tt])
      
      
      valid_time_points <- valid_time_points+1
    }, error=function(e){
      message("Warning: covariate matrix is computationally singular for time period ",tt)
    })
    
  }
  
  #---------------------------sensitivity analysis-----------------------
  message("Perform sensitivity analysis...\n")
  
  beta1 <- beta_arr[, 1:p]
  beta2 <- beta_arr[, (p + 1):(2 * p)]
  
  npoints <- nrow(eval_mat) 
  npoints_beta <- p-intercept
  
  # store results
  lambda_interval <- array(NA, c(npoints, 2, length(gamma_vals)))
  dimnames(lambda_interval) <- list(eval_values, c("lambda_lo", "lambda_hi"), gamma_vals)
  
  best_lambda <- array(NA, c(npoints, length(gamma_vals)))
  dimnames(best_lambda) <- list(eval_values, gamma_vals)
  
  best_slack <- array(NA, c(npoints, length(gamma_vals)))
  dimnames(best_slack) <- list(eval_values, gamma_vals)
  
  zero_attainable <- array(NA, c(npoints, length(gamma_vals)))
  dimnames(zero_attainable) <- list(eval_values, gamma_vals)
  
  lambda_interval_beta <- array(NA, c(npoints_beta, 2, length(gamma_vals)))
  dimnames(lambda_interval_beta) <- list(colMeans(beta2-beta1)[-1], c("lambda_lo", "lambda_hi"), gamma_vals)
  
  best_lambda_beta <- array(NA, c(npoints_beta, length(gamma_vals)))
  dimnames(best_lambda_beta) <- list(colMeans(beta2-beta1)[-1], gamma_vals)
  
  best_slack_beta <- array(NA, c(npoints_beta, length(gamma_vals)))
  dimnames(best_slack_beta) <- list(colMeans(beta2-beta1)[-1], gamma_vals)
  
  zero_attainable_beta <- array(NA, c(npoints_beta, length(gamma_vals)))
  dimnames(zero_attainable_beta) <- list(colMeans(beta2-beta1)[-1], gamma_vals)
  
  # optional: stop once overall zero-attainability occurs
  first_gamma_zero_attainable <- NA
  first_gamma_zero_attainable_beta <- NA
  
  #------------------------------------------------------------
  # helper: single-ratio min/max by Charnes-Cooper LP
  #------------------------------------------------------------
  solve_ratio_bounds <- function(obj_num, obj_den, lower_rho, upper_rho) {
    Tlen <- length(obj_num)
    D_mat <- diag(Tlen)
    
    mat <- rbind(D_mat, -D_mat)
    dir <- rep("<=", 2 * Tlen)
    rhs <- c(upper_rho, -lower_rho)
    
    new_obj <- c(obj_num, 0)
    
    new_mat <- cbind(mat, -rhs)
    new_dir <- rep("<=", 2 * Tlen)
    new_rhs <- rep(0, 2 * Tlen)
    
    new_mat <- rbind(new_mat, c(obj_den, 0))
    new_dir <- c(new_dir, "==")
    new_rhs <- c(new_rhs, 1)
    
    new_mat <- rbind(new_mat, c(rep(0, ncol(new_mat) - 1), 1))
    new_dir <- c(new_dir, ">=")
    new_rhs <- c(new_rhs, 0)
    
    lp_max <- Rglpk_solve_LP(
      obj = new_obj, mat = new_mat, dir = new_dir, rhs = new_rhs,
      max = TRUE, verbose = FALSE
    )
    lp_min <- Rglpk_solve_LP(
      obj = new_obj, mat = new_mat, dir = new_dir, rhs = new_rhs,
      max = FALSE, verbose = FALSE
    )
    
    c(min = lp_min$optimum, max = lp_max$optimum)
  }
  
  #------------------------------------------------------------
  # helper: slack LP for fixed lambda
  #------------------------------------------------------------
  solve_lambda_slack <- function(lambda, A, B, C, D, lower_rho, upper_rho, eps_denom = 1e-8) {
    Tlen <- length(A)
    
    # variables: rho_1,...,rho_T,s
    obj <- c(rep(0, Tlen), 1)
    
    v1 <- A - lambda * B
    v2 <- C - lambda * D
    
    mat <- rbind(
      c(v1, -1),
      c(-v1, -1),
      c(v2, -1),
      c(-v2, -1),
      cbind(diag(Tlen), rep(0, Tlen)),
      cbind(-diag(Tlen), rep(0, Tlen)),
      c(B, 0),
      c(D, 0),
      c(rep(0, Tlen), 1)
    )
    
    dir <- c(
      "<=", "<=", "<=", "<=",
      rep("<=", Tlen),
      rep("<=", Tlen),
      ">=", ">=",
      ">="
    )
    
    rhs <- c(
      0, 0, 0, 0,
      upper_rho,
      -lower_rho,
      eps_denom, eps_denom,
      0
    )
    
    sol <- Rglpk_solve_LP(
      obj = obj, mat = mat, dir = dir, rhs = rhs,
      max = FALSE, verbose = FALSE
    )
    
    # if solver does not end properly, return NA
    if (sol$status != 0) {
      return(list(
        status = sol$status,
        optimum = Inf,
        solution = NULL,
        rho = NULL,
        s = Inf,
        residual1 = NA,
        residual2 = NA,
        max_residual = NA
      ))
    }
    
    # get rho and s
    rho_sol <- sol$solution[1:Tlen]
    s_sol <- sol$solution[Tlen + 1]
    
    # check residual
    residual1 <- sum((A - lambda * B) * rho_sol)
    residual2 <- sum((C - lambda * D) * rho_sol)
    max_residual <- max(abs(residual1), abs(residual2))
    
    list(
      status = sol$status,
      optimum = sol$optimum,
      solution = sol$solution,
      rho = rho_sol,
      s = s_sol,
      residual1 = residual1,
      residual2 = residual2,
      max_residual = max_residual
    )
  }
  
  #------------------------------------------------------------
  # helper: adaptive search over lambda
  #------------------------------------------------------------
  adaptive_lambda_search <- function(lambda_lo, lambda_hi,
                                     A, B, C, D,
                                     lower_rho, upper_rho,
                                     tol_slack = 1e-6,
                                     grid_init = 15,
                                     max_refine = 5) {
    if (lambda_lo > lambda_hi) {
      return(list(attainable = FALSE, best_lambda = NA, best_slack = Inf))
    }
    
    current_lo <- lambda_lo
    current_hi <- lambda_hi
    best_lambda <- NA
    best_slack <- Inf
    
    for (iter in seq_len(max_refine)) {
      lambda_grid <- seq(current_lo, current_hi, length.out = grid_init)
      slack_vals <- rep(NA_real_, length(lambda_grid))
      
      best_sol <- NULL
      
      for (j in seq_along(lambda_grid)) {
        sol <- solve_lambda_slack(
          lambda = lambda_grid[j],
          A = A, B = B, C = C, D = D,
          lower_rho = lower_rho, upper_rho = upper_rho
        )
        slack_vals[j] <- sol$optimum
      }
      
      j_best <- which.min(slack_vals)
      
      if (slack_vals[j_best] < best_slack) {
        best_slack <- slack_vals[j_best]
        best_lambda <- lambda_grid[j_best]
        
        best_sol <- solve_lambda_slack(
          lambda = best_lambda,
          A = A, B = B, C = C, D = D,
          lower_rho = lower_rho, upper_rho = upper_rho
        )
      }
      
      if (best_slack <= tol_slack) {
        return(list(
          attainable = TRUE,
          best_lambda = best_lambda,
          best_slack = best_slack,
          residual1 = best_sol$residual1,
          residual2 = best_sol$residual2,
          max_residual = best_sol$max_residual,
          rho = best_sol$rho,
          status = best_sol$status
        ))
      }
      
      left_idx <- max(1, j_best - 1)
      right_idx <- min(length(lambda_grid), j_best + 1)
      
      current_lo <- lambda_grid[left_idx]
      current_hi <- lambda_grid[right_idx]
      
      if ((current_hi - current_lo) < 1e-10) break
    }
    
    list(
      attainable = (best_slack <= tol_slack),
      best_lambda = best_lambda,
      best_slack = best_slack
    )
  }
  
  #------------------------------------------------------------
  # main gamma loop for the CATE
  #------------------------------------------------------------
  for (gg in seq_along(gamma_vals)) {
    
    this_gamma <- gamma_vals[gg]
    cat("Start gamma:", this_gamma, "\n")
    
    lower_rho <- rep((1 / this_gamma)^lag, time_points)
    upper_rho <- rep(this_gamma^lag, time_points)
    
    point_attainable <- logical(npoints)
    
    for (point in seq_len(npoints)) {
      
      # ratio 1
      A <- c(beta1 %*% eval_mat[point, ])
      B <- weights[1, ] / mean(weights[1, ])
      
      # ratio 2
      C <- c(beta2 %*% eval_mat[point, ])
      D <- weights[2, ] / mean(weights[2, ])
      
      # Step 1: get lambda search interval from separate ratio bounds
      bounds1 <- solve_ratio_bounds(A, B, lower_rho, upper_rho)
      bounds2 <- solve_ratio_bounds(C, D, lower_rho, upper_rho)
      
      lambda_lo <- max(bounds1["min"], bounds2["min"])
      lambda_hi <- min(bounds1["max"], bounds2["max"])
      
      lambda_interval[point, "lambda_lo", gg] <- lambda_lo
      lambda_interval[point, "lambda_hi", gg] <- lambda_hi
      
      # no overlap => zero not attainable at this point
      if (lambda_lo > lambda_hi) {
        zero_attainable[point, gg] <- FALSE
        best_lambda[point, gg] <- NA
        best_slack[point, gg] <- Inf
        next
      }
      
      # Step 2: adaptive lambda search
      out <- adaptive_lambda_search(
        lambda_lo = lambda_lo,
        lambda_hi = lambda_hi,
        A = A, B = B, C = C, D = D,
        lower_rho = lower_rho,
        upper_rho = upper_rho,
        tol_slack = 1e-6,
        grid_init = 15,
        max_refine = 5
      )
      
      zero_attainable[point, gg] <- out$attainable
      best_lambda[point, gg] <- out$best_lambda
      best_slack[point, gg] <- out$best_slack
      
      point_attainable[point] <- isTRUE(out$attainable)
    }
    
    # stop early once zero is attainable for at least one point
    if (all(point_attainable)) {
      first_gamma_zero_attainable <- this_gamma
      cat("Stopping early: zero is attainable at gamma =", this_gamma, "\n")
      break
    }
  }
  
  #------------------------------------------------------------
  # main gamma loop for the beta
  #------------------------------------------------------------
  for (gg in seq_along(gamma_vals)) {
    
    this_gamma <- gamma_vals[gg]
    cat("Start gamma (beta):", this_gamma, "\n")
    
    lower_rho <- rep((1 / this_gamma)^lag, time_points)
    upper_rho <- rep(this_gamma^lag, time_points)
    
    point_attainable <- logical(npoints_beta)
    
    for (point in seq_len(npoints_beta)) {
      
      
      A <- c(beta1[,point + intercept])
      
      B <- weights[1, ] / mean(weights[1, ])
      
      C <- c(beta2[,point + intercept])
      
      D <- weights[2, ] / mean(weights[2, ])
      
      # Step 1: get lambda search interval from separate ratio bounds
      bounds1 <- solve_ratio_bounds(A, B, lower_rho, upper_rho)
      bounds2 <- solve_ratio_bounds(C, D, lower_rho, upper_rho)
      
      lambda_lo <- max(bounds1["min"], bounds2["min"])
      lambda_hi <- min(bounds1["max"], bounds2["max"])
      
      lambda_interval_beta[point, "lambda_lo", gg] <- lambda_lo
      lambda_interval_beta[point, "lambda_hi", gg] <- lambda_hi
      
      # no overlap => zero not attainable at this point
      if (lambda_lo > lambda_hi) {
        zero_attainable_beta[point, gg] <- FALSE
        best_lambda_beta[point, gg] <- NA
        best_slack_beta[point, gg] <- Inf
        next
      }
      
      # Step 2: adaptive lambda search
      out <- adaptive_lambda_search(
        lambda_lo = lambda_lo,
        lambda_hi = lambda_hi,
        A = A, B = B, C = C, D = D,
        lower_rho = lower_rho,
        upper_rho = upper_rho,
        tol_slack = 1e-6,
        grid_init = 15,
        max_refine = 5
      )
      
      zero_attainable_beta[point, gg] <- out$attainable
      best_lambda_beta[point, gg] <- out$best_lambda
      best_slack_beta[point, gg] <- out$best_slack
      
      point_attainable[point] <- isTRUE(out$attainable)
    }
    
    # stop early once zero is attainable for at least one point
    if (all(point_attainable)) {
      first_gamma_zero_attainable_beta <- this_gamma
      cat("Stopping early: zero is attainable at gamma =", this_gamma, "\n")
      break
    }
  }
  
  #------------------------------------------------------------
  # summary
  #------------------------------------------------------------
  find_first_gamma <- function(zero_attainable, gamma_vals) {
    n_x <- nrow(zero_attainable)
    
    pointwise_first_gamma <- sapply(seq_len(n_x), function(i) {
      idx <- which(zero_attainable[i, ] %in% TRUE)
      if (length(idx) > 0) gamma_vals[min(idx)] else NA
    })
    
    overall_first_gamma <- {
      idx <- which(apply(zero_attainable, 2, function(z) any(z %in% TRUE)))
      if (length(idx) > 0) gamma_vals[min(idx)] else NA
    }
    
    list(
      overall_first_gamma = overall_first_gamma,
      pointwise_first_gamma = pointwise_first_gamma
    )
  }
  
  gamma_summary <- find_first_gamma(zero_attainable, gamma_vals)
  gamma_summary_beta <- find_first_gamma(zero_attainable_beta, gamma_vals)
  
  res <- list(
    lambda_interval = lambda_interval,
    best_lambda = best_lambda,
    best_slack = best_slack,
    zero_attainable = zero_attainable,
    first_gamma_zero_attainable = first_gamma_zero_attainable,
    gamma_summary = gamma_summary,
    lambda_interval_beta = lambda_interval_beta,
    best_lambda_beta = best_lambda_beta,
    best_slack_beta = best_slack_beta,
    zero_attainable_beta = zero_attainable_beta,
    first_gamma_zero_attainable_beta = first_gamma_zero_attainable_beta,
    gamma_summary_beta = gamma_summary_beta,
    eval_values = eval_values
  )
  
  return(res)
  
}

