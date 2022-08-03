#' @useDynLib armbart
#' @importFrom Rcpp sourceCpp

# Calling the function in the R version
#' @export
r_bart <- function(x_train,
                   y,
                   x_test,
                   n_tree,
                   n_mcmc,
                   n_burn,
                   n_min_size,
                   tau, mu,
                   alpha, beta,
                   scale_boolean = TRUE,
                   K_bart = 2,
                   a_tau){


  # Saving a_min and b_max
  a_min <- NULL
  b_max <- NULL

  a_min <- min(y)
  b_max <- max(y)


  # Error of the matrix
  if(is.null(colnames(x_train)) || is.null(colnames(x_test)) ) {
    stop("Insert a valid NAMED matrix")
  }

  if(!is.vector(y)) {
    stop("Insert a valid y vector")
  }

  # Scale values
  if(scale_boolean) {
    # Normalizing y
    y_scale <- normalize_bart(y = y)

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu <- (4 * n_tree * (K_bart^2))
    nsigma <- naive_sigma(x = x_train,y = y_scale)
    d_tau <- rate_tau(x = x_train,y = y_scale,prob = 0.9,shape = a_tau)

  } else {

    # Not scaling the y
    y_scale <- y

    # Calculating \tau_{\mu} based on the scale of y
    # Need to change this value in case of non-scaling
    tau_mu <- (4 * n_tree * (K_bart^2))/((b_max-a_min)^2)

    nsigma <- naive_sigma(x = x_train,y = y_scale)
    d_tau <- rate_tau(x = x_train,y = y_scale,prob = 0.9,shape = a_tau)



  }

  # Getting the main function
  bart_obj <- bart(x_train = x_train,
                   y = y_scale,
                   x_test = x_test,
                   n_tree = n_tree,
                   n_mcmc = n_mcmc,
                   n_burn = n_burn,
                   n_min_size = n_min_size,
                   tau = tau, mu = mu,
                   tau_mu = tau_mu,
                   naive_sigma = nsigma, # naive sigma value
                   alpha = alpha, beta = beta,a_tau = a_tau, d_tau = d_tau)

  # Returning to the normal scale
  if(scale_boolean){
    bart_obj$y_train_hat_post<- unnormalize_bart(bart_obj$y_train_hat_post,
                                                 a = a_min, b = b_max)

    bart_obj$y_test_hat_post<- unnormalize_bart(bart_obj$y_test_hat_post,
                                                a = a_min, b = b_max)

    bart_obj$tau_post <- bart_obj$tau_post/((b_max-a_min)^2)
  }



  return(bart_obj)
}


zero_tau_prob_squared <- function(x, naive_tau_value, prob, shape) {

  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined
  return((stats::pgamma(naive_tau_value,
                        shape = shape,
                        rate = x) - (1 - prob))^2)
}

# Functions to find the zero for tau
zero_tau_prob <- function(x, naive_tau_value, prob, shape) {

  # Find the zero to the function P(tau < tau_ols) = 0.1, for a defined
  return(stats::pgamma(naive_tau_value,
                       shape = shape,
                       rate = x) - (1 - prob))
}


# Naive sigma_estimation
naive_sigma <- function(x,y){

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Adjusting the df
  df <- data.frame(x,y)
  colnames(df)<- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  # Getting sigma
  sigma <- stats::sigma(lm_mod)
  return(sigma)
}

# Naive tau_estimation
naive_tau <- function(x, y) {

  # Getting the valus from n and p
  n <- length(y)

  # Getting the value from p
  p <- ifelse(is.null(ncol(x)), 1, ncol(x))

  # Adjusting the df
  df <- data.frame(x,y)
  colnames(df)<- c(colnames(x),"y")

  # Naive lm_mod
  lm_mod <- stats::lm(formula = y ~ ., data =  df)

  # Getting sigma
  sigma <- stats::sigma(lm_mod)

  # Using naive tau
  # sigma <- sd(y)

  # Getting \tau back
  tau <- sigma^(-2)
  return(tau)
}


# Return rate parameter from the tau prior
rate_tau <- function(x, # X value
                     y, # Y value
                     prob = 0.9,
                     shape) {
  # Find the tau_ols
  tau_ols <- naive_tau(x = x,
                       y = y)

  # Getting the root
  min_root <-  try(stats::uniroot(f = zero_tau_prob, interval = c(1e-6, 10000),
                                  naive_tau_value = tau_ols,
                                  prob = prob, shape = shape)$root, silent = TRUE)

  if(inherits(min_root, "try-error")) {
    # Verifying the squared version
    min_root <- stats::optim(par = stats::runif(1), fn = zero_tau_prob_squared,
                             method = "L-BFGS-B", lower = 0,
                             naive_tau_value = tau_ols,
                             prob = prob, shape = shape)$par
  }
  return(min_root)
}
# Normalize BART function (Same way as theOdds code)
normalize_bart <- function(y) {

  # Defining the a and b
  a <- min(y)
  b <- max(y)

  # This will normalize y between -0.5 and 0.5
  y  <- (y - a)/(b - a) - 0.5
  return(y)
}

# Now a function to return everything back to the normal scale

unnormalize_bart <- function(z, a, b) {
  # Just getting back to the regular BART
  y <- (b - a) * (z + 0.5) + a
  return(y)
}

# Calculating a PI coverage
#' @export
pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){

  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)

  full_post_draw <- list()

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )

  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)

      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }

  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }

  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})

  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)

  return(pi_cov)
}
