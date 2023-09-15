# コンパイル
mod <- cmdstan_model("./stan/bsts-AR-tvc-smo-pred-pois-periphe-mat_MF3parm.stan")
mod <- cmdstan_model("./stan/bsts-AR-tvc-rd-smo-pred-pois-periphe-mat_MF3parm.stan", cpp_options = list(stan_threads = TRUE))
