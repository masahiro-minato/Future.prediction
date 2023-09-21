install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain()
#check_cmdstan_toolchain(fix = TRUE)

install_cmdstan(cores=4, overwrite = TRUE)
set_cmdstan_path()