# ライブラリー読込み
library(cmdstanr)
library(posterior)
library(rstan)
library(bayesplot)
library(ggthemes)
library(ggfortify)
library(ggmcmc)
library(gridExtra)
library(GGally)
library(ellipse)
library(extrafont)
library(Cairo)
library(KFAS)
library(scales)
library(openxlsx)
library(writexl)
library(remotes)
library(patchwork)
library(hexbin)
library(devtools)
library(parallel)
library(rstudioapi)
library(brms)
library(tidyverse)
library(Rcpp)
library(readxl)
library(cowplot)
library(ggsci)
library(gridGraphics)
library(magrittr)

# マルチコア対応
options(mc.cores = parallel::detectCores())
# Stanファイル自動保存
rstan_options(auto_write = TRUE)

# 関数の読込み
source( "./function/reg_coe.R",encoding = "utf8" )
source( "./function/ZIP_calc.R",encoding = "utf8" )
source( "./function/Dglm_poisson_multi.R",encoding = "utf8" )
source( "./function/plotSSM.CmdStanr.R",encoding = "utf8" )
source( "./function/ql_conv.R",encoding = "utf8" )
source( "./function/b_conv.R",encoding = "utf8" )
source( "./function/Date_conversion.R",encoding = "utf8" )
source( "./function/Create_Metis_MF3_SCname.R",encoding = "utf8" )
source( "./function/MF4Future.pred.R",encoding = "utf8" )
# source( ".././functions/functions_graph.R",encoding = "utf8" )
# source( "../Metis-MF3/function/geom_bar_SCname.R",encoding = "utf8" )
# source( "../Metis-MF3/function/Excel_Preprocessing.R",encoding = "utf8" )
# source( "../Metis-MF3/function/Prefecture.R",encoding = "utf8" )
# source( "../Metis-MF3/function/One_hot.R",encoding = "utf8" )

# install.packages("Cairo")
# install.packages("ggmcmc")
# install.packages("extrafont", dependencies=TRUE)
# install.packages("extrafont", dependencies = TRUE)
# install.packages("devtools")
# install.packages("rstan")
# install.packages("remotes")
# install.packages("ggthemes")
# install.packages("KFAS")
# install.packages("readxl")
# install.packages("ggsci")
