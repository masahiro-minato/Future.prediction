
fit<-readRDS("./Cmdstan_files/matrix/Metis-MF4-Peripheral.MF3parm.EM.pred.Peripheral.matrix.fit-29~120.rds")

# 表示パラメータの設定
parm <- "s_z3"
# 表示パラメータの抽出
parm.tib <- 
  as_tibble(fit$draws(parm) %>% as_draws_df)

# "s_z3","s_s3", "s_r3", "s_t3",

s_z3.tib <- 
  as_tibble(fit$draws("s_z3") %>% as_draws_df) %>% 
  select(1:4)
s_s3.tib <- 
  as_tibble(fit$draws("s_s3") %>% as_draws_df) %>% 
  select(1:4)
s_r3.tib <- 
  as_tibble(fit$draws("s_r3") %>% as_draws_df) %>% 
  select(1:4)
s_t3.tib <- 
  as_tibble(fit$draws("s_t3") %>% as_draws_df) %>% 
  select(1:4)
mu4.tib <- 
  as_tibble(fit$draws("mu4") %>% as_draws_df) %>% 
  select(starts_with("mu4"))
b_ar4.tib <- 
  as_tibble(fit$draws("b_ar4") %>% as_draws_df) %>% 
  select(starts_with("b_ar4"))
Intercept4.tib <- 
  as_tibble(fit$draws("Intercept4") %>% as_draws_df) %>% 
  select(starts_with("Intercept4"))
b_ope4.tib <- 
  as_tibble(fit$draws("b_ope4") %>% as_draws_df) %>% 
  select(starts_with("b_ope4"))
lambda4.tib <- 
  as_tibble(fit$draws("lambda4") %>% as_draws_df) %>% 
  select(starts_with("lambda4"))
sumgamma4.tib <- 
  as_tibble(fit$draws("sumgamma4") %>% as_draws_df) %>% 
  select(starts_with("sumgamma4"))
tvc4.tib <- 
  as_tibble(fit$draws("tvc4") %>% as_draws_df) %>% 
  select(starts_with("tvc4"))
r4.tib <- 
  as_tibble(fit$draws("r4") %>% as_draws_df) %>% 
  select(starts_with("r4"))
gamma4.tib <- 
  as_tibble(fit$draws("gamma4") %>% as_draws_df) %>% 
  select(starts_with("gamma4"))

# // ランダム効果成分 r
# r_pred[n,(T4+i)] = normal_rng(0, s_r3[n]);
# // 水準成分 μ
# mu_pred[n,(T4+i)] = normal_rng(Intercept4[n] + b_ar4[n] * (2 * mu_pred[n,(T4+i-1)] - mu_pred[n,(T4+i-2)]), s_z3[n]);
# // 稼働台数時変係数 b_ope
# b_ope_pred[n,(T4+i)] = normal_rng(b_ope_pred[n,(T4+i-1)], s_t3[n]);
# // 周期成分 γ
# sumgamma_pred[n,(T4-6+i)] = -sum(gamma_pred[n,(T4-6+i):(T4-6+i+5)]);
# gamma_pred[n,(T4+i)] = normal_rng(sumgamma_pred[n,(T4-6+i)], s_s3[n]);
# // 想定稼働台数
# Num_pred[n,T4+i] = Num_pred[n,T4+i-1] + Num_gain[n];
# // 時変係数成分 tvc
# tvc_pred[n,(T4+i)] = b_ope_pred[n,(T4+i)] * Num_pred[n,T4+i];
# // 状態推定値 λ
# lambda_pred[n,(T4+i)] = mu_pred[n,(T4+i)] + gamma_pred[n,(T4+i)] + tvc_pred[n,(T4+i)] + r_pred[n,(T4+i)];

# 標準偏差xに対する標準正規分布の乱数1個を返す関数
f <- function(x){
  return(rnorm(1,0,x))
}
# // ランダム効果成分 r
for (n in 1:4) {
  for (i in 93:96) {
    print(n)
    s_r3.tib <- 
      s_r3.tib %>% 
      mutate(
        !!str_c("r[",n,",",i,"]") := unlist(map(eval(parse(text = str_c("`s_r3[",n,"]`"))), f)))
  }
}
r.pred <- 
  s_r3.tib %>% 
  select(starts_with("r"))
  
#期待値x,標準偏差yに対する標準正規分布の乱数1個を返す関数
f2 <- function(x,y){
  return(rnorm(1,x,y))
}
# // 水準成分 μ
# mu_pred[n,(T4+i)] = normal_rng(Intercept4[n] + b_ar4[n] * (2 * mu_pred[n,(T4+i-1)] - mu_pred[n,(T4+i-2)]), s_z3[n]);
mu4.pred <- 
  mu4.tib %>% 
  cbind(s_z3.tib) %>% 
  cbind(Intercept4.tib) %>% 
  cbind(b_ar4.tib)
for (n in 1:4) {
  for (i in 93:96) {
    mu4.pred <- 
      mu4.pred %>% 
      mutate(
        !!str_c("mu4[",n,",",i,"]") := unlist(map2(eval(parse(text = str_c("`Intercept4[",n,"]` + ",
                                                                           "`b_ar4[",n,"]` *(2*`mu4[",n,",",(i-1),"]`- `mu4[",n,",",(i-2),"]`)"))),
                                                   eval(parse(text = str_c("`s_z3[",n,"]`"))), f2))
      )
  }
}
mu4.pred <- 
  mu4.pred %>% 
  select(starts_with("mu4"))

# // 稼働台数時変係数 b_ope
b_ope4.pred <- 
  b_ope4.tib %>% 
  cbind(s_t3.tib)
n <- 1
i <- 93
for (n in 1:4) {
  for (i in 93:96) {
    b_ope4.pred <- 
      b_ope4.pred %>% 
      mutate(
        !!str_c("b_ope4[",n,",",i,"]") := unlist(map2(eval(parse(text = str_c("`b_ope4[",n,",",(i-1),"]`"))),
                                                      eval(parse(text = str_c("`s_t3[",n,"]`"))), 
                                                      f2))
      )
  }
}
b_ope4.pred <- 
  b_ope4.pred %>% 
  select(starts_with("b_ope4"))
# // 周期成分 γ
# sumgamma_pred[n,(T4-6+i)] = -sum(gamma_pred[n,(T4-6+i):(T4-6+i+5)]);
# gamma_pred[n,(T4+i)] = normal_rng(sumgamma_pred[n,(T4-6+i)], s_s3[n]);
gamma4.pred <- 
  sumgamma4.tib %>% 
  cbind(s_s3.tib) %>% 
  cbind(gamma4.tib)

for (n in 1:4) {
  for (i in 93:96) {
    gamma4.pred <- 
      gamma4.pred %>% 
      mutate(
        !!str_c("sumgamma4[",n,",",i-6,"]") := eval(parse(text = str_c("apply(gamma4.pred[,c(
                                                                     \"gamma4[",n,",",(i-6),"]\",
                                                                     \"gamma4[",n,",",(i-5),"]\",
                                                                     \"gamma4[",n,",",(i-4),"]\",
                                                                     \"gamma4[",n,",",(i-3),"]\",
                                                                     \"gamma4[",n,",",(i-2),"]\",
                                                                     \"gamma4[",n,",",(i-1),"]\")],1,sum)"))),
        !!str_c("gamma4[",n,",",i,"]") := unlist(map2(eval(parse(text = str_c("`sumgamma4[",n,",",(i-6),"]`"))),
                                                      eval(parse(text = str_c("`s_s3[",n,"]`"))),
                                                      f2))
      )
  }
}
gamma4.pred <- 
  gamma4.pred %>% 
  select(starts_with("gamma4"))

# gamma4.pred %>% 
#   mutate(
#     `sumgamma4[1,93]` = apply(gamma4.pred[,c("gamma4[1,87]","gamma4[1,88]")],1,sum)
#   )

# // 想定稼働台数
# Num_pred[n,T4+i] = Num_pred[n,T4+i-1] + Num_gain[n];

Num.pred <- 
  tibble(`Num.pred[1,92]` = c(7500),
         `Num.pred[2,92]` = c(1500),
         `Num.pred[3,92]` = c(1000),
         `Num.pred[4,92]` = c(800))
Num.gain <- c(70,20,10,8)

for (n in 1:4) {
  for (i in 93:96) {
    print(i)
    Num.pred <- 
      Num.pred %>% 
      mutate(
        !!str_c("Num.pred[",n,",",i,"]") := eval(parse(text = str_c("`Num.pred[",n,",",i-1,"]` + Num.gain[",n,"]")))
      )
  }
}

# // 時変係数成分 tvc
# tvc_pred[n,(T4+i)] = b_ope_pred[n,(T4+i)] * Num_pred[n,T4+i];

tvc.pred <- 
  tvc4.tib %>% 
  cbind(b_ope4.pred) %>% 
  cbind(Num.pred)

for (n in 1:4) {
  for (i in 93:96) {
    print(i)
    tvc.pred <- 
      tvc.pred %>% 
      mutate(
        !!str_c("tvc4[",n,",",i,"]") := eval(parse(text = str_c("`b_ope4[",n,",",i,"]` * `Num.pred[",n,",",i,"]`")))
      )
  }
}
tvc.pred <- 
  tvc.pred %>% 
  select(starts_with("tvc4"))

# // 状態推定値 λ
# lambda_pred[n,(T4+i)] = mu_pred[n,(T4+i)] + gamma_pred[n,(T4+i)] + tvc_pred[n,(T4+i)] + r_pred[n,(T4+i)];

lambda.pred <- 
  lambda4.tib %>% 
  cbind(mu4.pred) %>% 
  cbind(gamma4.pred) %>% 
  cbind(tvc.pred) %>% 
  cbind(r.pred)
for (n in 1:4) {
  for (i in 93:96) {
    print(i)
    lambda.pred <- 
      lambda.pred %>% 
      mutate(
        !!str_c("lambda4[",n,",",i,"]") := eval(parse(text = str_c("`mu4[",n,",",i,"]` + `gamma4[",n,",",i,"]` + `tvc4[",n,",",i,"]` + `r[",n,",",i,"]`")))
      )
  }
}
lambda.pred <- 
  lambda.pred %>% 
  select(starts_with("lambda4"))
