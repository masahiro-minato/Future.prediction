#### EM件数の基本構造時系列モデルでの推定 ####

# ファイル読込
MF3_maintenance <- read_tsv("./tsv_data/Metis_MF3_2211.tsv")　    # 保守データ
MF3_machine <- read_tsv("./tsv_data/Metis_MIF_2211.tsv")          # 機器データ

# 現象項目
MF3_maintenance %>% distinct(Phenomenon) %>% dput()
# 周辺機名称
MF3_maintenance$Peripheral_name
# 列名
names(MF3_maintenance)
# 現象項目の中で上位15項目を抽出
Phe10.MF3 <- 
  head(levels(fct_infreq(MF3_maintenance$Phenomenon)), n=15)
Phe10.MF3
# 現象をインデックスで選択/すべての現象の場合は「index <- 0」で現象名は"全体"と表示される
index <- 7
if(index == 0){
  Phenom.MF3 <- "全体"
}else{
  Phenom.MF3 <- Phe10.MF3[index]
}
# 分析対象のEM現象表示
Phenom.MF3

# 保守データから日付ごとのEM件数を取得/"全体"の場合は現象列でのフィルターを使用しない
if(index != 0){
  EM.count.MF3 <- 
    MF3_maintenance %>% 
    dplyr::filter(is.na(Peripheral_name)) %>%
    dplyr::filter(Phenomenon == Phenom.MF3) %>%
    group_by(Maintenance_date) %>%
    summarise(
      EM.count = n()
    ) %>% 
    arrange(-EM.count) %>%
    ungroup()
}else{
  EM.count.MF3 <- 
    MF3_maintenance %>% 
    dplyr::filter(is.na(Peripheral_name)) %>%
    group_by(Maintenance_date) %>%
    summarise(
      EM.count = n()
    ) %>% 
    arrange(-EM.count) %>%
    ungroup()
}

# 市場機台数
MIF.count.MF3 <- 
  MF3_machine %>% 
  group_by(機種機番,納品年月日) %>% 
  summarise(
    num = n()
  ) %>% 
  group_by(納品年月日) %>% 
  summarise(
    N.date = n()
  ) %>% 
  ungroup()

# 結合
X.date.MF3 <- seq(as.POSIXct("2019-01-01"), as.POSIXct("2022-11-30"), by = "day")
MIF.date.MF3 <- tibble(X.date = X.date.MF3)
MF3.MIF_by.date <- 
  MIF.date.MF3 %>% 
  full_join(MIF.count.MF3, by=c("X.date" = "納品年月日")) %>% 
  full_join(EM.count.MF3, by=c("X.date" = "Maintenance_date")) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  mutate(
    MIF.cumsum = cumsum(N.date)
  )

# 基本構造時系列モデルの推定 --------------------------------------------------------

# データ長 #1430
length(MF3.MIF_by.date$X.date)
MF3.MIF_by.date$X.date[1066]

# データの準備
start_day = which(MF3.MIF_by.date$X.date == "2022-02-01 JST") 
end_day = which(MF3.MIF_by.date$X.date == "2022-11-30 JST")   #1430
pred_term = 20 # 予測期間

# データ
data_list <- list(
  y = MF3.MIF_by.date$EM.count[start_day:end_day],
  Num = MF3.MIF_by.date$MIF.cumsum[start_day:end_day]/max(MF3.MIF_by.date$MIF.cumsum[start_day:end_day]), # 正規化が必要
  T = end_day-start_day+1,
  N_max = max(MF3.MIF_by.date$MIF.cumsum[start_day:end_day]),
  grainsize = 250,
  pred_term = pred_term,
  Num_gain = 200/max(MF3.MIF_by.date$MIF.cumsum[start_day:end_day])
)
# サンプリング開始日/終了日の確認
MF3.MIF_by.date$X.date[start_day]
MF3.MIF_by.date$X.date[end_day]

# コンパイル
mod <- cmdstan_model("./stan/bsts-AR-tvc-rd-poisson-remod_predict.stan", cpp_options = list(stan_threads = TRUE))
# マルチコア対応
options(mc.cores = parallel::detectCores())
# path
output_dir = str_c("./csv/",Phenom.MF3)
output_basename = str_c(Phenom.MF3,".",format(Sys.time(), "%H-%M-%S"))
object.path = str_c("./Cmdstan_files/Metis-MF3.EM.predict.fit.",Phenom.MF3,"-",start_day,"~",end_day,".rds")
# csv保存ディレクトリーの作成
if(!dir.exists(output_dir)){
  dir.create(output_dir)}
# MCMCサンプリング
exTime <- system.time(
  fit <- mod$sample(data = data_list,
                    seed = 1234,
                    chains = 6,
                    parallel_chains = getOption("mc.cores", 24),
                    threads_per_chain = 2,
                    iter_warmup = 4000,
                    iter_sampling = 2000,
                    thin = 4,
                    # adapt_delta = 0.90,
                    max_treedepth = 15,
                    refresh = 500,
                    output_dir = output_dir,
                    output_basename = output_basename,
                    show_messages = FALSE)
)
# 実行時間
exTimeTable <- data.frame(user.self = exTime["user.self"], 
                          sys.self = exTime["sys.self"],
                          elapsed = exTime["elapsed"],
                          Phenomenon = Phenom.MF3,
                          row.names = "time")
exTimeTable
# 保存
write_tsv(exTimeTable, str_c("./time/exTimeTable.",Phenom.MF3,".",start_day,"~",end_day,".tsv"))

# 結果保存
fit$save_object(file = object.path)
# 読込
# fit <- read_rds(object.path)
# パラメータ表示
fit$print(c("s_w", "s_s", "s_r", "s_t", "b_ar", "b_ope[1]", "Intercept", "lp__"))
# rhatヒストグラム
fit %>% bayesplot::rhat() %>% hist
# b_ar 自己回帰係数
b_ar.mean <- round(mean((fit$draws("b_ar") %>% as_draws_df)$b_ar),3)

# 推定結果の図示 -----------------------------------------------------------------
# フォント設定
par(family="Noto Sans")
# x軸の年月日設定
MF3.MIF_by.date$X.date[start_day]
date.plot <- 
  seq(
    from = as.POSIXct(MF3.MIF_by.date$X.date[start_day]),
    by = "days",
    len = end_day-start_day+1+pred_term
  )
length(date.plot)
date_plot.MF3 <- tibble(X.date = date.plot)
date_plot.MF3.MIF_by.date.1 <- 
  date_plot.MF3 %>% 
  left_join(MIF.count.MF3, by=c("X.date" = "納品年月日")) %>% 
  left_join(EM.count.MF3, by=c("X.date" = "Maintenance_date")) %>% 
  dplyr::filter(X.date <= MF3.MIF_by.date$X.date[end_day]) %>% 
  mutate_if(where(is.numeric), ~replace(., is.na(.), 0))
date_plot.MF3.MIF_by.date.2 <- 
  date_plot.MF3 %>% 
  left_join(MIF.count.MF3, by=c("X.date" = "納品年月日")) %>% 
  left_join(EM.count.MF3, by=c("X.date" = "Maintenance_date")) %>% 
  dplyr::filter(X.date > MF3.MIF_by.date$X.date[end_day])
date_plot.MF3.MIF_by.date <- 
  bind_rows(date_plot.MF3.MIF_by.date.1,date_plot.MF3.MIF_by.date.2)

# グラフの表示期間の設定
limits = c(as.POSIXct("2022-04-01 JST"), as.POSIXct(MF3.MIF_by.date$X.date[end_day+pred_term]))
# ｙ軸目盛設定
breaks=seq(0,2000,10)

# すべての成分を含んだ状態推定値の図示
p_lambda_pois <- plotSSM.CmdStanr(fit = fit, 
                                  time_vec = date.plot,
                                  obs_vec = date_plot.MF3.MIF_by.date$EM.count,
                                  state_name = "lambda_exp", 
                                  graph_title = str_c("EM：",Phenom.MF3," λ(μ + γ + tvf + r)：すべての成分を含んだ状態推定値（自己回帰係数 b_ar=",b_ar.mean,"）"), 
                                  y_label = "件数",
                                  date_labels = "%Y/%m",
                                  date_breaks = "2 month") +
  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
  # 表示期間の設定
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") +
  scale_y_continuous(breaks=breaks) +
  annotate("text", x=as.POSIXct(MF3.MIF_by.date$X.date[end_day]), y=Inf, 
           hjust=-0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust=2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label="⇒ 【未来予測区間】", 
           size=5, colour="darkblue")
# 水準成分＋周期成分
p_mu_gamma_pois <- plotSSM.CmdStanr(fit = fit, 
                                    time_vec = date.plot,
                                    # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                                    state_name = "mu_gamma_exp", 
                                    graph_title = "μ + γ：水準成分＋周期成分", 
                                    y_label = "件数",
                                    date_labels = "%Y/%m",
                                    date_breaks = "2 month") +
  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
  # 表示期間の設定
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
  
# 水準成分＋ランダム成分
p_mu_r_pois <- plotSSM.CmdStanr(fit = fit, 
                                time_vec = date.plot,
                                # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                                state_name = "mu_r_exp", 
                                graph_title = "μ + r：水準成分＋ランダム成分", 
                                y_label = "件数",
                                date_labels = "%Y/%m",
                                date_breaks = "2 month") +
  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
  # 表示期間の設定
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 水準成分
p_mu_pois <- plotSSM.CmdStanr(fit = fit, 
                              time_vec = date.plot,
                              # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                              state_name = "mu_exp", 
                              graph_title = str_c("μ：水準成分（自己回帰係数 b_ar=",b_ar.mean,"）"), 
                              y_label = "件数",
                              date_labels = "%Y/%m",
                              date_breaks = "2 month") +
  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
  # 表示期間の設定
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 予測区間
p_pred_pois <- plotSSM.CmdStanr(fit = fit, 
                                time_vec = date.plot,
                                obs_vec = date_plot.MF3.MIF_by.date$EM.count,
                                state_name = "y_pred", 
                                graph_title = str_c("EM：",Phenom.MF3," 95%予測区間 （自己回帰係数 b_ar=",b_ar.mean,"）"), 
                                y_label = "件数",
                                date_labels = "%Y/%m",
                                date_breaks = "2 month",
                                fill.ribbon = "lightgreen") +
  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
  # 表示期間の設定
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") +
  scale_y_continuous(breaks=breaks) +
  annotate("text", x=as.POSIXct(MF3.MIF_by.date$X.date[end_day]), y=Inf, 
           hjust=-0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust=2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label="⇒ 【未来予測区間】", 
           size=5, colour="darkblue")
# ドリフト成分
p_drift_pois <- plotSSM.CmdStanr(fit = fit, 
                                 time_vec = date.plot[32:(end_day-start_day+1+pred_term)],
                                 # time_vec = date_plot[(start_day+31):end_day],
                                 state_name = "delta_exp",
                                 graph_title = "δ：ドリフト成分",
                                 y_label = "delta",
                                 date_labels = "%Y/%m",
                                 date_breaks = "2 month") +
                geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                # 表示期間の設定
                scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 周期成分
p_cycle_pois <- plotSSM.CmdStanr(fit = fit, 
                                 time_vec = date.plot,
                                 state_name = "gamma_exp", 
                                 graph_title = "γ：周期成分", 
                                 y_label = "gamma",
                                 date_labels = "%Y/%m",
                                 date_breaks = "2 month") +
                geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                # 表示期間の設定
                scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# ランダム成分
p_random_pois <- plotSSM.CmdStanr(fit = fit, 
                                  time_vec = date.plot,
                                  state_name = "r_exp", 
                                  graph_title = "r：ランダム成分", 
                                  y_label = "r",
                                  date_labels = "%Y/%m",
                                  date_breaks = "2 month") +
                  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                  # 表示期間の設定
                  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 水準成分＋時変係数成分
p_mu_tvc_pois <- plotSSM.CmdStanr(fit = fit, 
                                  time_vec = date.plot,
                                  # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                                  state_name = "mu_tvc_exp", 
                                  graph_title = "μ + tvf：水準成分＋時変係数成分", 
                                  y_label = "件数",
                                  date_labels = "%Y/%m",
                                  date_breaks = "2 month") +
                  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                  # 表示期間の設定
                  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 時変係数ｘ稼働台数
p_tvc_pois <- plotSSM.CmdStanr(fit = fit, 
                               time_vec = date.plot,
                               state_name = "tvc_exp", 
                               graph_title = "tvf：時変係数ｘ稼働台数", 
                               y_label = "件数",
                               date_labels = "%Y/%m",
                               date_breaks = "2 month") +
                  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                  # 表示期間の設定
                  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 水準成分+時変係数×1台稼働
p_mu_tvc_1_pois <- plotSSM.CmdStanr(fit = fit, 
                                    time_vec = date.plot,
                                    # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                                    state_name = "mu_tvc_1_exp", 
                                    graph_title = "μ + tvc：水準成分+時変係数ｘ1台稼働", 
                                    y_label = "件数",
                                    date_labels = "%Y/%m",
                                    date_breaks = "2 month") +
                    geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                    # 表示期間の設定
                    scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")
# 時変係数×1台稼働
p_tvc_1_pois <- plotSSM.CmdStanr(fit = fit, 
                                 time_vec = date.plot,
                                 # obs_vec = MF3.MIF_by.date$EM.count[start_day:end_day],
                                 state_name = "tvc_1", 
                                 graph_title = "tvc：時変係数", 
                                 y_label = "tvc",
                                 date_labels = "%Y/%m",
                                 date_breaks = "2 month") +
                  geom_vline(xintercept = as.POSIXct(MF3.MIF_by.date$X.date[end_day]), linetype = 2, color = "darkblue") +
                  # 表示期間の設定
                  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m")


# グラフ描画
plot <- plot_grid(p_lambda_pois, 
                  p_pred_pois,
                  p_mu_tvc_pois,
                  p_tvc_pois,
                  p_tvc_1_pois,
                  p_mu_pois,
                  p_drift_pois,
                  p_random_pois,
                  p_cycle_pois,
                  ncol = 1, 
                  align = "v")
# now add the title
title <- ggdraw() + 
  draw_label(
    str_c("Metis-MF3　　EM現象：", Phenom.MF3),
    fontface = 'bold',
    color = "darkblue",
    size = 30,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 100)
  )


plot.title <- 
  plot_grid(
    title, plot,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.03, 1)
  )

# OneDriveのパス
OneDrive.path <- "/mnt/c/Users/masah/OneDrive/ドキュメント/wsl/R_working/Future.prediction"
# 集合グラフ保存
ggsave(str_c("./PDF/Metis-MF3-EM_AR.tvc.period_poisson.予測区間",Phenom.MF3,"-",start_day,".pdf"), 
       plot = plot.title, device = cairo_pdf, dpi=300, width=40, height=30)
ggsave(str_c(OneDrive.path,"/PDF/Metis-MF3-EM_AR.tvc.period_poisson.予測区間",Phenom.MF3,"-",start_day,".pdf"), 
       plot = plot.title, device = cairo_pdf, dpi=300, width=40, height=30)

# 予測区間グラフ保存
ggsave(str_c("./PDF/Metis-MF3-EM-",Phenom.MF3,".予測区間-",start_day,".pdf"), 
       plot = p_pred_pois, device = cairo_pdf, dpi=300, width=20, height=5)
# 信頼区間グラフ保存
ggsave(str_c("./PDF/Metis-MF3-EM-",Phenom.MF3,".信頼区間-",start_day,".pdf"), 
       plot = p_lambda_pois, device = cairo_pdf, dpi=300, width=20, height=5)

# 取得EMデータのグラフ追記(仮値として乱数設定)
X.date.Feature <- seq(as.POSIXct(MF3.MIF_by.date$X.date[end_day]+3600*24), 
                      as.POSIXct(MF3.MIF_by.date$X.date[end_day]+3600*24*30), by = "day")
MIF.date.Feature <- tibble(X.date = X.date.Feature,
                           EM = floor(runif(pred_term, min=0, max=max(date_plot.MF3.MIF_by.date$EM.count,na.rm = TRUE))))
p_pred_pois.Feature <- 
  p_pred_pois +
  geom_point(alpha = 0.8, size = 1.5, shape=17, color="red",
             data = MIF.date.Feature, aes(x = X.date, y = EM))
# 予測区間+取得データグラフ保存
ggsave(str_c("./PDF/Metis-MF3-EM-",Phenom.MF3,".予測区間+取得データ-",start_day,".pdf"), 
       plot = p_pred_pois.Feature, device = cairo_pdf, dpi=300, width=20, height=5)

ggsave(str_c(OneDrive.path,"/PDF/Metis-MF3-EM-",Phenom.MF3,".予測区間+取得データ-",start_day,".pdf"), 
       plot = p_pred_pois.Feature, device = cairo_pdf, dpi=300, width=20, height=5)
