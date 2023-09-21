#### EM件数の基本構造時系列モデルでの推定 ####

# ファイル読込
MF3_maintenance <- read_tsv("./tsv_data/Metis_MF3_2303.tsv")　    # 保守データ
MF3_machine <- read_tsv("./tsv_data/Metis_MIF_2303.tsv")          # 機器データ

# 周辺機名称
MF3_maintenance %>% distinct(Peripheral_name) %>% dput()
# 周辺機名称
MF3_maintenance$Peripheral_name
# 列名
names(MF3_maintenance)
# 周辺機機種ごとの時系列でのEM件数の抽出
EM.count.MF3.Peripheral <- 
  MF3_maintenance %>% 
  dplyr::filter((Peripheral_name %in% c("COOK-C", "SINAI-H") & Treatment_location %in% c("ADF部")) | 
                (Peripheral_name %in% c("AMUR-C(HY)", "AMUR-C中綴じ","VOLGA-E") & Treatment_location %in% c("ﾊﾟﾝﾁ部","ﾌｨﾆｯｼｬｰ/ｿｰﾀｰ部"))) %>%
  group_by(Maintenance_date, Peripheral_name) %>%
  summarise(
    EM.count = n()
  ) %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(Peripheral_name),
              values_from = EM.count) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  # 列名の変更
  set_colnames(c("Maintenance_date", "EM.SINAI", "EM.VOLGA", "EM.COOK", "EM.AMUR", "EM.AMUR.HY")) %>%
  # 不要列の削除
  select("Maintenance_date", "EM.COOK", "EM.SINAI", "EM.VOLGA", "EM.AMUR", "EM.AMUR.HY")

# 市場機台数
MIF.count.MF3.Peripheral <- 
  MF3_machine %>% 
  group_by(納品年月日,ADF,フィニッシャー) %>% 
  summarise(
    N.date = n()
  ) %>% 
  # tidyr::drop_na() %>% # この行があると後処理がない場合のADFも削除される
  ungroup() %>% 
  pivot_wider(names_from = c(ADF,フィニッシャー),
              values_from = N.date) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
# 周辺機機種ごとの市場機台数を集計/機種名の含まれる列の行合計を計算する
# MIF.count.MF3.Peripheral <- 
#   MIF.count.MF3.Peripheral %>% 
  mutate(
    MIF.COOK = select(.,contains("COOK-C")) %>% rowSums(na.rm=TRUE),
    MIF.SINAI = select(.,contains("SINAI-H")) %>% rowSums(na.rm=TRUE),
    MIF.VOLGA = select(.,contains("VOLGA-E")) %>% rowSums(na.rm=TRUE),
    MIF.AMUR = select(.,contains("AMUR-C中綴じ")) %>% rowSums(na.rm=TRUE),
    MIF.AMUR.HY = select(.,contains("AMUR-C(HY)")) %>% rowSums(na.rm=TRUE)
  ) %>% 
  select(納品年月日, MIF.COOK, MIF.SINAI, MIF.VOLGA, MIF.AMUR, MIF.AMUR.HY)

# 結合
X.date.MF3 <- seq(as.POSIXct("2019-01-01"), as.POSIXct("2023-03-31"), by = "day")
MIF.date.MF3.Peripheral <- tibble(X.date = X.date.MF3)
MF3.MIF_by.date <- 
  MIF.date.MF3.Peripheral %>% 
  full_join(MIF.count.MF3.Peripheral, by=c("X.date" = "納品年月日")) %>% 
  full_join(EM.count.MF3.Peripheral, by=c("X.date" = "Maintenance_date")) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  # 市場機台数を累積和へ変換
  mutate_at(vars(contains("MIF")), cumsum)

# データの準備
start_day = which(MF3.MIF_by.date$X.date == "2019-02-01 JST")  # "2020-01-01 JST"
end_day = which(MF3.MIF_by.date$X.date == "2023-03-31 JST")
pred_term = 30 # 予測期間
data.num <- end_day - start_day + 1 + pred_term # 予測期間を含むデータ数
# EM観測値のみをマトリクスへ変換
names(MF3.MIF_by.date)
# 列番号の取得
X.date.colnum <- which(names(MF3.MIF_by.date) == "X.date")
# サンプリング期間におけるEM観測値をマトリクスへ変換と転置
y <- t(as.matrix(MF3.MIF_by.date[start_day:end_day,][,c(-X.date.colnum, -c(2:6))]))
# 稼働台数
Num.tib <- MF3.MIF_by.date[start_day:end_day,][,c(-X.date.colnum, -c(7:11))]
# 正規化のために最大値を算出
max_data <-
  Num.tib %>%
  summarize(across(where(is.numeric), \(x) max(x, na.rm = TRUE))) # R4.2.3以降
max_data$MIF.COOK
print(max_data)
# ベクトルへの変換
unlist(c(max_data))
# 正規化
for(col in names(max_data)){
  print(col)
  print(eval(parse(text = paste0("max_data$",col))))
  Num.tib <- Num.tib %>%
    mutate(
      !!col :=eval(parse(text = paste0("Num.tib$",col)))/eval(parse(text = paste0("max_data$",col)))
    )
}
# マトリクスへ変換と転置
Num <- t(as.matrix(Num.tib))
# データリストの作成
data_list <- list(
  y = y,
  N = 5,
  Num = Num, # 正規化が必要
  T = end_day-start_day+1,
  N_max = unlist(c(max_data)),
  grainsize = 250,
  pred_term = pred_term,
  Num_gain = c(200,20,100,50,20)/unlist(c(max_data)) # 月次での稼働台数増加予測数
)
# コンパイル
mod <- cmdstan_model("./stan/bsts-AR-tvc-rd-smooth-pred-poisson-periphe-matrix.stan", cpp_options = list(stan_threads = TRUE))
# マルチコア対応
options(mc.cores = parallel::detectCores())
# path
title <- str_c("Metis-MF3-Peripheral")
output_dir <- str_c("./csv/matrix")
output_basename <-str_c(title,".",format(Sys.time(), "%H-%M-%S"))
object.path <- str_c("./Cmdstan_files/matrix/",title,".EM.pred.Peripheral.matrix.fit-",start_day,"~",end_day,".rds")
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
                    iter_warmup = 4000, #10000,
                    iter_sampling = 2000, #5000,
                    thin = 4,
                    # adapt_delta = 0.90,
                    max_treedepth = 12,
                    refresh = 500,
                    init = 0, # エラー回避
                    output_dir = output_dir, # csvデータ保存ディレクトリ
                    output_basename = output_basename, # csvデータ保存名称
                    show_messages = FALSE
  )
)

# 実行時間
exTimeTable <- data.frame(user.self = exTime["user.self"], 
                          sys.self = exTime["sys.self"],
                          elapsed = exTime["elapsed"],
                          title = title,
                          row.names = "time")
exTimeTable
# 保存
write_tsv(exTimeTable, str_c("./time/exTimeTable.",title,".",start_day,"~",end_day,".tsv"))
# 結果保存
fit$save_object(file = object.path)
# サンプリング結果表示
fit$print(c( "s_z","s_s", "s_r", "s_t", "b_ar", "b_ope[1,1]", "b_ope[2,2]","b_ope[3,3]","Intercept", "lp__"), max_rows=100)
# rhatヒストグラム
fit %>% bayesplot::rhat() %>% hist(main=str_c("Histogram of ",title))
# 時系列グラフ用ｘ軸データ
time_vec <- 
  seq(
    from = as.POSIXct("2019-02-01 JST"),
    by = "days",
    len = data.num # end_day-start_day+1+pred_term
  )
length(time_vec)

names(MF3.MIF_by.date)[2:6]

# 表示パラメータの設定
parm <- "y_pred"
# 表示パラメータの抽出
parm.tib <- 
  as_tibble(fit$draws(parm) %>% as_draws_df)
# 分位数の算出
parm.quan <- 
  apply(parm.tib, 2, function(i){quantile(i,prob=c(0.025, 0.5, 0.975))})
parm.quan <- 
  as_tibble(parm.quan) %>% 
  # 不要列の削除
  select(-.chain, -.iteration, -.draw )
# 周辺機名称ベクトルの作成
Peripheral.MF3 <- substring(names(EM.count.MF3.Peripheral)[-1], 4,)
# グラフ表示用データフレームの作成
result_df <- 
  as_tibble(t(parm.quan),.name_repair = 'unique') %>% 
  mutate(
    Peripheral = rep(Peripheral.MF3, data.num)
  ) %>% 
  # 列名設定
  set_colnames(c("lwr", "fit", "upr", "Peripheral"))
# 時間軸項の追加
result_df$time <- 
  rep(time_vec, times = rep(length(Peripheral.MF3), length(time_vec)))
# レベル設定
result_df$Peripheral <- 
  factor(result_df$Peripheral, levels=Peripheral.MF3)
# DFとFINへ区分
result_df.ADF <- 
  result_df %>% 
  dplyr::filter(Peripheral %in% c("COOK","SINAI"))
result_df.FIN <- 
  result_df %>% 
  dplyr::filter(Peripheral %in% c("VOLGA","AMUR","AMUR.HY"))

# b_ar 自己回帰係数
b_ar <- 
  fit$draws("b_ar") %>% 
  as_draws_df
# b_ar平均値
b_ar.mean <- round(fit$draws("b_ar") %>% apply(3,mean),3)
print(b_ar.mean)

# 測定値データ
MF3.EM_by.date <- 
  MF3.MIF_by.date %>% 
  select(-starts_with("MIF.")) %>%
  set_colnames(c("X.date",Peripheral.MF3)) %>% 
  pivot_longer(cols = all_of(Peripheral.MF3),
               names_to = "Peripheral",
               values_to = "EM.num")
MF3.EM_by.date$Peripheral <- 
  factor(MF3.EM_by.date$Peripheral, levels = Peripheral.MF3)
# ADFのみ抽出
MF3.EM_by.date.ADF <- 
  MF3.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("COOK","SINAI"))
# FINのみ抽出
MF3.EM_by.date.FIN <- 
  MF3.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("VOLGA","AMUR","AMUR.HY"))

# フォント定義
par(family="Noto Sans")
# グラフタイトル
if (parm == "lambda_exp") {
  graph_title <- str_c("Metis-MF3 EM件数時系列分析 λ(μ + γ + tvf + r)：すべての成分を含んだ状態推定値")
  fill.graph <- "lightblue"
}else if (parm == "y_pred") {
  graph_title <- str_c("Metis-MF3 EM件数時系列分析 95%予測区間")
  fill.graph <- "lightgreen"
}
# 縦軸名称
y_label = "件数"
# 凡例設定
vec <- c()
for (i in 1:length(Peripheral.MF3)) {
  vec <- 
    vec %>% 
    append(str_c(i,".",Peripheral.MF3[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
vec.ADF <- c()
for (i in 1:2) {
  vec.ADF <- 
    vec.ADF %>% 
    append(str_c(Peripheral.MF3[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
vec.FIN <- c()
for (i in 3:5) {
  vec.FIN <- 
    vec.FIN %>% 
    append(str_c(Peripheral.MF3[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
# ADF+FIN
Att.labs <- vec 
names(Att.labs) <- Peripheral.MF3
# グラフの表示期間の設定
limits = c(as.POSIXct("2022-04-01 JST"), as.POSIXct(MF3.MIF_by.date$X.date[end_day] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,1000,50)
# 図示
p <- ggplot(data = result_df, aes(x = time)) + 
  theme_bw() + 
  labs(title = graph_title) +
  theme(plot.title = element_text(size = 18,  #font size and adjust
                                  hjust = 0.01,#adjust
  )) +
  ylab(y_label) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill=fill.graph) + 
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_point(alpha = 1.0, size = 0.9,
             data = MF3.EM_by.date, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  geom_vline(xintercept = as.POSIXct("2023-03-31 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  scale_y_continuous(breaks = breaks) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-03-31 JST"), y = Inf, 
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】", 
           size = 5, 
           colour = "darkblue")

plot(p)

# 予測区間+取得データグラフ保存
n <- 5
ggsave(str_c("./PDF/",title,".EM件数 予測区間-",start_day,"~",end_day,".pdf"),
       plot = p, device = cairo_pdf, dpi=300, width=20, height=(n*3))

# ADF
Att.labs <- vec.ADF
names(Att.labs) <- Peripheral.MF3[1:2]
# グラフの表示期間の設定
limits = c(as.POSIXct("2022-04-01 JST"), as.POSIXct(MF3.MIF_by.date$X.date[end_day] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,1000,50)
# 図示
p.ADF <- ggplot(data = result_df.ADF, aes(x = time)) + 
  theme_bw() + 
  labs(title = graph_title) +
  theme(plot.title = element_text(size = 18,  #font size and adjust
                                  hjust = 0.01,#adjust
  )) +
  ylab(y_label) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill=fill.graph) + 
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_point(alpha = 1.0, size = 0.9,
             data = MF3.EM_by.date.ADF, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  geom_vline(xintercept = as.POSIXct("2023-03-31 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  scale_y_continuous(breaks = breaks) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-03-31 JST"), y = Inf, 
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】", 
           size = 5, 
           colour = "darkblue")

plot(p.ADF)

# 予測区間+取得データグラフ保存
n <- 2
ggsave(str_c("./PDF/",title,".ADF_EM件数 予測区間-",start_day,"~",end_day,".pdf"), 
       plot = p.ADF, device = cairo_pdf, dpi=300, width=20, height=(n*5))

# FIN
Att.labs <- vec.FIN
names(Att.labs) <- Peripheral.MF3[3:5]
# グラフの表示期間の設定
limits = c(as.POSIXct("2022-04-01 JST"), as.POSIXct(MF3.MIF_by.date$X.date[end_day] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,100,2)
# 図示
p.FIN <- ggplot(data = result_df.FIN, aes(x = time)) + 
  theme_bw() + 
  labs(title = graph_title) +
  theme(plot.title = element_text(size = 18,  #font size and adjust
                                  hjust = 0.01,#adjust
  )) +
  ylab(y_label) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill=fill.graph) + 
  geom_line(aes(y = fit), linewidth = 1.2) +
  geom_point(alpha = 1.0, size = 0.9,
             data = MF3.EM_by.date.FIN, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  geom_vline(xintercept = as.POSIXct("2023-03-31 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  scale_y_continuous(breaks = breaks) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-03-31 JST"), y = Inf, 
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】", 
           size = 5, 
           colour = "darkblue")

plot(p.FIN)

# 予測区間+取得データグラフ保存
n <- 3
ggsave(str_c("./PDF/",title,".FIN_EM件数 予測区間-",start_day,"~",end_day,".pdf"), 
       plot = p.FIN, device = cairo_pdf, dpi=300, width=20, height=(n*5))

