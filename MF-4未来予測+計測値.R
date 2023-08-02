# ファイル読込
MF4_maintenance <- read_tsv("./tsv_data/Metis_MF4_2306.tsv")　    # 保守データ
MF4_machine <- read_tsv("./tsv_data/MF4_MIF_2306.tsv")            # 機器データ

# 周辺機名称
MF4_maintenance %>% distinct(Peripheral_name) %>% dput()
# 周辺機名称
MF4_maintenance$Peripheral_name
# 列名
names(MF4_maintenance)
# 周辺機機種ごとの時系列でのEM件数の抽出
EM.count.MF4.Peripheral <- 
  MF4_maintenance %>% 
  dplyr::filter((Peripheral_name %in% c("COOK-D", "CATHERINE") & Treatment_location %in% c("ADF部")) | 
                  (Peripheral_name %in% c("AMUR-D(HY)", "AMUR-D中綴じ","VOLGA-H") & Treatment_location %in% c("ﾊﾟﾝﾁ部","ﾌｨﾆｯｼｬｰ/ｿｰﾀｰ部"))) %>%
  group_by(Maintenance_date, Peripheral_name) %>%
  summarise(
    EM.count = n()
  ) %>% 
  ungroup() %>% 
  pivot_wider(names_from = c(Peripheral_name),
              values_from = EM.count) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  # 列名の変更
  set_colnames(c("Maintenance_date", "EM.CATHERINE", "EM.AMUR", "EM.COOK", "EM.VOLGA")) %>%
  # 不要列の削除
  select("Maintenance_date", "EM.COOK", "EM.CATHERINE", "EM.VOLGA", "EM.AMUR")

# 市場機台数
MIF.count.MF4.Peripheral <- 
  MF4_machine %>% 
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
  mutate(
    MIF.COOK = select(.,contains("COOK-D")) %>% rowSums(na.rm=TRUE),
    MIF.CATHERINE = select(.,contains("CATHERINE")) %>% rowSums(na.rm=TRUE),
    MIF.VOLGA = select(.,contains("VOLGA-H")) %>% rowSums(na.rm=TRUE),
    MIF.AMUR = select(.,contains("AMUR-D中綴じ")) %>% rowSums(na.rm=TRUE)
    # MIF.AMUR.HY = select(.,contains("AMUR-D(HY)")) %>% rowSums(na.rm=TRUE)
  ) %>% 
  select(納品年月日, MIF.COOK, MIF.CATHERINE, MIF.VOLGA, MIF.AMUR)

# 結合
X.date.MF4 <- seq(as.POSIXct("2023-02-01"), as.POSIXct("2023-06-30"), by = "day")
MIF.date.MF4.Peripheral <- tibble(X.date = X.date.MF4)
MF4.MIF_by.date <- 
  MIF.date.MF4.Peripheral %>% 
  full_join(MIF.count.MF4.Peripheral, by=c("X.date" = "納品年月日")) %>% 
  full_join(EM.count.MF4.Peripheral, by=c("X.date" = "Maintenance_date")) %>% 
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  # 市場機台数を累積和へ変換
  mutate_at(vars(contains("MIF")), cumsum)

# 周辺機名称ベクトルの作成
Peripheral.MF4 <- substring(names(EM.count.MF4.Peripheral)[-1], 4,)
# 測定値データ
MF4.EM_by.date <- 
  MF4.MIF_by.date %>% 
  select(-starts_with("MIF.")) %>%
  set_colnames(c("X.date",Peripheral.MF4)) %>% 
  pivot_longer(cols = all_of(Peripheral.MF4),
               names_to = "Peripheral",
               values_to = "EM.num")

MF4.EM_by.date$Peripheral <- 
  factor(MF4.EM_by.date$Peripheral, levels = Peripheral.MF4)

# 未来予測区間
MF4.6.EM_by.date <- 
  MF4.EM_by.date %>% 
  dplyr::filter(X.date > "2023-05-31")
# ADFのみ抽出
MF4.6.EM_by.date.ADF <- 
  MF4.6.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("COOK","CATHERINE"))
# FINのみ抽出
MF4.6.EM_by.date.FIN <- 
  MF4.6.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("VOLGA","AMUR"))
# 計測区間
MF4.EM_by.date <- 
  MF4.EM_by.date %>% 
  dplyr::filter(X.date <= "2023-05-31")
# ADFのみ抽出
MF4.EM_by.date.ADF <- 
  MF4.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("COOK","CATHERINE"))
# FINのみ抽出
MF4.EM_by.date.FIN <- 
  MF4.EM_by.date %>% 
  dplyr::filter(Peripheral %in% c("VOLGA","AMUR"))

start_day.MF4 = which(MF4.MIF_by.date$X.date == "2023-03-01 JST")
end_day.MF4 = which(MF4.MIF_by.date$X.date == "2023-05-31 JST")
pred_term = 30 # 予測期間
data.num.MF4 <- end_day.MF4 - start_day.MF4 + 1 + pred_term # 予測期間を含むデータ数
# path
title <- str_c("Metis-MF4-Peripheral.MF3parm")
output_dir <- str_c("./csv/matrix")
output_basename <-str_c(title,".",format(Sys.time(), "%H-%M-%S"))
object.path <- str_c("./Cmdstan_files/matrix/",title,".EM.pred.Peripheral.matrix.fit-",start_day.MF4,"~",end_day.MF4,".rds")

# サンプリングデータ読込
fit<-readRDS("./Cmdstan_files/matrix/Metis-MF4-Peripheral.MF3parm.EM.pred.Peripheral.matrix.fit-29~120.rds")
# サンプリング結果表示
fit$print(c( "s_z3","s_s3", "s_r3", "s_t3", "b_ar3", "b_ope4[1,1]", "b_ope4[2,2]","b_ope4[3,3]","Intercept4", "lp__"), max_rows=100)
# rhatヒストグラム
fit %>% bayesplot::rhat() %>% hist(main=str_c("Histogram of ",title))
# 時系列グラフ用ｘ軸データ
time_vec <- 
  seq(
    from = as.POSIXct("2023-03-01 JST"),
    by = "days",
    len = data.num.MF4 # end_day-start_day+1+pred_term
  )
length(time_vec)

names(MF4.MIF_by.date)[2:5]

# 表示パラメータの設定
parm <- "y_pred"
# 表示パラメータの抽出
parm.tib <- 
  as_tibble(fit$draws(parm) %>% as_draws_df)
# 分位数の算出
parm.quan <- 
  apply(parm.tib, 2, function(i){quantile(i,prob=c(0.025, 0.5, 0.975), na.rm=TRUE)})
parm.quan <- 
  as_tibble(parm.quan) %>% 
  # 不要列の削除
  select(-.chain, -.iteration, -.draw )
# 周辺機名称ベクトルの作成
Peripheral.MF4 <- substring(names(EM.count.MF4.Peripheral)[-1], 4,)
# グラフ表示用データフレームの作成
result_df <- 
  as_tibble(t(parm.quan),.name_repair = 'unique') %>% 
  mutate(
    Peripheral = rep(Peripheral.MF4, data.num.MF4)
  ) %>% 
  # 列名設定
  set_colnames(c("lwr", "fit", "upr", "Peripheral"))
# 時間軸項の追加
result_df$time <- 
  rep(time_vec, times = rep(length(Peripheral.MF4), length(time_vec)))
# レベル設定
result_df$Peripheral <- 
  factor(result_df$Peripheral, levels=Peripheral.MF4)
# DFとFINへ区分
result_df.ADF <- 
  result_df %>% 
  dplyr::filter(Peripheral %in% c("COOK","CATHERINE"))
result_df.FIN <- 
  result_df %>% 
  dplyr::filter(Peripheral %in% c("VOLGA","AMUR"))

# b_ar 自己回帰係数
b_ar <- 
  fit$draws("b_ar4") %>% 
  as_draws_df
# b_ar平均値
b_ar.mean <- round(fit$draws("b_ar4") %>% apply(3,mean),3)
print(b_ar.mean)

# # 測定値データ
# MF4.EM_by.date <- 
#   MF4.MIF_by.date %>% 
#   select(-starts_with("MIF.")) %>%
#   set_colnames(c("X.date",Peripheral.MF4)) %>% 
#   pivot_longer(cols = all_of(Peripheral.MF4),
#                names_to = "Peripheral",
#                values_to = "EM.num")
# MF4.EM_by.date$Peripheral <- 
#   factor(MF4.EM_by.date$Peripheral, levels = Peripheral.MF4)
# # ADFのみ抽出
# MF4.EM_by.date.ADF <- 
#   MF4.EM_by.date %>% 
#   dplyr::filter(Peripheral %in% c("COOK","CATHERINE"))
# # FINのみ抽出
# MF4.EM_by.date.FIN <- 
#   MF4.EM_by.date %>% 
#   dplyr::filter(Peripheral %in% c("VOLGA","AMUR"))

# フォント定義
par(family="Noto Sans")
# グラフタイトル
if (parm == "lambda_exp") {
  graph_title <- str_c("Metis-MF4 EM件数時系列分析 λ(μ + γ + tvf + r)：すべての成分を含んだ状態推定値")
  fill.graph <- "lightblue"
}else if (parm == "y_pred") {
  graph_title <- str_c("Metis-MF4 EM件数時系列分析 95%予測区間（MF3のパラメータでの予測）")
  fill.graph <- "lightgreen"
}
# 縦軸名称
y_label = "件数"
# 凡例設定
vec <- c()
for (i in 1:length(Peripheral.MF4)) {
  vec <- 
    vec %>% 
    append(str_c(i,".",Peripheral.MF4[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
vec.ADF <- c()
for (i in 1:2) {
  vec.ADF <- 
    vec.ADF %>% 
    append(str_c(Peripheral.MF4[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
vec.FIN <- c()
for (i in 3:4) {
  vec.FIN <- 
    vec.FIN %>% 
    append(str_c(Peripheral.MF4[i]," (自己回帰係数 b_ar = ",b_ar.mean[i],")"))
}
# ADF+FIN
Att.labs <- vec 
names(Att.labs) <- Peripheral.MF4
# グラフの表示期間の設定
limits = c(as.POSIXct("2023-03-01 JST"), as.POSIXct(MF4.MIF_by.date$X.date[end_day.MF4] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,30,5)

result_df$fit
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
  geom_point(alpha = 1.0, size = 2.0,
             data = MF4.EM_by.date, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  geom_point(alpha = 1.0, size = 2.0,
             data = MF4.6.EM_by.date, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  # geom_vline(xintercept = as.POSIXct("2023-05-01 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  # scale_y_continuous(breaks = breaks) +
  coord_cartesian(ylim = c(0,20)) +
  geom_vline(xintercept=as.POSIXct("2023-06-01 JST"), colour = "darkblue",linetype=2,alpha=0.7,linewidth=0.3) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-06-01 JST"), y = Inf,
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】",
           size = 3,
           colour = "darkblue")

plot(p)

# 予測区間+取得データグラフ保存
n <- 4
ggsave(str_c("./PDF/",title,".EM件数 予測区間-",start_day.MF4,"~",end_day.MF4," c(70,20,10,8)+計測値.pdf"),
       plot = p, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))

# ADF
Att.labs <- vec.ADF
names(Att.labs) <- Peripheral.MF4[1:2]
# グラフの表示期間の設定
limits = c(as.POSIXct("2023-03-01 JST"), as.POSIXct(MF4.MIF_by.date$X.date[end_day.MF4] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,4,1)
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
             data = MF4.EM_by.date.ADF, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  geom_point(alpha = 1.0, size = 0.9,
             data = MF4.6.EM_by.date.ADF, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  # geom_vline(xintercept = as.POSIXct("2023-03-31 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  # scale_y_continuous(breaks = breaks) +
  coord_cartesian(ylim = c(0,20)) +
  geom_vline(xintercept=as.POSIXct("2023-06-01 JST"), colour = "darkblue",linetype=2,alpha=0.7,linewidth=0.3) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-06-01 JST"), y = Inf,
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】",
           size = 3,
           colour = "darkblue")

plot(p.ADF)

# 予測区間+取得データグラフ保存
n <- 2
ggsave(str_c("./PDF/",title,".ADF_EM件数 予測区間-",start_day.MF4,"~",end_day.MF4," c(70,20,10,8)+計測値.pdf"), 
       plot = p.ADF, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))

# FIN
Att.labs <- vec.FIN
names(Att.labs) <- Peripheral.MF4[3:4]
# グラフの表示期間の設定
limits = c(as.POSIXct("2023-03-01 JST"), as.POSIXct(MF4.MIF_by.date$X.date[end_day.MF4] + 3600*24*pred_term))
# ｙ軸目盛設定
breaks=seq(0,4,1)
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
             data = MF4.EM_by.date.FIN, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  geom_point(alpha = 1.0, size = 0.9,
             data = MF4.6.EM_by.date.FIN, mapping =  aes(x = X.date, y = EM.num, color = Peripheral)) +
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 14),
        axis.text.y = element_text(size = 16)) +
  theme(legend.position = "none") +
  theme(strip.background = element_blank(), strip.text = element_text(size = 16, hjust = 0)) +
  facet_wrap(~Peripheral, ncol = 1, labeller = labeller(Peripheral = Att.labs)) +
  # geom_vline(xintercept = as.POSIXct("2023-03-31 JST"), linetype = 2, color = "darkblue") + # 垂直線
  scale_x_datetime(limits = limits, date_breaks = "1 month", date_labels = "%Y/%m") + # 表示期間の設定
  # scale_y_continuous(breaks = breaks) +
  coord_cartesian(ylim = c(0,5)) +
  geom_vline(xintercept=as.POSIXct("2023-06-01 JST"), colour = "darkblue",linetype=2,alpha=0.7,linewidth=0.3) +
  # 未来予測期間の文字表示
  annotate("text", x = as.POSIXct("2023-06-01 JST"), y = Inf,
           hjust = -0.01, # 文字のｘ方向位置調整 マイナスで右へ移動
           vjust = 2, 　　# 文字のｙ方向位置調整 プラスで下へ移動
           label = "⇒ 【未来予測区間】",
           size = 3,
           colour = "darkblue")

plot(p.FIN)

# 予測区間+取得データグラフ保存
n <- 2
ggsave(str_c("./PDF/",title,".FIN_EM件数 予測区間-",start_day.MF4,"~",end_day.MF4," c(70,20,10,8)+計測値.pdf"), 
       plot = p.FIN, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))
