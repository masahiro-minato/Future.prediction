fit <- MF4Future.pred(
  # 入力データ/tsv形式
  MF4_maintenance.tsv = "./tsv_data/Metis_MF4_2305.tsv",
  MF4_machine.tsv = "./tsv_data/MF4_MIF_2305.tsv",
  MF3_maintenance.tsv = "./tsv_data/Metis_MF3_2303.tsv",
  MF3_machine.tsv = "./tsv_data/Metis_MIF_2303.tsv",
  # MCMCサンプリング開始日と終了日
  start_day.MF4 = "2023-03-01 JST",
  end_day.MF4 = "2023-05-31 JST",
  # 予測期間
  pred_term = 30,
  # 日次での稼働台数増加予測数 c(COOK,SINAI,VOLGA,AMUR中綴じ)
  Num_gain = c(70,20,10,8),
  # csv保存ディレクトリ
  output_dir = "./csv/matrix",
  # 並列処理
  reduce_sum = FALSE,
  # MCMCサンプリング
  chains = 4,
  parallel_chains = getOption("mc.cores", 10),
  threads_per_chain = 2,
  iter_warmup = 1000, #6000,
  iter_sampling = 1000, #2000,
  thin = 4,
  max_treedepth = 10, #12
  # タイトル
  title = "Metis-MF4-Peripheral.MF3parm_tmp3",
  # サンプリング時間保存
  save_exTime = TRUE,
  # サンプリング結果保存
  save_fit = TRUE,
  # グラフ保存
  save_graph = TRUE)


p <- readRDS("./rds/Metis-MF4-Peripheral.MF3parm_tmp3.EM件数_予測区間-29~120_[70,20,10,8].rds")
p
