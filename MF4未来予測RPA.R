fit <- MF4Future.pred.QISS(
  path.excel = "./excel", 　　　　　　 # エクセルファイルのフォルダーパス
  pattern1 = "〇共通_保守*.+\\.xlsx$", # 保守データ
  pattern2 = "〇共通_周辺機装着情報*.+\\.csv$", # 周辺機装着情報
  pattern3 = "〇共通_機器*.+\\.xlsx$", # 機器データ
  EM = "Call",                         # "EM" or "Call"
  save = TRUE,        　　　　　　　　 # tsvファイル保存の有無
  path.tsv = "./tsv_data/", 　　　　　 # tsv保存フォルダーパス
  tsv.name.EM = "Metis_MF4_EM.tsv",    # 保守EMデータのtsvファイル名
  tsv.name.MIF = "Metis_MF4_MIF.tsv",  # MIFデータのtsvファイル名
  manuf.start.month = "2022/12/01",　　# 製造開始月
  # MF3データ/変更不要
  MF3_maintenance.tsv = "./tsv_data/Metis_MF3_2303.tsv",
  MF3_machine.tsv = "./tsv_data/Metis_MIF_2303.tsv",
  # MCMCサンプリング開始日
  start_day.MF4.JST = "2023-03-01 JST",
  # 予測期間
  pred_term = 30,
  # 日次での稼働台数増加予測数 c(COOK,SINAI,VOLGA,AMUR中綴じ)
  Num_gain = c(70,20,10,8),
  # csv保存ディレクトリ
  output_dir = "./csv/matrix",
  # MCMCサンプリング
  chains = 6,
  parallel_chains = getOption("mc.cores", 24),
  threads_per_chain = 2,
  iter_warmup = 2000, #6000,
  iter_sampling = 2000, #2000,
  thin = 4,
  max_treedepth = 12, #12
  refresh = 100,
  # タイトル
  title = "Metis-MF4-未来予測10月-2",
  # サンプリング時間保存
  save_exTime = TRUE,
  # サンプリング結果保存
  save_fit = TRUE,
  # グラフｙ軸目盛設定/NAにて自動設定するが、小数点表示を避ける場合は強制的に設定可能
  breaks = NA, 　　# seq(0,15,5)
  breaks.ADF = NA, # seq(0,15,5)
  breaks.FIN = NA, # seq(0,4,1)
  # グラフ保存
  save_graph = TRUE,
  save_graph.ADF = TRUE,
  save_graph.FIN = TRUE
)
