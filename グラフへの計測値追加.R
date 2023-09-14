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
# 周辺機名称ベクトルの作成
Peripheral.MF4 <- substring(names(EM.count.MF4.Peripheral)[-1], 4,)
# 結合
X.date.MF4 <- seq(as.POSIXct("2023-02-01"), as.POSIXct("2023-06-30"), by = "day")
Date.MF4.Peripheral <- tibble(X.date = X.date.MF4)
MF4.EM_by.date <- 
  Date.MF4.Peripheral %>% 
  full_join(EM.count.MF4.Peripheral, by=c("X.date" = "Maintenance_date")) %>%
  mutate_if(is.numeric, ~replace(., is.na(.), 0)) %>% 
  set_colnames(c("X.date",Peripheral.MF4)) %>% 
  pivot_longer(cols = all_of(Peripheral.MF4),
               names_to = "Peripheral",
               values_to = "EM.num")

# レベル設定
MF4.EM_by.date$Peripheral <- 
  factor(MF4.EM_by.date$Peripheral, levels = Peripheral.MF4)

# 表示期間
start_day.MF4 = which(Date.MF4.Peripheral$X.date == "2023-03-01 JST")
end_day.MF4 = which(Date.MF4.Peripheral$X.date == "2023-05-31 JST")

# 未来予測区間の実測値データ
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

# グラフデータの読込
p <- readRDS("./rds/Metis-MF4-Peripheral.MF3parm_tmp2.EM件数_予測区間-29~120_[70,20,10,8].rds")
p.ADF <- readRDS("./rds/Metis-MF4-Peripheral.MF3parm_tmp2.ADF_EM件数_予測区間-29~120_[70,20,10,8].rds")
p.FIN <- readRDS("./rds/Metis-MF4-Peripheral.MF3parm_tmp2.FIN_EM件数_予測区間-29~120_[70,20,10,8].rds")

# フォント定義
par(family="Noto Sans")
# 未来予測区間への実測値プロット
p.actualresults <- 
  p + geom_point(alpha = 1.0, size = 3.0,
             　　data = MF4.6.EM_by.date, 
             　　mapping =  aes(x = X.date, y = EM.num, color = Peripheral), shape = 18)

p.ADF.actualresults <- 
  p.ADF + geom_point(alpha = 1.0, size = 3.0,
                 data = MF4.6.EM_by.date.ADF, 
                 mapping =  aes(x = X.date, y = EM.num, color = Peripheral), shape = 18)

p.FIN.actualresults <- 
  p.FIN + geom_point(alpha = 1.0, size = 3.0,
                 data = MF4.6.EM_by.date.FIN, 
                 mapping =  aes(x = X.date, y = EM.num, color = Peripheral), shape = 18)

# グラフ保存
n <- 4
ggsave(str_c("./PDF/Metis-MF4-Peripheral.MF3parm_tmp2.EM件数_予測区間-29~120_[70,20,10,8]＆計測値.pdf"),
       plot = p.actualresults, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))
n <- 2
ggsave(str_c("./PDF/Metis-MF4-Peripheral.MF3parm_tmp2.ADF_EM件数_予測区間-29~120_[70,20,10,8]＆計測値.pdf"),
       plot = p.ADF.actualresults, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))
n <- 2
ggsave(str_c("./PDF/Metis-MF4-Peripheral.MF3parm_tmp2.FIN_EM件数_予測区間-29~120_[70,20,10,8]＆計測値.pdf"),
       plot = p.FIN.actualresults, device = cairo_pdf, dpi=300, width=10, height=(n*2.5))
