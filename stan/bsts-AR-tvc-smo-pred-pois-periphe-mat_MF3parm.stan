// MF3のパラメータでのMF4の未来予測
// DGLMによるポアソン分布時系列分析 
// 再パラメータ化による収束性の向上を図っている
// 複数データのサンプリングを可能にしている
// 複数データごとに説明変数である稼働台数を設定可能である

// functions {
//   // 並列処理用関数
//   real partial_sum_lpmf (array[] int slice_y, int start, int end, vector lambda){
//     return poisson_log_lpmf(slice_y | lambda[start:end]);
//   }
// }

data {
  // 共通
  int N;                               // データ数（項目数）
  int pred_term;                       // 予測期間の長さ
  array[N] real<lower=0> Num_gain;     // 増加台数予測値（1ステップ当たりの増加予測数）
  int<lower=1> grainsize;              // 粒度　250?
  // MF3
  int T3;                              // データ取得期間の長さ（時系列数）
  array[N,T3] int<lower=0> y3;         // 観測値 （Nは項目数,Tは時系列数）
  array[N] vector<lower=0>[T3] Num3;   // 稼働台数
  array[N] real<lower=0> N3_max;       // 稼働台数を正規化時の除算数
  // MF4
  int T4;                              // データ取得期間の長さ（時系列数）
  array[N,T4] int<lower=0> y4;         // 観測値 （Nは項目数,Tは時系列数）
  array[N] vector<lower=0>[T4] Num4;   // 稼働台数
  array[N] real<lower=0> N4_max;       // 稼働台数を正規化時の除算数
  
}

parameters {
  // MF3
  vector<lower=0>[N] s_z3;      // ドリフト成分の変動の大きさを表す標準偏差
  vector<lower=0>[N] s_s3;      // 季節変動の大きさを表す標準偏差
  vector<lower=0>[N] s_r3;      // ランダム効果の標準偏差
  vector<lower=0>[N] s_t3;      // 時変係数（稼働台数）の標準偏差
  vector[N] b_ar3;              // 自己回帰項の係数
  vector[N] Intercept3;         // 切片
  array[N] vector[T3] r_raw3;           // ランダム効果の再パラメータ化
  array[N] vector[T3] mu_err3;          // 水準成分の増減量
  array[N] vector[T3] b_ope_err3;       // 稼働台数項の係数の増減量
  array[N] vector[T3-6] gamma_err3;     // 季節成分の推定値
  // MF4
  vector<lower=0>[N] s_z4;      // ドリフト成分の変動の大きさを表す標準偏差
  vector<lower=0>[N] s_s4;      // 季節変動の大きさを表す標準偏差
  vector<lower=0>[N] s_r4;      // ランダム効果の標準偏差
  vector<lower=0>[N] s_t4;      // 時変係数（稼働台数）の標準偏差
  vector[N] b_ar4;              // 自己回帰項の係数
  vector[N] Intercept4;         // 切片
  array[N] vector[T4] r_raw4;           // ランダム効果の再パラメータ化
  array[N] vector[T4] mu_err4;          // 水準成分の増減量
  array[N] vector[T4] b_ope_err4;       // 稼働台数項の係数の増減量
  array[N] vector[T4-6] gamma_err4;     // 季節成分の推定値
}

transformed parameters {
  // MF3
  array[N] vector[T3] mu3;              // 水準+ドリフト成分の推定値
  array[N] vector[T3] b_ope3;           // 稼働台数項の係数
  array[N] vector[T3] lambda3;          // 観測値の期待値のlogをとった値
  // vector[T] lambda[N];             // 観測値の期待値のlogをとった値 旧表記
  array[N] vector[T3-6] sumgamma3;      // 周期成分
  array[N] vector[T3] tvc3;             // 時変係数*稼働台数成分
  array[N] vector[T3] r3;               // ランダム効果
  array[N] vector[T3] gamma3;           // 季節成分の推定値
  // MF4
  array[N] vector[T4] mu4;              // 水準+ドリフト成分の推定値
  array[N] vector[T4] b_ope4;           // 稼働台数項の係数
  array[N] vector[T4] lambda4;          // 観測値の期待値のlogをとった値
  // vector[T] lambda[N];             // 観測値の期待値のlogをとった値 旧表記
  array[N] vector[T4-6] sumgamma4;      // 周期成分
  array[N] vector[T4] tvc4;             // 時変係数*稼働台数成分
  array[N] vector[T4] r4;               // ランダム効果
  array[N] vector[T4] gamma4;           // 季節成分の推定値
  // MF#
  for (n in 1:N){
    // ランダム効果成分 r
    r3[n,] = s_r3[n] * r_raw3[n,];   // 再パラメータ化
    r4[n,] = s_r4[n] * r_raw4[n,];   // 再パラメータ化
    // 平準化トレンドモデル水準成分 μ, 稼働台数時変係数 b_ope
    // 1時点目
    mu3[n,1] = Intercept3[n] + s_z3[n] * mu_err3[n,1];
    b_ope3[n,1] = s_t3[n] * b_ope_err3[n,1];
    mu4[n,1] = Intercept4[n] + s_z4[n] * mu_err4[n,1];
    b_ope4[n,1] = s_t4[n] * b_ope_err4[n,1];
    // 2時点目
    mu3[n,2] = Intercept3[n] + b_ar3[n] * mu3[n,1] + s_z3[n] * mu_err3[n,2];
    mu4[n,2] = Intercept4[n] + b_ar4[n] * mu4[n,1] + s_z4[n] * mu_err4[n,2];
    // 2時点目以降
    for (i in 2:T3) {
      b_ope3[n,i] = b_ope3[n,(i-1)] + s_t3[n] * b_ope_err3[n,i];  // 再パラメータ化
    }
    for (i in 2:T4) {
      b_ope4[n,i] = b_ope4[n,(i-1)] + s_t4[n] * b_ope_err4[n,i];  // 再パラメータ化
    }
    // 3時点目以降
    for (i in 3:T3) {
      mu3[n,i] = Intercept3[n] + b_ar3[n] * (2 * mu3[n,(i-1)] - mu3[n,(i-2)]) + s_z3[n] * mu_err3[n,i];  // 再パラメータ化
    }
    for (i in 3:T4) {
      mu4[n,i] = Intercept4[n] + b_ar4[n] * (2 * mu4[n,(i-1)] - mu4[n,(i-2)]) + s_z4[n] * mu_err4[n,i];  // 再パラメータ化
    }
    // 周期成分 γ
    gamma3[n,1:6] = s_s3[n] * gamma_err3[n,1:6];
    for(i in 1:T3-6) {
      sumgamma3[n,i] = -sum(gamma3[n,(i):(i+5)]);
      gamma3[n,(i+6)] = sumgamma3[n,i] + s_s3[n] * gamma_err3[n,i]; // 再パラメータ化
    }
    gamma4[n,1:6] = s_s4[n] * gamma_err4[n,1:6];
    for(i in 1:T4-6) {
      sumgamma4[n,i] = -sum(gamma4[n,(i):(i+5)]);
      gamma4[n,(i+6)] = sumgamma4[n,i] + s_s4[n] * gamma_err4[n,i]; // 再パラメータ化
    }
    // 時変係数成分
    for(i in 1:T3) {
      tvc3[n,i] = b_ope3[n,i] * Num3[n,i];
    }
    for(i in 1:T4) {
      tvc4[n,i] = b_ope4[n,i] * Num4[n,i];
    }
    // 状態推定値 λ
    lambda3[n,1:T3] = mu3[n,1:T3] + gamma3[n,1:T3] + tvc3[n,1:T3] + r3[n,1:T3];
    lambda4[n,1:T4] = mu4[n,1:T4] + gamma4[n,1:T4] + tvc4[n,1:T4] + r4[n,1:T4];
  }
}

model {
  // 弱情報事前分布
  // MF3
  s_z3 ~ student_t(5, 0, 1);
  s_s3 ~ student_t(5, 0, 1);
  s_r3 ~ student_t(5, 0, 1);
  s_t3 ~ student_t(5, 0, 1);
  b_ar3 ~ student_t(3, 0, 1);
  Intercept3 ~ student_t(3, 0, 10);
  // MF4
  s_z4 ~ student_t(5, 0, 1);
  s_s4 ~ student_t(5, 0, 1);
  s_r4 ~ student_t(5, 0, 1);
  s_t4 ~ student_t(5, 0, 1);
  b_ar4 ~ student_t(3, 0, 1);
  Intercept4 ~ student_t(3, 0, 10);
  // データ数（N）ごとにサンプリングする
  for (n in 1:N) {
    // 時点ごとに加わるランダム効果
    // r ~ normal(0, s_r);
    r_raw3[n,] ~ normal(0,1);
    r_raw4[n,] ~ normal(0,1);
    // 稼働台数の時変係数
    // b_ope[2:T] ~ normal(b_ope[1:T-1], s_t);
    b_ope_err3[n,] ~ normal(0,1);
    b_ope_err4[n,] ~ normal(0,1);
    // 水準+自己回帰成分
    // mu[3:T] ~ normal(Intercept + b_ar * (2 * mu[2:T-1] - mu[1:T-2]), s_w);
    // mu_err[2]以降は、標準偏差1の正規乱数を得る
    mu_err3[n,] ~ normal(0,1);
    mu_err4[n,] ~ normal(0,1);
    // 季節成分
    // gamma[7:T] ~ normal(sumgamma[1:T-6], s_s);
    gamma_err3[n,1:T3-6] ~ normal(0,1);
    gamma_err4[n,1:T4-6] ~ normal(0,1);
  }
  // 観測方程式に従い観測値が得られる
  // 並列処理
  for (n in 1:N){
  //   target += reduce_sum(partial_sum_lpmf, y3[n,], grainsize, lambda3[n,]);
  //   target += reduce_sum(partial_sum_lpmf, y4[n,], grainsize, lambda4[n,]);
    // 並列処理を行わない場合
    y3[n,] ~ poisson_log(lambda3[n,]);
    y4[n,] ~ poisson_log(lambda4[n,]);
  }
}

generated quantities {
  // 予測範囲のパラメータ設定
  array[N] vector[T4+pred_term] mu_pred;                  // 水準+ドリフト成分の推定値
  array[N] vector[T4+pred_term] b_ope_pred;               // 稼働台数項の係数
  array[N] vector[T4+pred_term] lambda_pred;              // 観測値の期待値
  array[N] vector[T4-6+pred_term] sumgamma_pred;          // 周期成分
  array[N] vector[T4+pred_term] tvc_pred;                 // 時変係数*稼働台数成分
  array[N] vector[T4+pred_term] r_pred;                   // ランダム効果
  array[N] vector[T4+pred_term] gamma_pred;               // 季節成分の推定値
  array[N] vector[T4+pred_term] Num_pred;                          // 稼働台数

  // データ取得期間は元の変数と同値なので,そのまま代入する
  for (n in 1:N) {
    mu_pred[n,1:T4] = mu4[n,1:T4];
    b_ope_pred[n,1:T4] = b_ope4[n,1:T4];
    lambda_pred[n,1:T4] = lambda4[n,1:T4];
    sumgamma_pred[n,1:(T4-6)] = sumgamma4[n,1:(T4-6)];
    tvc_pred[n,1:T4] = tvc4[n,1:T4];
    r_pred[n,1:T4] = r4[n,1:T4];
    gamma_pred[n,1:T4] = gamma4[n,1:T4];
    Num_pred[n,1:T4] = Num4[n,1:T4];
  }

  // データ数（N）ごとに未来予測期間を算出
  for (n in 1:N) {
    for (i in 1:pred_term) {
      // ランダム効果成分 r
      r_pred[n,(T4+i)] = normal_rng(0, s_r3[n]);
      // 水準成分 μ
      mu_pred[n,(T4+i)] = normal_rng(Intercept4[n] + b_ar4[n] * (2 * mu_pred[n,(T4+i-1)] - mu_pred[n,(T4+i-2)]), s_z3[n]);
      // 稼働台数時変係数 b_ope
      b_ope_pred[n,(T4+i)] = normal_rng(b_ope_pred[n,(T4+i-1)], s_t3[n]);
      // 周期成分 γ
      sumgamma_pred[n,(T4-6+i)] = -sum(gamma_pred[n,(T4-6+i):(T4-6+i+5)]);
      gamma_pred[n,(T4+i)] = normal_rng(sumgamma_pred[n,(T4-6+i)], s_s3[n]);
      // 想定稼働台数
      Num_pred[n,T4+i] = Num_pred[n,T4+i-1] + Num_gain[n];
      // 時変係数成分 tvc
      tvc_pred[n,(T4+i)] = b_ope_pred[n,(T4+i)] * Num_pred[n,T4+i];
      // 状態推定値 λ
      lambda_pred[n,(T4+i)] = mu_pred[n,(T4+i)] + gamma_pred[n,(T4+i)] + tvc_pred[n,(T4+i)] + r_pred[n,(T4+i)];
    }
  }
  // パラメータ設定
  array[N] vector[T4+pred_term] lambda_exp;             // 状態推定値(EXP)
  array[N] vector[T4+pred_term] mu_gamma_exp;           // ランダム効果を除いた状態推定値
  array[N] vector[T4+pred_term] mu_r_exp;               // 周期成分を除いた状態推定値
  array[N] vector[T4+pred_term] mu_tvc_exp;             // 水準成分+時変係数成分
  array[N] vector[T4+pred_term] mu_exp;                 // 周期成分・ランダム効果を除いた状態推定値
  array[N] vector[T4+pred_term] gamma_exp;              // 周期成分(EXP)
  array[N] vector[T4+pred_term-31] delta_exp;           // ドリフト成分(EXP)
  array[N] vector[T4+pred_term] r_exp;                  // ランダム成分(EXP)
  array[N] vector[T4+pred_term] tvc_exp;                // 時変係数*稼働台数成分
  array[N] vector[T4+pred_term] mu_tvc_1_exp;           // 稼働台数を1台に設定(EXP)
  array[N] vector[T4+pred_term] mu_tvc_1;               // 稼働台数を1台に設定
  array[N] vector[T4+pred_term] tvc_1_exp;              // 稼働台数時変係数(EXP)
  array[N] vector[T4+pred_term] tvc_1;                  // 稼働台数時変係数,EM率相当
  array[N,T4+pred_term] int<lower=0> y_pred;            // 観測値の予測値

  for (n in 1:N) {
    // 指数関数を作用させる
    lambda_exp[n,] = exp(lambda_pred[n,]);
    mu_gamma_exp[n,] = exp(mu_pred[n,] + gamma_pred[n,]);
    mu_r_exp[n,] = exp(mu_pred[n,] + r_pred[n,]);
    mu_tvc_exp[n,] = exp(mu_pred[n,] + tvc_pred[n,]);
    mu_exp[n,] = exp(mu_pred[n,]);
    gamma_exp[n,] = exp(gamma_pred[n,]);
    delta_exp[n,1:(T4+pred_term-31)] = exp(mu_pred[n,32:(T4+pred_term)]-mu_pred[n,31:(T4-1+pred_term)]);
    r_exp[n,] = exp(r_pred[n,]);
    tvc_exp[n,] = exp(tvc_pred[n,]);
    mu_tvc_1_exp[n,] = exp(mu_pred[n,] + b_ope_pred[n,] / N4_max[n]);
    mu_tvc_1[n,] = mu_pred[n,] + b_ope_pred[n,] / N4_max[n];
    tvc_1_exp[n,] = exp(b_ope_pred[n,] / N4_max[n]);
    tvc_1[n,] = b_ope_pred[n,] / N4_max[n];
    // 観測値の予測値算出
    y_pred[n,] = poisson_log_rng(lambda_pred[n,]);
  }
}
