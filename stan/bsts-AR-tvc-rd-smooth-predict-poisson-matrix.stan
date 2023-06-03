// DGLMによるポアソン分布時系列分析 
// 再パラメータ化による収束性の向上を図っている
// 複数データのサンプリングを可能にしている
functions {
  // 並列処理用関数
  real partial_sum_lpmf (array[] int slice_y, int start, int end, vector lambda){
    return poisson_log_lpmf(slice_y | lambda[start:end]);
  }
}

data {
  int T;                      // データ取得期間の長さ（時系列数）
  int N;                      // データ数（項目数）
  int pred_term;              // 予測期間の長さ
  array[N,T] int<lower=0> y;  // 観測値 （Nは項目数,Tは時系列数）
  // int y[N,T];              // 観測値 旧表記
  vector[T] Num;              // 稼働台数
  real N_max;                 // 稼働台数を正規化時の除算数
  int Num_gain;               // 増加台数予測値（1ステップ当たりの増加予測数）
  int grainsize;              // 粒度　250?
}

parameters {
  vector<lower=0>[N] s_z;      // ドリフト成分の変動の大きさを表す標準偏差
  vector<lower=0>[N] s_s;      // 季節変動の大きさを表す標準偏差
  vector<lower=0>[N] s_r;      // ランダム効果の標準偏差
  vector<lower=0>[N] s_t;      // 時変係数（稼働台数）の標準偏差
  vector[N] b_ar;              // 自己回帰項の係数
  vector[N] Intercept;         // 切片
  array[N] vector[T] r_raw;           // ランダム効果の再パラメータ化
  array[N] vector[T] mu_err;          // 水準成分の増減量
  array[N] vector[T] b_ope_err;       // 稼働台数項の係数の増減量
  array[N] vector[T-6] gamma_err;     // 季節成分の推定値
}

transformed parameters {
  array[N] vector[T] mu;              // 水準+ドリフト成分の推定値
  array[N] vector[T] b_ope;           // 稼働台数項の係数
  array[N] vector[T] lambda;          // 観測値の期待値のlogをとった値
  // vector[T] lambda[N];             // 観測値の期待値のlogをとった値 旧表記
  array[N] vector[T-6] sumgamma;      // 周期成分
  array[N] vector[T] tvc;             // 時変係数*稼働台数成分
  array[N] vector[T] r;               // ランダム効果
  array[N] vector[T] gamma;           // 季節成分の推定値
  
  for (n in 1:N){
    // ランダム効果成分 r
    r[n,] = s_r[n] * r_raw[n,];   // 再パラメータ化
    // 平準化トレンドモデル水準成分 μ, 稼働台数時変係数 b_ope
    // 1時点目
    mu[n,1] = Intercept[n] + s_z[n] * mu_err[n,1];
    b_ope[n,1] = s_t[n] * b_ope_err[n,1];
    // 2時点目
    mu[n,2] = Intercept[n] + b_ar[n] * mu[n,1] + s_z[n] * mu_err[n,2];
    // 2時点目以降
    for (i in 2:T) {
      b_ope[n,i] = b_ope[n,(i-1)] + s_t[n] * b_ope_err[n,i];  // 再パラメータ化
    }
    // 3時点目以降
    for (i in 3:T) {
      mu[n,i] = Intercept[n] + b_ar[n] * (2 * mu[n,(i-1)] - mu[n,(i-2)]) + s_z[n] * mu_err[n,i];  // 再パラメータ化
    }
    // 周期成分 γ
    gamma[n,1:6] = s_s[n] * gamma_err[n,1:6];
    for(i in 1:T-6) {
      sumgamma[n,i] = -sum(gamma[n,(i):(i+5)]);
      gamma[n,(i+6)] = sumgamma[n,i] + s_s[n] * gamma_err[n,i]; // 再パラメータ化
    }
    // 時変係数成分
    for(i in 1:T) {
      tvc[n,i] = b_ope[n,i] * Num[i];
    }
    // 状態推定値 λ
    lambda[n,1:T] = mu[n,1:T] + gamma[n,1:T] + tvc[n,1:T] + r[n,1:T];
  }
}

model {
  // 弱情報事前分布
  s_z ~ student_t(5, 0, 1);
  s_s ~ student_t(5, 0, 1);
  s_r ~ student_t(5, 0, 1);
  s_t ~ student_t(5, 0, 1);
  b_ar ~ student_t(3, 0, 1);
  Intercept ~ student_t(3, 0, 10);
  // データ数（N）ごとにサンプリングする
  for (n in 1:N) {
    // 時点ごとに加わるランダム効果
    // r ~ normal(0, s_r);
    r_raw[n,] ~ normal(0,1);
    // 稼働台数の時変係数
    // b_ope[2:T] ~ normal(b_ope[1:T-1], s_t);
    b_ope_err[n,] ~ normal(0,1);
    // 水準+自己回帰成分
    // mu[3:T] ~ normal(Intercept + b_ar * (2 * mu[2:T-1] - mu[1:T-2]), s_w);
    // mu_err[2]以降は、標準偏差1の正規乱数を得る
    mu_err[n,] ~ normal(0,1);
    // 季節成分
    // gamma[7:T] ~ normal(sumgamma[1:T-6], s_s);
    gamma_err[n,1:T-6] ~ normal(0,1);
  }
  // 観測方程式に従い観測値が得られる
  // 並列処理
  for (n in 1:N){
    target += reduce_sum(partial_sum_lpmf, y[n,], grainsize, lambda[n,]);
    // 並列処理を行わない場合
    // y[n,] ~ poisson_log(lambda[n,]);
  }
}

generated quantities {
  // 予測範囲のパラメータ設定
  array[N] vector[T+pred_term] mu_pred;                  // 水準+ドリフト成分の推定値
  array[N] vector[T+pred_term] b_ope_pred;               // 稼働台数項の係数
  array[N] vector[T+pred_term] lambda_pred;              // 観測値の期待値
  array[N] vector[T-6+pred_term] sumgamma_pred;          // 周期成分
  array[N] vector[T+pred_term] tvc_pred;                 // 時変係数*稼働台数成分
  array[N] vector[T+pred_term] r_pred;                   // ランダム効果
  array[N] vector[T+pred_term] gamma_pred;               // 季節成分の推定値
  vector[T+pred_term] Num_pred;                          // 稼働台数

  // データ取得期間は元の変数と同値なので,そのまま代入する
  for (n in 1:N) {
    mu_pred[n,1:T] = mu[n,1:T];
    b_ope_pred[n,1:T] = b_ope[n,1:T];
    lambda_pred[n,1:T] = lambda[n,1:T];
    sumgamma_pred[n,1:(T-6)] = sumgamma[n,1:(T-6)];
    tvc_pred[n,1:T] = tvc[n,1:T];
    r_pred[n,1:T] = r[n,1:T];
    gamma_pred[n,1:T] = gamma[n,1:T];
  }
  Num_pred[1:T] = Num[1:T];
  // 想定稼働台数
  for (i in 1:pred_term) {
    Num_pred[T+i] = Num_pred[T+i-1] + Num_gain;
  }
  // データ数（N）ごとに未来予測期間を算出
  for (n in 1:N) {
    for (i in 1:pred_term) {
      // ランダム効果成分 r
      r_pred[n,(T+i)] = normal_rng(0, s_r[n]);
      // 水準成分 μ
      mu_pred[n,(T+i)] = normal_rng(Intercept[n] + b_ar[n] * (2 * mu_pred[n,(T+i-1)] - mu_pred[n,(T+i-2)]), s_z[n]);
      // 稼働台数時変係数 b_ope
      b_ope_pred[n,(T+i)] = normal_rng(b_ope_pred[n,(T+i-1)], s_t[n]);
      // 周期成分 γ
      sumgamma_pred[n,(T-6+i)] = -sum(gamma_pred[n,(T-6+i):(T-6+i+5)]);
      gamma_pred[n,(T+i)] = normal_rng(sumgamma_pred[n,(T-6+i)], s_s[n]);
      // 時変係数成分 tvc
      tvc_pred[n,(T+i)] = b_ope_pred[n,(T+i)] * Num_pred[T+i];
      // 状態推定値 λ
      lambda_pred[n,(T+i)] = mu_pred[n,(T+i)] + gamma_pred[n,(T+i)] + tvc_pred[n,(T+i)] + r_pred[n,(T+i)];
    }
  }
  // パラメータ設定
  array[N] vector[T+pred_term] lambda_exp;             // 状態推定値(EXP)
  array[N] vector[T+pred_term] mu_gamma_exp;           // ランダム効果を除いた状態推定値
  array[N] vector[T+pred_term] mu_r_exp;               // 周期成分を除いた状態推定値
  array[N] vector[T+pred_term] mu_tvc_exp;             // 水準成分+時変係数成分
  array[N] vector[T+pred_term] mu_exp;                 // 周期成分・ランダム効果を除いた状態推定値
  array[N] vector[T+pred_term] gamma_exp;              // 周期成分(EXP)
  array[N] vector[T+pred_term-31] delta_exp;           // ドリフト成分(EXP)
  array[N] vector[T+pred_term] r_exp;                  // ランダム成分(EXP)
  array[N] vector[T+pred_term] tvc_exp;                // 時変係数*稼働台数成分
  array[N] vector[T+pred_term] mu_tvc_1_exp;           // 稼働台数を1台に設定(EXP)
  array[N] vector[T+pred_term] mu_tvc_1;               // 稼働台数を1台に設定
  array[N] vector[T+pred_term] tvc_1_exp;              // 稼働台数時変係数(EXP)
  array[N] vector[T+pred_term] tvc_1;                  // 稼働台数時変係数,EM率相当
  array[N,T+pred_term] int<lower=0> y_pred;            // 観測値の予測値

  for (n in 1:N) {
    // 指数関数を作用させる
    lambda_exp[n,] = exp(lambda_pred[n,]);
    mu_gamma_exp[n,] = exp(mu_pred[n,] + gamma_pred[n,]);
    mu_r_exp[n,] = exp(mu_pred[n,] + r_pred[n,]);
    mu_tvc_exp[n,] = exp(mu_pred[n,] + tvc_pred[n,]);
    mu_exp[n,] = exp(mu_pred[n,]);
    gamma_exp[n,] = exp(gamma_pred[n,]);
    delta_exp[n,1:(T+pred_term-31)] = exp(mu_pred[n,32:(T+pred_term)]-mu_pred[n,31:(T-1+pred_term)]);
    r_exp[n,] = exp(r_pred[n,]);
    tvc_exp[n,] = exp(tvc_pred[n,]);
    mu_tvc_1_exp[n,] = exp(mu_pred[n,] + b_ope_pred[n,] / N_max);
    mu_tvc_1[n,] = mu_pred[n,] + b_ope_pred[n,] / N_max;
    tvc_1_exp[n,] = exp(b_ope_pred[n,] / N_max);
    tvc_1[n,] = b_ope_pred[n,] / N_max;
    // 観測値の予測値算出
    y_pred[n,] = poisson_log_rng(lambda_pred[n,]);
  }
}
