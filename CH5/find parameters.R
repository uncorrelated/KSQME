#
# p_H=0のときに、y_Lとpi_Lが所定の値にようになるように、s_Lとkapperの値を求める
#

# 共通の設定を読み込む
if(!any(ls()=="dparam")){
	source("common.R")
}

# ニュートン・ラフソン法で条件にあうパラメーターを探す
library(nleqslv)
objf <- function(p){
	with(dparam, {

		# 所定の値3つ
		p_H <- 0
		y_L <- -7.0
		pi_L <- -1.0/4

		# 自明な解（r_L=0）も探させる
		s_L <- p[1]
		kapper <- p[2]
		y_H <- p[3]
		pi_H <- p[4]
		r_H <- p[5]
		r_L <- p[6]

		r <- numeric(length(p))

		# 根を探すアルゴリズムなので、左辺=右辺を、右辺-左辺と書く
		r[1] <- (1 - p_H)*y_H + p_H*y_L - (r_H - ((1 - p_H)*pi_H + p_H*pi_L) - s_H) - y_H
		r[2] <- kapper*y_H + beta*((1 - p_H)*pi_H + p_H*pi_L) - pi_H
		r[3] <- (1 - p_L)*y_H + p_L*y_L - (r_L - ((1 - p_L)*pi_H + p_L*pi_L) - s_L) - y_L
		r[4] <- kapper*y_L + beta*((1 - p_L)*pi_H + p_L*pi_L) - pi_L
		r[5] <- r_star + phi*((1 - p_H)*pi_H + p_H*pi_L) - r_H
		r[6] <- 0 - r_L

		r # このベクトルの0を目指して最適化される
	})
}

init <- c(-2, 0.01, 1, 0.1, 0, 0) # 初期値重要
r_nleqslv <- nleqslv(init, objf, method = c("Newton"), control = list(allowSingular = TRUE))
if(1 < r_nleqslv$termcd){
	stop(r_nleqslv$message) # 収束しないときはエラーメッセージを出して止まる
}

dparam$s_L <- r_nleqslv$x[1]
dparam$kapper <- r_nleqslv$x[2]

print(sprintf("s_L=%.6f, kapper=%.6f", dparam$s_L, dparam$kapper))

