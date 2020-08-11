#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第7回「世代重複マクロモデル」のソースコードのRへの野良移植: 資産と消費のプロファイル, 図3
#
#==========================================================================
# MATLAB CODE WRITTEN FOR CH7 OF KEIZAI SEMINAR
# WRITTEN BY SAGIRI KITAO
# 09/14/2019
#
# ported to R by uncorrelated
#==========================================================================
#

source("common.R")

# ベースラインのパラメーター
param_baseline <- list(
	alpha = 0.40, 
	beta = 0.98,
	delta = 0.08,
	Nj = 61, # LIVE 60 YEARS MAX (80YRS OLD)
	NjW = 45, # WORKING YEARS (RETIRE AT NjW+1) : ENTER AT 21 AND RETIRE AT 65
	rho = 0.5 # SS REPLACEMENT RATE (0.0 OR 0.5)
)
param_baseline$meaJ <- rep(1/param_baseline$Nj, param_baseline$Nj) # AGE DISTRIBUTION (no population growth) 

r_baseline <- calcOLG(param_baseline)

ageA <- (1:(param_baseline$Nj + 1)) + 19
ageC <- (1:(param_baseline$Nj)) + 19

norm <- r_baseline$c[1]

# プロットする
xlim <- c(min(ageA), max(ageA))
par(oma=c(0, 0, 0, 0), mfrow=c(1,2), mar=c(4.5, 3, 2, 1), bg="white") # mfrow=c(1,2)は、2つのプロットを縦2分割して表示するためのオプション

plot(ageA, r_baseline$a/norm, type="l", main=expression( paste("Assets: ", a[j]) ), xlab="Age", ylab="", lwd=2)
plot(ageC, r_baseline$c/norm, type="l", main="Consumption", xlab="Age", ylab="", lwd=2)

# 図を保存する
dev.copy2eps(file="fig 3.eps", width=8, height=4)
# dev.copy(png, "fig 3.png", width=800, height=400, type="cairo", bg="white"); dev.off()

