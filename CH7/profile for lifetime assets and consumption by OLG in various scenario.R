#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第7回「世代重複マクロモデル」のソースコードのRへの野良移植: 複数シナリオの資産と消費のプロファイルの比較, 図4
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

# Scenario 1 のパラメーター
param_scenario_1 <- param_baseline
param_scenario_1$rho <- 0.25

# Scenario 2 のパラメーター
param_scenario_2 <- param_baseline
param_scenario_2$NjW <- 50

# Scenario 3 のパラメーター
param_scenario_3 <- param_baseline
param_scenario_3$Nj <- 66
param_scenario_3$meaJ <- rep(1/66, param_scenario_3$Nj)

# シナリオごとにそれぞれ計算
r_baseline <- calcOLG(param_baseline)
r_scenario_1 <- calcOLG(param_scenario_1)
r_scenario_2 <- calcOLG(param_scenario_2)
r_scenario_3 <- calcOLG(param_scenario_3)

# プロットする

# 資産aと消費cの計算期間が異なることに注意
ageA <- (1:(param_baseline$Nj + 1)) + 19
ageC <- (1:(param_baseline$Nj)) + 19

norm <- r_baseline$c[1]

xlim <- c(min(ageA), max(ageA))
lty <- c(1, 4, 3, 2)
lwd <- c(2, 2, 2, 2)
col <- c("black", "red", "purple", "blue")
ylim <- c(0, max(r_baseline$a, r_scenario_1$a, r_scenario_2$a, r_scenario_3$a))/norm
xlim <- c(19, param_scenario_3$Nj + 19)
par(oma=c(0, 0, 0, 0), mfrow=c(1,2), mar=c(4.5, 3, 2, 1), bg="white") # mfrow=c(1,2)は、2つのプロットを縦2分割して表示するためのオプション
plot(ageA, r_baseline$a/norm, type="l", main=expression( paste("Assets: ", a[j]) ), xlab="Age", ylab="", lty=lty[1], col=col[1], ylim=ylim, xlim=xlim, lwd=lwd[1])
lines(ageA, r_scenario_1$a/norm, lty=lty[2], col=col[2], lwd=lwd[2])
lines(ageA, r_scenario_2$a/norm, lty=lty[3], col=col[3], lwd=lwd[3])
lines((1:(param_scenario_3$Nj + 1)) + 19, r_scenario_3$a/norm, lty=lty[4], col=col[4], lwd=lwd[4])
legend("topleft", lty=lty, lwd=lwd, col=col, legend=c("Baseline", "Scenario 1", "Scenario 2", "Scenario 3"), bty="n", y.intersp=1.1, seg.len=2.5)

ylim <- c(min(r_baseline$c, r_scenario_1$c, r_scenario_2$c, r_scenario_3$c), max(r_baseline$c, r_scenario_1$c, r_scenario_2$c, r_scenario_3$c))
plot(ageC, r_baseline$c, type="l", main="Consumption", xlab="Age", ylab="", ylim=ylim, xlim=xlim, lwd=lwd[1])
lines(ageC, r_scenario_1$c, lty=lty[2], col=col[2], lwd=lwd[2])
lines(ageC, r_scenario_2$c, lty=lty[3], col=col[3], lwd=lwd[3])
lines((1:(param_scenario_3$Nj)) + 19, r_scenario_3$c, lty=lty[4], col=col[4], lwd=lwd[4])
legend("topleft", lty=lty, lwd=lwd, col=col, legend=c("Baseline", "Scenario 1", "Scenario 2", "Scenario 3"), bty="n", y.intersp=1.1, seg.len=2.5)

# 図を保存する
dev.copy2eps(file="fig 4.eps", width=8, height=4)
# dev.copy(png, "fig 4.png", width=800, height=400, type="cairo", bg="white"); dev.off()


