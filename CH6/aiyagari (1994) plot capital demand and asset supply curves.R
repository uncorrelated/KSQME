#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第6回「異質な個人を組み込んだマクロモデル」のソースコードのRへの野良移植: 需要供給曲線の計算とプロット, 図4
#
#==========================================================================
# Model of Aiyagari (1994)
# Written for Keizai Seminar #6
# By Sagiri Kitao
# 
# ported to R by uncorrelated
#==========================================================================
#

# 共通コードを読み込む
source("common.R")

#
# σに応じて遷移行列などを計算
#
tm_low <- calcTM(sig - 0.2)
tm_baseline <- calcTM(sig)
tm_high <- calcTM(sig + 0.2)

#==========================================================================
# COMPUTE INDIVIDUAL POLICY FUNCTION AND E(a)
#==========================================================================

NR = 20;
minR = -0.03;
maxR = (1-beta)/beta-0.001;
R = seq(minR, maxR, length.out=NR)

# 価値関数のタイプ（高精度/低精度）を定める
# (1) VALUE FUNCTION (USE THE SAME GRID FOR STATE AND CONTROL)
vfi_type <- "aiyagari_vfi1"

# (2) VALUE FUNCTION (USE FINER GRID FOR CONTROL)
# vfi_type <- "aiyagari_vfi2" # 時間がかかるが精度がアップ

# 演算ルーチンを指定する
# 1:R, 2:C, 4:R/MP, 8:C/MP
# 足し算で連続実行でき、ベンチマークを行える。例えば3を指定するとRとCのルーチンでそれぞれ計算（結果は同じ）
r_c_mp_switch <- 1

if(chkflag(r_c_mp_switch,1)){
	print("R Code:")
	print(system.time({
		source(paste(vfi_type, "R", sep=".")) # 価値関数を読み込む
		r_baseline <- sapply(R, function(r){
			with(tm_baseline, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
		r_high <- sapply(R, function(r){
			with(tm_high, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
		r_low <- sapply(R, function(r){
			with(tm_low, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
	}))
}

if(chkflag(r_c_mp_switch,2)){
	print("C Code:")
	print(system.time({
		vfi <- loadCModule()
		r_baseline <- sapply(R, function(r){
			with(tm_baseline, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
		r_high <- sapply(R, function(r){
			with(tm_high, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
		r_low <- sapply(R, function(r){
			with(tm_low, {
				vfi(r, alpha, beta, delta, mu, b, s, prob)$a
			})
		})
	}))
}

if(chkflag(r_c_mp_switch,3)){
	print("R & doParallel:")
	print(system.time({
		source(paste(vfi_type, "R", sep="."))

		# マルチスレッド実行に必要なパッケージのインストールを確認
		chkPkg(c("foreach", "iterators", "parallel", "doParallel")) 

		# 並列処理のためのパッケージをロード
		library(doParallel)
		library(foreach)

		# 利用可能コア数を調べて、並列処理のためのクラスターを構成
		cl <- makeCluster(detectCores())
		registerDoParallel(cl)

		r_baseline <- foreach(i = 1:length(R), .combine='c') %dopar% {
			with(tm_baseline, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}
		r_high <- foreach(i = 1:length(R), .combine='c') %dopar% {
			with(tm_high, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}
		r_low <- foreach(i = 1:length(R), .combine='c') %dopar% {
			with(tm_low, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}

		# 構成クラスターを破棄
		stopCluster(cl)
	}))
}

if(chkflag(r_c_mp_switch,4)){
	print("C & doParallel:")
	print(system.time({
		# マルチスレッド実行に必要なパッケージのインストールを確認
		chkPkg(c("foreach", "iterators", "parallel", "doParallel")) 

		# 並列処理のためのパッケージをロード
		library(doParallel)
		library(foreach)

		# 利用可能コア数を調べて、並列処理のためのクラスターを構成
		cl <- makeCluster(detectCores())
		registerDoParallel(cl)

		r_baseline <- foreach(i = 1:length(R), .combine='c') %dopar% {
			vfi <- loadCModule()
			with(tm_baseline, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}
		r_high <- foreach(i = 1:length(R), .combine='c') %dopar% {
			vfi <- loadCModule()
			with(tm_high, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}
		r_low <- foreach(i = 1:length(R), .combine='c') %dopar% {
			vfi <- loadCModule()
			with(tm_low, {
				# VALUE FUNCTION ITERATION
				vfi(R[i], alpha, beta, delta, mu, b, s, prob)$a
			})
		}

		# 構成クラスターを破棄
		stopCluster(cl)
	}))
}

# 推定結果をまとめておく
r_aiyagari <- data.frame(R = R, high = r_high, baseline = r_baseline, low = r_low)

# タブ区切りのファイルに保存（グラフを描き直すときに読み込めば、上の時間のかかる計算を省略できる） 
write.table(r_aiyagari, "r_aiyagari.txt", row.names=FALSE, quote=FALSE, sep="\t")
# r_aiyagari <- read.table("r_aiyagari.txt", header=TRUE, sep="\t")

#==========================================================================
# COMPUTE K
#==========================================================================

R_K = seq(0, 0.05, 0.005)
K = with(tm_baseline, labor*(alpha/(R_K+delta))^(1/(1-alpha)))

# プロットする
with(r_aiyagari, {
	par(bg="white")

	col <- c("red", "red", "red", "black")
	lty <- c(3, 2, 4, 1)
	lwd <- c(2, 2, 2, 2)

	plot(baseline, R, type="l", main="Capital Market Equilibrium", xlab=expression(paste(E, group("[", a, "]")," and ", "K")), ylab="Interest rate", xlim=c(0, 10), lty=lty[1], col=col[1], lwd=lwd[1])
	lines(high, R, lty=lty[2], col=col[2], lwd=lwd[2])
	lines(low, R, lty=lty[3], col=col[3], lwd=lwd[3])
	lines(K, R_K, lty=lty[4], col=col[4], lwd=lwd[4])

	sge <- c(3, 1, 2, 4) # ラベルの表示順番を入れ替えるためcol lty lwdを並べ替えるための配列
	legend("bottomright", lty=lty[sge], lwd=lwd[sge], col=col[sge], legend=c(expression(sigma==0.2), expression(sigma==0.4), expression(sigma==0.6), expression(K)), bty="n", y.intersp=1.5, seg.len=3, inset=0.05)
	abline(h=0, lty=3, col="gray")
})

# 図を保存する
dev.copy2eps(file="fig 4.eps", width=6, height=6)
# dev.copy(png, "fig 4.png", width=480, height=480, type="cairo", bg="white"); dev.off()

