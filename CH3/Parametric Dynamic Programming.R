#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第3回「動的計画法」のソースコードのRへの野良移植: 第5節無限期間モデルの解法, 図4
#

source("common.R")

param <- list(
	### カリブレーション ###
	beta  = 0.96, # 割引因子
	gamma = 1.0,  # 相対的危険回避度(異時点間の代替の弾力性の逆数)
	alpha = 0.40, # 資本分配率
	delta = 1.00, # 固定資本減耗

	### 離散化用のパラメータ ###
	nk   = 21,    # グリッドの数
	kmax = 0.5,   # 資本グリッドの最大値
	kmin = 0.05  # 資本グリッドの最小値 (0にすると生産が出来なくなる)
)

# 固定資本減耗を変えるときは以下2行をアンコメント
# param$delta <- 0.08
# param$kmax <- 10.0

### 収束の基準 ###
maxit = 1000;    # 繰り返し計算の最大値
tol  = 1.0e-005; # 許容誤差(STEP 2)
dif1 = 1.0;      # 価値関数の繰り返し誤差
dif2 = 1.0;      # 政策関数の繰り返し誤差

# パラメーターを変えて比較するかも知れないので、数値計算部分を関数化する
calcPDP <- function(param){
	with(param, {
		cat("\n\n-+- Solve a neoclassical growth model -+-\n\n");
		cat("\n-+- PARAMETER VALUES -+-\n");
		cat(sprintf("beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f\n", beta, gamma, alpha, delta))
		cat(sprintf("kmin=%5.2f, kmax=%5.2f, #grid=%d\n", kmin, kmax, nk))

		## STEP 1(a): グリッド生成
		kgrid = seq(kmin, kmax, length.out=nk);

		## STEP 1(b): 価値関数・政策関数の初期値を当て推量
		pfcn0 = matrix(0, nk, 1);
		vfcn0 = CRRA(kgrid^alpha + (1-delta)*kgrid, gamma);

		pfcn1 = matrix(0, nk, 1);
		vfcn1 = matrix(0, nk, 1);

		## 価値関数・政策関数の経路を記録(なくても可)
		# 行列サイズの自動拡張でパフォーマンスが落ちないようにあらかじめ、必要サイズを取得しておく
		vpath <- matrix(NA, nk, maxit)
		ppath <- matrix(NA, nk, maxit)
		vpath[, 1] = vfcn0;
		ppath[, 1] = pfcn0;

		# 収束途中の繰り返し計算誤差を保存
		# 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき) ← 速度をお求めの場合はCかFortranでモジュールを…
		dif <- matrix(NA, 2, maxit)

		for(it in 1:maxit){ # itがループのカウンターになる

			cat(sprintf('iteration index: %i \n', it));
			cat(sprintf('value function iteration error: %e\n', dif1));
			cat(sprintf('policy function iteration error: %e\n', dif2));

			# 次期の価値関数をスプライン関数で近似する
			vnext <- splinefun(kgrid, vfcn0)

			for(i in 1:nk){
				# 価値関数が最大値となる来季資本（kprime）を初期値初期値0.01から探す
				# nlmが最小値を求める関数なので、目的関数内でベルマン方程式に-1をかけ、計算された最小値に-1をかけて最大値を求める
				r_nlm <- nlm(function(p, k, vnext, alpha, beta, gamma){
						# 変数のように関数を定義できると言うか、関数も変数のR
						-1 * BellmanEq(p, k, vnext, alpha, beta, gamma)
					}, 0.01, k=kgrid[i], vnext=vnext, alpha=alpha, beta=beta, gamma=gamma)
				pfcn1[i] <- r_nlm$estimate
				vfcn1[i] <- r_nlm$minimum * -1
			}

			# 繰り返し計算誤差を計算
			dif1 = max(abs((vfcn1-vfcn0)/vfcn0));
			dif2 = max(abs((pfcn1-pfcn0)/pfcn0));

			# 繰り返し計算誤差を保存
			dif[1, it] = dif1;
			dif[2, it] = dif2;

			# 価値関数の経路を記録(なくても問題なし)
			vpath[, it] = vfcn0;
			ppath[, it] = pfcn0;

			# 価値関数・政策関数をアップデート
			vfcn0 = vfcn1;
			pfcn0 = pfcn1;

			# 繰り返し計算誤差がtol未満になったら、収束したとして終了
			if(dif1 < tol){
				cat(sprintf("The iteration finished successfully: dif1=%f < %f", dif1, tol))
				break
			}
		}

		# 繰り返し上限に達していたらエラー
		if(maxit <= it){
			stop("iteration limit exceeded.")
		}

		# 最終的な政策関数が得られてから消費関数を計算
		cfcn = kgrid^alpha + (1-delta)*kgrid - pfcn0[,1]

		## 代数的解の計算
		AA = (1.0-beta)^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
		BB = alpha/(1.0-alpha*beta);
		v_true = AA + BB*log(kgrid);
		p_true = beta*alpha*(kgrid^alpha);

		## オイラー方程式から誤差を測定
		cons = kgrid^alpha + (1-delta)*kgrid - pfcn0[, 1];
		LHS  = dCRRA(cons, gamma);
		kp   = pfcn0[,1];
		interp1 <- splinefun(kgrid, pfcn0[,1])
		kpp  = interp1(kp);
		cons = kp^alpha + (1-delta)*kp - kpp;
		rent = alpha*kp^(alpha-1.0) - delta;
		RHS  = beta*(1+rent)*dCRRA(cons, gamma);
		err  = RHS/LHS-1.0;

		# 計算結果はリストにまとめておく
		list(kgrid=kgrid, v_true=v_true, p_true=p_true, vfcn=vfcn0, pfcn=pfcn0, cfcn=cfcn, dif=dif[,1:it], err=err)
	})
}

# 実行時間の計測と表示
print(system.time({
	r_pdp <- calcPDP(param)
}))


# 他で比較するために、誤差を保存
write.table(data.frame(
	kgrid = r_pdp$kgrid, 
	error = r_pdp$err
), file="err_pdp.csv", col.names=TRUE, row.names=FALSE, sep=",")

if(!file.exists(file_err_ddp)){
	stop(sprintf("Couldn't find the file:'%s'. You need to run 'Discretized Dynamic Programming.R' beforehand."))
}
err_ddp <- read.table(file_err_ddp, sep=",", header=TRUE)

with(r_pdp, {
	dev.new(width=10, height=5)
	par(oma=c(0, 0, 3, 0), mfrow=c(1, 2), mar=c(4.5, 4.5, 1, 1), bg="white")

	lty <- c(1, 2)
	lwd <- c(2, 2)
	col <- c("black", "blue")
	legends <- c("error of the value function", "error of the policy function")
	ylim <- c(0, 1)
	dif[is.infinite(dif)] <- ylim[2] # Infはプロットされないので、最大値を入れておく
	plot(1:ncol(dif), dif[1, ], main="", xlab="The number of iteration", ylab="Iteration error", type="l", lty=lty[1], lwd=lwd[1], col=col[1], ylim=ylim)
	lines(1:ncol(dif), dif[2, ], lty=lty[2], lwd=lwd[2], col=col[2])
	legend("topright", lty=lty, lwd=lwd, col=col, legend=legends, bty="n", y.intersp=1.1, seg.len=2.5, inset=0.02)

	lty <- c(1, 2)
	lwd <- c(2, 2)
	col <- c("black", "red")
	legends <- c("Discretized", "Prameteric (Continuous)")
	plot(err_ddp$kgrid, err_ddp$error, main="", xlab=expression(paste("Current Priod Assets: ", k[t])), ylab="Error of the Euler equation", type="l", lty=lty[1], lwd=lwd[1], col=col[1])
	lines(kgrid, err, lty=lty[2], lwd=lwd[2], col=col[2])
	legend("bottomright", lty=lty, lwd=lwd, col=col, legend=legends, bty="n", y.intersp=1.1, seg.len=2.5, inset=0.02)

	mtext(side=3, line=1, outer=T, text="Fig.4 the errors of Prameteric Dynamic Programming (Continuous method)", cex=1.5)

	dev.copy2eps(file="fig 4.eps", width=8, height=4) # 図を保存する
	# dev.copy(png, "fig 3.png", width=800, height=400, type="cairo", bg="white"); dev.off()
})

