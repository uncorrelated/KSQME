#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第5回「ニューケインジアン・モデルの新展開」のソースコードのRへの野良移植: 準線形線形モデル
#

# ディープパラメーター
dparam <- list(
	beta = 0.998951, # 主観的割引因子
	gamma = 0.005187, # 全要素生産性A_tの上昇率; γ = bar{γ} + z, z_t = ρ_z * z_{t-1} + ϵ_z
	nu = 0.166667, # 需要の価格弾力性の逆数
	tau = 2.83, # リスク回避度に関するパラメーター
	phi = 17.845363, # 調整費用の大きさに関するパラメーター
	kappa = 0.780000, # （本文中に言及なしで、線形モデルに利用）
	psi_1 = 1.800000, # インフレ率に対するR_{n,t}の反応係数
	psi_2 = 0.630000, # 産出ギャップに対するR_{n,t}の反応係数
	rho_R = 0.770000, # 金融政策ショックの慣性を表すパラメーター
	rho_g   = 0.98, # 財政政策ショックの慣性を表すパラメーター
	rho_z   = 0.88, # 全要素生産性ショックの慣性を表すパラメーター
	sigma_R = 0.22, # ϵ_Rの標準偏差; なぜか非線形モデルの100倍
	sigma_g = 0.71, # ϵ_gの標準偏差; なぜか非線形モデルの100倍
	sigma_z = 0.31, # ϵ_zの標準偏差; なぜか非線形モデルの100倍
	g_bar = log(1.250000), # 自然産出量Y^*の計算に使うパラメーターg 〜 AR(1)の定常状態の値
	G_over_Y = 0.200000, # 政府支出の産出量比
	pi_target = 1+3.30/400 # 目標インフレ率; Juliaのソースコード中で、piA=3.30; 1+piA/400となっている
)

# 定常状態（steady-state）の所得、消費、インフレ率、名目利子率を求める
# 式(20)〜(23)をR=R_n/Π==R_{n-1}/Π, Π^*=Π=…, Y^*=Y_t=…, G = (1 - g^-1)Yと置いて整理
# R_nが名目金利を粗利で表したもので、_nなしRは実質値なことに注意
# Hirose and Sunakawa (2019) p.98にマトメあり/expした値なのかどうか混乱しないように注意

attach(dparam)
ss <- list(
	c  = (1.0-nu)^(1.0/tau), # 式(21)(23)から
	pi = pi_target # 目標インフレ率
	) 
ss$y <- exp(g_bar)*ss$c # 式(22)から
ss$R <- beta^(-1)*exp(gamma)*ss$pi # 式(20)から
ss$R <- 1.450000 # なぜか準線形モデルでは…100*(ss$R - 1)?
detach(dparam)

# set bounds
bound <- with(dparam,{
	# 境界の大きさの程度を示す変数
	m.pctrn = 0.05*100 # なぜか非線形モデルの100倍
	m.mg = 2.0
	m.mz = 2.0
	m.mr = 2.0
	list(
		Rmin = -m.pctrn, # deviation from rstar
		Rmax = m.pctrn,
		gmin = -m.mg*sigma_g/sqrt(1-rho_g^2),
		gmax = m.mg*sigma_g/sqrt(1-rho_g^2),
		zmin = -m.mz*sigma_z/sqrt(1-rho_z^2),
		zmax = m.mz*sigma_z/sqrt(1-rho_z^2),
		rmin = -m.mr*sigma_R,
		rmax = m.mr*sigma_R
	)
})

slopecon <- with(bound, {
	m <- matrix(c(2/(Rmax-Rmin), -(Rmax+Rmin)/(Rmax-Rmin), 2/(gmax-gmin), -(gmax+gmin)/(gmax-gmin), 2/(zmax-zmin), -(zmax+zmin)/(zmax-zmin), 2/(rmax-rmin), -(rmax+rmin)/(rmax-rmin)), 4, 2, byrow=TRUE)
	rownames(m) <- c("R", "g", "z", "r")
	colnames(m) <- c("2/(max-min)", "-(max+min)/(max-min)")
	m
})

# randall-romero/CompEconRパッケージの一部
source("randall-romero/csize.R")
source("randall-romero/ckron.R")
source("randall-romero/qnwnorm.R")

# 謎のグリッド拡張
poly2s <- function(x){
	if(is.matrix(x)){
		nc <- ncol(x)
		nr <- nrow(x)
	} else {
		nc <- length(x)
		nr <- 1
	}
	r <- matrix(NA, nr, 1 + 2*nc)
	r[, 1] <- 1.0
	r[, 2*(1:nc)] <- x
	r[, 2*(1:nc) + 1] <- 2*x^2 - 1
	r
}

calcC <- function(fc0, fp0, rnpast, znow, rnow, ZLB=FALSE){
	with(dparam, {
		if(ZLB){
			c0 <- fc0 - (1/tau)*(-ss$R - fp0 - rho_z*znow)
		} else {
			invtau <- 1/tau
			c0 <- fc0 + invtau*(1-(1-rho_R)*psi_1*beta)*fp0 - invtau*rho_R*rnpast + invtau*rho_z*znow - invtau*rnow
			c0 <- c0/(1+invtau*(1-rho_R)*(kappa*psi_1+psi_2))
		}
		c0
	})
}

calcPi <- function(c, fp0){
	with(dparam, {
		kappa*c + beta*fp0
		})
}

calcY <- function(c, gnow){
	c + gnow
}

calcR <- function(rnpast, c, pi, rnow){
	with(dparam, {
            rho_R*rnpast + (1-rho_R)*(psi_1*pi + psi_2*c) + rnow
	})
}

r_solve <- with(dparam, {
	# initial values
	# NOTE: We use an index-function approach with a pair of policy functions.
	# One assumes the ZLB always binds and the other assumes the ZLB never
	# binds. The next period's policy function is given by a weighted average
	# of the policy function in the ZLB regime and the policy function in the
	# non-ZLB regime with an indicator function. The value of the indicator
	# function is one when the notional rate is greater than the ZLB,
	# otherwise zero.

	initPolicyF <- function(c, pi, R, y, c_ZLB, pi_ZLB, R_ZLB, y_ZLB, np=2, ns=4){
		nv <- 1 + np*ns # グリッドの長さ
		m <- matrix(rep(c(c, pi, R, y, c_ZLB, pi_ZLB, R_ZLB, y_ZLB), each=nv), nv, 2*ns)
		colnames(m) <- c("c", "pi", "R", "y", "c_ZLB", "pi_ZLB", "R_ZLB", "y_ZLB")
		m
	}

	# 謎のグリッド生成コード
	makegrid <- function(np=2, ns=4){
		nv <- 1 + np*ns # グリッドの長さ
		# ある1変数が-1/1のときは、他の3変数は必ず0の模様
		m <- matrix(c(0, rep(c(-1, 1, rep(0, nv)), ns - 1), -1, 1), nv, ns)
		colnames(m) <- c("R", "g", "z", "r")
		m
	}

	xgrid <- makegrid()
	bbt <- poly2s(xgrid) # 何の略語なのかなど不明
	bbtinv <- solve(bbt)

	# gh nodes and weights
	m.ngh <- 3 # 変数ごとのグリッドポイントの数
	m.nexog <- 3 # 変数の数; g,z,rの3変数の組み合わせで、3以外は動かない
	r_qnwnorm <- qnwnorm(m.ngh, 0.0, 1.0)

	# 組み合わせを示す行列をつくる
	# args: nc:変数の種類の数, nk:それぞれの変数が取る値の種類の数
	cmatrix <- function(nc, nk){
		l <- nk^nc
		m <- matrix(NA, l, nc)
		v <- 1:nk
		for(j in 1:nc){
			m[, j] <- v
			v <- rep(v, each=nk)
		}
		m
	}

	m <- cmatrix(m.nexog, m.ngh)

	ghweights <- with(r_qnwnorm, {
		w <- 1
		for(j in 1:m.nexog){
			w <- w * weights[m[, j]]
		}
		w
	})

	ghnodes <- with(r_qnwnorm, {
		m <- matrix(r_qnwnorm$xpoints[m], m.ngh^m.nexog, m.nexog)
		colnames(m) <- c("g", "z", "r")
		m
	})

	rm(m)

	policy_0 <- with(ss, { initPolicyF(0, 0, 0, 0, 0, 0, 0, 0) }) # 非線形モデルと初期値が異なる
	policy_1 <- with(ss, { initPolicyF(0, 0, 0, 0, 0, 0, 0, 0) })

	# Juliaのコードのfcvec0n, fpvec0n, fcvec0b, fpvec0bとfcvec1n, fpvec1n, fcvec1b, fpvec1bを2つの行列にまとめる
	fvec0 <- {
		fcss <- 0 # 非線形モデルと初期値が異なる
		fpss <- 0 # Juliaのソースコードはfpss = bet*phi*css^(-tau)*yss*(piss-piss)*pissで必ずゼロ
		nv <- ncol(bbt)
		ns <- 4 # 求める変数の数？
		m <- matrix(rep(c(fcss, fpss, fcss, fpss), each=nv), nv, ns)
		# ZLBにn:Not bounded, b:Bounded
		colnames(m) <- c("c0n", "pi0n", "c0b", "pi0b")
		m
	}

	fvec1 <- {
		m <- fvec0
		colnames(m) <- c("c1n", "pi1n", "c1b", "pi1b")
		m
	}


	# variables rnot0, c0, pi0, y0, rn0, fc0, fp0, gnow, znow, rnow
	Rpast <- with(bound, (Rmax-Rmin)/2*xgrid[,1] + (Rmax+Rmin)/2)
	gnow <- with(bound, (gmax-gmin)/2*xgrid[,2] + (gmax+gmin)/2)
	znow <- with(bound, (zmax-zmin)/2*xgrid[,3] + (zmax+zmin)/2)
	rnow <- with(bound, (rmax-rmin)/2*xgrid[,4] + (rmax+rmin)/2)

	tol <- 1e-05 # 最小更新幅
	maxit <- 1000 # 最大ループ回数
	msteps <- numeric(maxit) # 更新幅を保存するベクトル

	for(i in 1:maxit){
		# current period's c and pi (obtained by current period's fc and fp)
		# Juliaのコードのfc0/1とfp0/1はfvec0/1にまとめている

		# in the non-ZLB regime
		c0n <- calcC(fvec0[, "c0n"], fvec0[, "pi0n"], Rpast, znow, rnow)

		pi0n <- calcPi(c0n, fvec0[, "pi0n"])
		y0n = calcY(c0n, gnow)
		R0n = calcR(Rpast, c0n, pi0n, rnow)

		# in the ZLB regime
		c0b <- calcC(fvec0[, "c0b"], fvec0[, "pi0b"], Rpast, znow, rnow)
		pi0b <- calcPi(c0b, fvec0[, "pi0b"])
		y0b = calcY(c0b, gnow)
		R0b = calcR(Rpast, c0b, pi0b, rnow)

		# fitting polynomials
		# fvec0 = bbt %*% coefになる模様
		coef <- bbtinv %*% fvec0

		# g,z,Rへのショックをまとめて計算しておく
		shocks <- t(ghnodes) * c(sigma_g, sigma_z, sigma_R)

		# 政策関数のグリッドはm.nv=nrow(fvec1)個
		for(iv in 1:nrow(fvec1)){

			# 確率変数g,z,Rの下での期待値を計算
			#  3変数をまとめて扱うグリッドポイントごとに、その区間の累積密度を乗じた値を返す
			#  合計すれば、期待値

			# g,z,Rへのショックの組み合わせはm.nz=ncol(shocks)パターン
			r_sapply <- sapply(1:ncol(shocks), function(iz){

				gp  <- rho_g*gnow + shocks["g", iz]
				zp  <- rho_z*znow + shocks["z", iz]
				rp  <- rep(shocks["r", iz], length(rnow))

				# non-ZLBとZLBのfc1とfp1を同時に計算
				X <- matrix(c(R0n, gp, zp, rp), 4, byrow=TRUE) * slopecon[, 1] + slopecon[, 2]
				rownames(X) <- c("R", "gp", "zp", "rp")
				p2s <- poly2s(t(X))
				P <- p2s %*% coef

				# first assume the ZLB is not binding, and use coeffcn and coeffpn
				fc1n <- P[iv, "c0n"] # coeffcn
				fp1n <- P[iv, "pi0n"] # coeffpn

				# next period's c and pi (obtained by next period's fc and fp)
				c1n <- calcC(fc1n, fp1n, R0n[iv], zp[iv], rp[iv])
				pi1n <- calcPi(c1n, fp1n)
				y1n = calcY(c1n, gp[iv])
				R1n = calcR(R0n[iv], c1n, pi1n, rp[iv])

				if (is.na(R1n) || R1n < -1*ss$R){
					c1n <- calcC(fc1n, fp1n, R0n[iv], zp[iv], rp[iv], TRUE)
					pi1n <- calcPi(c1n, fp1n)
					y1n = calcY(c1n, gp[iv])
					R1n = -1*ss$R
				}

				fcxn <- c1n
				fpxn <- pi1n

				# in the ZLB regime
				fc1b <- P[iv, "c0b"] # coeffcb
				fp1b <- P[iv, "pi0b"] # coeffpb

				# next period's c and pi (obtained by next period's fc and fp)
				c1b <- calcC(fc1b, fp1b, R0b[iv], zp[iv], rp[iv])
				pi1b <- calcPi(c1b, fp1b)
				y1b = calcY(c1b, gp[iv])
				R1b <- calcR(R0n[iv], c1n, pi1n, rp[iv])

				if (is.na(R1b) || R1b < -1*ss$R){
					c1b <- calcC(fc1b, fp1b, R0n[iv], zp[iv], rp[iv], TRUE)
					pi1b <- calcPi(c1b, fp1b)
					y1b = calcY(c1b, gp[iv])
					R1b = -1*ss$R
				}

				fcxb <- c1b
				fpxb <- pi1b

				r <- c(
					ghweights[iz]*fcxn,
					ghweights[iz]*fpxn,
					ghweights[iz]*fcxb,
					ghweights[iz]*fpxb
				)

				names(r) <- c(
					"fc0n", "fp0n", "fc0b", "fp0b"
				)

				r
			})

			 # ラベルはJuliaのコードとの対応; prefixのfはfuture?
			fvec1[iv, "c1n"] <- sum(r_sapply["fc0n", ])
			fvec1[iv, "pi1n"] <- sum(r_sapply["fp0n", ])

			fvec1[iv, "c1b"] <- sum(r_sapply["fc0b", ])
			fvec1[iv, "pi1b"] <- sum(r_sapply["fp0b", ])
		}



		policy_1[, "c"] <- c0n
		policy_1[, "pi"] <- pi0n
		policy_1[, "R"] <- R0n
		policy_1[, "y"] <- y0n

		policy_1[, "c_ZLB"] <- c0b
		policy_1[, "pi_ZLB"] <- pi0b
		policy_1[, "R_ZLB"] <- R0b
		policy_1[, "y_ZLB"] <- y0b

	        # calculate the norm between the old and new policy functions
	        msteps[i] <- max(abs(policy_1 - policy_0))

		print(sprintf("%d %f %f", i, max(abs(policy_1[,1:4] - policy_0[,1:4])), max(abs(policy_1[,5:8] - policy_0[,5:8]))))

		# update the policy functions
		# なぜか新旧の政策関数のウェイト付き平均をとって、政策関数を更新
		m.damp <- 0.7
		policy_0 <- m.damp*policy_0 + (1 - m.damp)*policy_1
		fvec0 <- m.damp*fvec0 + (1 - m.damp)*fvec1

		# 新旧の政策関数の差がtol未満であれば終了
		if(tol > msteps[i]){
			break;
		}
	}

	list(policy=policy_0, coef=coef, xgrid=xgrid)
})

