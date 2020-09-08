#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第8回「定量的マクロ経済学のフロンティア」のソースコードのRへの野良移植: 第2節 クルセル＝スミスモデル（Krusell and Smith Model）
#
#==========================================================================
# The original Fortran code was written by T. Yamada. at 09/14/2019
# ported to R by uncorrelated
#==========================================================================
#

source("common.R")

# ディレクトリの存在チェック
if(!dir.exists(rdir)){
	if(!dir.create(rdir)){
		stop(sprintf("Couldn't create the directory: %s", rdir))
	}
}


attach(dparam)
attach(gparam)
attach(sparam)

###### regression coefficients ######
lreg <- list(
	icept = numeric(nz),
	icept1 = numeric(nz),
	R2 = numeric(nz),
	slope = rep(1.0, nz),
	slope1 = numeric(nz))

# lreg$icept <- c(0.2, 0.1)
# lreg$slope <-c(1.0, 0.9)

lreg$icept <- c(0.0, 0.0)
lreg$slope <-c(1.0, 1.0)

# transition probability of (z,z'): business cycle
# call gen_prob_matrix(durgd, 0.5d0, nz, pi_z, tran_zz)

###### transition probability ######
trans_prob <- with(new.env(), {
	# Transition probability of idiosyncratic and aggregate shocks.

	unemp <- 0.5
	duration <- durgd

	pi_dist = numeric(nz) # グローバル変数のpi_distとは関係が無い
	pi_dist[1] = 1.0 - unemp
	pi_dist[2] = unemp

	tran_zz = array(0, c(nz, nz)) # pr(z,z')
	tran_zz[2,2] = (duration-1.0) / duration
	tran_zz[2,1] = 1.0 - tran_zz[2, 2]
	tran_zz[1,1] = 1.0 - (1.0 - tran_zz[2, 2])*(pi_dist[2]/pi_dist[1])
	tran_zz[1,2] = 1.0 - tran_zz[1, 1]

        # copy-n-paste from A.Smith's code をだいたいそのまま使う
        pgg00 = (durug-1.0)/durug # tran_ee(2,2)|good
        pbb00 = (durub-1.0)/durub # tran_ee(2,2)|bad
        pbg00 = 1.25*pbb00
        pgb00 = 0.75*pgg00
        pgg01 = (unempg - unempg*pgg00)/(1.0-unempg)
        pbb01 = (unempb - unempb*pbb00)/(1.0-unempb)
        pbg01 = (unempb - unempg*pbg00)/(1.0-unempg)
        pgb01 = (unempg - unempb*pgb00)/(1.0-unempb)
        pgg = (durgd-1.0)/durgd       # tran_zz(1,1)
        pgb = 1.0 - (durbd-1.0)/durbd # tran_zz(1,2)

        pgg10 = 1.0 - (durug-1.0)/durug # tran_ee(1,2)|good
        pbb10 = 1.0 - (durub-1.0)/durub # tran_ee(1,2)|bad
        pbg10 = 1.0 - 1.25*pbb00
        pgb10 = 1.0 - 0.75*pgg00
        pgg11 = 1.0 - (unempg - unempg*pgg00)/(1.0-unempg)
        pbb11 = 1.0 - (unempb - unempb*pbb00)/(1.0-unempb)
        pbg11 = 1.0 - (unempb - unempg*pbg00)/(1.0-unempg)
        pgb11 = 1.0 - (unempg - unempb*pgb00)/(1.0-unempb)
        pbg = 1.0 - (durgd-1.0)/durgd # tran_zz(2,1)
        pbb = (durbd-1.0)/durbd       # tran_zz(2,2)

	prob = array(0, c(ne, nz, ne, nz)) # pr(e,z,e',z')

        prob[1, 1, 1, 1] = pgg*pgg11
        prob[1, 1, 1, 2] = pbg*pbg11
        prob[1, 1, 2, 1] = pgg*pgg01
        prob[1, 1, 2, 2] = pbg*pbg01
        prob[1, 2, 1, 1] = pgb*pgb11
        prob[1, 2, 1, 2] = pbb*pbb11
        prob[1, 2, 2, 1] = pgb*pgb01
        prob[1, 2, 2, 2] = pbb*pbb01
        prob[2, 1, 1, 1] = pgg*pgg10
        prob[2, 1, 1, 2] = pbg*pbg10
        prob[2, 1, 2, 1] = pgg*pgg00
        prob[2, 1, 2, 2] = pbg*pbg00
        prob[2, 2, 1, 1] = pgb*pgb10
        prob[2, 2, 1, 2] = pbb*pbb10
        prob[2, 2, 2, 1] = pgb*pgb00
        prob[2, 2, 2, 2] = pbb*pbb00

	#   prob:     transition probability of (e,z)
	#   tranz_zz: transition probability of z
	list(ez=prob, zz=tran_zz)
})

# 第1列が好景気のとき、第2列が不景気のとき…遷移行列であれば縦横が逆
pi_dist <- matrix(c(1.0-unempg, unempg, 1.0-unempb, unempb), ne, nz)

###### macro variables ######
L_agg = with(new.env(), {
	L_agg = numeric(nz)
	for(i in 1:nz){
		L_agg[i] = sum(pi_dist[, i] * endow) # なぜか内積をとる
	}
	L_agg
})

detach(sparam)
detach(gparam)
detach(dparam)

writeLines("-+-+-+- Solving Krusell and Smith model -+-+-+-", con=stderr())
writeLines("--- calibrated parameters ---", con=stderr())
with(dparam, {
	writeLines(sprintf("discount factor: %10.4f", beta), con=stderr())
	writeLines(sprintf("relative risk aversion: %10.4f", gamma), con=stderr())
	writeLines(sprintf("capital share: %10.4f", alpha), con=stderr())
	writeLines(sprintf("depreciation rate: %10.4f", delta), con=stderr())
	writeLines(sprintf("TFP level (g): %10.4f", tfp[1]), con=stderr())
	writeLines(sprintf("TFP level (b): %10.4f", tfp[2]), con=stderr())
	writeLines(sprintf("unemployment rate (g): %10.4f", unempg), con=stderr())
	writeLines(sprintf("unemployment rate (b): %10.4f", unempb), con=stderr())
	writeLines(sprintf("duration of unemp period (g): %10.4f", durug), con=stderr())
	writeLines(sprintf("duration of unemp period (b): %10.4f", durub), con=stderr())
	writeLines("", con=stderr())
})

dll <- dyn.load(paste("KnS1998", .Platform$dynlib.ext, sep = ""))

markov <- function(tran, ng, ns){
	.Call("Cmarkov", tran, as.integer(ng), as.integer(ns))
}

end_grid_method <- function(lreg){
	with(dparam,{
		with(gparam, {
			# Fortranのfactor_price()は2つに分割
			wage_f <- function(agg_cap, agg_lab, tfp, alpha, delta){
				tfp*(1.0-alpha)*agg_cap^alpha * agg_lab^(-alpha)
			}

			rent_f <- function(agg_cap, agg_lab, tfp, alpha, delta){
				tfp* alpha *agg_cap^(alpha-1.0)*agg_lab^(1.0-alpha) - delta
			}

			approx_agg <- function(agg_cap, state, icept, slope){
				exp(icept[state] + slope[state]*log(agg_cap))
			}

			RHS_Euler <- function(asset, e_state, K_agg, z_state, consf, xgrid){
				.Call("RHS_Euler", asset, as.integer(e_state), K_agg, L_agg, as.integer(z_state), consf, xgrid, kgrid, trans_prob$ez, as.integer(na), as.integer(ne), as.integer(nk), as.integer(nz), alpha, beta, gamma, delta, endow, tfp)
			}

			#### 4重ループを展開した数列をつくる ####
			i4 <- with(new.env(), {
				m <- na*ne*nk*nz
				i <- 1:m
				n <- i - 1
				
				n <- n - 0*nz
				m <- m/nz
				iz <- as.integer(n / m) # current TFP level
				
				n <- n - iz*m
				m <- m/nk
				ik <- as.integer(n / m) # current aggregate capital
				
				n <- n - ik*m
				m <- m/ne
				ie <- as.integer(n / m) # current employment status

				n <- n - ie*m
				m <- m/na
				ina <- as.integer(n / m) # current saving

				list(z=iz+1, k=ik+1, e=ie+1, na=ina+1)				
			})

			#### 2重ループを展開した数列をつくる ####
			i2 <- with(new.env(), {
				m <- nk*nz
				i <- 1:m
				n <- i - 1
				
				n <- n - 0*nz
				m <- m/nz
				iz <- as.integer(n / m)
				
				n <- n - iz*m
				m <- m/nk
				ik <- as.integer(n / m)

				list(z=iz+1, k=ik+1)			
			})

			# initial guess (hand-to-mouth)
			rent <- rent_f(kgrid[i4$k], L_agg[i4$z], tfp[i4$z], alpha, delta)
			wage <- wage_f(kgrid[i4$k], L_agg[i4$z], tfp[i4$z], alpha, delta)
			con0 <- coh0 <- array(wage*endow[i4$e] + (1.0+rent)*aprime, c(na, ne, nk, nz))

			# given current capital and approximate aggregation, compute the next period's capital
			agg_cap_next = array(approx_agg(kgrid[i2$k], i2$z, lreg$icept, lreg$slope), c(nk, nz))

			maxit <- 1000
			toler <- 1.0e-4

			debug_c <- 0

			for(it in 1:maxit){
				# factor prices in the next period
				K_agg <- c(agg_cap_next)[i4$k + (i4$z - 1)*nk] # agg_cap_next[i4$k, i4$z] # 配列からのデータ取り出し方法を要確認

				# main part of EGM
				rhs = RHS_Euler(aprime[i4$na], i4$e, K_agg, i4$z, con0, coh0)

				con1 = array(rhs^(-1.0/gamma), c(na, ne, nk, nz))
				coh1 = con1 + aprime

				# check converegence
				metric = max(abs(con0 - con1)/con0)

				# update policy funciton
				con0 = con1
				coh0 = coh1

				if (metric < toler) {
					break
				}
			}

			###### policy function ######
			policy = array(0, c(nd, ne, nk, nz))

			for(z in 1:nz){
				for(k in 1:nk){
					rent <- rent_f(kgrid[k], L_agg[z], tfp[z], alpha, delta)
					wage <- wage_f(kgrid[k], L_agg[z], tfp[z], alpha, delta)
#					writeLines(sprintf("z:%d, k:%d, rent:%f, wage:%f", z, k, rent, wage), con=stderr())
					for(e in 1:ne){
						asset = (coh0[, e, k, z] - wage*endow[e]) / (1.0 + rent)
						sp <- splinefun(asset, aprime)
#						ap <- approxfun(asset, aprime)
						for(i in 1:nd){
							if (asset[1] > grid[i]){
								policy[i, e, k, z] = amin
#							} else if (asset[gparam$na] < grid[i]){
#								policy[i, e, k, z] = max(amin, ap(grid[i]))
							} else {
								policy[i, e, k, z] = max(amin, sp(grid[i]))
							}

							if(is.na(policy[i, e, k, z])){
								writeLines(sprintf("na: %d %d %d %d %f <> %f %f", i, e, k, z, asset[gparam$na], grid[i], amin), con=stderr())
								writeLines(sprintf("ap:%f sp:%f", ap(grid[i]), sp(grid[i])), con=stderr())
							}
						}
					}
				}
			}

			policy
		})
	})
}

set.seed(225)
law_of_motion_sim <- function(policy){
	# パラメーターによってはPathにNAが大量に入るので、本当はエラーチェックがいる
	.Call("law_of_motion_sim", policy, as.integer(sparam$numi), as.integer(sparam$nums), trans_prob$ez, trans_prob$zz, kgrid, grid, as.integer(dparam$ne))
}

regress <- function(k_path, z_path, lreg){
	# discard first 1000 periods
	z_path <- z_path[-(1:1000)]
	k_path <- k_path[-(1:1000)]
	
	z_path <- z_path[k_path > 0]
	k_path <- k_path[k_path > 0]

	# find regression parameters
	for(z in 1:dparam$nz){

		index_x <- z_path == z
		index_x[length(index_x)] = FALSE
		index_y <- c(FALSE, index_x[-length(index_x)])

		x = log(k_path[index_x])
		y = log(k_path[index_y])

		# 単回帰は簡単に書ける
		lreg$slope1[z] = cov(x, y)/var(x)
		lreg$icept1[z] = mean(y) - lreg$slope1[z]*mean(x)
		tss = length(y)*var(y)
		rss = sum((y - lreg$icept1[z] - lreg$slope1[z]*x)^2)
		lreg$R2[z] = 1 - rss/tss
	}

	lreg
}

###### variables for main loop ######
maxit <- 300
toler <- 1e-4
adj <- 0.5

dnorm <- numeric(maxit)
beta0 <- matrix(0, maxit, 2)
beta1 <- matrix(0, maxit, 2)

for(i in 1:maxit){
	writeLines(sprintf("iteration counter: %d", i), con=stderr())

	policy <- end_grid_method(lreg)
	path <- law_of_motion_sim(policy)
	lreg <- regress(path$k, path$z, lreg)

	attach(lreg)

	# metric of error
	metric1 = pmax(abs((icept-icept1)/icept))
	metric2 = pmax(abs((slope-slope1)/slope))
	metric  = max(metric1, metric2)

	beta0[i, ] <- icept1
	beta1[i, ] <- slope1

	j <- ifelse(10 < i, i - 9, 1)
	dnorm[i] <- sqrt(sum(var(beta0[j:i, ])^2 + var(beta1[j:i, ])^2))
 
	# update coefficients
	lreg$icept = adj*icept + (1.0-adj)*icept1
	lreg$slope = adj*slope + (1.0-adj)*slope1

	writeLines("", con=stderr())
	writeLines(sprintf("error:     %10.4f%%", metric*100.0), con=stderr())
	writeLines(sprintf("intercept: %10.4f", lreg$icept), con=stderr())
	writeLines(sprintf("slope:     %10.4f", lreg$slope), con=stderr())
	writeLines(sprintf("R^2:       %10.4f", R2), con=stderr())
	writeLines(sprintf("norm(var): %10.4f", dnorm[i]), con=stderr())
	writeLines("", con=stderr())

	detach(lreg)

	# （収束しないので）過去10回の回帰係数の分散がtolerance未満であれば、収束したと見なす ← 反則
	if(10 < i && toler > dnorm[i]){
		writeLines(sprintf("The iteration maybe finished successfully: norm(var(%d:%d))=%f", i - 9, i, dnorm[i]), con=stderr())
		break;
	}

	if (!is.na(metric) && metric < toler){
		writeLines(sprintf("The iteration finished successfully: metric=%f < %f", metric, toler), con=stderr())
		break;
	}
}

if(maxit <= i){
	stop("iteration limit exceeded.")
}

# output simulation results
setwd(rdir)
saveRDS(policy, file=objf)

con_z <- file("z_path.txt", "w")
con_k <- file("k_path.txt", "w")
con_u <- file("u_path.txt", "w")
writeLines(paste(path$z, sep="\n"), con = con_z)
writeLines(paste(path$k, sep="\n"), con = con_k)
writeLines(paste(path$u, sep="\n"), con = con_u)
close(con_z)
close(con_k)
close(con_u)

setwd("..")

