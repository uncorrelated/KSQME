#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第8回「定量的マクロ経済学のフロンティア」のソースコードのRへの野良移植: 第3節HANKモデルとフォワード・ガイダンス・パズル, 図2
#

library(MNS2016) # 速度にうるさいJuliaファンが目立つので、Rcppでパッケージ化したC++/Eigen/Spectraのコードを呼び出すことにした。

# parameters
param <- list(
	gamma = 2.0,
	psi = 2.0,
	bl = 0.0,
	mu = 1.2,

	tau_bar = c(0, 0, 2.031292672482304),
	b4yss = 5.5/4.0,

	max_a = 75.0,
	min_a = 0.0,
	max_a_ns = 75.0,
	min_a_ns = 0.0,
	n_a = 200,
	n_a_ns = 1000,
	grid_a_shift = 0.354983310608304, # why this number?
	bl = 0,

	rho_z = 0.96566,
	var_z_inn = 0.01695^(1/2), # variance^(1/2)
	n_z = 3,

	max_beta_iter = 50,
	max_egm_iter = 1000,
	max_egm_const_iter = 1000,
	beta_err_tol = 1e-4, #1e-12,
	egm_err_tol = 1e-8, #1e-15,
	egm_const_err_tol = 1e-12, 
	ns_err_tol = .Machine$double.eps
)

param$sub_mon_com = with(param, (mu-1.0)/mu)

W_SS = with(param, 1.0/(mu*(1.0 - sub_mon_com))) # = 1
R_SS = 1.005;
PI_SS = 1.0;

# Discretize the endogenous state space for assets
grid_a_preshift = with(param, seq(log(min_a_ns + grid_a_shift), log(max_a_ns + grid_a_shift), length.out=n_a))
grid_a = exp(grid_a_preshift) - param$grid_a_shift;

grid_a_preshift_NS = with(param, seq(log(min_a + grid_a_shift),log(max_a + grid_a_shift), length.out=n_a_ns))
grid_a_NS = exp(grid_a_preshift_NS) - param$grid_a_shift;

# Discretize the exogenous stochastic process
r_rouwenhorst = with(param, Rouwenhorst(rho_z, var_z_inn, n_z)) # Rのコードだが、MNS2016パッケージに定義を置いてある
grid_z = exp(r_rouwenhorst$state_values);
egn = eigen(t(r_rouwenhorst$p));
p <- which.min(egn$values) # 結果的にJuliaのコードが最小の固有値のベクターを選んでいるので、数字があうようにする
s_z = egn$vectors[, p]/sum(egn$vectors[, p]) # この固有ベクトルでいい理由は謎

# Set interest rates and wages to the steady-state levels.
RToday = R_SS;
wToday = W_SS;

# bisection
bisec_min = 0.75;
bisec_max = 0.995;

meanY  = 1.0; # The initial guess of Yss.
beta = (bisec_min + bisec_max)/2;

# Initial guess for the policy functions for consumption and labor supply

# repのeach/timesに慣れてくると、以下のようにループなしで処理できる事が多い
pf_c = with(param, matrix(0.3 + 0.1*grid_a, n_a, n_z))
pf_n = with(param, matrix( (wToday*grid_z*t(pf_c)^(-gamma))^(1.0/psi), n_a, n_z, byrow=TRUE))

# ifで分岐するとベクトルで処理できないのでifelseを使うことになるが、真のときの式を偽のときに実行してもNAやInfにならないのであれば、以下のように論理式をかけて足してしまう方が高速なことが多い。論理式は真のときは1、偽のときは0として評価される。
values <- with(param, {
	grid_a_z <- with(param, rep(grid_a, each=3) + grid_z)
	matrix( 0 + (gamma == 1.0)*log(grid_a_z)/(1.0-beta) + (gamma != 1.0)*((grid_a_z)^(1.0-gamma)/(1.0-gamma))/(1.0-beta), n_a, n_z, byrow=TRUE)
})

pf_sav <- with(param, matrix(0, n_a, n_z))

print(system.time({
	for(beta_iter in 1:param$max_beta_iter){

		beta = (bisec_min + bisec_max)/2

		tauToday = with(param, 1.0/s_z[3]*b4yss*4.0*meanY*(RToday-1.0)/RToday/tau_bar[3]);
		dToday   = meanY*(1.0 - wToday);

		cat("HH_opt_EGM: ")
		etime <- system.time({
			r_HH_opt_EGM <- HH_opt_EGM(beta, param, grid_a, grid_z, r_rouwenhorst$p, wToday, RToday, tauToday, dToday, pf_c, pf_n, pf_sav);
		})
		cat(sprintf("%f seconds elapsed.\n", etime["elapsed"]))

		pf_c <- r_HH_opt_EGM$pf_c
		pf_n <- r_HH_opt_EGM$pf_n
		pf_sav <- r_HH_opt_EGM$pf_sav

		cat("HH_dist: ")
		etime <- system.time({
			r_HH_dist <- HH_dist(beta, param, grid_a, grid_a_NS, grid_z, r_rouwenhorst$p, wToday, RToday, tauToday, dToday, pf_c, pf_n, pf_sav);
		})
		cat(sprintf("%f seconds elapsed.\n", etime["elapsed"]))

		meanC = r_HH_dist$meanC
		meanN = r_HH_dist$meanN
		meanA = r_HH_dist$meanA

		cat(sprintf("[%d %f %f %f]\n", beta_iter, meanC, meanN, meanA))

		meanY = meanC;

		# Evaluate convergence
		if(abs(meanA/(4.0*meanY) - param$b4yss) < param$beta_err_tol) break; 

		# Update the guess using the bisection method
		if ( meanA/(4.0*meanY) > param$b4yss){
			bisec_max = beta;
		} else {
			bisec_min = beta;
		}
	}
}))

# 計算結果をプロット
with(r_HH_dist, {
	par(oma=c(0, 0, 0, 0), mfrow=c(3,1), mar=c(4.5, 4.5, 2, 1), bg="white") # mfrow=c(3,1)は、3つのプロットを横3分割して表示するためのオプション
	xlim <- c(0, 40)
	drawPlot <- function(title, y){
		plot(grid_a_NS, y, main=title, xlab="Assets", ylab="Density", xlim=xlim, type="l")
	}
	drawPlot(expression(paste("Productivity: ", S[L])), dist[,1])
	drawPlot(expression(paste("Productivity: ", S[M])), dist[,2])
	drawPlot(expression(paste("Productivity: ", S[H])), dist[,3])
})

# プロットをepsファイルに保存
dev.copy2eps(file="fig 2.eps", width=8, height=6)
# dev.copy(png, "fig 2.png", width=800, height=600, type="cairo", bg="white"); dev.off()

