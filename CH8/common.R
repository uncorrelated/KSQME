#
# 計算とプロットの共通部分
#
rdir <- "output"
objf <- "policy.obj"

grid_exp <- function(min, max, n, times=0){
	# 何回も対数や指数をとるよりは、対数の底を変える方がわかりやすい。値は同じにはならないが、
	# (exp(1)^times)^seq(log(min + 1, base=exp(1)^times), log(max + 1, base=exp(1)^times), length.out=n) - 1
	# ではいけないのか？
	l <- max
	for(i in 1:times){
		l <- log(l + 1.0)
	}
	s <- seq(min, l, length.out=n)
	for(i in 1:times){
		s <- exp(s) - 1.0
	}
	s
}

dparam <- list(
	####### preference ######
	beta  = 0.99,  # discount factor
	gamma = 1.00,  # relative risk aversion

	####### production ######
	alpha = 0.36,  # capital share
	delta = 0.025, # depreciation rate

	###### idiosyncratic unemployment shock ######
	ne = 2, # employed, unemployed
	endow = c(1.0, 0.05), # ui = 5%

	####### aggregate productivity shock and unemployment dynamics ######
	nz = 2, #aggregate state = 2
	unempg = 0.04, # unemployment rate when good
	unempb = 0.1,  # unemployment rate when bad
	durug  = 1.5,  # duration of unemployment when good
	durub  = 2.5,  # duration of unemployment when bad
	durgd  = 8.0,  # duration of boom (8 quarters)
	durbd  = 8.0,  # duration of recession (8 quarters)
	tfp = c(1.01, 0.99) # TFP shock
)

gparam <- list(
	###### grid parameters ######
	na = 21, # 101
	nd = 501, # 1001
	nk = 6,
	amax  = 300.0,
	amin  = 0.0,
	kmax = 40.0,
	kmin = 30.0
)

sparam <- list(
	###### simulation part ######
	nums = 6000, #6000 # 11000
	numi = 1000 #1000 # 5000
)

attach(dparam)
attach(gparam)

###### grid ######
# aprime = numeric(na) # 来季資産
# grid = numeric(nd)
# kgrid = numeric(nk)
kgrid = seq(kmin, kmax, length.out=nk)
aprime = grid_exp(amin, amax, na, 3)
grid =  grid_exp(amin, amax, nd, 3)

detach(gparam)
detach(dparam)


