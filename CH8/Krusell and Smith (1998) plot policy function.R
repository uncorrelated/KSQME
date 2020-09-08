#
# 経済セミナー連載「定量的マクロ経済学と数値計算」
# 第8回「定量的マクロ経済学のフロンティア」のソースコードのRへの野良移植: 第2節 クルセル＝スミスモデル（Krusell and Smith Model）, 図1
#
source("common.R")

policy <- with(list(path_to_obj = paste(rdir, objf, sep="/")),{
	if(!file.exists(path_to_obj)){
		stop(paste("Couldn't open the file: ", path_to_obj, sep=""))
	}
	readRDS(path_to_obj)
})

par(mar=c(4, 4.5, 3, 1), bg="white")

len <- 375
y <- c(policy[1:len, 1, 3, 1], policy[1:len, 2, 3, 1])
ylim <- c(min(y), max(y))
lwd <- c(2, 2, 2)
lty <- c(1, 2, 3)
col <- c("black", "blue", "gray")
plot(grid[1:len], policy[1:len, 1, 3, 1], type="l", lwd=lwd[1], lty=lty[1], col=col[1], ylim=ylim, xlab=expression( paste("Current Period Assets: ", a[t]) ), ylab=expression( paste("Next Period Assets: ", a[t + 1]) ), main="Policy Function near the Borrowing Constraint")
lines(grid[1:len], policy[1:len, 2, 3, 1], lwd=lwd[2], lty=lty[2], col=col[2])
lines(grid[1:len], grid[1:len], lwd=lwd[2], lty=lty[3], col=col[3])
legend("topleft", lty=lty, lwd=lwd, col=col, legend=c("employed", "unemployed", expression(paste("45-degree line: ", a[t]==a[t+1]))), bty="n", y.intersp=1.1, seg.len=2.5)

dev.copy2eps(file="fig 1.eps", width=4, height=4)
# dev.copy(png, "fig 1.png", width=400, height=400, type="cairo", bg="white"); dev.off()
