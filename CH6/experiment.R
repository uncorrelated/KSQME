#
# 実際には使われていない実験コード
#

labor0 = 1.125709056856582;  # LABOR IN THE BASELINE (USED ONLY IN EXPERIMENTS)

# ADJUST PRODUCTIVITY GRID S.T. labor REMAINS UNCHANGED (IN EXPERIMENTS)
if(FALSE){
	adj=labor0/labor;
	s=s*adj;
	labor <- c(s %*% invdist)
}

