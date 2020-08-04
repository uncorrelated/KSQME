[経済セミナー連載の「定量的マクロ経済学と数値計算」のサポート・リポジトリー](https://github.com/keizai-seminar-quant-macro)にR用のコードが部分的にしか用意されていなかったので，ちょっと書いてみたものです。

実行方法はgit cloneなどをしてファイルを展開した後に、Rに読み込ませるだけです。例えば、第3節の離散化による2期間モデルの解法は、[two-period discretized method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20discretized%20method.R)を以下のように実行すると計算します。

	source("two-period discretized method.R")

`getwd()`で確認できるワーキングディレクトリに`two-period discretized method.R`が無いと呼び出せないので注意してください。また、他のファイルやパッケージに依存しているファイルもあります。ファイルはすべて展開し、ソースコード中のコメントやエラーメッセージに応じてパッケージをインストールしてください。


#### 第2回 2期間モデルと数値計算の概観
- [CH2/two-period discretized method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20discretized%20method.R) 第3節 離散化による2期間モデルの解法
- [CH2/two-period continuous method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20continuous%20method.R) 第4節 操作変数を連続変数にする: 最適化
- [CH2/two-period root-finding method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20root-finding%20method.R) 第5節1 非線形方程式のゼロ点を探す
- [CH2/two-period projection method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20projection%20method.R) 第5節2 射影法

#### 第4回 オイラー方程式と多項式近似
- [CH4/time iteration method.R](https://github.com/uncorrelated/KSQME/blob/master/CH4/time%20iteration%20method.R) 第2節 分権経済と時間反復法
- [CH4/approximation to the Runge function.R](https://github.com/uncorrelated/KSQME/blob/master/CH4/approximation%20to%20the%20Runge%20function.R) 第3節 多項式近似
- [CH4/time iteration chebyshev approximation method.R](https://github.com/uncorrelated/KSQME/blob/master/CH4/time%20iteration%20chebyshev%20approximation%20method.R) 第3節3 時間反復法への多項式近似の応用

#### 第5回 ニューケインジアン・モデルの新展開
- [CH5/two-state NK model closed-form.R](https://github.com/uncorrelated/KSQME/blob/master/CH5/two-state%20NK%20model%20closed-form.R) 第2節 2状態モデル
- [CH5/N-state NK model times iteration method.R](https://github.com/uncorrelated/KSQME/blob/master/CH5/N-state%20NK%20model%20times%20iteration%20method.R) 第3節 N状態モデル
- [CH5/compare a non-linear NK model to linear one.R](https://github.com/uncorrelated/KSQME/blob/master/CH5/compare%20a%20non-linear%20NK%20model%20to%20linear%20one.R) 第4節 非線形モデル

