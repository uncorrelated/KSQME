[経済セミナー連載の「定量的マクロ経済学と数値計算」のサポート・リポジトリー](https://github.com/keizai-seminar-quant-macro)にR用のコードが部分的にしか用意されていなかったので、ちょっと書いてみたものです。

実行方法は`git clone`などをしてファイルを展開した後に、Rに読み込ませるだけです。例えば、第3節の離散化による2期間モデルの解法は、[two-period discretized method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20discretized%20method.R)を以下のように読み込ませて計算します。

	source("two-period discretized method.R")

`getwd()`で確認できるワーキングディレクトリに[two-period discretized method.R](https://github.com/uncorrelated/KSQME/blob/master/CH2/two-period%20discretized%20method.R)が無いと呼び出せないので注意してください。また、他のファイルやパッケージに依存しているファイルもあります。ファイルはすべて展開し、ソースコード中のコメントやエラーメッセージに応じてパッケージをインストールしてください。


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
- [CH5/compare a quasi-linear NK model to linear one.R](https://github.com/uncorrelated/KSQME/blob/master/CH5/compare%20a%20quasi-linear%20NK%20model%20to%20linear%20one.R) 第4節 非線形モデル（の中で言及された準線形モデル; Juliaのコードとプロット結果が微妙に異なるので要精査）

#### 第6回 異質な個人を組み込んだマクロモデル
- [CH6/aiyagari (1994) compute eq K and r.R](https://github.com/uncorrelated/KSQME/blob/master/CH6/aiyagari%20(1994)%20compute%20eq%20K%20and%20r.R) 第2節 ビューリーモデル（図3のプロット）
- [CH6/aiyagari (1994) plot capital demand and asset supply curves.R](https://github.com/uncorrelated/KSQME/blob/master/CH6/aiyagari%20(1994)%20plot%20capital%20demand%20and%20asset%20supply%20curves.R) 第2節 ビューリーモデル（図4のプロット）

第6回の処理はこれまでと比較してかなり重たい上に、Dynamic Programmingであるためか、元のMatlabのソースコードがベクトルや行列演算の形に直すのは困難に思えるものだったので、RからCで書かれた関数の拡張を呼び出すことにより、処理の高速化が行えるオプションを用意しました。LinuxやMacOSでは、以下のようにCのソースコードをコンパイルした後に、Rのソースコード中で指定するオプション`r_c_mp_switch`を変更することで、数十倍の高速化が図れます。

	R CMD SHLIB aiyagari_vfi1.c
	R CMD SHLIB aiyagari_vfi2.c

Windowsでも同様にdllを作れ実行できるはずですが、まだ試していません。

図4のプロットでは、これもオプションとして、さらにマルチコア対応となっています。Cの拡張とマルチコアを利用することで、純粋なRのコードと比較してかなりの処理時間の短縮が可能になります。手元の古めの計算機では70〜140倍といった差異になりました。

#### 第7回 世代重複マクロモデル
- [CH7/profile for lifetime assets and consumption by OLG.R](https://github.com/uncorrelated/KSQME/blob/master/CH7/profile%20for%20lifetime%20assets%20and%20consumption%20by%20OLG.R) 第2節 世代重複モデル（図3のプロット）
- [CH7/profile for lifetime assets and consumption by OLG in various scenario.R](https://github.com/uncorrelated/KSQME/blob/master/CH7/profile%20for%20lifetime%20assets%20and%20consumption%20by%20OLG%20in%20various%20scenario.R) 第2節 世代重複モデル（図4のプロット）

