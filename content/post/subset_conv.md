+++
title = '如何加速子集卷积'
date = 2024-10-06T00:40:13+08:00
draft = true
mathJax = true
summary = '简单复习一下子集卷积和一些简单的优化方法'
+++

# 子集卷积

我们定义在集合 $\mathcal{N}=\left\lbrace 1,\dots ,n\right\rbrace$ 的[幂集](https://en.wikipedia.org/wiki/Power_set)上的函数 $f,g:2^{\mathcal{N}}\to R$，那么对于所有 $\mathcal{S}\subseteq \mathcal{N}$ 来说 $f$ 和 $g$ 的子集卷积就是

$$
\left(f*g\right)(\mathcal{S}):=\sum_{\mathcal{T}\subseteq \mathcal{S}}f(\mathcal{S})g\left(\mathcal{S}\setminus \mathcal{T}\right)
$$

且计算子集卷积的时间为 $O\left(n^22^n\right)$。

## 问题转化

如果直接从 Björklund et. al. 的论文来解释，对于我来说太过困难，如果我们可以将问题转换为我们已经熟悉的东西，那么就简单了。

如果我们引入 $n$ 个未定元 $x_1,\dots,x_n$，显然我们可以将函数 $f,g$ 和多项式 $F\in R\left\lbrack x_1,\dots,x_n\right\rbrack$ 的系数联系起来。只需要令

$$
\left\lbrack x_1^{\left\lbrack 1\in \mathcal{S}\right\rbrack}\dots x_n^{\left\lbrack n\in \mathcal{S}\right\rbrack}\right\rbrack F=f(\mathcal{S})
$$

其中 $\left\lbrack 1\in \mathcal{S}\right\rbrack$ 等为 [Iverson 括号](https://en.wikipedia.org/wiki/Iverson_bracket)。$G$ 的系数也同样对应函数 $g$。

我们观察子集卷积的某个子集的贡献，对应到多项式乘法中如果两个集合有交集则会贡献到 $x_k^2$ 的系数，子集卷积所需要的信息还是可以得到保留。

这样的转化有一个好处：我们可以将系数截断至 $\mathbf{x}^2$ 也就是 $R\left\lbrack \mathbf{x}^2\right\rbrack /\mathbf{x}^2$ 上的乘法。后面我们用 $\mathbf{x}$ 来替代 $x_1,\dots,x_n$。

后面我们考虑计算 $FG\bmod{\mathbf{x}^2}$。

## FFT

首先我们知道环 $R\left\lbrack x_1,\dots,x_n\right\rbrack=R\left\lbrack x_1,\dots,x_{n-1}\right\rbrack \left\lbrack x_n\right\rbrack$ 这样迭代定义的，那么对于每一元而言，我们只需计算其在 $\pm 1$ 处的值，那么就可以计算出

$$
FG\bmod{\left(\mathbf{x}^2-1\right)}
$$

这要求 $2_R$ 存在乘法逆元。但是最后我们仍然无法得到 $FG\bmod{\mathbf{x}^2}$，因为原先贡献给 $\mathbf{x}^2$ 的系数现在被贡献到 $\mathbf{x}^0$ 的系数中，我们没有办法分离出来。

这里我们考虑引入一元 $z$ 然后计算 $FG\bmod{\left(\mathbf{x}^2-z\right)}\in R\left\lbrack z\right\rbrack\left\lbrack \mathbf{x}\right\rbrack$，显然

$$
FG\bmod{\mathbf{x}^2}=\left\lbrack z^0\right\rbrack \left(FG\bmod{\left(\mathbf{x}^2-z\right)}\right)
$$

对于结果的每一项我们定义 $\deg t:=\deg_zt+\sum_{j=1}^n\deg_{x_j}t$，考虑乘法之后带来的影响：对于每一项而言，如果因为对 $\mathbf{x}^2-z$ 取模导致 $\deg t$ 减少 $k$ 那么说明 $\deg_zt$ 增加 $k$ 且 $\sum_{j=1}^n\deg_{x_j}t$ 减少 $2k$。例如考虑 $(x+y)^k$ 和 $(x+1)^k$，因为前者是[齐次多项式](https://en.wikipedia.org/wiki/Homogeneous_polynomial)，所以我们只需计算后者也可以知道前者中每一项对应 $y$ 的度数。

这已经几乎给出了算法：因为我们最后相当于需要令 $z\gets 0$，所以我们需要使得不含有 $z$ 的结果的项可以被分离出来，也就是说含有 $z$ 的项和不含有 $z$ 的项的系数是不能分离的。那么我们将 $F$ 和 $G$ 分成 $n+1$ 个多项式 $F=\sum_{j=0}^nF_j$ 其中 $F_j$ 中的每一项都是齐次多项式也就是度数都为 $j$，$G$ 同理。然后计算 $\left(\sum_{j=0}^nF_j\right)\left(\sum_{j=0}^nG_j\right)\bmod{\left(\mathbf{x}^2-z\right)}$，我们在模 $\mathbf{x}^2-1$ 下进行 $O\left(n^2\right)$ 次的齐次多项式的乘法，并且将次数减少的项丢掉（根据上述方法我们可以还原出 $z$ 的次数，但是这不是必要的）。

## Zeta/Moebius 变换

上面的限制中要求了 $2_R$ 存在乘法逆元，而因为可以任取两个点，如果我们选取 $0,1$ 的话，就是转而计算

$$
FG\bmod{\left(\mathbf{x}^2-z\mathbf{x}\right)}\in R\left\lbrack z\right\rbrack\left\lbrack \mathbf{x}\right\rbrack
$$

事实上这一选择是有很大的好处的，这也正是 Björklund et. al. 给出的算法。Zeta 和 Moebius 变换分别对应 FFT 和 iFFT 变换。

## 优化

TODO

## 附录

TODO

# 参考文献

1. A. Björklund, T. Husfeldt, P. Kaski, and M. Koivisto. [Fourier meets Möbius: fast subset convolution](https://arxiv.org/abs/cs/0611101).
