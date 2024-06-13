+++
title = 'Chirp Z 变换'
date = 2024-06-12T22:39:13+08:00
draft = false
mathJax = true
summary = '学习一下逆 Chirp Z 变换'
+++

# Chirp Z 变换 (CZT)

这个东西很多人都很熟悉了，我们换一个不需要翻转多项式的方法。

**问题**：给出 $f(x)=\sum_{i=0}^{m-1}f_ix^i\in\mathbb{C}\left\lbrack x\right\rbrack$ 和 $q\in\mathbb{C}$，求出 $f(1),f(q),\dots ,f(q^{n-1})$。

当 $q=0$ 时显然，我们考虑 $q\neq 0$ 时。考虑

$$
ij=\binom{i}{2}+\binom{-j}{2}-\binom{i-j}{2}
$$

其中 $\dbinom{a}{2}=\dfrac{a(a-1)}{2}$。那么我们构造

$$
\begin{aligned}
G(x)&:=\sum_{i=-(m-1)}^{n-1}q^{-\binom{i}{2}}x^i, \\\\
F(x)&:=\sum_{i=0}^{m-1}f_iq^{\binom{-i}{2}}x^i
\end{aligned}
$$

那么对于 $i=0,\dots ,n-1$ 有

$$
\begin{aligned}
\left\lbrack x^i\right\rbrack\left(G(x)F(x)\right)&=\sum_{j=0}^{m-1}\left(\left\lbrack x^{i-j}\right\rbrack G(x)\right)\left(\left\lbrack x^j\right\rbrack F(x)\right) \\\\
&=\sum_{j=0}^{m-1}f_jq^{\binom{-j}{2}-\binom{i-j}{2}} \\\\
&=q^{-\binom{i}{2}}\cdot f\left(q^i\right)
\end{aligned}
$$

一次卷积，利用循环卷积优化。时间为 $O\left(\mathsf{M}\left(n+m\right)\right)$。

**备注**：$q^{\binom{i+1}{2}}=q^{\binom{i}{2}}\cdot q^i$，$\dbinom{-i}{2}=\dbinom{i+1}{2}$。

# 逆 Chirp Z 变换 (ICZT)

**问题**：给出 $f(1),f(q),\dots ,f(q^{n-1})$ 和 $q\in\mathbb{C}\setminus\left\lbrace 0\right\rbrace$ 且 $q^i\neq q^j,\forall i\neq j$ 求出 $f(x)=\sum_{i=0}^{n-1}f_ix^i\in\mathbb{C}\left\lbrack x\right\rbrack_{\lt n}$ 的系数。

回顾 Lagrange 插值公式：对于 $x_0,\dots ,x_{n-1}\in\mathbb{C}$ 有 $x_i\neq x_j,\forall i\neq j$，此时

$$
f(x)=\sum_{i=0}^{n-1}f\left(x_i\right)\prod_{0\leq j\lt n\atop j\neq i} \frac{x-x_j}{x_i-x_j}
$$

> L'Hôpital 法则：
> 对于多项式的导数我们定义为 $f'=\sum_{i}if_ix^{i-1}$，对于 $g\in\mathbb{C}\left\lbrack x\right\rbrack$ 有乘法法则：$\left(fg\right)'=f'g+fg'$。若对于 $\alpha\in\mathbb{C}$ 有 $f(\alpha)=0$ 那么 $\left(fg\right)'(\alpha)=f'(\alpha)g(\alpha)$，所以若 $f'(\alpha)\neq 0$ 那么 $(fg/f)(\alpha)=g(\alpha)=(fg)'(\alpha)/f'(\alpha)$。

令

$$
M(x):=\prod_{i=0}^{n-1}\left(x-x_i\right)
$$

在这里 $M(x_i)=0\neq M'(x_i)$ 因为对任意 $i=0,\dots ,n-1$ 都可以写成 $M(x)=(x-x_i)\prod_{0\leq j\lt n\atop j\neq i}(x-x_j)$，所以

$$
M'(x_i)=\prod_{0\leq j\lt n\atop j\neq i}\left(x_i-x_j\right)
$$

此时 Lagrange 插值公式被改为

$$
f(x)=M(x)\cdot \sum_{i=0}^{n-1}\frac{f\left(x_i\right)}{M'(x_i)}\frac{1}{x-x_i}
$$

再代入我们的和 $q$ 相关的情况就是

$$
f(x)=M(x)\cdot\sum_{i=0}^{n-1}\frac{f\left(q^i\right)}{M'(q^i)}\frac{1}{x-q^i}
$$

在这里 $M(x)=\prod_{j=0}^{n-1}\left(x-q^j\right)$，不妨设 $n-1=2k$ 和

$$
H(x):=\prod_{j=0}^{k-1}\left(x-q^j\right)
$$

那么

$$
M(x)=H(x)\cdot q^{k^2}\cdot H\left(x/q^k\right)
$$

因为 $T(n)=T(n/2)+O\left(\mathsf{M}(n)\right)+O(n)\implies T(n)=O(\mathsf{M}(n))$ 可以在 $O\left(\mathsf{M}\left(n\right)\right)$ 求出 $M(x)$，然后利用 CZT 求出 $M'(1),\dots ,M'(q^{n-1})$，设 $c_i=f\left(q^i\right)/M'\left(q^i\right)$ 那么

$$
f(x)=M(x)\cdot\sum_{i=0}^{n-1}\frac{c_i}{x-q^i}
$$

我们知道 $f(x)\in\mathbb{C}\left\lbrack x\right\rbrack_{\lt n}$ 所以只需计算出 $\sum_{i=0}^{n-1}\frac{c_i}{x-q^i}\bmod{x^n}$ 即可：

$$
\begin{aligned}
\sum_{i=0}^{n-1}\frac{c_i}{x-q^i}\bmod x^n&=-\sum_{i=0}^{n-1}\left(\sum_{j=0}^{n-1}c_iq^{-i(j+1)}x^j\right) \\\\
&=-\sum_{j=0}^{n-1}C\left(q^{-j-1}\right)x^j
\end{aligned}
$$

其中

$$
C(x)=\sum_{i=0}^{n-1}c_ix^i
$$

可以用 CZT 求出 $C\left(q^{-1}\right),\dots ,C\left(q^{-n}\right)$。

noshi91 指出，$M(x)$ 和 $M'(1),\dots ,M'\left(q^{n-1}\right)$ 可以用更快的方法计算。

我们重新考虑 $\prod_{0\leq j\lt n\atop j\neq i}\left(q^i-q^j\right)$，遵循 noshi91 的定义，令 $s_i:=\prod_{j=1}^i\left(1-q^j\right)$ 且 $s_0:=1$，那么

$$
\begin{aligned}
\prod_{0\leq j\lt n\atop j\neq i}\left(q^i-q^j\right)&=\left(\prod_{j=0}^{i-1}\left(q^i-q^j\right)\right)\left(\prod_{j=i+1}^{n-1}\left(q^i-q^j\right)\right) \\\\
&=\left(\prod_{j=0}^{i-1}q^j\left(q^{i-j}-1\right)\right)\left(\prod_{j=i+1}^{n-1}q^i\left(1-q^{j-i}\right)\right) \\\\
&=(-1)^iq^{\sum_{j=0}^{i-1}j}\left(\prod_{j=0}^{i-1}\left(1-q^{i-j}\right)\right)\cdot q^{i(n-i-1)}\left(\prod_{j=i+1}^{n-1}\left(1-q^{j-i}\right)\right) \\\\
&=(-1)^iq^{\binom{i}{2}}\left(\prod_{k=1}^{i}\left(1-q^k\right)\right)\cdot q^{i(n-i-1)}\left(\prod_{k=1}^{n-i-1}\left(1-q^k\right)\right) \\\\
&=(-1)^iq^{\binom{i}{2}}s_i\cdot q^{i(n-i-1)}s_{n-i-1}
\end{aligned}
$$

注意到 $q^{\binom{i+1}{2}}\cdot q^{(i+1)(n-(i+1)-1)}=\left(q^{\binom{i}{2}}\cdot q^{i}\right)\left(q^{i(n-i-1)}\cdot q^{n-2i-2}\right)$ 可递推计算。

## $q$-模拟 ($q$-analog)

我们定义 $n\in\mathbb{N}$ 的 $q$-模拟为

$$
\left\lbrack n\right\rbrack_q:=\begin{cases}
0,&\text{if }n=0, \\\\
1+q+\cdots +q^{n-1},&\text{otherwise}.
\end{cases}
$$

注意到 $q=1$ 时，有 $\left\lbrack n\right\rbrack_q=n$。并且 $q\left\lbrack n\right\rbrack_q =q+q^2+\cdots +q^n$ 所以

$$
\left\lbrack n\right\rbrack_q=\begin{cases}
n,&\text{if }q=1, \\\\
\dfrac{1-q^n}{1-q},&\text{otherwise}.
\end{cases}
$$

我们也可以发现

$$
\begin{aligned}
\left\lbrack -n\right\rbrack_q&=\frac{1-q^{-n}}{1-q} \\\\
&=q^{-n}\frac{q^n-1}{1-q} \\\\
&=-q^{-n}\left\lbrack n\right\rbrack_q
\end{aligned}
$$

然后我们可以定义 $q$-阶乘

$$
n!_q:=\begin{cases}
1,&\text{if }n=0, \\\\
\left\lbrack 1\right\rbrack_q\left\lbrack 2\right\rbrack_q\cdots \left\lbrack n\right\rbrack_q,&\text{otherwise}.
\end{cases}
$$

最后我们定义 $q$-二项式系数

$$
\binom{n}{k}_q:=\begin{cases}
\dfrac{n!_q}{k!_q(n-k)!_q},&\text{if }0\leq k\leq n, \\\\
0,&\text{otherwise}.
\end{cases}
$$

为了证明下面的定理，我们先引入一些常见的性质。

若 $q\neq 1$ 我们有

$$
\begin{aligned}
\left\lbrack n\right\rbrack_q&=\frac{1-q^n}{1-q} \\\\
&=\frac{1-q^k+q^k-q^n}{1-q} \\\\
&=\frac{1-q^k}{1-q}+q^k\frac{1-q^{n-k}}{1-q} \\\\
&=\left\lbrack k\right\rbrack_q +q^k\left\lbrack n-k\right\rbrack_q
\end{aligned}
$$

而当 $q=1$ 时我们有 $n=k+(n-k)$ 所以上式仍然成立。因此

$$
\begin{aligned}
\binom{n+1}{k}_q&=\frac{(n+1)!_q}{k!_q(n+1-k)!_q} \\\\
&=\frac{n!_q}{k!_q(n+1-k)!_q}\cdot \left(\left\lbrack k\right\rbrack_q +q^k\left\lbrack n+1-k\right\rbrack_q\right) \\\\
&=\frac{n!_q}{(k-1)!_q(n+1-k)!_q}\cdot \left\lbrack k\right\rbrack_q^{-1}\cdot \left(\left\lbrack k\right\rbrack_q +q^k\left\lbrack n+1-k\right\rbrack_q\right) \\\\
&=\frac{n!_q}{(k-1)!_q(n+1-k)!_q}+q^k\frac{n!_q\left\lbrack n+1-k\right\rbrack_q}{k!_q(n+1-k)!_q} \\\\
&=\binom{n}{k-1}_q+q^k\binom{n}{k}_q
\end{aligned}
$$

因为 $q$-二项式系数也有对称性

$$
\binom{n}{k}_q=\binom{n}{n-k}_q
$$

我们将上式 $k$ 用 $n+1-k$ 替换得到了

$$
\begin{aligned}
\binom{n+1}{k}_q&=\binom{n+1}{n+1-k}_q \\\\
&=\binom{n}{n-k}_q+q^{n+1-k}\binom{n}{n+1-k}_q \\\\
&=\binom{n}{k}_q+q^{n-k+1}\binom{n}{k-1}_q
\end{aligned}
$$

这被称为 $q$-Pascal 递推式。

### Rothe 的 $q$-二项式定理

对于变量 $q,a,x$ 有

$$
\prod_{i=0}^{n-1}\left(a+q^ix\right)=\sum_{k=0}^n\binom{n}{k}_qq^{\binom{k}{2}}a^{n-k}x^k
$$

特别的，当 $n=0$ 时左式定义为 $1$。

**证明**：设 $r_n(x,a)$ 等于右式，应用上述 $q$-Pascal 递推式，我们有

$$
\begin{aligned}
r_{n+1}(x,a)&=\sum\_{k=0}^{n+1}\binom{n+1}{k}_qq^{\binom{k}{2}}x^ka^{n+1-k} \\\\
&=\sum\_{k=0}^{n+1}\left(\binom{n}{k-1}_q+q^k\binom{n}{k}_q\right)q^{\binom{k}{2}}x^ka^{n+1-k} \\\\
&=\sum\_{k=1}^{n+1}\binom{n}{k-1}_qq^{\binom{k}{2}}x^ka^{n+1-k}+\sum\_{k=0}^{n}\binom{n}{k}_qq^{\binom{k}{2}}\left(qx\right)^ka^{n+1-k} \\\\
&=\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j+1}{2}}x^{j+1}a^{n-j}+\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}\left(qx\right)^ja^{n+1-j} 
\end{aligned}
$$

第三行中我们使用 $\dbinom{n}{-1}_q=\dbinom{n}{n+1}_q=0$，前面我们提到对于 $j\in\mathbb{N}$ 有 $\dbinom{j+1}{2}=\dbinom{j}{2}+j$，那么

$$
\begin{aligned}
r_{n+1}(x,a)&=\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}q^jx^{j+1}a^{n-j}+\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}\left(qx\right)^ja^{n+1-j} \\\\
&=\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}\left(qx\right)^ja^{n-j}x+\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}\left(qx\right)^ja^{n+1-j} \\\\
&=\sum\_{j=0}^{n}\binom{n}{j}_qq^{\binom{j}{2}}\left(qx\right)^ja^{n-j}(x+a) \\\\
&=(a+x)r_n(qx,a)
\end{aligned}
$$

重复下去我们有

$$
\begin{aligned}
r_{n+1}(x,a)&=(a+x)r_n(qx,a) \\\\
&=(a+x)(a+qx)r_{n-1}\left(q^2x,a\right) \\\\
&=\cdots \\\\
&=(a+x)(a+qx)\cdots \left(a+q^nx\right)r_0\left(q^{n+1}x,a\right)
\end{aligned}
$$

而我们已经定义了 $r_0(u,v)=1$。

回到我们的问题，那么

$$
\begin{aligned}
M(x)&=\prod_{j=0}^{n-1}\left(x-q^j\right) \\\\
&=\sum_{k=0}^n\binom{n}{k}_qq^{\binom{k}{2}}(-1)^kx^{n-k}
\end{aligned}
$$

而 $\dbinom{n}{k}_q=\dfrac{\left(1-q\right)\cdots \left(1-q^n\right)}{\left(1-q\right)\cdots \left(1-q^k\right)\cdot \left(1-q\right)\cdots \left(1-q^{n-k}\right)}=\dfrac{s_n}{s_ks\_{n-k}}$。

**备注**：若 $q^n=1$ 不能再使用上式，我们按照定义计算。

# 参考文献

1. 37zigen. [多項式補間：アルゴリズム](https://37zigen.com/lagrange-interpolation/).
2. noshi91. [標本点が等比数列を成す場合に補間多項式を計算するアルゴリズム](https://noshi91.github.io/algorithm-encyclopedia/polynomial-interpolation-geometric).
3. Bostan, A. (2010). [Fast algorithms for polynomials and matrices. JNCF 2010. Algorithms Project, INRIA](https://specfun.inria.fr/bostan/publications/exposeJNCF.pdf).
4. Warren P. Johnson. An Introduction to $q$-analysis.
