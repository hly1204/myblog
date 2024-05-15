+++
title = '再谈 Bostan–Mori 算法'
date = 2024-05-13T19:11:17+08:00
draft = false
mathJax = true
summary = '复习一下 Bostan–Mori 算法'
+++

# 再谈 Bostan–Mori 算法

Bostan–Mori 算法可以快速完成下面的求算。如果没有特殊指定，我们一般讨论有理函数和形式幂级数或形式 Laurent 级数。

反形式 Laurent 级数在这里可以帮助我们简化一些理解。

Elegia（李白天）指出论文中的 MSB-优先算法本质为 LSB-优先算法的转置，如果熟悉转置原理，也可以以此理解并且将优化后的 LSB-优先算法机械化的转置得到优化后的 MSB-优先算法，这比直接考虑优化更简单。

我们从头介绍这个算法。

## $\frac{1}{Q(x)}\bmod{x^n}$

先从形式幂级数的乘法逆元开始，其中 $Q(x)\in\mathbb{C}\left\lbrack\left\lbrack x\right\rbrack\right\rbrack$，那么

$$
\frac{1}{Q(x)}\bmod{x^n}=\frac{1}{Q(x)Q(-x)}\cdot Q(-x)\bmod{x^n}
$$

令 $V(x^2)=Q(x)Q(-x)$ 那么只需求出

$$
\frac{1}{V(x)}\bmod{x^{\left\lceil n/2\right\rceil}}
$$

问题变为原来的一半，简单分析所需要的 FFT 次数为递归最外层的两倍：加上优化之后可以做到 $12\mathsf{E}(n)$。传统的 Newton 法简单优化之后可以做到 $10\mathsf{E}(n)$。我不清楚能否做的至少一样好。

### $Q(x)Q(-x)Q(\mathrm{i}x)Q(-\mathrm{i}x)$

这样做可以将问题进一步缩小，但是付出的代价反而变得更高，仔细思考这其实就是我们在下一次中做的事情。虽然这没有意义，但是仍然告诉我们存在 Radix-3 FFT 也方便的写法。

## $\left\lbrack x^k\right\rbrack\frac{P(x)}{Q(x)}$

为了求算 C-recursive 数列的一项，我们有其递推式和初值，那么就可以写出数列的有理函数表示 $P(x),Q(x)\in\mathbb{C}\left\lbrack x\right\rbrack$ 且 $Q(x)\neq 0$，然后使用 Bostan–Mori 算法，考虑

$$
\frac{P(x)}{Q(x)}=\frac{P(x)Q(-x)}{Q(x)Q(-x)}=\frac{U_0(x^2)+xU_1(x^2)}{V(x^2)}
$$

那么

$$
\left\lbrack x^k\right\rbrack\dfrac{P(x)}{Q(x)}=\left\lbrack x^{\left\lfloor k/2\right\rfloor}\right\rbrack \frac{U_{k\bmod{2}}(x)}{V(x)}
$$

我们付出两次多项式乘法的代价将 $k$ 变为原先的一半。

**优化**：考虑倍长 DFT 点值和取奇偶项的方法（这在原论文中有详细描述，如果熟悉 FFT 算法，会发现倍长 DFT 点值的算法其实就是 Twisted FFT 的一步：对 $f(\zeta x)$ 求 DFT 点值；取奇偶则是考虑（Radix-2 FFT 求出的）DFT 的点值一定是 $f(1),f(-1),f(\mathrm{i}),f(-\mathrm{i}),...$ 这样 **成对** 出现）。

## $\left\lbrack x^{\left\lbrack L,R\right)}\right\rbrack\frac{P(x)}{Q(x)}$

我们先考虑更简单的问题：

$$
\left\lbrack x^{\left\lbrack L,R\right)}\right\rbrack\frac{1}{Q(x)}=\left\lbrack x^{\left\lbrack L,R\right)}\right\rbrack\frac{1}{Q(x)Q(-x)}\cdot Q(-x)
$$

我们需要求出 $\left\lbrack x^{\left\lbrack L-\deg{Q},R\right)}\right\rbrack\dfrac{1}{Q(x)Q(-x)}$ 然后作一次乘法并取出 $x^L,\dots ,x^{R-1}$ 的系数。令 $V(x^2)=Q(x)Q(-x)$ 那么我们只需求出

$$
\left\lbrack x^{\left\lbrack \left\lceil\frac{L-\deg{Q}}{2}\right\rceil,\left\lceil\frac{R}{2}\right\rceil\right)}\right\rbrack\frac{1}{V(x)}
$$

就可以还原出 $\left\lbrack x^{\left\lbrack L-\deg{Q},R\right)}\right\rbrack\dfrac{1}{Q(x)Q(-x)}$。那么我们只需求出 $\left\lbrack x^{\left\lbrack L-\deg{P},R\right)}\right\rbrack\dfrac{1}{Q(x)}$ 再和 $P(x)$ 作一次乘法即可求出 $\left\lbrack x^{\left\lbrack L,R\right)}\right\rbrack\frac{P(x)}{Q(x)}$。

上面的算法虽然已经可以工作，但是每一次的递归的时间复杂度与 $R-L$ 相关，我们希望能至少在递归求算时摆脱 $R-L$，更具体的，我们先考虑求算 $\left\lbrack x^{\left\lbrack L,L+\deg Q+1\right)}\right\rbrack \dfrac{1}{Q(x)}$，考虑

$$
\left\lbrack x^{\left\lbrack L,L+\deg Q+1\right)}\right\rbrack \frac{1}{Q(x)}=\left\lbrack x^{\left\lbrack L,L+\deg Q+1\right)}\right\rbrack \dfrac{1}{Q(x)Q(-x)}\cdot Q(-x)
$$

那么我们需要求出

$$
\left\lbrack x^{\left\lbrack L-\deg Q,L+\deg Q+1\right)}\right\rbrack \dfrac{1}{Q(x)Q(-x)}
$$

那么对于 $V(x^2)=Q(x)Q(-x)$ 而言，我们只需求出

$$
\left\lbrack x^{\left\lbrack \left\lceil \frac{L-\deg Q}{2}\right\rceil,\left\lceil \frac{L+\deg Q+1}{2}\right\rceil\right)}\right\rbrack \frac{1}{V(x)}
$$

这是因为

$$
\left\lbrack x^{k}\right\rbrack\dfrac{1}{Q(x)Q(-x)}=
\begin{cases}
\left\lbrack x^{k/2}\right\rbrack\dfrac{1}{V(x)},&\text{if }k\equiv 0\pmod{2}, \\\\
0,&\text{otherwise}.
\end{cases}
$$

我们知道 $L+\deg Q$ 和 $L-\deg Q$ 的奇偶性是一样的，所以

$$
\left\lceil \frac{L+\deg Q+1}{2}\right\rceil -\left\lceil \frac{L-\deg Q}{2}\right\rceil =
\begin{cases}
\deg Q+1,&\text{if }L+\deg Q\equiv 0\pmod{2}, \\\\
\deg Q,&\text{otherwise}.
\end{cases}
$$

这样我们可以写出伪代码

$$
\begin{array}{ll}
&\textbf{Algorithm }\operatorname{Slice-CoefficientsA}(Q,L)\text{:} \\\\
&\textbf{Input}\text{: }Q(x)\in\mathbb{C}\left\lbrack x\right\rbrack,L\in\mathbb{Z}\text{.} \\\\
&\textbf{Output}\text{: }\left\lbrack x^{\left\lbrack L,L+\deg Q+1\right)}\right\rbrack Q(x)^{-1}\text{.} \\\\
1&\textbf{if }L\leq 0\textbf{ then return }\left\lbrack x^{\left\lbrack L,L+\deg Q+1\right)}\right\rbrack Q(x)^{-1} \\\\
&\text{Use other algorithm to compute }Q(x)^{-1} \\\\
2&V(x^2)\gets Q(x)Q(-x) \\\\
3&k\gets \left\lceil \frac{L-\deg Q}{2}\right\rceil \\\\
4&(t_k,\dots ,t_{k+\deg Q})\gets \operatorname{Slice-CoefficientsA}\left(V,k\right) \\\\
5&T(x)\gets x^{(L-\deg Q)\bmod{2}}\sum_{j=0}^{\deg Q}t_{j+k}x^{2j} \\\\
6&\textbf{return }\left\lbrack x^{\left\lbrack \deg Q,2\deg Q+1\right)}\right\rbrack T(x)Q(-x)
\end{array}
$$

注意当 $\deg Q=0$ 时这个算法会死循环，我们特殊处理即可（实际上我们也不会采用 $\deg Q+1$ 作为切片的长度，而是用 $2^\mathbb{N}$）。

**优化**：对于 $V(x)$ 我们可以保持其 DFT 点值的形式传递，然后在下一次使用倍长点值的算法，对于 $T(x)$ 因为它是奇/偶函数，所以 DFT 点值只需要求一半，乘法可以利用循环卷积以及返回值可以只取高位的 DFT 点值。

## 扩展 C-recursive 递推数列

这个问题是我们知道 $Q(x)$ 本身和 $Q(x)^{-1}$ 的一部分连续的系数比如 $\left\lbrack x^{\left\lbrack L,L+\deg Q\right)}\right\rbrack Q(x)^{-1}$ 和 $L\geq 0$，我们希望求出 $\left\lbrack x^{\left\lbrack L+\deg Q,L+2\deg Q\right)}\right\rbrack Q(x)^{-1}$，这等价于我们要求某个 $P(x)$ 且 $\deg P\lt \deg Q$ 使得 $\dfrac{P(x)}{Q(x)}$ 的前 $\deg Q$ 项与 $\left\lbrack x^{\left\lbrack L,L+\deg Q\right)}\right\rbrack Q(x)^{-1}$ 相同。简单来说：递推关系（有理函数的分母）是不变的，我们所做的只是更换初值（有理函数的分子）。

具体的，考虑

$$
\frac{P(x)}{Q(x)}=\sum_{j\geq 0}a_jx^j
$$

我们现在希望将递推前进 $d$ 项，那么就是

$$
\frac{U(x)}{Q(x)}=\sum_{j\geq d}a_jx^{j-d}=\frac{P(x)}{Q(x)x^d}-\frac{Q(x)\sum_{j=0}^{d-1}a_jx^j}{Q(x)x^d}
$$

### $x^k\bmod{Q(x)}$

我们可以用上面的算法来求算 $x^k\bmod{Q(x)}$，这需要我们意识到 Euclid 除法（带余除法）是在反形式 Laurent 级数环 $\mathbb{C}((x^{-1}))$ 上计算乘法逆元和乘法。例如我们熟悉的 $a(x)$ 除以 $b(x)$（设 $\deg a\geq \deg b$ 且 $b(x)\neq 0$）在 $\mathbb{C}((x^{-1}))$ 上有

$$
\frac{a(x)}{b(x)}=q(x)+\frac{r(x)}{b(x)}
$$

其中 $b(x)^{-1}=\sum_{j\geq\deg b}p_jx^{-j}$ 此时 $\deg \left(b^{-1}\right)=-\deg b$ 是对于多项式环上 **度数** 定义的自然扩展。而多项式 $q(x)$ 是 $a(x)$ 除以 $b(x)$ 的商且 $\deg q=\deg a-\deg b$ 然后多项式 $r(x)$ 就是 $a(x)$ 除以 $b(x)$ 的余数。

不妨认为 $q(x)$ 是『**整数**』部分，而 $\dfrac{r(x)}{b(x)}$ 就是『**小数**』部分（类比于整数除法），并且多项式没有令人头疼的『**进位**』。

那么 $x^k\bmod{Q(x)}$ 就是在反形式 Laurent 级数环上求出 $\left\lbrack x^{\left\lbrack -\deg Q,0\right)}\right\rbrack\dfrac{x^k}{Q(x)}=\left\lbrack x^{\left\lbrack -\deg Q-k,-k\right)}\right\rbrack\dfrac{1}{Q(x)}$ 后再作一次乘法。

# 参考文献

1. Alin Bostan, Ryuhei Mori. [A Simple and Fast Algorithm for Computing the $N$-th Term of a Linearly Recurrent Sequence](https://arxiv.org/abs/2008.08822).
2. Yasunori Kinoshita, Baitian Li. [Power Series Composition in Near-Linear Time](https://arxiv.org/abs/2404.05177).
