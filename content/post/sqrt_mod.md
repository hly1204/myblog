+++
title = '模平方根'
date = 2024-09-21T13:19:48+08:00
draft = false
mathJax = true
summary = '简单复习一下 Tonelli–Shanks 算法'
+++

# 有限域平方根的 Tonelli–Shanks 算法

Cipolla 的算法相信大家都会，我们就来看看 Tonelli–Shanks 算法，我们要明确 Tonelli–Shanks 算法是基于『**离散对数**』的。

如果需要扩展到模奇素数幂次 $p^k$ 的平方根，我们可以照搬『**形式幂级数环**』$\mathbb{F}_p\lbrack\lbrack x\rbrack\rbrack$ 上对于平方根的 Newton 法处理，所以后文仅讨论简单的模奇素数 $p$ 下的求算。

## 记号

- $p$ 为奇素数；

- $\operatorname{ord}a$ 为 $a$ 在模 $p$ 意义下的阶，我们省略 $p$ 因为上下文中 $p$ 是明确的。即 $\operatorname{ord}a=\min\left\lbrace k:k\gt 0,a^k\equiv 1\pmod{p}\right\rbrace$；

- $\operatorname{ind}_{g}a=k$ 为选取在模 $p$ 下的原根 $g$ 后 $g^k\equiv a\pmod{p}$。

## 算法

假设 $a$ 是二次剩余，不妨设 $p-1=2^n\cdot m$ 其中 $m$ 为奇数。$r$ 为二次非剩余，令 $c\equiv r^m\pmod{p}$ 和 $b\equiv a^{(m-1)/2}\pmod{p}$，存在整数 $e\in\left\lbrace 0,\dots,2^n-1\right\rbrace$ 满足 $ab^2\equiv c^e\pmod{p}$ 且

$$
\left(abc^{-e/2}\right)^2\equiv a\pmod{p}
$$

**证明**：考虑 $\operatorname{ord}c$，因为 $c^{2^{n-1}}\equiv r^{\frac{p-1}{2}}\equiv -1\pmod{p}$ 所以 $c^{2^{\lt n}}\not\equiv 1\pmod{p}$ 且 $c^{2^n}\equiv 1\pmod{p}$ 所以 $\operatorname{ord}c=2^n$。

又因为 $ab^2\equiv a^m\pmod{p}$ 是 $x^{2^n}\equiv 1\pmod{p}$ 的解，所以 $a^m$ 是 $c$ 的幂次，记为 $a^m\equiv c^e\pmod{p}$。

因为 $a$ 是二次剩余，所以

$$
\begin{aligned}
c^{2^{n-1}\cdot e}&\equiv \left(-1\right)^e &\pmod{p} \\\\
&\equiv a^{2^{n-1}\cdot m} &\pmod{p} \\\\
&\equiv a^{\frac{p-1}{2}} &\pmod{p} \\\\
&\equiv 1 &\pmod{p}
\end{aligned}
$$

$e$ 是偶数，那么

$$
\begin{aligned}
\left(abc^{-e/2}\right)^2 &\equiv a^2b^2c^{-e} &\pmod{p} \\\\
&\equiv a^{m+1}c^{-e} &\pmod{p} \\\\
&\equiv a &\pmod{p}
\end{aligned}
$$

Tonelli 和 Shanks 指出可以逐比特计算 $e$，令 $e=e_0+2e_1+4e^2+\cdots$ 其中 $e_k\in\left\lbrace 0,1\right\rbrace$。然后开始时 $e_0=0$，紧接着计算 $e_1,e_2,\dots$，由下式给出

$$
\left(c^ec^{-\left(e\bmod{2^k}\right)}\right)^{2^{n-1-k}}=c^{2^{n-1}\cdot e_k}=
\begin{cases}
1, &\text{if }e_k=0, \\\\
-1, &\text{if }e_k=1.
\end{cases}
$$

因为 $c^e$ 和 $e\bmod{2}$ 已知，可以递推出 $e_1,\dots$。不能直接写作 $\left(c^{\left\lfloor e/2^k\right\rfloor}\right)^{2^{n-1}}$ 因为我们只知道 $e$ 的最低位。

### 扩展

Adleman，Manders 和 Miller 将 Tonelli–Shanks 算法扩展到计算立方根，五次方根等，我们来看看他们是怎么做的。

#### 立方根

我们先要判断三次非剩余是否存在，使用 Euler 判别准则即可，若三次非剩余不存在，意味着 $\gcd(p-1,3)=1$ 此时可以直接计算出立方根有且只有一个。后文假设三次非剩余存在。

设 $a$ 是三次剩余，同样设 $p-1=3^n\cdot m$ 其中 $m$ 不是 $3$ 的倍数。$r$ 为三次非剩余，$c\equiv r^m\pmod{p}$。然后我们如法炮制计算出 $e$ 使得 $a^m\equiv c^e\pmod{p}$。这时候我们记 $m=3j+k$ 且 $k\in\left\lbrace 1,2\right\rbrace$，那么

$$
a^{3j}a^k\equiv c^e\pmod{p}
$$

令

$$
t=\begin{cases}
c^{e/3}a^{-j}, &\text{if }k=1, \\\\
c^{-e/3}a^{j+1}, &\text{if }k=2.
\end{cases}
$$

满足 $t^3\equiv a\pmod{p}$。当 $k=2$ 时我们也有更通用的处理方式：考虑 $a^{2m}\equiv c^{2e}\pmod{p}$，那么令 $2m=3j'+1$ 即可，这一方式应该可以推广至五次方根等。

## 附录

### Euler 判别准则

设 $n\in\mathbb{Z}_{\gt 0}$，$\gcd(a,p)=1$，那么 $a$ 是一个 $n$ 次剩余当且仅当

$$
a^{(p-1)/d}\equiv 1\pmod{p}
$$

其中 $d=\gcd(p-1,n)$。$n$ 次剩余的数量为 $\frac{p-1}{d}$，每个都是 $d$ 个整数模 $p$ 下的 $n$ 次方。

**证明**：考虑选取一原根 $g$，那么 $\mathbb{F}_p^{\times}$ 的元素和 $g^0,\dots,g^{p-2}$ 一一对应。对于方程 $x^n\equiv a\pmod{p}$ 两边取索引，我们有

$$
n\operatorname{ind} _ {g}x\equiv \operatorname{ind} _ {g}a\pmod{(p-1)}
$$

方程有解当且仅当 $d\mid \operatorname{ind} _ {g}a$。因为我们要求解 $\operatorname{ind} _ {g}x$，显然只有 $n\operatorname{ind} _ {g}x\equiv 0\equiv \operatorname{ind} _ {g}a\pmod{d}$ 才有解，而 $\operatorname{ind} _ {g}a\equiv 0\pmod{d}$ 意味着 $a\equiv g^{kd}\pmod{p}$ 对于某个 $k\in\mathbb{Z}$ 成立。

最后，显然有 $g^d,\dots,g^{\frac{p-1}{d}d\bmod{(p-1)}}$ 为 $\frac{p-1}{d}$ 个不同的 $n$ 次剩余。

# 参考文献

1. Daniel. J. Bernstein. [Faster Square Roots in Annoying Finite Fields](https://cr.yp.to/papers.html#sqroot).
2. Leonard Adleman, Kenneth Manders, Gary Miller. On taking roots in finite fields.
