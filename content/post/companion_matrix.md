+++
title = '友矩阵'
date = 2024-06-23T11:53:52+08:00
draft = false
mathJax = true
summary = '简单学习一下友矩阵'
+++

# 友矩阵 (companion matrix)

对于首一多项式

$$
g(x)=g_0+g_1x+\cdots +g_{d-1}x^{d-1}+x^d\in\mathbb{K}\left\lbrack x\right\rbrack
$$

的友矩阵式是

$$
\mathbf{C} _ g:=
\begin{bmatrix}
&&&-g_0 \\\\
1&&&-g_1 \\\\
&\ddots &&\vdots \\\\
&&1&-g_{d-1}
\end{bmatrix}
\in\mathbb{K}^{d\times d}
$$

如果我们将 $\mathbf{C} _ g$ 乘以一个向量 $B_b$，其中

$$
b(x):=b_0+b_1x+\cdots +b_{d-1}x^{d-1}\in\mathbb{K}\left\lbrack x\right\rbrack_{\lt d}
$$

和

$$
B_b:=\begin{bmatrix}
b_0 \\\\
b_1 \\\\
\vdots \\\\
b_{d-1}
\end{bmatrix}
$$

那么

$$
\underbrace{\begin{bmatrix}
&&&-g_0 \\\\
1&&&-g_1 \\\\
&\ddots &&\vdots \\\\
&&1&-g_{d-1}
\end{bmatrix}} _ {\mathbf{C} _ g}
\underbrace{\begin{bmatrix}
b_0 \\\\
b_1 \\\\
\vdots \\\\
b_{d-1}
\end{bmatrix}} _ {B_b}=
\underbrace{\begin{bmatrix}
-g_0b_{d-1} \\\\
b_0-g_1b_{d-1} \\\\
\vdots \\\\
b_{d-2}-g_{d-1}b_{d-1}
\end{bmatrix}} _ {B_{xb\bmod{g}}}
$$

显然有 $\left(\mathbf{C} _ {g}\right)^kB_b=B_{x^kb\bmod{g}}$。我们也可以观察 $\left(\mathbf{C}_g\right)^k$ 的性质：$\left(\mathbf{C}_g\right)^{k+1}$ 是 $\left(\mathbf{C}_g\right)^k$ 删掉第一列然后在末尾加入一列生成的，比如

$$
\begin{aligned}
\mathbf{C} _ g&=\begin{bmatrix}B_{x\bmod{g}}&B_{x^2\bmod{g}}&\cdots &B_{x^d\bmod{g}}\end{bmatrix}, \\\\
\left(\mathbf{C} _ g\right)^2&=\begin{bmatrix}B_{x^2\bmod{g}}&B_{x^3\bmod{g}}&\cdots &B_{x^{d+1}\bmod{g}}\end{bmatrix}, \\\\
\vdots \\\\
\left(\mathbf{C} _ g\right)^k&=\begin{bmatrix}B_{x^k\bmod{g}}&B_{x^{k+1}\bmod{g}}&\cdots &B_{x^{k+d}\bmod{g}}\end{bmatrix}
\end{aligned}
$$

也就是说我们要求 $\left(\mathbf{C} _ {g}\right)^k$ 的一列可以通过『**多项式**』相关的算法解决。

## 常系数齐次线性递推数列

对于数列 $\left(a_j\right)_{j\geq 0}$ 和其递推式

$$
a_n:=\sum_{j=1}^{d}c_ja_{n-j},\quad (n\geq d)
$$

其中 $c_j$ 不全为零，我们的目标是在给出初值 $a_0,\dots ,a_{d-1}$ 和递推式中的 $c_1,\dots ,c_d$ 后求出 $a_k$。

这里 $\left(a_j\right)_{j\geq 0}$ 被称为 $d$ 阶的常系数齐次线性递推数列。

我们不妨尝试用矩阵描述：令

$$
C(x):=x^d-\sum_{j=1}^dc_{d-j+1}x^{j-1}
$$

那么对于 $t\in\mathbb{N}$ 有

$$
\begin{bmatrix}
a_{t+1} \\\\
a_{t+2} \\\\
\vdots \\\\
a_{t+d}
\end{bmatrix}=\underbrace{\begin{bmatrix}
&1&& \\\\
&&\ddots & \\\\
&&&1 \\\\
c_d&c_{d-1}&\cdots &c_1
\end{bmatrix}} _ {\left(\mathbf{C} _ C\right)^{\intercal}}
\begin{bmatrix}
a_t \\\\
a_{t+1} \\\\
\vdots \\\\
a_{t+d-1}
\end{bmatrix}
$$

且

$$
\begin{bmatrix}
a_{k} \\\\
a_{k+1} \\\\
\vdots \\\\
a_{k+d-1}
\end{bmatrix}=\underbrace{\begin{bmatrix}
&1&& \\\\
&&\ddots & \\\\
&&&1 \\\\
c_d&c_{d-1}&\cdots &c_1
\end{bmatrix}^k} _ {\left(\left(\mathbf{C} _ C\right)^{\intercal}\right)^k=\left(\left(\mathbf{C} _ C\right)^{k}\right)^{\intercal}}
\begin{bmatrix}
a_0 \\\\
a_{1} \\\\
\vdots \\\\
a_{d-1}
\end{bmatrix}
$$

我们知道 $\left(\left(\mathbf{C} _ C\right)^{k}\right)^{\intercal}$ 的第一行为 $B_{x^k\bmod{C}}$，这就是求解常系数齐次线性递推数列的 Fiduccia 算法。
