---
title: Empirical Processes
---

# Foreword 

Now that we are comfortable with the essential probability, we cover the basic theory of empirical processes. This theory is very important moving forward. 

## More Convergence

I will now cover more theorems in weak convergence that are going to be of use. 

**Theorem 1:** Let $(\Omega_1, \cal{F}_1, P)$ be a probability space. Also let $(\Omega_2, \cal{F}_2)$ and $(\Omega_3, \cal{F}_3)$ be measurable topological spaces. Assume that $X_n:\Omega_1 \to \Omega_2$ is a sequence of random variables which converges to $X$ in law. If $g: \Omega_2 \to \Omega_3$ is a continuous function, then $g(X_n)$ converges to $g(X)$ in law. 

**Proof:** If $f: \Omega_3 \to R$ is a bounded and continuous function, then so is $f \circ g$.

$$
\int f\circ g d\mu_n \to \int f \circ g d\mu
$$

Thus, 
$$
\int f dg_*\mu_n \to \int f dg_*\mu
$$
where $g_*\mu$ is the pushforward of the pushforward. This proves the theorem. 

**Theorem 2:** $X_n$ and $Y_n$ are sequences of ranodm variables taking values on the euclidean space. 

1) If $X_n$ weakly converges to $0$, then it converges in probability (probability distribution of $0$ is the dirac delta function). 

2) If $X_n$ and $Y_n$ respectively converge to $X$ and $0$ weakly, then $X_n + Y_n$ weakly converges to $X$. 


**Theorem 3:** Assume that a sequence of random variables ${X_n}$ converges to $X$ weakly. If f is a continuous function and if 
$$
\sup_n E[f(X_n)] < C
$$
then the following hold: 
1) The integration $E[|f(X)|]$ takes a finite value. 
2) $E[|f(X)|] \leq \limsup_{n \to \infty} E[|F(X_n)|]$

**Definition (Asymptotically Uniformly Integrable)** A sequence of real valued random variables $X_n$ is said to be AUI if 
$$
\lim_{M \to \infty} \lim_{n \to \infty} \sup_{N \geq n} E[|X_N]_{\{|X_N| \geq M\}} = 0
$$

In Statistical Learning Theory, there is this case: 
$$
Z_n = f(\xi_n) + a_nX_n
$$
We prove the following conditions. 
1) $Z_n$ is AUI. 
2) $\xi_n$ weakly converges to $\xi$. 
3) The function $f$ is continuous. 
4) $a_n$ converges to $0$. 
5) $X_n$ weakly converges.

Then, we have that $Z_n \to f(\xi_n)$ weakly, and thus, their expectations converge. 
