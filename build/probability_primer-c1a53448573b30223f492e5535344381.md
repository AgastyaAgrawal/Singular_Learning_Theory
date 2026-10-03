---
title: A Primer on Probability
---

# Foreword 

I will now comprehensively cover the essential probability theory that you need to be equipped with to be able to handle the rest of the material. I will start with convergence in probability. 

# Modes of Convergence 

Let us fix our triple $(\Omega, \mathcal{A}, \mu)$. Random variables are measurable functions from this space to $\mathbb{R}$. The key intuition to understanding these definitions is to try coming up with them. We want to have some convergence definitions of random variables. But we know real analysis, and we want to involve the sample space in the definition. Hence it is natural to define convergence of random variables in terms of the measure of events which describe their convergence. 
1) Almost Sure convergence: An event happens almost surely if the measure of the set in which the event does not happen is 0. Let us use this for convergence. 
$$
X_n \to X \text{ a.s} \iff \mu(\{\omega: X_n(\omega) \not\to X(\omega)\}) = 0
$$

This convergence is in relation to pointwise convergence event. 

2) Convergence in Measure: A sequence $(X_n)$ converges to $X$ if $\mu(|X_n - X| > \epsilon) \to 0$ for every $\epsilon > 0$. Think about this definition a bit. We are calculating the measure of an event, and we need to have some convergence property here. Hence the necessity of such a "contrived" definition. 

I have to be clear about the notation here. By $\mu(|X_n - X| > \epsilon)$, I mean $\mu(\omega: |X_n(w) - X(w)| > \epsilon)$, since that is where our measure is defined. I will not explain this later and continue on with the abuse of notation. I would also like to point out that this makes sense because our triple is fixed. More generally, we can have the random variables codomain be different spaces, but our preimage must be the same measure space for each random variable. 

This convergence is in relation to the function distance. 
(To talk about function convergence, one can look at pointwise convergence, or one can think about the function difference and consider some norm there. Since we are dealing with real numbers here, we took the standard norm on the real numbers). 

3) $L^p$ convergence: $X_n \to X$ in $L^p$, if $\int |X_n - X|^p d\mu \to 0$. 

Framed another way, this is $E[|X_n - X|^p] \to 0$. 

This is using the $L^p$ metric and considering the expectation instead. We will always be considering $p \geq 1$, as that is when it is a metric space. 

4) There is another convergence, called convergence in distribution/convergence in law/weak convergence, which we will cover shortly. This will require a much more detailed conversation. 

I will now refer to the convergences by their indices. 

$$
1, 3 \implies 2 \implies 4
$$

Moreover, convergence in measure implies that there is a subsequence that converges almost surely (this is the version of the converse of $1\implies 2$), and convergence in $L^p$ implies convergence in $L^r$ if $p \geq r$. 

Let us state them precisely and prove them (except $4$). 

**Proposition 1:** On a probability space, $X_n \to X$ almost surely, then $X_n \to X$ in probability. The converse need not be true. 

**Proof:** Fix $\epsilon > 0$. We need to show that $\mu(|X_n - X| > \epsilon) \to 0$. Let $A_n = \{\omega: |X_m(\omega) - X(\omega)|>\epsilon\ \text{for some } m \geq n\}$

For the coounterexample, consider $\mathbb{R}$ with lebesgue measure and $X_n = I_{[n, \infty)]}$, and $X = 0$. 

**Proposition 2:** If $X_n \to X$ in measure, then there is a subsequence of $(X_n)$ which converges to $X$ a.e. 

For each $m$, $\lim_n \mu(|X_n - X| > \dfrac{1}{2^m}) = 0$. Fix $n_m$ such that the measure is less than $\dfrac{1}{2^m}$. Let $Y_i = X_{n_i}$. Now we are done by Borel Cantelli Lemma. 

**Proposition 3:**  Let $p \geq 1$. If $X_n \to X$ in $L^p$, then $X_n \to X$ in measure. 

This follows directly from Chebyshev.

$$
\mu(|X_n - X| > \epsilon) \leq \dfrac{||X_n - X||^p}{\epsilon^p} \to 0
$$
where the norm is the $L^p$ norm. 

**Propostition 4:** Let $p \geq 1$. If $X_n \to X$ in $L^p$, and $1 \leq r \leq p$, then $X_n \to X$ in $L^r$, 


4. Weak Convergence

We will discuss this now. 

### Convergence in Expectation 

We discussed three convergence modes, and will discuss one more now. But when does convergence of random variables imply that their expectations will converge as well? Turns out, this only happens in $L^p$ convergence. 

One counterexample kills all the three other cases. On $([0, 1], \lambda)$, take $X_n = n1_{(0, 1/n)}$. For each $\omega >0$, $X_n(\omega)$ is eventually $0$, hence it converges to $0$ almost surely, hence in probability as well as distribution, yet $E[X_n] = 1$ for every $n$. 

(For the last case, $L^p$ convergence implies $L^1$ convergence which is convergence of expectations).

The random variables have the property of their expectations converging with one additional required property: Uniform Integratibility. This result is called Vitali theorem. 

**Vitali Theorem for convergence in probability (hence for almost surely):** 

If $X_n \to X$ in probability (and $X_n$ are $L^1$ integrable) then 
$$
\{X_n\} \text{ uniform integrable } \iff X \in L^1 \text{ and }  X_n \to X \text{ in } L^1.
$$
Note that if we have the same statement for $|X_n|^p$, then we have $L^p$ convergence. 

**Vitali Theorem for weak convergence:** 

Let $X_n \implies X$. If $\{X_n\}$ are uniformly integrable, then $X \in L^1$ and $E_{X_n}[X_n] \to E_X[X]$. The converse holds if the means converge for $|X_n|$. 

Why did I not just repeat the same theorem as the previous case? Because here, you possibly have different sample spaces as preimages for each $X_n$. 

Let me now state the uniform integrability definition. 

**Definition:** A family $\{X_i\}$ of real valued random variables is uniformly integrable if 
$$
\lim_{M \to \infty} \sup_{i \in I} \mathbb{E}[|X_i|1_{|X_i| > M}]= 0
$$
Notice that $M$ is uniform over all $i$. This condition says that the tails uniformly carry neglible mass. Also see that this condition is on the laws alone. 

Here is an important theorem (de la Valle Poussin):

$$
\{X_n\} \text{ uniformly integrable } \iff \sup_iE[|X_i|^{1 + \epsilon}] < \infty
$$

**Definition:** A sequence $X_i$ of nonnegative random variables is said to satisfy the asymptotic expectation condition with index $k$ if there exists $\epsilon_0 > 0$ such that
$$
\sup_n \mathbb{E}[X_n^{k + \epsilon}] < \infty
$$



## Weak Convergence

Consider a probability space $(\omega, \cal{A}, P)$. 

Almost sure convergence requires knowing almost all sample points. This notion is called convergence with probability one. Strong Law of Large Numbers and Kolmogorov Three Series (not necessary for now) Theorem are results of this convergence. 

Convergence in probability does not require knowledge of the values of the variables at individual sample points. It needs the knowledge of certain probability of events. The Weak Law of Large Numbers is a result of this type of convergence. 

Convergence in $L^p$ requires the functions in $L^p$. 

In practice, we do not have minute knowledge about random variables. We only perceive expectations, or sum total of certain characteristics. In practice, we often perceive the moments of a random variable, which are meaninful only when they exist. However, if you take any bounded continuous function on $\mathbb{R}$, $E(f(X))$ is always meaningful and can be regarded as a statistic. Here is a nice result which makes this significant as a statistic. 

**Result:** If we have two random variables $X$ and $Y$, and if $E(f(X)) = E(f(Y))$ for every bounded continuous function $f: R \to R$, then $X$ and $Y$ have the same distribution. 

This leads to the following definition: 

$X_n \to X$ in distribution if $E(f(X_n)) \to E(f(X))$ for every bounded continuous function $f$. This is also called convergence in law, denoted $X_n \implies X$

---
**Remark:** The first three convergence modes, the measure $\mu$ in the preimage (sample space) matters. For the convergnce in distribtion, only the pushforwarsd of $\mu$ matter by the random variables. These are the measures on $\mathbb{R}$, also called the law. The sample space is irrelevent. That is why you can have convergence in distribution even if the random variables are wildly different, hence it is a weak form of convergence. This distinction is important when it comes to convergence, and is underemphasized generally. One can do this precisely because of the change of measure (also called the law of unconcious statistician, now you know why):
$$
\int_{\Omega}g(X)d\mu = \int_{\mathbb{R}}g d\mu_X
$$
where $\mu_X$ is the pushforward measure on $\mathbb{R}$.

---

A sequence of probabilities $(\mu_n)$ on $\mathbb{R}$ converges weakly to a probability $\mu$ if $\int f d\mu_m \to \int f d\mu$ for every bounded continunous function $f$ This is denoted $\mu_n \implies \mu$. 


We now prove the result. We first use change of measure on the condition, to get two probability distributions (which are pushforwards) $P$ and $Q$ such that 
$$
\int f dP = \int g dQ
$$
for every bounded continuous function $f$ on $R$. 

We need to show $P(B) = Q(B)$ for every borel set $B$. The proof is quite simple. Clearly, indicator function of $B$ as the candidate of $f$ is the best, unfortunately it is not continuous. So what? Just take a continuous sequence of functions converging to the indicator function and then use Dominated Convergence theorem to take limit inside. We are done. 

Alternatively, another simple proof: Notice that $\sin (x)$ and $\cos(x)$ are bounded continuous functions. Now invoke the uniqueness of characteristic functions to conclude. 

By the way, according to the first proof, we can generalize the result to restrict the set of functions to $C^\infty$ functions on a compact support.


## Exercises
This is a good time to take a break from theory. I now suggest you to go through the plethora of examples (or exercises) on this topic. They are in another section for this topic available in the notebook. 

## Theorems in Weak Convergence 

We need to develop more theory in weak convergence, since looking at each bounded continuous function seems like quite an impractical job to do. Here is a nice toolkit. 

**Theorem 1:** Let $(X_n)$ and $X$ be random variables. $X_n \sim \mu_n$ and $X \sim \mu$. Distribution functions of $\mu_n$ and $\mu$ are $F_n$ and $F$ respectively. Following statements are equivalent: 

1) $X_n \implies X$.
2) $\mu_n \implies \mu$.
3) $F_n(x) \to F(x)$ for every continuity point $x$ of $F$.
4) $F_n(x) \to F(x)$ for a dense set of points $x$ in $\mathbb{R}$. 

The second statement reduces our job by transferring the entire action to the real line. Now we only need to calculate integrals on $\mathbb{R}$. The third statement reduces our job further. We need not look at every bounded continuous function, just evaluate distribution at a lot points and show convergence of thos numbers. The final statement reduces it even furhter. We do not need evaluation at uncountably many points, just a nice dense set (such as $\mathbb{Q}$, or some other set) works. 

I will hold off on proving this in this notebook, but the reader can look at other resources which prove this. It is just heavy work in analysis. 

**Theorem 2:** The above four conditions are also equivalent to the folowing. 

5) $\int f d\mu_n \to \int f d\mu$ for all bounded $C^\infty$ functions which have bounded derivatives. 
6) $\int f d\mu_n \to \int f d\mu$ for all real bounded unifromly continuous $f$. 

**Theorem 3:** The above six conditions are also equivalent to the following: 



For the next theorem, recall: 

$$
\varphi(t) = E(e^{itX}) = E(\cos tx) + iE(\sin tx); t\in \mathbb{R}
$$
This is the characteristic function of $X$, or its distribution function $\mu$. It is also the fourier transform of $\mu$. This is a complex valued function whose domain is $\mathbb{R}$. 

**Theorem 4 (Levy's Continuity Theorem):** The above nine conditions are equivalent to the following: 

10) $\varphi_n(t) \to \varphi(t)$ for each $t \in \mathbb{R}$. 
11) $\varphi_n(t) \to \varphi(t)$ uniformly over each bounded interval. 

This theorem is pretty peak, you can use it to prove the Central Limit Theorem taught in the introductory courses. 

**Theorem 5 (Skorokhod Theorem):** The above eleven conditions are equivalent to the following: 

12) There is a sequence of random variables $Z_n$ and $Z$ on $([0, 1], B, \lambda)$ such that $Z_n \sim \mu_n$ for each $n$, $Z \sim \mu$, and $Z_n \to Z$ a.e wrt $\lambda$. 

If the above conditions hold, then $f(Z_n) \to f(Z)$ a.e wrt $\lambda$, and thus the expectations converge (by DCT). But $E[f(X_n)] = [f(Z_n)]$ and $E[f(X) = f(Z)]$, thus we are done with the forward direction. For the converse, one needs to look at the quantile function / inverse distribution function. 

**Theorem:** Weak Convergence is robust. $X_n \implies X$ and $\varphi: R \to R$ is continuous, and $Y = \varphi(X)$, then $Y_n \implies Y$. 
