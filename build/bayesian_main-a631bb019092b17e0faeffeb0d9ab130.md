# Basic Bayesian Theory

## Foreword 

We will continue where we left off. In the last post, I introduced a few important observables along with the setup. We are interested in their asymptotics, and the behaviors of these random variables. But before we dive into them, we first see some motivation to their definitions, and how they are all connected to each other (WAIC definition was pretty non intuitive to me the first time I saw it). Somewhere along the middle, I will slightly detach from Watanabe's notation to do this and approach this topic differently, before converging back to his notation. This is primarily done to make the relationships between the observables more clearly and a give a more detailed introduction. Once this is done, we will start with the asymptotics. 

How do we see the relationships between those observables though, and even try to begin their asymptotics? To help you get settled with the mindset, I will give you the spoiler beforehand: Generating Functions. 


For an arbitrary triple of a true distribution, a statistical model, and a prior, the behaviors of the free energy, the losses and WAIC are derived by the following procedure: 

1) Firstly, we define the formal relation between a true distribution and a statistical model. 

2) Secondly, definitions of Bayesian Observables and their normalized ones are introduced. 

3) Thirdly, the Cumulant Generating Functions of the Bayesian prediction is defined. 

4) The Basic Theory of Bayesian Statistics is proved using the Cumulant Generating Functions. 

## The Formal Relation 

We make 4 main definitions right now. 

> Definition 1 (Realizability): Let $W \subset \mathbb{R}^d$ be the set of all parameters. If there exists $w_0 \in W$ such that $q(x) = p(x|w_0)$ almost surely, then $q(x)$ is said to be realizable by $p(x|w)$. For a given pair of true distribuiton and the statistical model, the set of true parameters is defined by 
>$$
W_{00} = \{w \in W: q(x) = p(x|w) \text{for arbitrary }x \text{ s.t }q(x) > 0\}
$$

Thus, realizability holds iff the set of true parameters is not an empty set (one can see this by using the fact that the integral of both the functions is 1 and then using the almost sure condition). 


We can alternatively write 
$$
W_{00} = \{w \in W: \int q(x)\log\dfrac{q(x)}{p(x)}dx = 0\}
$$

The average loss function is defined by 
$$
L(w) = -\int q(x) \log p(x|w) dx
$$
It follows that 
$$
L(w) = S + K(q(x)||p(x|w))
$$
Hence, if $q(x)$ is realizable by a statistical model, then the average log loss function is minimized iff $w \in W_{00}$, and the minimum values is the entropy of the true distribution. 

Some motivation: As we continue the learning process, our statistical model gets closer to the true distribution. The realizability assumption makes it simple for us and tells us that our statistical model can represent the true distribution, which need not be the case in general. 

> Definition 2 (Regularity): For a given pair of $q(x)$ and $p(x|w)$, let $W_0 = \{w \in W: L(w) = \min_{w'}L(w')\}$, which is called the set of optimal parameters for the minimum average logg loss . $q(x)$ is said to be regular for $p(x|w)$ if the following three conditions hold: 
>
> 1) $W_0$ is a singleton set. 
> 2) $w_0 \in W_0$ is in the interior of $W$ ($W$ is equipped with the subspace topology from $\mathbb{R^n}$). 
> 3) The Hessian Matrix of the average log loss function at $w_0$ is positive definite. 

Let us see a simple consequence of these definitions. If $W$ is a compact set and if $L(w)$ is a continuous function, then $L(w)$ has a minimum point, hence $W_0$ is not an empty set. In this case, $q(x)$ is realizable by $p(x|w)$ iff $W_{00} = W_0$ (reason this by yourself). 

Some Motivation: The regularity condition tells us that we can approximate the the loss landscape near its minima by a parabola. This follows by doing a taylor approximation, doing a singular value decomposition on the hessian, and then using the fact that all the eigenvalues are positive (do the exact calculations yourself). One of the main claims is that this is a terrible approximation in many cases and do not hold for the loss landscapes of the general neural networks. 

> Definition 3 (Essential Uniqueness): Assume that $W_0$ is not an empty set. If there exists a unique probability density function $p_0(x)$ such that for arbitrary $w_0 \in W_0$,
$p(x|w_0) = p_0(x)$, then it is said that the optimal probability density function is essentially unique. 

If $q(x)$ is realizable by $p(x|w)$, then the optimal probability density function is essentially unique, because $p_0(x) = q(x)$. Thus, realizability $\implies$ Essential Uniqueness. 

> Definition 4 (Relarively Finite Variance of Log Density Ratio Function): For a given pair $w_0 \in W_0$ and $w \in W$, the log density ratio function is defined by 
>
>$$
f(x, w_0, w) = \log \dfrac{p(x|w_0)}{p(x|w)}
$$
>
>If there exists $c_0 > 0$ such that for an arbitrary pair $w_0, w$, 
>
>$$
\mathbb{E}_X[f(x, w_0, w)] \geq c_0 \mathbb{E}_X[f(x, w_0, w)^2]
$$
>
> then the function $f(x, w_0, w)$ has a relatively finite variance. 

Remark: The function $f$ has a relatively finite variance iff 
$$
\sup_{w \notin W_0} \dfrac{E_X[f^2]}{E_X[f]} < \infty
$$

Some motivation: What is the purpose of the theory that we are developing? We want to study the model's distance from the truth, and how that develops over time. We want to abstract out the truth's randomness from our process. Let us see what this condition is really about. 

Let us look at the average log loss function (do you prefer to look at the asymptotics of a sum or a product? Take log.) whose definition is derived naturally from the average loss function. 
$$
\mathbb{E}_X[\log p(X|w)] = \mathbb{E}_X[\log q(x)] - \mathbb{E}_X[\log \dfrac{q(x)}{p(x|w)}]
$$
$$
 = -S - K(w)
$$
Removing the entropy considerations from $S$, we get the "excess log loss over the truth per data point", which is what the definition is in the realizable case (but look at the next lemma to understand the point in the general case). Thus the importance. This remark also says the mathematical statement to note, that the expectation of the log density ratio is the KL distance between the true distribution and the posterior in the realizable case. I have explained the ratio, I will explain the condition in a short while. 

Lemma: Assume that $w_0 \in W_0$ and $w \in W$. If $f(x, w_0, w)$ has a relatively finite variance, then the optimal probability density is essentially unique. 

Assume that $w_1$ and $w_2$ are arbitrary elements of $W_0$. By the definition of $W_0$, $0 = L(w_1) - L(w_2) = \int q(x)f(x, w_1, w_2)dx \geq c_0 \int q(x)f(x, w_1, w_2) dx$ and as the integrand is non negative, then $p(x|w_1) = p(x|w_2)$ almost surely. Do observe that this is not a clever proof, rather just making use of all the given information. 

Thus, we may denote $f(x, w_0, w) = f(x, w)$. It thus follows that $p(x|w) = p_0(x)\exp(-f(x, w))$. 

Thus, relatively finite variance $\implies$ essential uniqueness. However, the converse need not be true. It does hold with a few added conditions. 

Lemma: Assume that $W$ is a compact set and that $q(x)$ is realizable by $p(x|w)$ and that the log desity ratio function $f(x, w) = \log \dfrac{q(x)}{p(x|w)}$ is a continuous function of $(x, w)$. If there exists $c_1, c_2 > 0$ such that for an arbitrary $w \in W$, 

$$
\int_{|x| > c_1} q(x)f(x, w)^2dx \leq c_2\int_{|x| \leq c_1} q(x)f(x, w)dx
$$
then $f(x, w)$ has a relatively finite variance. 

I refer to the reader the book for the proof. It is not important to the rest of the theory, just some technical details. 

Lemma: Assume that $W$ is a compact set and that for an arbitrary pair $w_0 \in W_0$ and $w \in W$, the second derivatives of $E_X[f(x, w_0, w)]$ and $E_X[f(x, w_0, w)^2]$ are continuous functions. If $q(x)$ is regular for $p(x|w)$, then the log density ratio function has a relatively finite variance. 

The proof of this lemma is right there with the proof to the previous lemma. 

Assume that $W$ is compact. We thus have the following relations: 
1) $\{$Regular$\} \subset \{$ Relatiely Finite Variance $\}$. 
2) $\{$Realizable$\} \subset \{$ Relatively Finite Variance $\}$. 
3) $\{$Relatively Finite Variance $\} \subset \{$ Essentially unique $\}$. 

# Normalized Variables

We have the average and empirical loss functions:

$$L(w) = -E_X [\log p(X|w)]$$

$$L_n(w) = -\frac{1}{n} \sum_{i=1}^n \log p(X_i|w)$$

Now, $f(x, w_0, w) = \log \frac{p(x|w_0)}{p(x|w)}$ has a relatively finite variance.

Now, since $f(x, w_0, w)$ has a relatively finite variance, it is essentially unique, we thus denote it by $f(x, w)$.

See that $L(w) = E_X \left[ \log \frac{1}{p(X|w)} \right] = E_X \left[ \log \frac{p(X|w_0)}{p(X|w)} \right] - E_X [\log p(X|w_0)]$

$$= E_X [f(x, w_0, w)] - E_X [\log p(X|w_0)]$$

$$= E_X [f(x, w_0, w)] + L(w_0)$$

We define the normalized average loss function as $E_X[f(x, w_0, w)] = K(w)$
A key thing to note in all the definitions is that we will be rewriting them in terms of $f(x, w)$. The formal advantage is explained in a short while, but the understanding is that we want to remove the entropy of the unique probability density and consider only the KL divergence between $p_0(X)$ and $p(X|w)$ (Posterior), as that is what we want to analyse (how that changes over time and what are the asymptotics).

(For this to make sense, note that $E_X [f(x, w_0, w)] = K(q(x) || p(X|w)) = K(w)$)

Similarly, we have the normalized empirical log loss function
$$K_n(w) = \frac{1}{n} \sum_{i=1}^n f(X_i, w)$$

---

Furthermore, a key property of our normalized fn. is that
$K(w) \ge 0$, and $K(w) = 0 \iff w \in W_0$.

The partition function is defined as $Z_n = \int \prod_{i=1}^n p(X_i|w) \varphi(w) dw$

Thus, $Z_n = \int \exp \left( \sum_{i=1}^n \log p(X_i|w) \right) \varphi(w) dw$

$$= \int \exp \left( - \sum_{i=1}^n \log \frac{1}{p(X_i|w)} \right) \varphi(w) dw$$

Thus, $\frac{Z_n}{\prod_{i=1}^n p_0(X_i)} = \int \exp \left( - \sum_{i=1}^n \log \frac{1}{p(X_i|w)} \right) \cdot \exp \left( - \sum_{i=1}^n \log p_0(X_i) \right) \varphi(w) dw$

$$\Rightarrow \frac{Z_n}{\prod_{i=1}^n p_0(X_i)} = \int \exp \left( -n \cdot \frac{1}{n} \left( \sum_{i=1}^n \log \frac{p_0(X_i)}{p(X_i|w)} \right) \right) \varphi(w) dw$$

$$= \int \exp ( -n K_n(w) ) \varphi(w) dw$$

We define this as the normalized marginal likelihood $Z_n^{(0)}$.

Thus, $Z_n = \prod_{i=1}^n p_0(X) Z_n^{(0)} = \exp(-n L_n(w_0)) Z_n^{(0)}$.

Let us recall that the posterior is

$$p(w|X^n) = \frac{1}{Z_n} \prod_{i=1}^n p(X_i|w) \varphi(w)$$

$$= \frac{\exp(-n L_n(w)) \cdot \varphi(w)}{\exp(-n L_n(w_0)) Z_n^{(0)}}$$

$$= \frac{\exp(-n K_n(w)) \cdot \varphi(w)}{Z_n^{(0)}}$$

---

Free energy is the negative log of the partition function. The normalized free energy is defined similarly.

$$F_n^{(0)} = -\log \int \exp(-n K_n(w)) \varphi(w) dw$$

Here is the table for the generalization, cross validation, and training losses and the WAIC, and also their normalisation.

| | | Normalisation |
| :--- | :--- | :--- |
| Generalization loss: | $-E_X [\log p(X|X^n)]$ | $-E_X [\log E_w [\exp(-f(X,w))]]$ |
| Training loss: | $-\frac{1}{n} \sum_{i=1}^n \log p(X_i|X^n)$ | $-\frac{1}{n} \sum_{i=1}^n \log E_w [\exp(-f(X_i,w))]$ |
| Cross Validation loss: | $\frac{1}{n} \sum_{i=1}^n \log E_w \left[ \frac{1}{p(X_i|w)} \right]$ | $\frac{1}{n} \sum_{i=1}^n \log E_w [\exp(f(X_i,w))]$ |
| WAIC: | $T_n + \frac{1}{n} \sum_{i=1}^n V_w [\log p(X_i|w)]$ | $T_n^{(0)} + \frac{1}{n} \sum_{i=1}^n V_w [f(X_i,w)]$ |

Let us connect all these observables.

---

Defn: For $\alpha \neq 0$, we define the $\alpha$-predictive at $x$ as

$$p^{(\alpha)}(x) = \left( E_w \left[ p(x|w)^\alpha \right] \right)^{1/\alpha}$$

Note that this is the $f$-power mean of $p(x|w)$ with respect to the predictive posterior distribution.
Essentially, the observables come from this family.

Recall that the moment generating function of a random variable $X$ with PDF $p(x)$ is given by $M_X(\alpha) = E[e^{\alpha X}]$
and the cumulant generating function is $G_X(\alpha) = \log E[e^{\alpha X}]$

We will look at the cumulant generating function of the log likelihood at $x$ (with a little notational change) wrt the posterior, which is exactly the log of the ($\alpha$-predictive)$^\alpha$.

$$g_x(\alpha) = \log E_w [p(x|w)^\alpha] = \log E_w [e^{\alpha \log p(x|w)}]$$

We can now write the observables in terms of the CGF.

$$G_n = -E_X [\log E_w [p(x|w)]] = -E_X [g_x(1)]$$
$$T_n = -\frac{1}{n} \sum_{i=1}^n g_{X_i}(1)$$
$$C_n = \frac{1}{n} \sum_{i=1}^n \log \frac{1}{E_w(p(X_i|w)^{-1})} = \frac{1}{n} \sum_{i=1}^n g_{X_i}(-1)$$

Thus they are all derived from the CGF. We will discuss WAIC shortly, and also note how it comes about.

---

Let us discuss the normalization first.

We define $F_X(\alpha) = \log E_w [e^{\alpha f(x,w)}]$ in the normalized coordinates.

The corresponding definitions in terms of $F_X(\alpha)$ give the normalized variables.

Thm: We also get $F_X(\alpha) = g_x(-\alpha) + \alpha \log q(x)$.

The point of introducing the CGF is to be able to do a Taylor expansion here.

Doing a Taylor approximation near 0,

$$g(\alpha) = g(0) + \alpha g'(0) + \frac{\alpha^2}{2!} g''(0) + \dots$$

where $g'(0) = E_w [\log p(x|w)] = C_1$
$g''(0) = V_w [\log p(x|w)] = C_2$

So, $g_x(1) + g_x(-1) = C_2 + \frac{g^{(4)}(\theta)}{24} + \frac{g^{(4)}(-\theta)}{24}$

Summing over the training points, we get

That is, $T_n + C_n + \frac{1}{n} \sum_{i=1}^n C_2(X_i) = \dots$

$$-T_n + C_n = \frac{1}{n} \sum_{i=1}^n V_w [\log p(X_i|w)] + R_n$$
$$\Rightarrow C_n = T_n + \frac{1}{n} \sum_{i=1}^n V_w [\log p(X_i|w)] + R_n$$

This gives rise to WAIC, and using the normalized CGF gives rise to the normalized WAIC.

---

Caveat: Watanabe does an abuse of notation (harmless) and defines it the following way. We will continue with this.

$$G_n(\alpha) = E_X \left[ \log E_w [ e^{\alpha \log p(x|w)} ] \right]$$
$$T_n(\alpha) = \frac{1}{n} \sum_{i=1}^n \log E_w [ e^{\alpha \log p(X_i|w)} ]$$

Thus, his definition is just the mean of our definitions. However, do note that the two means are different, since $E_X$ and $E_w$ do not commute.

I will now state the result which explains why we need to look at the normalised observables.

$$F_n/n = L_n(w_0) + O_p(\log n / n)$$
$$G_n = L(w_0) + O_p(1/n)$$
$$C_n = L_n(w_0) + O_p(1/n)$$
$$T_n = L_n(w_0) + O_p(1/n)$$
$$W_n = L_n(w_0) + O_p(1/n)$$

All these asymptotics abstract away from the prior. In the realizable case, they abstract away from the statistical model as well.
We need to clarify the behavior of these random variables, look at higher order terms. The normalised versions of these all converge to 0 in probability, removing the effects of $L_n(w_0)$ and $L(w_0)$ and focusing on what is necessary.

---

We will now explain the asymptotics. We will be considering the normalized observables (and drop the (0) superscript):

Theorem: Let $G'(0), T'(0), G''(0), T''(0)$ be random variables, and assume that

$$\sup_{|\alpha| \le 1} \left| \left( \frac{d}{d\alpha} \right)^3 G_n(\alpha) \right| = O_p\left(\frac{1}{n}\right)$$

$$\sup_{|\alpha| \le 1} \left| \left( \frac{d}{d\alpha} \right)^3 T_n(\alpha) \right| = O_p\left(\frac{1}{n}\right)$$

(That is inside radius 1, the third derivatives of the taylor expansion are bounded).

Then
$$G_n = -G_n(1) = -G_n'(0) - \frac{1}{2} G_n''(0) + O_p(1/n)$$
$$T_n = -T_n(1) = -T_n'(0) - \frac{1}{2} T_n''(0) + O_p(1/n)$$
$$C_n = T_n(-1) = -T_n'(0) + \frac{1}{2} T_n''(0) + O_p(1/n)$$
$$W_n = -T_n(1) + T_n''(0) = -T_n'(0) + \frac{1}{2} T_n''(0) + O_p(1/n)$$

Proof: By the mean value theorem, given $\alpha$, there exists $\alpha^*$ s.t $|\alpha^*| \le |\alpha|$ and that
$$G_n(\alpha) = G_n(0) + \alpha G_n'(0) + \frac{1}{2} \alpha^2 G_n''(0) + \frac{1}{6} \alpha^3 G_n'''(\alpha^*)$$
$\alpha = 1$ gives the theorem (similar method for dealing with $T_n$)

---

We thus get the theoretical behaviors of the free energy, generalization loss, and the rest by the following procedure.

Recipe For Bayesian Theory Construction

1. An arbitrary triple of $(q(x), p(x|w), \varphi(w))$ is fixed. The set of parameters is denoted $W$, and $X^n \sim q(X)$.

2. The empirical and average loss functions are defined by
$$L_n(w) = -\frac{1}{n} \sum_{i=1}^n \log p(X_i|w)$$
$$L(w) = -\int q(x) \log p(x|w) dw$$

Find the optimal parameters minimizing $L(w)$.
$$W_0 = \{ w \in W : L(w) = \min_{w' \in W} L(w') \}$$

3. Check that the log density ratio fn. has a relatively finite variance. Then we have essential uniqueness.

4. Define
$$K(w) = \int f(x,w) q(x) dx$$
$$K_n(w) = \frac{1}{n} \sum_{i=1}^n f(X_i,w)$$

The normalised partition function is given by
$$Z_n^{(0)} = \int \exp(-n K_n(w)) \varphi(w) dw$$

Then the free energy $F_n = n L_n(w_0) - \log Z_n^{(0)}$.

---

5. The average by the posterior is equal to

$$E_w[ \cdot ] = \frac{\int (\cdot) \exp(-n K_n(w)) \varphi(w) dw}{\int \exp(-n K_n(w)) \varphi(w) dw}$$

Calculate $E_w[f(x,w)]$, $V_w[f(x,w)]$

Then $E_w[K(w)] = E_X E_w[f(X,w)]$
$$E_w[K_n(w)] = \frac{1}{n} \sum_{i=1}^n E_w[f(X_i,w)]$$
$$E_X V_w[f(X,w)] = \frac{1}{n} \sum_{i=1}^n V_w[f(X_i,w)] \text{ are obtained.}$$

Finally, based on the basic theorem,

$$G_n = L(w_0) + E_w[K(w)] - \frac{1}{2} E_X V_w[f(X,w)] + O_p(1/n)$$
$$C_n = L_n(w_0) + E_w[K_n(w)] + \frac{1}{2n} \sum_{i=1}^n V_w[f(X_i,w)] + O_p(1/n)$$
$$T_n = L_n(w_0) + E_w[K_n(w)] - \frac{1}{2n} \sum_{i=1}^n V_w[f(X_i,w)] + O_p(1/n)$$

and $W_n$ having the same expansion of $C_n$
(which by the way proves that $C_n = W_n + O_p(1/n)$)

## Further Steps 

We have already setup the framework and the basic theory, and also discovered important results (which were not restricted to regular models). We are now ready to find the behavior of the regular posterior distribution, and then generalize that behavior to
the general case. 


