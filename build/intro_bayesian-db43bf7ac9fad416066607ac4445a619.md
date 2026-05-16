We start with Bayesian Statistics. Watanabe's theory is fundamentally based on generalizing classical results in Bayesian Statistics, so it is important to get a strong grip and understand this classical theory well before moving on. It also gives us the complete understanding of the framework we are working in, and is the first essential thing to master. 

## Connection with Machine Learning and Setup 

Machine Learning Models are primarliy consisting of two frameworks (or a combination of them): Frequentist and Bayesian.

The setup is that we have a true data generating distribution $p_{data}(x)$, and consider a set of arbitrary samples $X = \{x_1, ..., x_n\}$ . We take a statistical model $p_{model}(x, \theta)$ (which is a parametric family of probability distributions) which aims to estimate the true distribution. 

The likelihood function of our statistical model is defined as 

$$
P_{model}(X; \theta) = \prod_{i = 1}^n p_{model}(x_i; \theta)
$$

The frequentist approach is to find the optimal $\theta$ which maximizes this likelihood function. 

The KL divergence from probability distribution $p$ to $q$ is defined as 
$$
D(p||q) = \int_{-\infty}^{\infty} p(x) \log \dfrac{p(x)}{q(x)}dx
$$

This is the main measure that we will use to associate similarity between probability distributions (even though it is not really a metric, it is clear that it is not even symmetric). 

It can be easily seen that finding the optimal $\theta$ (called the maximum likelihood estimator) is equivalent to minimizing the KL divergence from the empirical true distribution to our statistical model, which is a function of $\theta$. An approximation to the local optimal parameter is often approached via (stochastic) gradient descent. This is also the case in neural networks, which are essentially function approximators. We use SGD to approximate to the local optimal parameter vector. 

We will not delve into the frequentist approach more here (you may refer to Goodfellow et al). We will move on to the Bayesian approach here. Thus, when we refer to neural networks here, an important distinction is that now this is not the standard neural networks where SGD is used. Still, we gain many insights from this approach that also carry to the standard networks. 

In the Bayesian approach, instead of considering just the optimal parameter, we consider a probability distributin over the space of parameters itself. Initially, this is called the prior function, and as we observe the data from the true distribution, we update this prior function to successively obtain a posterior functin, which is an estimate over the entire parameter space to what generates the true distribution function. 

Specifically, we consider an appropriate prior function $\varphi(w)$ and a statistical model $p(x|w)$. These are chosen by us, and this choice often determines what estimate our bayesian method will given us. We assume there is a true date generating distribution $q(x)$, from which we draw $N$ samples *independently*, $\{x_1, ..., x_n\}$. This sample induces a function $p(w|x_1, ..., x_n)$, which is the update of our prior function. This further induces $p(x|x_1, ..., x_n)$, which is our estimate of the true distribution. This process goes on as we make more samples. As can be seen, this is more computatinally intensive. However, this approach is superior in many cases, we will specifically see an example later on. We will now define everything mathematically. To summarie, here is the procedure: 

1) Construct the universe and the mathematical laws between bayesian observables which hold for any arbitrary: true distribution, statistical model, and a prior. 

2) Evaluate how appropriate the statistical model and the prior is using these laws.

3) Employ the most suitable pair.

## Introduction to Bayesian Statistics 

The posterior function is obtained through Bayes' rule. 

$$
p(w|x_1, ..., x_n) = \dfrac{p(x_1, ..., x_n|w) \varphi(w)}{\int p(x_1, ..., x_n|w)\varphi(w)dw}
$$

But neither do we know the statistical model, nor do we know the prior. Thus a meaningful approach is to just start with something, evaluate how good it is, and then update it. The evaluation is done through the mathematical laws described above.

This gives rise to the estimated pdf of $x$, called the predictive distribution: 
$$
\hat{p}(x)  =\int p(x|w) P(w|x_1, ..., x_n) dw
$$

Expected: $\hat{p}(x) \approx q(x)$ if $(p(x|w), \varphi(x))$ is appropriate for $q(x)$. We want to evaluate the tuple appropriateness without information about $q(x)$. We develop the machinery for that. 

## True Distribution

A realized value of $X^n$ in a trial is denoted $x^n = (x_1, \ldots, x_n)$. In practical applications, while we may not know $q(x)$, we assume its existence.

Let us just revise the basics first as they will be important in the calculations that we make. 

Let $f : X^n \rightarrow f(X^n) \in \mathbb{R}$

$$E[f(X^n)] = \int \cdots \int f(x^n) \prod_{i=1}^n q(x_i) dx_i$$

Do observe that we are able to take the product here because of the independent sampling. 

$$\mathbb{V}[f(X^n)] = E[f(X^n)^2] - E[f(X^n)]^2$$

The average entropy of the true distribution is defined as:
$$S = - \int q(x) \log q(x) dx$$

The empirical entropy is defined as:
$$S_n = -\frac{1}{n} \sum_{i=1}^n \log q(X_i)$$

By definition, one can see that $E[S_n] = S$. I outline it nonetheless, so that you can get comfortable with the calculations. 

$$E[S_n] = -\frac{1}{n} \sum_{i=1}^n E[\log(q(X_i))]$$
$$= -\frac{1}{n} \sum_{i=1}^n \int \log q(x) \cdot q(x) dx$$
$$= -\int \log q(x) \cdot q(x) dx = S$$

Similarly, one can see that the variance of the empirical entropy is:
$$\mathbb{V}[S_n] = \frac{1}{n} \left[ \int q(x)(\log q(x))^2 dx - S^2 \right]$$

The average and empirical entropies of the true distribution which is a conditional distribution is defined similarly:

$$S = - \int q(x,y) \log q(y|x) dx dy$$
$$S_n = - \frac{1}{n} \sum_{i=1}^n \log q(y_i|x_i)$$

## Model, Prior and Posterior

Let $W \subset \mathbb{R}^d$. Let $X^n$ be independent real random values subjected to $q(x)$. For an arbitrary pair $(p(x|w), \varphi(w))$, the posterior probability density is defined by

$$p(w|X^n) = \frac{1}{Z(X^n)} \varphi(w) \prod_{i=1}^n p(X_i|w)$$

where $Z(X^n)$ is defined by 
$$
Z(X^n) = \int \varphi(w) \prod_{i=1}^n p(X_i|w) dw
$$

which is called the partition function/marginal likelihood/evidence. 

Expected value over the posterior distribution is denoted $E_w[\cdot]$.
Do note that $E_w[f(w)] = \frac{1}{Z(X^n)} \int f(w) \varphi(w) \prod_{i=1}^n p(X_i|w) dw$

This expected value is a random variable as it depends on $X^n$.
(Better to say, it is the expected value over a conditional probability density and hence is a random variable)

The posterior gives rise to the predictive density function:
$p(x|X^n) \stackrel{def}{=}$
$$= E_w[p(x|w)] = \int p(x|w) p(w|X^n) dw$$

(estimate $w$ from $X^n$, estimate $x$ from $w$, vary over all $w$)

If $\int \varphi(w) dw < \infty$, the prior is called proper, because it is normalized so that $\int \varphi(w) dw = 1$. Even for an improper prior, posterior and predictive probability densities can be defined if $Z(X^n)$ is finite and well defined.

## An Important Example - The Exponential Family

In many simple statistical models, the posterior converges to the normal distribution as $n \to \infty$. We see such a case in the example referred to below. However, even in some simple cases, and many others, this fails. This is the key problem resulting in the new theory. 

At this point, I highly recommend referring to this example: 

We are now going to prove the formulae given in the example. 

If the statistical model is of the form
$$
p(x|\theta) =  u(x)\exp (v(\theta)^T w(x))
$$
where $u$ is a real valued function (and the other two are vector valued), then this distribution is said to belong to the exponential family. Furthermore, if the distribution of the parameter $\theta \in \Theta$ depends on some hyperparameter $\phi$, and can be written as
$$
\varphi(\theta|\phi) = \frac{\exp (v(\theta)^T \phi)}{z(\phi)}
$$
where $z(\phi)$ is the normalizing factor, then $\varphi(\theta|\phi)$ is said to be a conjugate prior distribution. 

In the case when the distribution is of the form
$$
p(x|\mu, \sigma^2) = \frac{1}{\sqrt{2\pi\sigma^2}}\exp \{-{\frac{(x-\mu)^2}{2\sigma^2}}\}
$$, 

we can take $u(x) = \dfrac{1}{\sqrt{2\pi}}$, $v(\theta)^T = [\dfrac{1}{\sigma^2}, \dfrac{\mu}{\sigma^2}, \dfrac{\mu^2}{\sigma^2}, \log\sigma^2]$

and $w(x)^T = [\dfrac{-x^2}{2}, x, \dfrac{-1}{2}, \dfrac{-1}{2}]$. Thus it is of the exponential family. 

Now, as we know $p(\theta|x_1, ..., x_n) = \dfrac{\prod_{i=1}^np(x_i|\theta)\varphi(\theta|\phi)}{Z(x_1, ..., x_n)}$. Let us calculate the numerator.
$$
\prod_{i=1}^np(x_i|\theta)\varphi(\theta|\phi) = \frac{\exp (v(\theta)^T\phi)}{z(\phi)}\prod^n_{i=1}u(x_i)\exp (v(\theta)^Tw(x_i))
$$

$$
=   \exp(v(\theta)^T\phi)\cdot \exp(v(\theta)^T\sum_{i = 1}^{n}w(x_i))\cdot \frac{1}{z(\phi)}\prod_{i = 1}^nu(x_i)
$$

$$
=\exp(v(\theta)^T(\phi + \sum_{i=1}^{n}w(x_i)))\frac{1}{z(\phi)}\prod_{i=1}^{n}u(x_i)
$$

Let us denote $\phi_n = \phi + \sum_{i=1}^{n}w(x_i)$. Then we have that the numerator is: 
$$
= \frac{\exp(v(\theta)^T\phi_n)}{z(\phi_n)}\cdot\frac{z(\phi_n)}{z(\phi)}\prod_{i=1}^{n}u(x_i)
$$

Let us get the $Z(x_1, ..., x_n)$ for which we need to integrate the numerator with respect to $\theta$. and here we use a nice hack. We know that the integral of the prior with respect to $\theta$ is 1, regardless of what $\phi$ is. So set $\phi = \phi_n$ and we see that the first term integrates out to 1, while the second term is a scalar number independent of $\theta$. 

Hence, 
$$
Z(x_1, ..., x_n) = \frac{z(\phi_n)}{z(\phi)}\prod_{i=1}^{n}u(x_i)
$$
and thus, the posterior is 
$$
p(\theta|x_1, ..., x_n) = \frac{\exp(v(\theta)^T\phi_n)}{z(\phi_n)}
$$
which is also from the exponential family!

Finally, the predictive probability is given by 
$$
p(x|x_1, .., x_n) = \frac{Z(x_1, ..., x_n, x)}{Z(x_1, ..., x_n)} = \frac{u(x)z(\phi_n + w(x))}{z(\phi_n)}
$$
One may notice that we are using a different formula for the predictive density, bypassing the integral definition. This comes directly from using the bayes rule in the given definition (check it yourself), and it is computationally more useful in some cases to use this instead. 

For the example given at the start of the section, it is just a matter of inputting numbers into the formulae. 

## Estimation and Generalization 

We need an objective measure which indicates the difference between true and estimated probability density to evaluate how accurate the predictive density is. 

Let $X^n$ be a sample taken independently from $q(x)$ and $p(x|X^n)$ be a predictive density using a statistical model $p(x|w)$ and a prior $\varphi(w)$. We are going to make two definitions:

$$T_n = -\frac{1}{n} \sum_{i=1}^n \log p(X_i | X^n)$$

$$G_n = - \int q(x) \log p(x | X^n) dx$$

Notice how both of these quantities are random variables. 

Thus 
$$G_n - S = - \int q(x) \log p(x|X^n) dx + \int q(x) \log q(x) dx$$
$$= \int q(x) \log \frac{q(x)}{p(x|X^n)} dx$$
$$= K(q(x) || p(x | X^n))$$

Thus $G_n \geq S$, with equality iff $q(x) = p(x | X^n)$.
That is, the smaller $G_n$ is, the more precise our estimate is according to the KL divergence. 

$G_n - S$ and $T_n - S_n$ are called generalization errors and training errors resp.

An observation: As entropy does not depend on either a model / prior, smaller generalization error is equivalent to lower KL divergence. 

Definition: Assume $n \geq 2$. Let $X^n \setminus X_i$ be a set of random variables (leave one out).

Cross validation loss is defined by $C_n = -\frac{1}{n} \sum_{i=1}^n \log p(X_i | X^n \setminus X_i)$ and $C_n - S_n$ is called the cross validation error.

We now prove an important theorem, which has three statements regarding the definitions that we made. 

> Theorem: Assume that $X^n$ is independent. Then the following holds.
>
> 1) Assume that $E[G_n], E[C_n]$ are finite values. Then 
>
>$$
>E[C_n] = E[G_{n-1}]
>$$
>
>2) The cross validation loss satisfies the following.
>
>$$C_n = \frac{1}{n} \sum_{i=1}^n \log E_w \left[ \frac{1}{p(X_i|w)} \right]$$
>
>3) For an arbitrary set of $X^n$, $C_n \geq T_n$. $C_n = T_n$ iff $p(X_i|w)$ is a const function of $w$ on $\{w \in W, p(w|X^n) > 0\}$.

\
*Note: $E_w[f(w)] = \int f(w) p(w|X^n) dw$ this is just the integral with respect to the posterior distribution*

(1)
Here is the proof of the first statement.
$$E[C_n] = -\frac{1}{n} \sum_{i=1}^n E \left[ \log p(X_i | X^n \setminus X_i) \right]$$

$$
= -\frac{1}{n} \sum_{i=1}^n \int q(x) \log p(x | X^n \setminus X_i) dx
$$

$$
= -\frac{1}{n} \int q(x) \sum_{i=1}^n \log p(x | X^n \setminus X_i) dx
$$

$$
= -\frac{1}{n} \sum_{i=1}^n \underset{X^n \setminus X_i}{E} \left[ \underset{X_i}{\underset{\uparrow \text{Fubini}}{E}} [\log p(X_i | X^n \setminus X_i)] \right]
$$

$$
= -\frac{1}{n} \sum_{i=1}^n E_{X^n \setminus X_i} \left[ \int q(x) \log p(x | X^n \setminus X_i) dx \right]
$$

$$
= -\frac{1}{n} \sum_{i=1}^n E_{X^n \setminus X_i} [ - G_{n-1}(X^n \setminus X_i) ] = E_{X^{n-1}} [G_{n-1}]
$$

While it was not mentioned what the expectation is being taken over in the statement, the proof clarifies it. In any case, the answer to the clarification is the canonical and the most standard answer. 

(2) We now prove the second statement:

$$C_n = -\frac{1}{n} \sum_{i=1}^n \log p(X_i | X^n \setminus X_i)$$

Thus, 
$$C_n = -\frac{1}{n} \sum_{i=1}^n \log \int p(X_i | w) p(w | X^n \setminus X_i) dw$$

$$= -\frac{1}{n} \sum_{i=1}^n \log \frac{\int p(X_i | w) \prod_{j \neq i} p(X_j | w) \cdot \varphi(w) dw}{Z(X^n \setminus X_i)}$$

$$= -\frac{1}{n} \sum_{i=1}^n \log \left( \frac{\int \prod_{i=1}^n p(X_i | w) \cdot \varphi(w) dw}{\int \prod_{j \neq i} p(X_j | w) \cdot \varphi(w) dw} \right)$$

$$= \frac{1}{n} \sum_{i=1}^n \log \left( \frac{\int \prod_{j \neq i} p(X_j | w) \cdot \varphi(w) dw}{\int \prod_{i=1}^n p(X_i | w) \cdot \varphi(w) dw} \right)$$

Call the integrand in the denominator $A$.

$$= \frac{1}{n} \sum_{i=1}^n \log \left( \frac{\int \frac{A}{p(X_i|w)} dw}{\int A dw} \right)$$

$$= \frac{1}{n} \sum_{i=1}^n \log E_w \left[ \frac{1}{p(X_i|w)} \right]$$

\
We introduced cross validation as a measure to evaluate the accuracy of our estimation. However, rhere are two issues with cross validation:

(1) Although the averages of $C_n, G_n$ are equal the variances need not be equal. However, we do have this relation:
$$
\sigma(G_n - S) = \sigma(C_n - S_n) + O(\frac{1}{n})
$$.

(2) In the second statement, if the average by the posterior is numerically approximated, then
$$ISCV = \frac{1}{n} \sum_{i=1}^n \log \hat{E}_w \left[ \frac{1}{p(X_i|w)} \right]$$
is called the importance sampling cross validation loss.

$$CV = -\frac{1}{n} \sum_{i=1}^n \log E_w^{(-i)} \left[ p(X_i|w) \right]$$
is fundamentally different from the former ($E_w^{(-i)}$ is expectation with respect to $X^n \setminus X_i$).

(3) Let us prove the third statement now. 
$$C_n - T_n = \frac{1}{n} \sum_{i=1}^n \log \left[ E_w[p(X_i|w)] E_w \left[ 1/p(X_i|w) \right] \right] \geq 0$$
By Cauchy Schwarz. Equality holds iff $p(X_i|w)^{1/2} \propto p(X_i|w)^{-1/2}$ as a funciton of $w$ $\Rightarrow p(X_i|w)$ is a const fn of $w$.

---
We introduce another measure now, and it is often better than the cross validaion loss. There are also many cases where WAIC can be employed whereas cross validation cannot. 

Def: Assume $n \geq 1$. Let $X^n$ be a set of random variables. The widely applicable information criterion (WAIC) is defined by
$$W_n = T_n + \frac{1}{n} \sum_{i=1}^n \mathbb{V}_w \left[ \log p(X_i|w) \right]$$
$W_n - S_n$ is called the WAIC error.

Result: If $X^n$ is independent, WAIC is asymptotically equal to cross validation loss.
$$W_n = C_n + \mathcal{O}(1/n^2)$$
$$E[W_n] = E[C_n] + \mathcal{O}(1/n^2)$$

Remark: $C_n$ and WAIC can be employed to evaluate a stats model and a prior even if the prior is improper.

Just to summarize, we have introduced three instruments of measure: 

1) Generalization Error: $G_n - S$
2) Cross Validation Error: $C_n - S_n$
3) WAIC Error: $W_n - S_n$

In numerical experiments, we often care about minimizing errors instead of the loss itself due to the lower variance. 


## Marginal likelihood or Partition Fn


If a prior satisfies $\int \varphi(w) dw = 1$, then the marginal likelihood (partition function) satisfies

$$\int Z(x_1, \dots, x_n) dx_1 \dots dx_n = \int \varphi(w) dw \int \prod p(x_i | w)dx_1 \dots dx_n  = 1$$
We have slyly used Fubini Theorem above. 
Thus $E_{X^n}[Z(X^n)] = 1$ if the prior is proper.

$Z(X^n)$ can be thus be understood as an estimated pdf of $X^n$ using a statistical model $p(x|w)$ and a prior $\varphi(w)$. Thus it is sometimes written as $p(X^n)$. Thus this is an estimate of $q(X^n)$ by looking at the KL divergence. 

Definition: The Free energy, or the minus log marginal likelihood is defined by
$$F_n = - \log Z(X^n)$$

We look at this quantity also as an estimate of $q(X^n)$.

Using the notation $q(X^n) = \prod q(x_i)$ and $p(X^n) = Z(X^n)$, firstly note that
$$H(X^n) = - \int q(X^n) \log q(X^n) dx^n = - \sum_{i=1}^n \int q(X^n) \log q(x_i) dx^n$$
$$= - \sum_{i=1}^n \int q(x_i) \log q(x_i) dx_i = nS$$

Thus, 
$$E[F_n] - nS = E[-\log p(X^n)] + \int q(X^n) \log q(X^n) dx^n$$
$$= - \int q(X^n) \log p(X^n) dx^n + \int q(X^n) \log q(X^n) dx^n$$


Then $E[F_n] - nS = \int q(X^n) \log \frac{q(X^n)}{p(X^n)} dx^n = K(q(X^n) || p(X^n))$

Smaller $E[F_n]$ means KL divergence decreases ($\Rightarrow KL \geq 0$) hence $p(X^n)$ is a better estimate for $q(X^n)$.

Thus $E[G_n] - S$ is the average KL divergence from $q(x)$ to $p(x|X^n)$ whereas $E[F_n] - nS$ is their sum.

We now prove yet another important theorem. 

>Theorem: Let $n \geq 1$. The average generalization loss is equal to the increase in free energy.
>$$E[G_n] = E[F_{n+1}] - E[F_n]$$
>Thus $E[F_n] = \sum_{i=1}^{n-1} E[G_i] + E[F_1]$

Proof:
For an arbitrary fn $f(x)$, $\int q(x) f(x) = E_{X_{n+1}} [f(X_{n+1})]$

Now, $G_n = - \int q(x) \log p(x|X^n) dx$
$$= - E_{X_{n+1}} [\log p(X_{n+1} | X^n)]$$
$$= - E_{X_{n+1}} \left[ \log \frac{\int p(X_{n+1}|w) \varphi(w) \prod p(X_i|w) dw}{\int \varphi(w) \prod p(X_i|w) dw} \right]$$

$$= - E_{X_{n+1}} [\log Z(X^{n+1})] + \log Z(X^n)$$
Thus, $E[G_n] = E[F_{n+1}] - E[F_n]$.

Remark: As $F_n = - \log Z(X^n)$, the correspondence between free energy and marginal likelihood is one to one.

However, in general, asymptotic order of the marginal likelihood as a random variable is not equal to its average, whereas for free energy that is the case. We have not proved this yet. Thus, for asymptotic statistics, free energy is a more convenient random variable. 

We can illustrate the failure for the former: Let the marginal likelihood ratio for $r(X^n) = \dfrac{Z(X^n)}{q(X^n)}$ Then $E[r(X^n)] = 1$. However, this is a result: $r(X^n) \rightarrow 0$ in probability.


## Meaning of marginal likelihood

Assume that $p_0(p,\varphi)$ is the prior distributions of a model $p(x|w)$ and a prior $\varphi(w)$.
Then $p(X^n | p,\varphi) = \int \prod p(X_i|w) \varphi(w) dw = Z(X^n)$

By Bayes thm, $p(p,\varphi | X^n) = \dfrac{p(X^n | p,\varphi) p_0(p,\varphi)}{p(X^n)}$

Thus if $n$ is sufficiently large, maximizing $Z(X^n)$ is maximizing $p(p,\varphi | X^n)$

## Conditional independent Cases

We will make the definitions for the conditionally independent case. They are quite similar. 

Let us assume $X^n$ is dependent but $Y^n$ is conditionally independent.

For an arbitrary function $f : (X^n, Y^n) \rightarrow f(X^n, Y^n) \in \mathbb{R}$,
$$E[f(X^n, Y^n)] = \int \dots \int f(X^n,Y^n) \prod q(y_i|x_i) dy_1 \dots dy_n$$
which is a function of $X^n$.

$$S = -\frac{1}{n} \sum_{i=1}^n \int q(y_i|x_i) \log q(y|x_i) dy$$
$$S_n = -\frac{1}{n} \sum_{i=1}^n \log q(Y_i|X_i)$$

$$p(w | X^n, Y^n) = \frac{1}{Z(X^n, Y^n)} \varphi(w) \prod p(Y_i | X_i, w)$$
$$p(Y | X, X^n, Y^n) = \int p(Y | X, w) p(w | X^n, Y^n) dw$$

Everything else is defined similarly. But in this case, $E[C_n] \neq E[G_{n-1}]$
And $E[W_n] = E[G_n] + O(1/n)$

For example: 
1) Regression problem $f(y_i|x_i)$ for a fixed set $\{x_i\}$ are studied. Cross validation cannot be employed.

2) Consider the time series expressed by the relation
 $$
 Z_t = a_1 Z_{t-1} + a_2 Z_{t-2} + a_3 Z_{t-3} + \text{Gaussian Noise}
 $$.
This can be understood as a regression problem
$$(Z_{t-1}, Z_{t-2}, Z_{t-3}) = X_t \rightarrow Y_t = Z_t.$$

Thus $X_t$ is dependent, so cross validation loss cannot be employed. WAIC, however, can be. It is superior. 

## Exercises 

I now refer you to the first set of exercises. 

