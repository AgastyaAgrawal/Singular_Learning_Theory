These are the solutions to the selected exercises from Watanabe's green book. Majority of them are from chapter 1. Please refer to the exercises from the pdf, which is freely available online. 

# Problems

## Problem 1
(a) Let $w_0 = (1, 1, \dots, 1) \in \mathbb{R}^{10}$, and let $W$ be a random variable on $\mathbb{R}^{10}$ which is subject to
$$p(w) = c \{ \exp(- ||w||^2) + 100 \exp(-10 ||w - w_0||^2) \}$$


Let $w = x w_0 + v$. Then, 

$$
||w||^2 = ||w_0||^2 x^2 + ||v||^2 = 10 x^2 + ||v||^2
$$
and 

$$
||w - w_0||^2 = ||(x-1) w_0 + v||^2 = 10(x-1)^2 + ||v||^2
$$

Thus, 
$$
p(x, v) = c \{ \exp( - (10 x^2 + ||v||^2) ) + 100 \exp( -10 (10(x-1)^2 + ||v||^2) ) \}

$$
To maximize $p(x,v)$ we can do so separately.
We minimize $||v||^2 \Rightarrow v=0$.

$\text{argmax } p(x, v) = (\text{argmax } p(x, 0), 0)$

$$p(x, 0) = c \{ \exp(-10 x^2) + 100 \exp(-100(x-1)^2) \}$$

By checking the derivative we have $x \in (0, 1)$
As the function is dominated by the second term, $x \approx 1$.
Then $w \approx w_0$.

(b) $E[W] \approx 0$

$$E[W] = \int_{\mathbb{R}^{10}} w p(w) dw$$

We know that $\int_{\mathbb{R}^{10}} p(w) dw = 1$
So
$$c \left[ \int_{\mathbb{R}^{10}} e^{-||w||^2} dw + 100 \int_{\mathbb{R}^{10}} e^{-10 ||w - w_0||^2} dw \right] = 1$$

Remember: $\int_{\mathbb{R}^d} \exp(-\alpha ||x - \mu||^2) dx = \left( \frac{\pi}{\alpha} \right)^{d/2}$

Then $c \left[ \pi^5 + 100 \cdot \left( \frac{\pi}{10} \right)^5 \right] = 1$.

$$E[W] = c \left[ \underset{\downarrow \text{odd}}{\int_{\mathbb{R}^{10}} w \exp(-||w||^2) dw} + \int_{\mathbb{R}^{10}} w \exp(-10 ||w - w_0||^2) dw \right]$$

$$= c \left[ \int_{\mathbb{R}^{10}} w \exp(-10 ||w - w_0||^2) dw \right]$$

$$= c (10^{-5} \pi^5 w_0) \approx 0.000999 w_0 \approx 0.$$

Even though the pdf peaks at $w_0$, the vol is very small as the var is $1/10$ in all $10$ dims, Ah the curse of dimensionality. 

The takeaway from this is that MAP can be a pretty bad estimator for the right parameter. Notice that as we consider higher dimensions, the distance between the expected value and the MAP estimate increases to arbitrarily large values. 

## Problem 2 -  Fluctuation Dissipation Theorem

Let $\beta > 0$, and $H(x) : \mathbb{R}^n \rightarrow \mathbb{R}$. Say $X \in \mathbb{R}^n$ is subject to a pdf.

$$p(x | \beta) = \frac{1}{Z(\beta)} \exp(-\beta H(x)) \quad \text{where } Z(\beta) = \int \exp(-\beta H(x)) dx$$

We need to prove that $\dfrac{\partial E[H(X)]}{\partial \beta} = - \mathbb{V}[H(X)]$

$$E[H(X)] = \int \frac{1}{Z(\beta)} \exp(-\beta H(x)) H(x) dx$$

$$\frac{\partial Z(\beta)}{\partial \beta} = \int \frac{\partial}{\partial \beta} \exp(-\beta H(x)) dx = \int -H(x) \exp(-\beta H(x)) dx$$
$$= - Z(\beta) E[H(X)]$$

$$\frac{\partial f(\beta)}{\partial \beta} = \int -H^2(x) \exp(-\beta H(x)) dx = - Z(\beta) E[H^2(X)]$$

Thus,
$$\frac{\partial E[H(X)]}{\partial \beta} = \frac{- Z^2(\beta) E[H^2(X)] + Z^2(\beta) E^2[H(X)]}{Z^2(\beta)}$$

$$= E^2[H(X)] - E[H^2(X)]$$
$$= - \mathbb{V}[H(X)]$$

Thus, $\dfrac{\partial E[H(X)]}{\partial \beta} = - \mathbb{V}[H(X)]$

## Problem 3
Let $p(x|a)$ be a statistical model of $x \in \{0, 1\}$ defined by

$$p(x|a) = a^x (1-a)^{1-x} \quad \text{where } 0 \leq a \leq 1, \quad \varphi(a) = 1 \text{ is the prior}.$$

Let $X^n$ be independently subject to $p(x|a_0)$.

Let $n_1 = \sum_{i=1}^n x_i, \quad n_2 = n - n_1$

1) Find the MLE.

$$\log p(x_i | a) = x_i \log a + (1-x_i) \log(1-a)$$

$$f(a) = \sum \log p(x_i | a) = n_1 \log a + n_2 \log(1-a)$$

Better way: Find $\text{argmax } \prod a^{x_i} (1-a)^{1-x_i}$

$$= a^{n_1} (1-a)^{n-n_1}$$

As $0 \leq a \leq 1$, we can use AM-GM to maximize. The argmax satisfies:

$$\frac{a}{n_1} = \frac{1-a}{n-n_1}$$
$$n a - n_1 a = n_1 - n_1 a$$
$$\Rightarrow a = \frac{n_1}{n}$$


2) Estimated probability distribution $p(x | \hat{a})$ (Frequentist estimation)

$$p(1 | \hat{a}) = \hat{a} = \frac{n_1}{n}, \quad p(0 | \hat{a}) = 1 - \hat{a} = \frac{n - n_1}{n} = \frac{n_2}{n}$$

3) Bayesian Predictive distribution $p(x | X^n)$

Let us calculate the posterior first.

$$p(a | X^n) = \frac{p(X^n | a) \varphi(a)}{\int_0^1 p(X^n | a) \varphi(a) da}$$

$$= \frac{p(X^n | a)}{\int_0^1 p(X^n | a) da} \quad (\text{this can be calculated by multiplication})$$

So, $p(x | X^n) = \int_0^1 p(x | a) p(a | X^n) da$

$$= \frac{\int_0^1 p(x | a) p(X^n | a) da}{\int_0^1 p(X^n | a) da}$$

Now, $p(X^n | a) = \prod p(x_i | a) = a^{\sum x_i} (1-a)^{n - \sum x_i} = a^{n_1} (1-a)^{n_2}$

Now, $\int_0^1 p(X^n | a) da = \int_0^1 a^{n_1} (1-a)^{n_2} da = \beta(n_1 + 1, n_2 + 1)$


$$p(1 | X^n) = \frac{\int_0^1 a^{n_1+1} (1-a)^{n_2}}{\beta(n_1+1, n_2+1)} = \frac{\beta(n_1+2, n_2+1)}{\beta(n_1+1, n_2+1)}$$

$$p(0 | X^n) = \frac{\beta(n_1+1, n_2+2)}{\beta(n_1+1, n_2+1)}$$

Now, $\beta(\alpha, \beta) = \frac{\Gamma(\alpha) \Gamma(\beta)}{\Gamma(\alpha+\beta)}$

Thus, $p(1 | X^n) = \frac{n_1+1}{n+2}$

$$p(0 | X^n) = 1 - \frac{n_1+1}{n+2} = \frac{n_2+1}{n+2}$$