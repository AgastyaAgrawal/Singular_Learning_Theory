---
title: The Laplace Method
---

## Foreword 

We are going to discuss an important technique of estimating asymptotics of integrals. It is very relevant to the kind of functions we deal in regular statistical theory as well. This method is called the Laplace Method. This read is quite nice, it involves some low cortisol exercises and theory, which includes some nice techniques and facts. 

# The Laplace Method for Integrals 

We want to approximate 
$$
F(t) = \int_{-\infty}^{\infty}\varphi(x, t)dx \qquad(t \to \infty)
$$

The idea is simple: If by chance the function has a sharp peak around some $x_0$ for all large $t$, pehaps we can just consider the integral in that neighborhood and estimate the entire integral by that value. 

WLOG, say $x_0 = 0$. What we are saying is that it would be quite nice if: 
$$
\int_{-\infty}^{\infty} - \int_{-\delta}^{\delta} << \int_{-\delta}^{\delta}
$$
It would be nicer if in $|x| < \delta$, we can approximate $\varphi$ by some simpler function, perhaps by taylor approximation. 

**Example**
$$
\int_{-\infty}^{\infty}e^{-tx^2}\log(1 + x + x^2)dx \approx \int_{-1/2}^{1/2}e^{-tx^2}(x + \frac{x^2}{2} - \frac{2}{3}x^3)dx
$$

## Some Nice Facts 

The integrands of the type $e^{-tx^2}x^k$ will come often, so let us deal with them beforehand. 

For k odd, since the integrand is an odd function, the integral will of course be $0$. 
$$
\int_{-\infty}^{\infty}e^{-tx^2}x^kdx = 0
$$
For k even, we do a simple substitution and then use the gamma integral (we thus require $Re(t) > 0$).
$$
\int_{-\infty}^{\infty}e^{-tx^2}x^kdx = 2\int_{0}^{\infty}e^{-tx^2}x^kdx
$$
Set $tx^2 = y$ and substitute that, $2txdx = dy$. We get
$$
t^{-(k + 1)/2}\int_0^{\infty}e^{-y}y^{(k+1)/2}dy = t^{-(k + 1)/2}\Gamma(\frac{k + 1}{2})
$$
We get these three estimates/equations as well. 
$$
\int_{-\infty}^{\infty} |e^{-tx^2}x^k| dx = \mathrm{O}\{(\operatorname{Re} t)^{-\frac{1}{2}(k+1)}\} \qquad (\operatorname{Re} t > 0).
$$

$$
\int_{0}^{\infty} e^{-tx}x^k dx = t^{-k-1}k! \qquad (\operatorname{Re} t > 0),
$$

$$
\int_{0}^{\infty} |e^{-tx}x^k| dx = \mathrm{O}\{(\operatorname{Re} t)^{-k-1}\} \qquad (\operatorname{Re} t > 0),
$$

## General Case 

Let us tackle the general case. Just observe that we are making exactly the same assumptions that regularity imposes on $-K_n(\omega)$. 

We will consider the function 
$$
F(t) = \int_{-\infty}^{\infty}e^{th(x)}dx
$$
with the following assumptions on $h(x)$:
1) $h(x)$ is a real and continuous function. It attains its absolute maximum, and it does so only at a single point. WLOG assume that is at $x = 0$, and WLOG, assume that the maximum value is $0$. This means that $h$ is negative everywhere else. 

2) There exist real numbers $b, c$ such that
$$
h(x) \leq -b \quad \text{if } \quad |x| \geq c
$$
This is somewhat more general. This happens whenever $h(x) \to -\infty$ whenver $|x| \to \infty$. 

3) The integral should converge for sufficiently large values of $t$. We will assume for simplicity it happens at $t = 1$. 

4) Assume $h'(x)$ exists in some neighborhood of $0$, and that $h''(0)$ exists and is less than $0$. Of course, by maximality, $h'(0) = 0$. 

It follows from our assumptions that given any positive $\delta$, there exists $\eta(\delta$) such that $h(x) \leq -\eta(\delta)$ in $|x| \geq \delta$.

Draw a graph and make an argument, this is pretty easy. 

Now since $th(x) = (t - 1)h(x) + h(x)$, 
$$
\int_{-\infty}^{-\delta} + \int_{\delta}^{\infty} < e^{-(t-1)\eta(\delta)} \int_{-\infty}^{\infty} e^{h(x)} dx \qquad (t > 1).
$$

Why did we do this? We cannot directly bound $th(x)$, as then the upper bound will diverge to $\infty$. The integral on the RHS is just a number now. As $t \to \infty$,RHS goes to $0$.

Thus, outside a $\delta$ neighborhood of 0, the integral is pretty small. Let us check the integral size in the neighborhood of $0$. 


We do a taylor approximation, as illustrated in the first example. We justify the use formally. 
$$
|h(x) - \frac{x^2h''(0)}{2}| = O(x^2)
$$
Now, this is not a one liner since the third derivative may not exist. So you can show this by considering $g(x) = h(x) - \dfrac{x^2h''(0)}{2}$, and then using mean value theorem on it. Easy stuff. 

Choosing a small enough $\delta$ it is reasonable to replace $h(x)$ by the second order taylor term then. 

$$
\int_{-\delta}^{\delta} e^{\frac{1}{2}tx^2(h''(0)-2\varepsilon)} dx < \int_{-\delta}^{\delta} e^{th(x)} dx < \int_{-\delta}^{\delta} e^{\frac{1}{2}tx^2(h''(0)+2\varepsilon)} dx.
$$

Now the upper bound and lower bound can both be calculated, we already have show how. 
$$
\int_{-\infty}^{\infty} e^{th(x)} dx \sim (2\pi)^{\frac{1}{2}}(-th''(0))^{-\frac{1}{2}} \qquad (t \to \infty).
$$

## Asymptotic Expansions

For simplicity, we will assume that in some interval $|x| \leq \delta$,
$$
h(x) = a_2x^2 + a_3x^3 + ...
$$
and that $a_2 < 0$. 


We will now discuss the more general case 
$$
F(t) = \int_{-\infty}^{\infty}g(x)e^{th(x)}dx
$$

where $g(x)$ is an integrable function and in the interval $|x| \leq \delta$, equal to a convergent power series
$$
g(x) = b_0 + b_1x + b_2x^2 + ...
$$

We need some rough estimate saying that the contributions of the intervals +-$(\delta, \infty)$. We assume that for each positive integer $M$, we have 
$$
\int_{-\infty}^{-\delta}g(x)e^{th(x)}dx = O(t^{-m}) \quad \int_{\delta}^{\infty}g(x)e^{th(x)}dx = O(t^{-m})
$$
if $t \to \infty$. 
We will also assume $h(x) = -O(x^2)$. 

We consider $\exp(ta_2x^2)$ as the main factor of the integrand. The remaining factor 
$$
g(x)\exp\{tx^3(a_3 + a_4x + a_5x^2 + ...)\}
$$
We expand this as a double power series 
$$
P(tx^3, x) = \sum_{m = 0}^{\infty}\sum_{n = 0}^{\infty}c_{mn}(tx^3)^mx^n
$$
We want to approximate $P$ uniformly by its partial sums, and thus restrict $tx^3$ to some finite interval. That is, we will use this power series only if $|x| \leq t^{-1/3}$. We abbreviate $t^{-1/3} = \tau$. We consider a delta interval $0 \leq \tau \leq \delta$. 
We first notice that +-$(\tau, \delta)$ can be neglected.




