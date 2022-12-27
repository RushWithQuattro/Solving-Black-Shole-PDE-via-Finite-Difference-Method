# Solving-Black-Shole-PDE-via-Finite-Difference-Method

This project focus on pricing different styles of options under the Black-Scholes-Merton framwork by numerically solving the BSM PDE using 3 finite difference methods. 

In this project, we provide the python implmentation of two finite difference method(Implicit and Crank-Nicolson Method) to solve the PDEs of 3 styles of options(European, American, and Bermuda with path dependent payoff). We present and compare the accuracy of solutions and the efficiency of excuting the algorithms.

### **Types of Options**

Under the Black-Scholes-Merton framework, the stock price without dividend yield is modeled by a exponential martingale following SDE:

> $dS_t=𝑢S_tdt+σS_tdW_t$

Under risk neutral measure $Q$, This exponential martinagle can be represented as
> $dS_t=rS_tdt+σS_tdW_t$

An typical vanilla call option has the pay-off function $Max(S_T-K,0)$ at maturity. While a typical put option has the pay-off function $Max(K-S_T,0)$. In this projects, we only deal with options with the above form of pay-off function, but with different conditions on when can the option be excercised. 

One can show that, under non-arbitrage condition, a delta hedged portoflio consists of $∂C/∂S$ shares of underlying stock and 1 european style option $C$ can be formed in a self-financing way. Then, the following holds true for this delta hedged portfolio.
> $r(C-\frac{∂C}{∂S} S_t)dt = (\frac{∂C}{∂t}+\frac{1}{2}\frac{∂^2C}{∂^2S}S_tσ^2)dt$

Re-arrange the terms we obtain the Black-Schole PDE for european option with non-dividend paying underlying,

> $\frac{∂C}{∂t}+ rS_t\frac{∂C}{∂t}+\frac{1}{2}\frac{∂^2C}{∂^2S}S_tσ^2=rC$


In the case of the American options, it has the same pay-off function as its european syle counterpart. However, unlike the european options, which can be only be excercised at the agreed maturity date, the holder of american option can excercise it at any time during the contract period, as long as she feels beneficial to do so. 

Then, we can show that the American option has the following dynamics,

> $min(-\frac{∂C}{∂t}-rS_t\frac{∂C}{∂t}-\frac{1}{2}\frac{∂^2C}{∂^2S}S_tσ^2-rC, C-w(S-K))=0$

A Bermuda option is smilar to american option in the way that they both allow excercise before the contract maturity. The difference is that bermuda option can only be excercised at a fixed set of date before maturity. 

