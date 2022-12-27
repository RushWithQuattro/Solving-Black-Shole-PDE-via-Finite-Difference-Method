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


### **Pricing by solving the PDEs**

Now with the PDEs describing the dynamics of the options, the problem of option pricing now becomes solving the PDEs numerically given a set of boundaries(typically involves the terminal pay-offs at maturity). For this project we employs the finite difference method. The main idea behind the finite difference method is to discretize the continues evolution of the option value on a grid consists of discrete points. By specifying the proper boundaries of the grid one can iteratively compute the values at any points on the grid. There are 3 following ways in which such method can be implemented.

### **The Implicit Method** 
Under this method, the components of the Black-Schole PDE can be approximated as follows
> *Forward Approximation* of $\frac{∂C}{∂t}$ as $\frac{C_{t+1,j}-C_{t,j}}{δt}$

> *Central Approximation* of $\frac{∂C}{∂S}$ as $\frac{C_{t,j+1}-C_{t,j-1}}{2δS}$

> *Standard Approximation* of $\frac{∂^2C}{∂^2S}$ as $\frac{C_{t,j+1}+C_{t,j-1}-2C_{t,j}}{(δS)^2}$

Plug in the above terms to the Black_Schole PDE, then for each time step we can obtain the following matrix formuation,

> $BC_t=C_{t+1}+K_t$

We solve the above equation using LU decompsition at each time step. 


### **The Crank-Nicolson Method** 
Under this method, the components of the Black-Schole PDE can be approximated as follows
> *Central Approximation* of $\frac{∂C_{t-(1/2),j}}{∂t}$ as $\frac{C_{t,j}-C_{t-1,j}}{δt}$

> *Central Approximation* of $\frac{∂C_{t-(1/2),j}}{∂S}$ as $\frac{1}{2}(\frac{C_{t-1,j+1}-C_{t-1,j-1}}{2δS}+\frac{C_{t,j+1}-C_{t,j-1}}{2δS})$

> *Standard Approximation* of $\frac{∂^2C}{∂^2S}$ as $\frac{1}{2}(\frac{C_{t-1,j+1}-2C_{t-1,j}+C_{t-1,j-1}}{δS^2}+\frac{C_{t,j+1}-2*C_{t,j}+C_{t,j-1}}{δS^2})$

Plug in the above terms to the Black_Schole PDE, then for each time step we can obtain the following matrix formuation,

> $AC_{t}=DC_{t+1}+K_{t-1}+K_{t}$

We solve the above equation using LU decompsition at each time step. 




###**Stability and Rate of Convergence**

One concern about finite difference method is that some incorrectly chosen parameters(in this case, $σ,r,δt$) might cause the algorithm to be unstable(i.e fail to converge to a finite solution). The main source of this problem lies at the coefficient matrix ( A, D and B). For implicit method, it is stable if the infinite norm of $B^{-1} <=1$; for Crank Nicolson method, it is stable if the infinite norm of $C^{-1}D <=1$. Luckily, these conditions are true for any $σ,r,δt$ used for the above methods. However, for other methods like the explicit method, stability is not guaranteed. This is why in this project we mainly focus on the above two methods. 

Furthermore, one can show theoretically that Crank Nicolson Method has the fastest rate of convergence. In this project we will verify this statement.


###**Algorithm**
We implement the following algorithm to interatively solving for the value of european option. The general form is
```
# implicit method, European Call
for 1:T 
   A*C_t=C_t+1 
```

Considering the early excercise feature of the American and Bermudan option,we compare the intrinsic value of the option against the pay-off of early-excercise to determine the value at that time step. Thus, the modified algorithm takes the general form


```
# implicit method, American Call
for 1:T:
   A*C_t=C_t+1
   for 1:K:
    if C_tk<max(S_k-K): 
      C_tk=S_k-K
  
```
```
# implicit method, Bermudan Call
for 1:T:
   A*C_t=C_t+1
   if t in agreed_excercise_date:
    for 1:K:
      if C_tk<max(S_k-K): 
        C_tk=S_k-K
  
```



