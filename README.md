# Approximating Roots with Computation

by Rushil A.

[toc]

## Outline

In this mini-project I'll be using numerical methods and computation in order to approximate roots of an equation.

The specific method used here is called the *Newton-Raphson method* (alternatively known as just ‘*Newton’s method*’), named after Isaac Newton and Joseph Raphson.



## Theory

### Newton-Raphson method

> For a single-variable function <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/> defined for a real variable <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> with derivative <img src="/tex/06e7cc81ea7c4442d159c33723c273db.svg?invert_in_darkmode&sanitize=true" align=middle width=13.60737509999999pt height=24.7161288pt/> and an initial guess <img src="/tex/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> for a root of <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/>. If the function satisfies assumptions, and the initial guess is close, then
> <p align="center"><img src="/tex/def96efa0840a1e9096d89e564d4a52e.svg?invert_in_darkmode&sanitize=true" align=middle width=158.18011769999998pt height=38.83491479999999pt/></p>
>
> where <img src="/tex/277fbbae7d4bc65b6aa601ea481bebcc.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> is a better approximation of the root than <img src="/tex/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/>. Geometrically, <img src="/tex/d3d2154b26f6e94a6e22311613d14210.svg?invert_in_darkmode&sanitize=true" align=middle width=45.079973399999986pt height=24.65753399999998pt/> is the intersection of the x-axis and the tangent of the graph of f at <img src="/tex/75f733ac6a2c0cc4c6633d8c13968b79.svg?invert_in_darkmode&sanitize=true" align=middle width=77.07780795pt height=24.65753399999998pt/>: that is, the improved guess is the unique root of the linear approximation at the initial point. The process is repeated as 
> <p align="center"><img src="/tex/71efb9a688dda717db93231d01e5a15c.svg?invert_in_darkmode&sanitize=true" align=middle width=179.54446289999998pt height=38.83491479999999pt/></p>
> until a sufficiently precise value is reached. This algorithm is first in the class of Householder's methods, succeeded by Halley's method. The method can also be extended to complex functions and to systems of equations. 
>
> ([*Source*](https://en.wikipedia.org/wiki/Newton's_method))



### Derivation of the equation

The iterative equation <img src="/tex/a66c3ed9cfd3a3dbb6c3594ee19f5c3b.svg?invert_in_darkmode&sanitize=true" align=middle width=135.24718514999998pt height=33.20539859999999pt/> is the form in the A-level specification too, and we’ll use the same in our program here.
However, we weren’t taught the derivation in class, so I’ll attempt one here.

The geometric description above can be represented in a graph diagram, such as the one below:

> ![img](https://ds055uzetaobb.cloudfront.net/brioche/uploads/7KrMvNiT7l-newtons-method.png?width=1200)
>
> ([*Source*](https://brilliant.org/wiki/newton-raphson-method/))

So we now have a tangent at <img src="/tex/d7084ce258ffe96f77e4f3647b250bbf.svg?invert_in_darkmode&sanitize=true" align=middle width=17.521011749999992pt height=14.15524440000002pt/> that has the following characteristics:

- It goes through the point <img src="/tex/72abfde4985a4ee3584cb93adb7d518c.svg?invert_in_darkmode&sanitize=true" align=middle width=79.38001499999999pt height=24.65753399999998pt/>
- It has gradient <img src="/tex/61f945ed1db6043d4a22807c1730260b.svg?invert_in_darkmode&sanitize=true" align=middle width=45.55764674999998pt height=24.7161288pt/>

So now to get the equation of the tangent, we can put it into point slope form <img src="/tex/09179b5af2eefcc301cfae01e22d66c0.svg?invert_in_darkmode&sanitize=true" align=middle width=139.56607994999996pt height=24.65753399999998pt/>:
<p align="center"><img src="/tex/ee6948fd0a250a1b17a969ddc568479f.svg?invert_in_darkmode&sanitize=true" align=middle width=197.77598609999998pt height=17.2895712pt/></p>
We want our value of <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> when <img src="/tex/a42b1c71ca6ab3bfc0e416ac9b587993.svg?invert_in_darkmode&sanitize=true" align=middle width=38.78604674999999pt height=21.18721440000001pt/> (or the root of the equation):
<p align="center"><img src="/tex/b1d85b4c3dae61a96d539571c6e06769.svg?invert_in_darkmode&sanitize=true" align=middle width=197.34598785pt height=17.2895712pt/></p>
Hence if we solve for <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>, which will be the x-intercept and closer to the root (our <img src="/tex/14e12a1273c346610e9daaf5e3aee29a.svg?invert_in_darkmode&sanitize=true" align=middle width=34.16493134999999pt height=14.15524440000002pt/>):
<p align="center"><img src="/tex/ffc076d223f13ddfcd341c625b6b277a.svg?invert_in_darkmode&sanitize=true" align=middle width=431.13298469999995pt height=38.83491479999999pt/></p>
we obtain our desired form.



### Limitations and Practical Considerations

> #### Difficulty in calculating the derivative of a function
>
> Newton's method requires that the derivative can be calculated directly. An analytical expression for the derivative may not be easily obtainable or could be expensive to evaluate. In these situations, it may be appropriate to approximate the derivative by using the slope of a line through two nearby points on the function. Using this approximation would result in something like the *secant method* whose convergence is slower than that of Newton's method.
>
> #### Overshoot
>
> If the first derivative is not well behaved in the neighbourhood of a particular root, the method may overshoot, and diverge from that root. This includes points of inflection, local maxima or minima around <img src="/tex/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> or the root.
>
> #### Features of graph around root
>
> If a stationary point of the function is encountered, the derivative is zero and the method will terminate due to division by zero. 
>
> #### Poor initial estimate
>
> A large error in the initial estimate can contribute to non-convergence of the algorithm. 
>
> #### Discontinuity
>
> This will only work for functions defined for all real numbers. Asymptotes will also yield undefined values.
>
> ([*Source*](https://en.wikipedia.org/wiki/Newton%27s_method#Failure_of_the_method_to_converge_to_the_root))



## Implementation

### Newton-Raphson Method

Firstly, we’ll define our function. Let’s say for the sake of our example we’re trying to find the roots of <img src="/tex/190ff5d661c908c19efa91c6c5979e78.svg?invert_in_darkmode&sanitize=true" align=middle width=112.92207629999997pt height=26.76175259999998pt/>.

```python
def f(x):
    return (x**2-4*x-7)
```

Next, we’ll define the derivative of our function (so we can use it in place of the <img src="/tex/61f945ed1db6043d4a22807c1730260b.svg?invert_in_darkmode&sanitize=true" align=middle width=45.55764674999998pt height=24.7161288pt/> in the equation).

```python
def f_prime(x):
    return 2*x-4
```

Now we’ve got all our setup done, let’s construct a function that will actually carry out the Newton-Raphson method.

```python
def nrf(guess, n):
    print("x_0 = " + str(guess))
    for i in range(n+1):
        next_guess = guess - f(guess) / f_prime(guess)
        print("x_" + str(i+1) + " = " + str(next_guess))
        guess = next_guess
```

Now, let’s call it. Say our starting value is 8, and we want it to iterate 10 times.

```python
nrf(8, 10)
```

This provides the following output:

```s
x_0 = 8
x_1 = 5.916666666666666
x_2 = 5.36258865248227
x_3 = 5.316938934730458
x_4 = 5.316624805231569
x_5 = 5.3166247903554
x_6 = 5.3166247903554
x_7 = 5.3166247903554
x_8 = 5.3166247903554
x_9 = 5.3166247903554
x_10 = 5.3166247903554
```

So our final value for the root closest to 8 is approximately 5.32 (to 3 significant figures.)
The true value of that positive root is <img src="/tex/0df30062e1f5cd381d16442cf5a8ae41.svg?invert_in_darkmode&sanitize=true" align=middle width=58.44749129999999pt height=28.511366399999982pt/> and our approximation of 5.3166… is actually very close, and didn’t take that many iterations to reach there. This is because our starting value of 8 is very close to the actual root.

We can also find the other root of this quadratic equation by starting at another value, say -4, iterating the same number of times.

```python
nrf(-4, 10)
```

This provides the following output:

```
x_0 = -4
x_1 = -1.9166666666666665
x_2 = -1.3625886524822695
x_3 = -1.3169389347304572
x_4 = -1.3166248052315688
x_5 = -1.3166247903553998
x_6 = -1.3166247903553998
x_7 = -1.3166247903553998
x_8 = -1.3166247903553998
x_9 = -1.3166247903553998
x_10 = -1.3166247903553998
```

The other root is <img src="/tex/694c374714055c86ddfca9eaf6925810.svg?invert_in_darkmode&sanitize=true" align=middle width=58.44749129999999pt height=28.511366399999982pt/> and our approximate of -1.3166… is also close. It reached a static value in just 5 iterations.



### Secant Method - Tackling the derivative problem

#### Theory

In the [Limitations and Practical Considerations](#Limitations and Practical Considerations), we considered the possibility that we could be working with a function that may not have an easily obtainable derivative. In this case, I’ll work around this problem by using an alternative method called the *Secant Method*. Although different in name, the only main difference between the *Newton-Raphson method* and the *Secant method* is the way that <img src="/tex/61f945ed1db6043d4a22807c1730260b.svg?invert_in_darkmode&sanitize=true" align=middle width=45.55764674999998pt height=24.7161288pt/> is calculated.

Instead of using an expression to programatically write in the derivative of the function, eg. `x**2-4*x-7`‘s derivative being `2*x-4` we will use a finite difference method to calculate the derivative at a point of <img src="/tex/7997339883ac20f551e7f35efff0a2b9.svg?invert_in_darkmode&sanitize=true" align=middle width=31.99783454999999pt height=24.65753399999998pt/>. The ‘quotient difference’, on which this method relies on, is nothing more than the fundamental idea that differentiation relies on.

We take differentiation from first principles to be as such:
<p align="center"><img src="/tex/bb89ca4ff05639cd0c8423f4dea9f8c6.svg?invert_in_darkmode&sanitize=true" align=middle width=204.00071175pt height=34.7253258pt/></p>
So, if we use a small value of <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> and just use our function for <img src="/tex/7997339883ac20f551e7f35efff0a2b9.svg?invert_in_darkmode&sanitize=true" align=middle width=31.99783454999999pt height=24.65753399999998pt/> we can get an approximation for the value of <img src="/tex/9ce3fa8c71f5905e328dcae5b1d69e2d.svg?invert_in_darkmode&sanitize=true" align=middle width=36.60970994999999pt height=24.7161288pt/> for any <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/>. This means we don’t actually have to evaluate the expression for the derivative, providing a workaround for the limitation discussed above.

#### Implementation

All we have to do now, is to change our `f_prime` function to compute the quotient difference at some <img src="/tex/332cc365a4987aacce0ead01b8bdcc0b.svg?invert_in_darkmode&sanitize=true" align=middle width=9.39498779999999pt height=14.15524440000002pt/> with some difference value <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/>.

```python
def f_prime(x, h):
    return (f(x + h) - f(x)) / h
```

Of course, the smaller our <img src="/tex/2ad9d098b937e46f9f58968551adac57.svg?invert_in_darkmode&sanitize=true" align=middle width=9.47111549999999pt height=22.831056599999986pt/> is, the better the approximation for the derivative, but the more expensive it is computationally.



## Final Code

```python
def f(x):
    # Define your function here
    return (x**2-4*x-7)

def f_prime(x):
    # Define your function's first derivative here
    return 2*x-4

def f_prime_sec(x, h):
    """Calculates the value of the derivative of a function 

    Args:
        x (float): The value for which we are calculating the derivative value
        h (float): Quotient difference

    Returns:
        float: Approximate derivative of f(x) at x, with quotient difference h
    """
    return (f(x + h) - f(x)) / h
    
def approximate_root(guess, n, sec=False, h=1*10**-3):
    """Iteratively estimates the roots of a function by using the Newton-Raphson Method

    Args:
        guess (float): Initial value to iterate from
        n (int): Number of iterations
        sec (bool, optional): Whether or not to use the secant method. Defaults to False.
        h (float, optional): Quotient difference interval. Defaults to 1x10^-3.
    """
    print("x_0 = " + str(guess))
    for i in range(n):
        f_prime_x = f_prime_sec(guess, h) if (sec == True) else f_prime(guess)
        next_guess = guess - f(guess) / f_prime_x # Evaluating x_n+1
        print("x_" + str(i+1) + " = " + str(next_guess))
        guess = next_guess
```

Function can be invoked as follows:

```python
approximate_root(8, 10) # Uses the derivative manually calculated
approximate_root(-4, 10, sec=True, h=1*10**-4) # Uses secant method to approximate derivative, setting the value of h
```



## Points for discussion

- To see how the Secant method compares to the Newton-Raphson method depending on the type of graph considered or other factors.

- To look into other numerical methods that have other use-cases, as well as how they compare to the two considered above.

- I wonder if there’s a more efficient way to write up my algorithm. Also, this is a simple enough algorithm, and doesn’t use any libraries so isn’t difficult to port to other languages. Maybe try porting it to something like C, where it’d actually be somewhat efficient - who knows, maybe even in Assembly?

- It’d be cool to see both of these methods applied to:

  - More than one dimension.
  - Complex functions
  - Systems of equations

- If it’s possible (in some, or even all cases) to mathematically prove whether or not the iterative method will converge given a particular starting value and the equation of the graph.

- I’ve heard about the ‘Laplace Transform’ when talking about numerical methods, I wonder what that is.

- On the Wikipedia page for the Secant method, it discusses the “order of convergence” and the golden ratio. It’d be interesting to see where the golden ratio fits into it all.

  > The iterates <img src="/tex/d7084ce258ffe96f77e4f3647b250bbf.svg?invert_in_darkmode&sanitize=true" align=middle width=17.521011749999992pt height=14.15524440000002pt/> of the secant method converge to a root of <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/>, if the initial values <img src="/tex/e714a3139958da04b41e3e607a544455.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> and <img src="/tex/277fbbae7d4bc65b6aa601ea481bebcc.svg?invert_in_darkmode&sanitize=true" align=middle width=15.94753544999999pt height=14.15524440000002pt/> are sufficiently close to the root. The order of convergence is φ, where
  > <p align="center"><img src="/tex/544cdafa7c09d5606c335f5fc85a3f53.svg?invert_in_darkmode&sanitize=true" align=middle width=180.90827864999997pt height=36.65224035pt/></p>
  > is the golden ratio. In particular, the convergence is superlinear, but not quite quadratic.
  >
  > This result only holds under some technical conditions, namely that <img src="/tex/190083ef7a1625fbc75f243cffb9c96d.svg?invert_in_darkmode&sanitize=true" align=middle width=9.81741584999999pt height=22.831056599999986pt/> be twice continuously differentiable and the root in question be simple (i.e., with multiplicity 1).
  >
  > ([*Source*](https://en.wikipedia.org/wiki/Secant_method#Convergence))