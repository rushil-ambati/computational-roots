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

approximate_root(-4, 10, True)