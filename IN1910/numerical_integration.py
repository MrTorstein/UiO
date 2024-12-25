from typing import Callable
from numpy import linspace

def left_riemann_sum(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    h = (b - a) / n
    
    x = linspace(a, b, n)
    
    Sum = 0
    for i in range(n - 1):
        Sum += f(x[i])
    
    return h * Sum

def midpoint(f: Callable[[float], float], a: float, b: float, n: int) -> float:
    h = (b - a) / n
    
    x = linspace(a, b, n)
    
    Sum = 0
    for i in range(n - 1):
        Sum += (f(x[i]) + f(x[i + 1]))
    
    return 0.5 * h * Sum

def integrate(f: Callable[[float], float], a: float, b: float, n: int, method: str = "midpoint") -> float:
    if method == "left_riemann_sum":
        return left_riemann_sum(f, a, b, n)
    elif method == "midpoint":
        return midpoint(f, a, b, n)
    else:
        raise ValueError("Integration method given is not supported. Only left riemann sum and midpoint is supported.")