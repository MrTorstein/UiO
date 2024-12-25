def add(x: float, y: float) -> float:
    return x + y

def divide(x: float, y: float) -> float:
    return x / y

def factorial(x: int) -> int:
    if x < 0:
        print("Factorial not defined for negative numbers!")
        raise ValueError
    elif x == 0:
        return 1
    else:
        product = 1
        
        for factor in range(1, x + 1):
            product *= factor
        
        return product

def sin(x: float, N: int = 20) -> float:
    Sum = 0
    for n in range(0, N):
        Sum += (-1) ** n * x ** (2 * n + 1) / (factorial(2 * n + 1))
    return Sum

def sqrt(y: float, x0: float = 1, tol: float = 1e-12) -> float:
    x_n = 0.5 * (y / x0 + x0)
    
    if abs(x_n - x0) >= tol:
        print(abs(x_n - x0))
        x_n = sqrt(y, x_n)
    
    return x_n