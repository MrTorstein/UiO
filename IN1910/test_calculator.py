from calculator import add, divide, factorial, sin, sqrt
from pytest import mark, raises
from numpy import pi, sqrt as Sqrt

@mark.parametrize("addend, sum", [((1, 2), 3), ((0.1, 0.2), 0.3)])
def test_add(addend, sum):
    """ Testing addition function """
    tol = 1e-13 # Tolerance for error
    
    assert abs(add(*addend) - sum) <= tol

@mark.parametrize("fraction, quotient", [((2, 1), 2), ((1, 2), 0.5), ((0.1, 0.2), 0.5)])
def test_divide(fraction, quotient):
    """ Testing division function """
    
    assert divide(*fraction) == quotient

@mark.parametrize("number, result", [(0, 1), (1, 1), (5, 120),(10, 3628800)])
def test_factorial(number, result):
    """ Testing factorial function """
    
    assert factorial(number) == result

@mark.parametrize("angle, result", [(0, 0), (pi / 4, 1  / Sqrt(2)), (pi / 2, 1), (3 * pi / 2, -1)])
def test_sin(angle, result):
    """ Testing sinus function """
    
    tol = 1e-13 # Tolerance for error
    
    assert abs(sin(angle) - result) <= tol

@mark.parametrize("number, root", [(0, 0), (1, 1), (4, 2), (5, 2.23606797749979)])
def test_sqrt(number, root):
    """ Testing sqrt function """
    
    tol = 1e-10 # Tolerance for error
    
    assert abs(sqrt(number) - root) <= tol

def test_factorial_raises_ValueError_for_negatives():
    """ Test for unwanted negativ input in factorial function """
    from pytest import raises
    
    with raises(ValueError):
        factorial(-1)

def test_factorial_raises_TypeError_for_decimals():
    """ Test for unwanted decimal input in factorial function """
    from pytest import raises
    
    with raises(TypeError):
        factorial(0.5)

def test_divide_raises_ZeroDivisionError():
    """ Test for handeling of division by zero """
    from pytest import raises
    
    with raises(ZeroDivisionError):
        divide(1, 0)