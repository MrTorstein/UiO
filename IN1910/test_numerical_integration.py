from numerical_integration import left_riemann_sum, midpoint, integrate
from numpy import e, log
from pytest import mark, raises

def test_left_riemann_sum():
    """ Test for the left riemann sum integration function """
    
    tol = 0.001 # Tolerance for error
    
    assert abs(left_riemann_sum(lambda x: 3 * x ** 2 * e ** (x ** 3), 0, 1, 9000) - (e - 1)) <= tol

def test_left_midpoint():
    """ Test for the midpoint integration function """
    
    tol = 0.001 # Tolerance for error
    
    assert abs(midpoint(lambda x: 3 * x ** 2 * e ** (x ** 3), 0, 1, 7000) - (e - 1)) <= tol

@mark.parametrize("method, steps", [("left_riemann_sum", 9000), ("midpoint", 7000)])
def test_integrate(method, steps):
    """ Parametrized tests for integrate method """
    
    tol = 0.001 # Tolerance for error
    
    assert abs(integrate(lambda x: 3 * x ** 2 * e ** (x ** 3), 0, 1, steps, method) - (e - 1)) <= tol

def test_integrate_raises_ValueError_for_invalid_method():
    """ Test for handeling wrong input for method in integrate """
    
    with raises(ValueError):
        integrate(lambda x: 3 * x ** 2 * e ** (x ** 3), 0, 1, 2, "feil")