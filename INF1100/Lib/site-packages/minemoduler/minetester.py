from temp import *
def test_F2C():
		a = 32; b = 0; c = 0;
		computed = F2C(a)
		expected = 0
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%.14f er ikke lik %.14f" %(computed, expected)
		assert success, msg

def test_F2K():
		a = -459.67; b = 0; c = 0;
		computed = F2K(a)
		expected = 0
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%.14f er ikke lik %.14f" %(computed, expected)
		assert success, msg

def test_C2F():
		a = 0; b = 0; c = 0;
		computed = C2F(a)
		expected = 32
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%.14f er ikke lik %.14f" %(computed, expected)
		assert success, msg

def test_C2K():
		a = 0; b = 0; c = 0;
		computed = C2K(a)
		expected = 273.15
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%.14f er ikke lik %.14f" %(computed, expected)
		assert success, msg

def test_K2C():
		a = 100; b = 0; c = 0;
		computed = K2C(a)
		expected = -173.14999999999998
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%.14f er ikke lik %.14f" %(computed, expected)
		assert success, msg

def test_K2F():
		a = 100; b = 0; c = 0;
		computed = K2F(a)
		expected = -279.67
		tol = 1.e-14
		success = abs(computed - expected) < tol
		msg = "%g er ikke lik %g" %(computed, expected)
		assert success, msg