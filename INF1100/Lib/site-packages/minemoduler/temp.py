__all__ = ["F2C", "F2K", "C2F", "C2K", "K2C", "K2F"]

def F2C(F):
	return (F - 32)/1.8



def F2K(F):
	return (F + 459.67)* 5./9.



def C2F(C):
	return C*1.8 + 32



def C2K(C):
	return C + 273.15



def K2C(K):
	return K - 273.15



def K2F(K):
	return K*9./5.-459.67



if __name__ == "__main__":
	print "run as program"
	from minetester import test_F2C, test_F2K, test_C2F, test_C2K, test_K2C, test_K2F
	test_F2C(); test_F2K(); test_C2F(); test_C2K(); test_K2C(); test_K2F()