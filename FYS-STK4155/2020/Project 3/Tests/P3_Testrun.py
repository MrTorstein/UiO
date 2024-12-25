import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Class.P3_Class import NN_Reg, Data_reg

"""
Test of the NN_Reg and Data_reg for a 1 dimensional 
"""

test_array = np.linspace(0, 1, 2)
test_res = np.array([[1, 1]])
success = 0
tol = 1.

# Setting up regression classes
A = NN_Reg(test_array, test_array, epochs = 500)
B = Data_reg(test_array, test_array)


A.Set_up_data(deg = 1) # Setting up Neural network
A.Train(ratio = 0.01) # Training network
Y_pred_1 = A.Pred_data(X = test_res) # Predict probabilities for test data

MSE_NN = A.EstPredErr(Y_pred_1, test_res) # Calculating mean square error
if MSE_NN > tol:
    print("NN Test Fail")
    print(MSE_NN)
else:
    success += 1

B.Set_up_data(deg = 1, lam = 1e-5) # Setting up linear regression methods

# Performing Ordinary least squares
B.PlaneOLSReg(ratio = 0.01)
Y_pred_2 = B.Pred_data(test_res)
MSE_OLS = B.EstPredErr(Y_pred_2, test_res)
if MSE_OLS > tol:
    print("OLS Test Fail")
    print(MSE_OLS)
else:
    success += 1

# Performing Ridge regression
B.PlaneRidgeReg(ratio = 0.01)
Y_pred_3 = B.Pred_data(test_res)
MSE_Ridge = B.EstPredErr(Y_pred_3, test_res)

if MSE_Ridge > tol:
    print("Ridge Test Fail")
    print(MSE_Ridge)
else:
    success += 1

if success == 3:
    print("Test is a success")