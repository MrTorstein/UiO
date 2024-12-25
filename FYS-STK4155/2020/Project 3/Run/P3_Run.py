import os, sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import QTable
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Class.P3_Class import NN_Reg, Data_reg


# Extracting data
# Redshift
dr7_spcobj = QTable.read("Data/gal_info_dr7_v5_2.fit")
redshift = dr7_spcobj["Z"]
# Star Formation rate
dr7_sfr = QTable.read("Data/gal_totsfr_dr7_v5_2.fits")
SFR = np.abs(dr7_sfr["AVG"])

# Choosing number of steps between datapoints used
# There are over 700 000 points and this is to many
steps = 100

# Setting up regression classes
A = NN_Reg(redshift[::steps], SFR[::steps] // 50, epochs = 500)
B = Data_reg(redshift[::steps], SFR[::steps])

# Defining mean square estimate arrays
MSE_NN = np.zeros(30)
MSE_OLS = np.zeros(30)
MSE_Ridge = np.zeros(30)
# Looping complexity degrees
for i in range(1, 31):
    A.Set_up_data(deg = i) # Setting up Neural network
    A.Train() # Training network

    test_pred = A.Pred_probabilities() # Predict probabilities for test data

    # Using probabilities to produce prediction
    Y_pred_1 = np.zeros(len(test_pred))
    for j in range(len(test_pred)):
        Y_pred_1[j] = np.random.choice([1, 100], 1, p = (1 - test_pred[j]))
    MSE_NN[i - 1] = A.EstPredErr(Y_pred_1, A.Y_test) # Calculating mean square error


    B.Set_up_data(deg = i, lam = 1e-5) # Setting up linear regression methods

    # Performing Ordinary least squares
    B.PlaneOLSReg()
    Y_pred_2 = B.Pred_data(B.X_test)
    MSE_OLS[i - 1] = B.EstPredErr(Y_pred_2, B.Y_test)
    
    # Performing Ridge regression
    B.PlaneRidgeReg()
    Y_pred_3 = B.Pred_data(B.X_test)
    MSE_Ridge[i - 1] = B.EstPredErr(Y_pred_3, B.Y_test)

# Plot MSE
plt.plot(range(1, 31), MSE_NN)
plt.plot(range(1, 31), MSE_OLS)
plt.plot(range(1, 31), MSE_Ridge)
plt.legend(["NN", "OLS", "Ridge"])
plt.xlabel("Complexity of model")
plt.ylabel("MSE")
plt.savefig("Run/Figure1.png")


# Using a good fit to make a prediction of the whole dataset and ploting result
degree = 25
A.Set_up_data(deg = degree) # Setting up Neural network
A.Train() # Training network

pred = A.Pred_probabilities(X = A.X) # Predict probabilities for test data

# Using probabilities to produce predictions
Y_pred_1 = np.zeros(len(pred))
for j in range(len(pred)):
    Y_pred_1[j] = np.random.choice([1, 100], 1, p = (1 - pred[j]))
MSE_NN = A.EstPredErr(Y_pred_1, A.y) # Calculating mean square error


B.Set_up_data(deg = degree, lam = 1e-5) # Setting up linear regression methods

# Performing Ordinary least squares
B.PlaneOLSReg()
Y_pred_2 = B.Pred_data(B.X)
MSE_OLS = B.EstPredErr(Y_pred_2, B.y)

# Performing Ridge regression
B.PlaneRidgeReg()
Y_pred_3 = B.Pred_data(B.X)
MSE_Ridge = B.EstPredErr(Y_pred_3, B.y)

# Plot Predictions
plt.figure()
plt.plot(A.X[:, 1], Y_pred_1, "ro")
plt.plot(B.X[:, 1], Y_pred_2, "go")
plt.plot(B.X[:, 1], Y_pred_3, "yo")
plt.plot(redshift[::steps], SFR[::steps], "b.")
plt.title("Complexity of models = %i"%degree)
plt.xlabel("Readshift z")
plt.ylabel("Star Formation Rate [Solar mass per year]")
plt.legend(["NN Predictions", "OLS Prediction", "Ridge Prediction", "Data Points"])
plt.savefig("Run/Figure2.png")

plt.show()