"""
Some testruns of the Data reg class
"""
from P1_Class import Data_reg
import numpy as np
from random import random
from imageio import imread
import matplotlib.pyplot as plt


# Make data.
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)
x, y = np.meshgrid(x, y)

def FrankeFunction(x, y):
    term1 = 0.75 * np.exp(-(0.25 * (9 * x - 2) ** 2) - 0.25 * ((9 * y - 2) ** 2))
    term2 = 0.75 * np.exp(-((9 * x + 1) ** 2) / 49.0 - 0.1 * (9 * y + 1))
    term3 = 0.5 * np.exp(-(9 * x - 7) ** 2 / 4.0 - 0.25 * ((9 * y - 3) ** 2))
    term4 = -0.2 * np.exp(-(9 * x - 4) ** 2 - (9 * y - 7) ** 2)
    return term1 + term2 + term3 + term4

z = FrankeFunction(x, y) + np.random.normal(0, 10 ** (-4), x.shape)
"""
A = Data_reg(x, y, z)
A.Set_up_data()
A.PlaneOLSReg()
A.Plot_fitted_and_data()
print("R^2 =", A.R2(A.Pred_data(A.X_test), A.z_test))


# OLS #
deg_list = range(3)

EPE = np.zeros((len(deg_list), 2))
err = np.zeros(len(deg_list))
bias = np.zeros(len(deg_list))
var = np.zeros(len(deg_list))
MSE_est = np.zeros(len(deg_list))
for i in deg_list:
    A.Set_up_data(deg = i)
    A.PlaneOLSReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap #
    err[i], bias[i], var[i] = A.Bootstrapper()

    # Cross Validation #
    MSE_est[i] = A.Cross_Validationer()

plt.figure()
plt.plot(deg_list, EPE[:, 0])
plt.plot(deg_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel("Polynomial degree")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("OLS: Error plot")


plt.figure()
plt.plot(deg_list, err)
plt.plot(deg_list, bias)
plt.plot(deg_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel("Polynomial degree")
plt.title("OLS: Bootstrap")


plt.figure()
plt.plot(deg_list, MSE_est)
plt.plot(deg_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel("Polynomial degree")
plt.ylabel("MSE")
plt.title("OLS: Cross Validation")

plt.show()


# Ridge #
lam_list = np.linspace(-2, 5, 10)

# EPE #
EPE = np.zeros((len(lam_list), 2))
err = np.zeros(len(lam_list))
bias = np.zeros(len(lam_list))
var = np.zeros(len(lam_list))
MSE_est = np.zeros(len(lam_list))
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    A.PlaneRidgeReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap #
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Ridge")

    # Cross Validation #
    MSE_est[i] = A.Cross_Validationer(method = "Ridge")

plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Ridge: Error plot")


plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.title("Ridge: Bootstrap")


plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("MSE")
plt.title("Ridge: Cross Validation")

plt.show()


# Lasso #

# EPE #
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    A.PlaneLassoReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train, "Lasso"), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test, "Lasso"), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap #
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Lasso")

    # Cross Validation #
    MSE_est[i] = A.Cross_Validationer(method = "Lasso")


plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Lasso: Error plot")


plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.title("Lasso: Bootstrap")


plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("MSE")
plt.title("Lasso: Cross Validation")

plt.show()

# Load Data of Oslo fjord
terrain = imread("n59_e010_1arc_v3.tif")
# Can only look at n X n data so must remove some values in y direction
z = terrain[0:int(terrain.shape[1] / 10), 0:int(terrain.shape[1] / 10)]

# Show data
plt.figure()
plt.title("Terrain of Finnemarka")
plt.imshow(z, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")
plt.show()

x = np.linspace(0, z.shape[1], z.shape[1])
y = np.linspace(0, z.shape[0], z.shape[0])
x, y = np.meshgrid(x, y)


A = Data_reg(x, y, z)

# OLS #
deg_list = range(3)

# Defining arrays
EPE = np.zeros((len(deg_list), 2))
MSE_est = np.zeros(len(deg_list))
err = np.zeros(len(deg_list))
bias = np.zeros(len(deg_list))
var = np.zeros(len(deg_list))

# Performing analysis
for i in deg_list:
    # EPE 
    A.Set_up_data(deg = i)
    A.PlaneOLSReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap 
    err[i], bias[i], var[i] = A.Bootstrapper()

    # Cross Validation 
    MSE_est[i] = A.Cross_Validationer()

# Plotting results 
# EPE 
plt.figure()
plt.plot(deg_list, EPE[:, 0])
plt.plot(deg_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel("Polynomial degree")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("OLS: Error plot")

# Bootstrap 
plt.figure()
plt.plot(deg_list, err)
plt.plot(deg_list, bias)
plt.plot(deg_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel("Polynomial degree")
plt.title("OLS: Bootstrap")

# Cross Validation 
plt.figure()
plt.plot(deg_list, MSE_est)
plt.plot(deg_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel("Polynomial degree")
plt.ylabel("MSE")
plt.title("OLS: Cross Validation")


plt.show()
"""

terrain = imread("n59_e010_1arc_v3.tif")
# Reducing the area we look at to make the analysis manageable
z = terrain[0:int(terrain.shape[1] / 5), 0:int(terrain.shape[1] / 5)]

# Show data image
plt.figure()
plt.title("Terrain of Finnemarka")
plt.imshow(z, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")

x = np.linspace(0, z.shape[1], z.shape[1])
y = np.linspace(0, z.shape[0], z.shape[0])
x, y = np.meshgrid(x, y)


A = Data_reg(x, y, z) # Initialising the class
A.Set_up_data(deg = 20)
A.PlaneOLSReg()
z_pred = A.Pred_data(A.X)

plt.figure()
plt.title("Pred of Finnemarka")
plt.imshow(z_pred, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")

plt.show()