""" This takes a long time to run, but suit yourselves """
import sys
sys.path.append("../")

import numpy as np
from random import random
from imageio import imread
import matplotlib.pyplot as plt
from Class.P1_Class import Data_reg

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

z = FrankeFunction(x, y) + np.random.normal(0, 1, x.shape) # Add noise

# Perform OLS for 1st-5th degree and plotting the result
deg_list = range(6)

A = Data_reg(x, y, z)
R_list = []
for i in deg_list:
    A.Set_up_data(deg = i)
    A.PlaneOLSReg()
    R_list.append(A.R2(A.Pred_data(A.X_test), A.z_test)) # Calculating the R^2 value
A.Plot_fitted_and_data(figname = "Figure10.png", Block = False)

plt.figure()
plt.title(r"$R^2$ for OLS")
plt.xlabel("Polynomial degree")
plt.ylabel(r"$R^2$")
plt.plot(deg_list, R_list)
plt.savefig("Figure00.png")
plt.show()


# OLS analysis for Franke function #
deg_list = range(30)

# EPE #
EPE = np.zeros((len(deg_list), 2))
for i in deg_list: # Loop over degrees and calculating MSE for training and test dataset
    A.Set_up_data(deg = i)
    A.PlaneOLSReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

# Plotting result
plt.figure()
plt.plot(deg_list, EPE[:, 0])
plt.plot(deg_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel("Polynomial degree")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("OLS: Error plot")
plt.savefig("Figure01.png")


# Bootstrap #
err = np.zeros(len(deg_list))
bias = np.zeros(len(deg_list))
var = np.zeros(len(deg_list))
for i in deg_list: # Performing Bootstrap over degrees for OLS
    A.Set_up_data(deg = i)
    err[i], bias[i], var[i] = A.Bootstrapper()

# Plotting result
plt.figure()
plt.plot(deg_list, err)
plt.plot(deg_list, bias)
plt.plot(deg_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel("Polynomial degree")
plt.title("OLS: Bootstrap")
plt.savefig("Figure02.png")


# Cross Validation #
MSE_est = np.zeros(len(deg_list))

for i in deg_list: # CV over degree for OLS
    A.Set_up_data(deg = i)
    MSE_est[i] = A.Cross_Validationer()

# Plotting
plt.figure()
plt.plot(deg_list, MSE_est)
plt.plot(deg_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel("Polynomial degree")
plt.ylabel("MSE")
plt.title("OLS: Cross Validation")
plt.savefig("Figure03.png")

plt.show()


# Ridge analysis. The same as OLS just with lambda instead of polydegree #
lam_list = np.linspace(-10, 5, 100)

# EPE #
EPE = np.zeros((len(lam_list), 2))
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    A.PlaneRidgeReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])
plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Ridge: Error plot")
plt.savefig("Figure04.png")


# Bootstrap #
err = np.zeros(len(lam_list))
bias = np.zeros(len(lam_list))
var = np.zeros(len(lam_list))
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Ridge")

plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.title("Ridge: Bootstrap")
plt.savefig("Figure05.png")


# Cross Validation #
MSE_est = np.zeros(len(lam_list))

for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    MSE_est[i] = A.Cross_Validationer(method = "Ridge")


plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("MSE")
plt.title("Ridge: Cross Validation")
plt.savefig("Figure06.png")

plt.show()


# Lasso analysis like what we did Ridge above #

# EPE #
EPE = np.zeros((len(lam_list), 2))
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    A.PlaneLassoReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train, "Lasso"), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test, "Lasso"), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])
plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Lasso: Error plot")
plt.savefig("Figure07.png")


# Bootstrap #
err = np.zeros(len(lam_list))
bias = np.zeros(len(lam_list))
var = np.zeros(len(lam_list))
for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Lasso")

plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.title("Lasso: Bootstrap")
plt.savefig("Figure08.png")


# Cross Validation #
MSE_est = np.zeros(len(lam_list))

for i in range(len(lam_list)):
    A.Set_up_data(lam = lam_list[i])
    MSE_est[i] = A.Cross_Validationer(method = "Lasso")


plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("MSE")
plt.title("Lasso: Cross Validation")
plt.savefig("Figure09.png")

plt.show()


# Performing the analysis with real data
# Load Data of Oslo fjord
terrain = imread("../Data/n59_e010_1arc_v3.tif")
# Reducing the area we look at to make the analysis manageable
z = terrain[0:int(terrain.shape[1] / 5), 0:int(terrain.shape[1] / 5)]

# Show data image
plt.figure()
plt.title("Terrain of Finnemarka")
plt.imshow(z, cmap = "gist_earth")
plt.colorbar().set_label("meters above sea level [m]")
plt.xlabel("x")
plt.ylabel("y")
plt.savefig("Figure20.png")
plt.show()

x = np.linspace(0, z.shape[1], z.shape[1])
y = np.linspace(0, z.shape[0], z.shape[0])
x, y = np.meshgrid(x, y)


A = Data_reg(x, y, z) # Initialising the class

# OLS Analysis like for the franke function #
deg_list = range(30)

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
plt.savefig("Figure11.png")

# Bootstrap 
plt.figure()
plt.plot(deg_list, err)
plt.plot(deg_list, bias)
plt.plot(deg_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel("Polynomial degree")
plt.title("OLS: Bootstrap")
plt.savefig("Figure12.png")

# Cross Validation 
plt.figure()
plt.plot(deg_list, MSE_est)
plt.plot(deg_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel("Polynomial degree")
plt.ylabel("MSE")
plt.title("OLS: Cross Validation")
plt.savefig("Figure13.png")


plt.show()

# Ridge analysis #

# Defining arrays
EPE = np.zeros((len(lam_list), 2))
err = np.zeros(len(lam_list))
bias = np.zeros(len(lam_list))
var = np.zeros(len(lam_list))
MSE_est = np.zeros(len(lam_list))

# Performing analysis
for i in range(len(lam_list)):
    # EPE
    A.Set_up_data(lam = lam_list[i])
    A.PlaneRidgeReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Ridge")

    # Cross Validation
    MSE_est[i] = A.Cross_Validationer(method = "Ridge")


# Plotting results 
# EPE 
plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Ridge: Error plot")
plt.savefig("Figure14.png")

# Bootstrap
plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.title("Ridge: Bootstrap")
plt.savefig("Figure15.png")


# Cross Validation
plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Ridge")
plt.ylabel("MSE")
plt.title("Ridge: Cross Validation")
plt.savefig("Figure16.png")

plt.show()

# Lasso analysis #

# EPE #
# Performing analysis
for i in range(len(lam_list)):
    # EPE
    A.Set_up_data(lam = lam_list[i])
    A.PlaneLassoReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train, "Lasso"), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test, "Lasso"), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Lasso")

    # Cross Validation
    MSE_est[i] = A.Cross_Validationer(method = "Lasso")

# Plotting results 
# EPE 
plt.figure()
plt.plot(lam_list, EPE[:, 0])
plt.plot(lam_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Lasso: Error plot")
plt.savefig("Figure17.png")

# Bootstrap #
plt.figure()
plt.plot(lam_list, err)
plt.plot(lam_list, bias)
plt.plot(lam_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.title("Lasso: Bootstrap")
plt.savefig("Figure18.png")


# Cross Validation #
plt.figure()
plt.plot(lam_list, MSE_est)
plt.plot(lam_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel(r"$\lambda$ from Lasso")
plt.ylabel("MSE")
plt.title("Lasso: Cross Validation")
plt.savefig("Figure19.png")


plt.show()

# Performing analysis of best fit models

# OLS Best fit polynomial#
# Defining arrays
EPE = np.zeros((len(deg_list), 2))
MSE_est = np.zeros(len(deg_list))
err = np.zeros(len(deg_list))
bias = np.zeros(len(deg_list))
var = np.zeros(len(deg_list))

# Performing analysis
# EPE 
A.Set_up_data(deg = 20)
A.PlaneOLSReg()
temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
EPE_OLS = np.array([temp_MSE_train, temp_MSE_test])

# Bootstrap 
err_OLS, _, _ = A.Bootstrapper()

# Cross Validation 
MSE_est_OLS = A.Cross_Validationer()


# Ridge degree analysis #
# Performing analysis
for i in deg_list:
    # EPE
    A.Set_up_data(deg = i, lam = -0.5)
    A.PlaneRidgeReg()
    temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
    temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
    EPE[i] = np.array([temp_MSE_train, temp_MSE_test])

    # Bootstrap
    err[i], bias[i], var[i] = A.Bootstrapper(method = "Ridge")

    # Cross Validation
    MSE_est[i] = A.Cross_Validationer(method = "Ridge")


# Plotting results 
# EPE 
plt.figure()
plt.plot(deg_list, EPE[:, 0])
plt.plot(deg_list, EPE[:, 1])

plt.legend(["train", "test"])
plt.xlabel("Polynomial degree")
plt.ylabel("Estimated prediction error")
plt.yscale("log")
plt.title("Ridge: Error plot")
plt.savefig("Figure21.png")

# Bootstrap
plt.figure()
plt.plot(deg_list, err)
plt.plot(deg_list, bias)
plt.plot(deg_list, var)

plt.yscale("log")
plt.legend(["err", "bias", "var"])
plt.xlabel("Polynomial degree")
plt.title("Ridge: Bootstrap")
plt.savefig("Figure22.png")


# Cross Validation
plt.figure()
plt.plot(deg_list, MSE_est)
plt.plot(deg_list, err)

plt.yscale("log")
plt.legend(["MSE", "Bootstrap_MSE"])
plt.xlabel("Polynomial degree")
plt.ylabel("MSE")
plt.title("Ridge: Cross Validation")
plt.savefig("Figure23.png")

plt.show()


# Best fit Ridge model
# EPE
A.Set_up_data(deg = 8, lam = -0.5)
A.PlaneRidgeReg()
temp_MSE_train = A.EstPredErr(A.Pred_data(A.X_train), A.z_train)
temp_MSE_test = A.EstPredErr(A.Pred_data(A.X_test), A.z_test)
EPE_Ridge = np.array([temp_MSE_train, temp_MSE_test])

# Bootstrap
err_Ridge, _, _ = A.Bootstrapper(method = "Ridge")

# Cross Validation
MSE_est_Ridge = A.Cross_Validationer(method = "Ridge")

# Printing the resulting table for best fit models
print("        OLS   Ridge")
print("MSE =", EPE_OLS[0], EPE_Ridge[0])
print("BS MSE =", err_OLS, err_Ridge)
print("CV MSE =", MSE_est_OLS, MSE_est_Ridge)