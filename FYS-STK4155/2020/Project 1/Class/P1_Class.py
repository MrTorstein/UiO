# Imports
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
from sklearn.utils import resample
from sklearn.linear_model import Lasso
from sklearn.preprocessing import scale
from sklearn.model_selection import train_test_split, KFold

class Data_reg():
    """
    Takes x, y, z and can calculate Linear regression of 2D plane dataset
    """

    def __init__(self, x, y, z):
        # Initiate data and scale #
        self.x = scale(x, axis = 1)
        self.y = scale(y, axis = 0)
        self.z = z
    
    def Set_up_data(self, deg = 5, lam = 0):
        # Set up data matrix X #
        self.Lambda = lam
        self.X = np.zeros((len(self.z), np.sum(range(deg + 2))))
        for i in range(deg + 1):
            for j in range(i + 1):
                self.X[:, int(np.sum(range(i + 1))) + j] = self.x[0] ** (i - j) * self.y[:, 0] ** j

    def _Split_data(self, ratio = 0.2):
        # Split the Data into test and training groups #
        self.X_train, self.X_test, self.z_train, self.z_test = train_test_split(self.X, self.z, test_size = ratio)

    def PlaneOLSReg(self, X_train = None, z_train = None, ratio = 0.2):
        # Train the OLS model #
        self._Split_data(ratio)

        if X_train is None:
            X_train = self.X_train
        if z_train is None:
            z_train = self.z_train

        self.beta = np.linalg.pinv(X_train.T.dot(X_train)).dot(X_train.T).dot(z_train)
        
        return self.beta
    
    def PlaneRidgeReg(self, X_train = None, z_train = None, ratio = 0.2):
        # Train the Ridge model #
        self._Split_data(ratio)

        if X_train is None:
            X_train = self.X_train
        if z_train is None:
            z_train = self.z_train

        self.beta = np.linalg.pinv(X_train.T.dot(X_train) + self.Lambda * np.identity(len(X_train[0]))).dot(X_train.T).dot(z_train)
        
        return self.beta
    
    def PlaneLassoReg(self, X_train = None, z_train = None, ratio = 0.2):
        # Train the Lasso model #
        self._Split_data(ratio)

        if X_train is None:
            X_train = self.X_train
        if z_train is None:
            z_train = self.z_train

        self.beta = Lasso(alpha = self.Lambda).fit(X_train, z_train)
        
        return self.beta

    def Pred_data(self, X, method = None):
        # Use a specified method to predict data #
        if method == "Lasso":
            z_pred = self.beta.predict(X)
        else:
            z_pred = X @ self.beta
        return z_pred
    
    def Plot_fitted_and_data(self, figname = "Figure99.png"):
        # Plot the resulting fit beside the original surface
        fig = plt.figure()

        ax = fig.add_subplot(1, 2, 1, projection = '3d')
        surf = ax.plot_surface(self.x, self.y, self.z, cmap = cm.viridis, linewidth = 0, antialiased = False)
        fig.colorbar(surf, shrink = 0.5, aspect = 5)
        plt.title('Function')
        plt.xlabel("x")
        plt.ylabel("y")

        ax = fig.add_subplot(1, 2, 2, projection = '3d')
        surf = ax.plot_surface(self.x, self.y, self.Pred_data(self.X), cmap = cm.viridis, linewidth = 0, antialiased = False)
        fig.colorbar(surf, shrink = 0.5, aspect = 5)
        plt.title('Fitted Function')
        plt.xlabel("x")
        plt.ylabel("y")

        plt.savefig(figname)
        plt.show()

    def EstPredErr(self, z_pred, z):
        # Calculate MSE as the estimated prediction error
        if len(z_pred.shape) > 2:
            liste = np.zeros(len(z_pred[:, :]))
            for i in range(len(liste)):
                liste[i] = np.mean((z - z_pred[:, :, i]) ** 2)
            return liste
        else:
            return np.mean((z - z_pred) ** 2)

    def R2(self, z_pred, z):
        # Calcualte the R^2 value for the preditions
        if len(z_pred.shape) > 2:
            liste = np.zeros(len(z_pred[:, :]))
            for i in range(len(liste)):
                liste[i] = 1 - np.sum((z - z_pred[:, :, i]) ** 2) / np.sum((z - np.mean(z_pred[:, :, i])) ** 2)
            return liste
        else:
            return 1 - np.sum((z - z_pred) ** 2) / np.sum((z - np.mean(z_pred)) ** 2)

    def Bootstrapper(self, n_boots = 10 ** 2, method = "OLS"):
        # Perform the bootstrap analysis #
        if method == "OLS":
            reg_func = self.PlaneOLSReg
        elif method == "Ridge":
            reg_func = self.PlaneRidgeReg
        elif method == "Lasso":
            reg_func = self.PlaneLassoReg
        else:
            print("method has to be specified")
        
        self._Split_data()
        n = len(self.X_test)
        pred = np.zeros((self.z_test.shape[0], self.z_test.shape[1], n_boots))
        
        for i in range(n_boots):
            X_, z_ = resample(self.X_train, self.z_train)
            reg_func(X_train = X_, z_train = z_)
            pred[:, :, i] = self.Pred_data(self.X_test, method)
        
        err = np.mean(np.mean((pred - self.z_test[:, :, np.newaxis]) ** 2, axis = 2))
        bias = np.mean((self.z_test - np.mean(pred, axis = 2)) ** 2)
        var = np.mean(np.var(pred, axis = 2))

        return err, bias, var
    
    def Cross_Validationer(self, k = 5, method = "OLS"):
        # Performing the Cross Validation #
        if method == "OLS":
            reg_func = self.PlaneOLSReg
        elif method == "Ridge":
            reg_func = self.PlaneRidgeReg
        elif method == "Lasso":
            reg_func = self.PlaneLassoReg
        else:
            print("method has to be specified")
        
        kf = KFold(n_splits = k, shuffle = True)
        e_K = np.zeros(k)

        j = 0
        for train_ind, test_ind in kf.split(self.X):
            X_train, z_train, X_test, z_test = self.X[train_ind], self.z[train_ind], self.X[test_ind], self.z[test_ind]
            reg_func(X_train = X_train, z_train = z_train)
            z_test_pred = self.Pred_data(X_test, method)

            e_K[j] = np.mean((z_test - z_test_pred) ** 2)
            j += 1
        
        MSE_est = np.mean(e_K)
        return MSE_est
