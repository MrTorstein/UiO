import numpy as np
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split




class Data_reg():
    """
    Modified version of the class from project 1
    
    Takes x and y and can calculate Linear regression 1D dataset
    """

    def __init__(self, x, y):
        # Initialising datavariables #
        self.x = x
        self.y = y

        self.X = False
    
    def Set_up_data(self, deg = 5, lam = 0):
        # Setting up initial data matrix #
        self.Lambda = lam
        self.n_beta = deg + 1

        self.X = np.zeros((len(self.y), self.n_beta))
        for i in range(deg + 1):
            self.X[:, i] = self.x ** i

    def _Split_data(self, ratio = 0.2):
        # Spliting data in train and test #
        self.X_train, self.X_test, self.Y_train, self.Y_test = train_test_split(self.X, self.y, test_size = ratio)

    def PlaneOLSReg(self, X_train = None, Y_train = None, ratio = 0.2):
        # Performing ordinary least squares regression #
        if type(self.X) != type(False):
            self._Split_data(ratio)

        if X_train is None:
            X_train = self.X_train
        if Y_train is None:
            Y_train = self.Y_train

        self.beta = np.linalg.pinv(X_train.T.dot(X_train)).dot(X_train.T).dot(Y_train)
        
        return self.beta
    
    def PlaneRidgeReg(self, X_train = None, Y_train = None, ratio = 0.2):
        # Performing Ridge regression #
        if type(self.X) != type(False):
            self._Split_data(ratio)

        if X_train is None:
            X_train = self.X_train
        if Y_train is None:
            Y_train = self.Y_train

        self.beta = np.linalg.pinv(X_train.T.dot(X_train) + self.Lambda * np.identity(len(X_train[0]))).dot(X_train.T).dot(Y_train)
        
        return self.beta

    def Pred_data(self, X):
        # Predicting data with the regression model #
        Y_pred = X @ self.beta
        return Y_pred
    
    def EstPredErr(self, y_pred, y):
        # Calculating mean squares estimate #
        if len(y_pred.shape) > 1:
            liste = np.zeros(len(y_pred[:]))
            for i in range(len(liste)):
                liste[i] = np.mean((y - y_pred[:, i]) ** 2)
            return liste
        else:
            return np.mean((y - y_pred) ** 2)



class NN_Reg:
    """
    Neural network class to perform two category classification
    """

    def __init__(self, x, y, n_h_nodes = 10, n_c = 2, epochs = 10, batch_size = 100, eta = 0.01):
        # Initialising some parameters #

        self.x = x
        self.y = y

        self.n_h_nodes = n_h_nodes
        self.n_c = n_c

        self.epochs = epochs
        self.batch_size = batch_size
        self.eta = eta

        self.b_i = 0.01

        self.a_func = False

    def Set_up_data(self, deg = 5):
        # Setting up initial data #
        self.n_beta = deg + 1

        self.X = np.zeros((len(self.y), self.n_beta))
        for i in range(deg + 1):
            self.X[:, i] = self.x ** i
    
    def _Split_data(self, ratio = 0.2):
        # Spliting data in train and test #
        self.X_train_i, self.X_test, self.Y_train_i, self.Y_test = train_test_split(self.X, self.y, test_size = ratio)

        self.n_X = len(self.Y_train_i)
        self.iterations = self.n_X // self.batch_size

        onehot_vector = np.zeros((self.n_X, self.n_c))
        onehot_vector[range(self.n_X), self.Y_train_i.astype(int)] = 1
        self.Y_train_o = onehot_vector

        onehot_vector = np.zeros((len(self.Y_test), self.n_c))
        onehot_vector[range(len(self.Y_test)), self.Y_test.astype(int)] = 1
        self.Y_test_o = onehot_vector

    def Set_activation_function(self, func = lambda x: 1/(1 + np.exp(-x))):
        # Setting up an activation function #
        self.a_func = func

    def _Find_b_and_w(self):
        # Find initial bias and weight #
        self.w_h = np.random.randn(self.n_beta, self.n_h_nodes)
        self.b_h = np.zeros(self.n_h_nodes) + self.b_i

        self.w_o = np.random.randn(self.n_h_nodes, self.n_c)
        self.b_o = np.zeros(self.n_c) + self.b_i

    def _FF(self):
        # feed-forward function for training #
        self.Y_h = np.matmul(self.X_train, self.w_h) + self.b_h
        self.a_h = self.a_func(self.Y_h)
        self.Y_o = np.matmul(self.a_h, self.w_o) + self.b_o
        
        exp_term = np.exp(self.Y_o)
        self.P = exp_term / np.sum(exp_term, axis = 1, keepdims = True)

    def _FF_o(self, X):
        # feed-forward function for output #
        Y_h = np.matmul(X, self.w_h) + self.b_h
        a_h = self.a_func(Y_h)
        Y_o = np.matmul(a_h, self.w_o) + self.b_o

        exp_term = np.exp(Y_o)
        return 1 - exp_term / np.sum(exp_term, axis = 1, keepdims = True)

    def _Backpropagation(self):
        # perform the backpropagation after the feed forward #
        e_o = np.squeeze(self.P) - self.Y_train
        if len(np.shape(e_o)) < 2:
            e_o = e_o[:, np.newaxis]
        e_h = np.matmul(e_o, self.w_o.T) * self.a_h * (1 - self.a_h)

        self.w_o_gradient = np.matmul(self.a_h.T, e_o)
        self.b_o_gradient = np.sum(e_o, axis = 0)

        self.w_h_gradient = np.matmul(self.X_train.T, e_h)
        self.b_h_gradient = np.sum(e_h, axis = 0)

        self.w_o = self.w_o - self.eta * self.w_o_gradient
        self.b_o = self.b_o - self.eta * self.b_o_gradient
        self.w_h = self.w_h - self.eta * self.w_h_gradient
        self.b_h = self.b_h - self.eta * self.b_h_gradient

    def Pred_data(self, X = False):
        if type(X) == type(False):
            X = self.X_test

        # Uses the model to predict values #
        return np.argmax(self._FF_o(X), axis = 1)

    def Pred_probabilities(self, X = False):
        if type(X) == type(False):
            X = self.X_test
        
        # Predict values and return their probabilities #
        return self._FF_o(X)

    def Train(self, ratio = 0.2):
        # Train the model using train data #
        self._Split_data(ratio)

        if self.a_func == False:
            print("Setting activation function to the Sigmoid function")
            self.Set_activation_function()
        
        self._Find_b_and_w()

        data_indices = np.arange(self.n_X)
        
        for i in range(self.epochs):
            for j in range(self.iterations):
                # Pick datapoints with replacement
                r_d = np.random.choice(data_indices, size = self.batch_size, replace = False)

                # minibatch training data
                self.X_train = self.X_train_i[r_d]
                self.Y_train = self.Y_train_o[r_d]
                
                self._FF()
                self._Backpropagation()

    def EstPredErr(self, y_pred, y):
        if len(y_pred.shape) > 1:
            liste = np.zeros(len(y_pred[:]))
            for i in range(len(liste)):
                liste[i] = np.mean((y - y_pred[:, i]) ** 2)
            return liste
        else:
            return np.mean((y - y_pred) ** 2)