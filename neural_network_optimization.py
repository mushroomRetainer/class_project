
import numpy as np
from pyswarm import pso
from sklearn.cross_validation import train_test_split
from sklearn.linear_model import LinearRegression

class Neural_netork:
    '''
    status options: 'on', 'off', 'learn'
        on = returns a trained neural network set of biases
        off = returns all zeros
        learn = returns random biases within given bounds
    '''
    def __init__(self, status, num_fans, upper_bound = 1, lower_bound = -1, random_seed = 12345):
        self.status = status
        self.num_fans = num_fans
        self.upper_bound = upper_bound
        self.lower_bound = lower_bound
        if random_seed is not None:
            np.random.seed(random_seed)
        
    def get_biases(self, basin_temp__F = None, T_amb__F = None, RH_percent = None, T_hotbed__F = None):
        if self.status == 'off':
            return np.zeros(self.num_fans)
        elif self.status == 'learn':
            return np.random.rand(self.num_fans) * (self.upper_bound - self.lower_bound) + self.lower_bound
        else:
            
            def obj_fun(biases):
                x = []
                x.append(basin_temp__F)
                for bias in biases:
                    x.append(bias)
                x.append(T_amb__F)
#                x.append(RH_percent) # TODO need to make sure the training set has all the data
                x.append(T_hotbed__F)
                x = np.array(x).reshape(1, -1)
                return self.linear_model.predict(x)[0]
            
            optimal_biases, optimal_power = optimize(obj_fun, [-1]*12, [1]*12)
            return optimal_biases
            
    # tutorial: https://bigdata-madesimple.com/how-to-run-linear-regression-in-python-scikit-learn/
    def train(self, X_data, Y_data):
        print('Training Model')
        X_train, X_test, Y_train, Y_test = train_test_split(X_data, Y_data)
        self.linear_model = LinearRegression()
        self.linear_model.fit(X_train, Y_train)
        print('Train score:',self.linear_model.score(X_train, Y_train))
        print('Test score:',self.linear_model.score(X_test, Y_test))
#        pass

def optimize(obj_fun, lb, ub):
    xopt, fopt = pso(obj_fun, lb, ub)
#    print('found PSO solution:',xopt,fopt)
    return xopt, fopt




# PSO example: https://pythonhosted.org/pyswarm/

#def banana(x):
#    x1 = x[0]
#    x2 = x[1]
#    return x1**4 - 2*x2*x1**2 + x2**2 + x1**2 - 2*x1 + 5
#
#def con(x):
#    x1 = x[0]
#    x2 = x[1]
#    return [-(x1 + 0.25)**2 + 0.75*x2]
#
#lb = [-3, -1]
#ub = [2, 6]
#
#xopt, fopt = pso(banana, lb, ub, f_ieqcons=con)
    