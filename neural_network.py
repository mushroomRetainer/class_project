
import numpy as np

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
        
    def get_biases(self):
        if self.status == 'off':
            return np.zeros(self.num_fans)
        elif self.status == 'learn':
            return np.random.rand(self.num_fans) * (self.upper_bound - self.lower_bound) + self.lower_bound
        else:
            return None # TODO: add this in once network is trained
        
