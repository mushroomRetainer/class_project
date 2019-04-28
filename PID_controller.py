
class PID:
    '''
    equations: http://apmonitor.com/pdc/index.php/Main/ProportionalIntegralDerivative
    '''
    
    def __init__(self, SP, K_c, tau_I = None, tau_D = None, u_bias = 0, u_lb = None, u_ub = None): 
        
        # store setpoint and bias
        self.SP = SP
        self.u_bias = u_bias
        
        # store tuning variables
        self.K_c = K_c
        self.tau_I = tau_I
        self.tau_D = tau_D
        
        # store clamping limits
        self.u_lb = u_lb
        self.u_ub = u_ub
        
        # initalize the error for integration
        self.error_integral = 0
        
        # initialize the previous_PV for the derivative term
        self.prev_PV = None
        
    def update_SP(self, new_SP):
        self.SP = new_SP
        
    def output(self, PV, dt):
        
        # calculate/update all error values
        error = self.SP - PV
        self.error_integral += error * dt
        # can't calculate derivative until second timestep
        if self.prev_PV is not None:
            dPV_dt = (PV - self.prev_PV) / dt
        else:
            dPV_dt = 0
            
        # calculate u
        u = self.u_bias  + self.K_c * error 
        if self.tau_I is not None:
            u += self.K_c / self.tau_I * self.error_integral
        if self.tau_D is not None:
            u += -self.K_c * self.tau_D * dPV_dt
        
        # store previous PV for next derivative calculation
        self.prev_PV = PV
        
        # anti-reset windup
        if self.u_ub is not None and u > self.u_ub:
            self.error_integral -= error * dt
            u = self.u_ub
        elif self.u_lb is not None and u < self.u_lb:
            self.error_integral -= error * dt
            u = self.u_lb
        
        return u
    
    
    
'''
simple test of the PID controller
system roughly equates to a tank with a hole in it that you can pump water into or out of (plus a disturbance that randomly adds or removes more)
controls the water level
'''
    
if __name__ == "__main__":
    
    import numpy as np
    
    np.random.seed(1234)
    
    def process(u, PV, dt, noise = False):
        if noise:
            noise_dt = (2*np.random.rand() - 1) * 0.25 # adds a random value between [-0.25,0.25]
            return max(0, PV + (u - 0.1 * PV**(1/2) + noise_dt) * dt) # can't go negative
        else:
            return max(0, PV + (u - 0.1 * PV**(1/2)) * dt) # can't go negative
    
    K_c = 0.75
    tau_I = 2
    tau_D = None
    u_lb = -1
    u_ub = 3
    noise = True
    
    PV_0 = 5
    
    SP = np.ones(101)
    SP[:] = 10
    SP[25:] = 5
    SP[50:] = 0
    SP[75:] = 20
    
    time = np.linspace(0,100,101)
    u = np.zeros(101)
    PV = np.zeros(101)
    PV[0] = PV_0
    
    pid = PID(SP[0], K_c, tau_I, tau_D, u_lb=u_lb, u_ub=u_ub)
    
    for i in range(1,len(time)):
        pid.update_SP(SP[i])
        dt = time[i] - time[i-1]
        PV[i] = process(u[i-1], PV[i-1], dt, noise=noise)
        u[i] = pid.output(PV[i], dt)
    
    import matplotlib.pyplot as plt
    
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(time, PV, label='PV')
    plt.plot(time, SP, label='SP')
    plt.legend()
    
    plt.subplot(2,1,2)
    plt.plot(time, u, label='u')
    plt.legend()
    plt.xlabel('time')