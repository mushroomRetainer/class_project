
import numpy as np


class Util:
    
    def __init__(self):
        # unit conversions
        self.Pa_per_bar = 100000
        self.Pa_per_kPa = 1000
        self.m_per_ft = 0.3048
        self.m3_s_per_CFM = 1/2118.88 #CFM = ft^3/min
        self.m3_s_per_gal_min = 1/15850.3
        self.CtoK = 273.15
        self.J_per_kJ = 1000
        self.P_standard__Pa = 101325 # = 1 atm
        self.T_standard__K = 298
        self.s_per_hr = 60**2
        self.s_per_min = 60
        self.min_per_hr = 60
        self.kw_per_hp = 0.7457
        
        # constants, material properties:
        self.MW_water__kg_mol = 18.02 / 1000 # kg/mol
        self.MW_air__kg_mol = 28.97 / 1000 # kg/mol
        self.R = 8.31451 #Pa m^3 / mol K
        self.T_ref__K = 25 + self.CtoK # 25 degC
        self.D_AB_vapor_air_298K__m2_s = .26e-4 # given at standard T and P, source: 7th edition, Principles of Heat and Mass Transfer, Incropera et. al
        
        
    def degF_to_degK(self,T_F):
        return (T_F - 32)*5/9 + self.CtoK
    
    def degC_to_degF(self,T_C):
        return T_C * 9/5 + 32
    
    def degK_to_degF(self,T_K):
        return (T_K - self.CtoK) * 9/5 + 32
    
    def deltaF_to_deltaK(deltaT_F):
        return deltaT_F*5/9
    
    
    # convert to an index from grid points in 2D
    def get_index(self, x_i, y_j, n_x, n_y):
        return y_j*n_x + x_i
        
    # return the indexes along the horizontal for a given row (default is first row)
    def get_horizontal_indexes(self, n_x, n_y, y_j = 0):
        return np.arange(y_j*n_x, (y_j+1)*n_x)
    
    # return the indexes along the vertical for a given column (default is first column)
    def get_vertical_indexes(self, n_x, n_y, x_i = 0):
        return np.arange(x_i, n_x*n_y, n_x)
 
    # convert an index to grid points in 2D
    def get_xi_yi(self,index, n_x, n_y):
        x_i = index % n_x
        y_j = int(index / n_x)
        return x_i, y_j
    
    def get_sphere_SA(self,radius):
        return 4*np.pi*(radius)**2
    
    def get_sphere_Volume(self,radius):
        return 4/3*np.pi*(radius)**3
    
    def get_circle_area(self,radius):
        return np.pi*(radius)**2
    
    def D_AB_vapor_air__m2_s(self,Temperature__K, Pressure__Pa):
        '''
        D_AB is proportional to p^(-1)*T^(3/2), assuming ideal gas behavior
        the given value is at standard T and P
        '''
        return self.D_AB_vapor_air_298K__m2_s * (Pressure__Pa / self.P_standard__Pa)**(-1) * (Temperature__K / self.T_standard__K)**(3/2)
    
    def mass_fraction_vapor(self,C_vapor_i__kg_m3, T_air_i__K, P_amb__Pa):
        C_vapor_i__mol_m3 = C_vapor_i__kg_m3 / self.MW_water__kg_mol
        C_total_i__mol_m3 = P_amb__Pa / (self.R * T_air_i__K)
        C_air_i__mol_m3 = C_total_i__mol_m3 - C_vapor_i__mol_m3
        C_air_i__kg_m3 = C_air_i__mol_m3 * self.MW_air__kg_mol
        C_total_i__kg_m3 = C_vapor_i__kg_m3 + C_air_i__kg_m3
        return C_vapor_i__kg_m3 / C_total_i__kg_m3
    
    def Re(self,velocity_inf, l_characteristic, kinematic_viscocity):
        return velocity_inf * l_characteristic / kinematic_viscocity
    
    def Sc(self, kinematic_viscocity, D_AB):
        return kinematic_viscocity / D_AB
    
    def Pr(self, C_p, dynamic_viscocity, k_thermal_conductivity):
        return C_p * dynamic_viscocity / k_thermal_conductivity
    
    def Nu_correlation(self, Re_D, Pr):
        '''
        heat transfer correlation for falling drop, evaluate properties at T_inifity
        returns Nu when given Re and Pr
        returns Sh when given Re and Sc
        '''
        return 2 + 0.6 * Re_D**(1/2) * Pr**(1/3)
    
    def h_mass_transfer(self, Sh, diameter, D_AB):
        '''
        convective mass transfer coefficient for falling drop
        note that the sphere diameter is the characteristic length
        D_AB is the binary mass transfer coefficient
        '''
        return Sh * D_AB/diameter
    
    def h_heat_transfer(self, Nu, diameter, k_fluid):
        '''
        convective heat transfer coefficient for falling drop
        note that the sphere diameter is the characteristic length
        k_fluid is the heat transfer coefficient for the fluid
        '''
        return Nu * k_fluid/diameter
    
    def flowrate_out__CFM(self, power_draw_single_fan__HP, max_airflow__CFM = 1510000, max_power_draw__HP = 182.9528*10/12, efficiency = 1): # from Jake's data as well as Mark's spec sheet (using 10/12 of the the maxpower)):
        '''
        cube rule for fans : airflow_1 / airflow_2 = power_1^3 / power_2^3
        max power is taken across a year of data, when fans are operating at 100% capacity (this was averaged in an excel spreadsheet, so it's hardcoded as an optional input above)
        max airflow is the design condition flow taken from spec sheet (also hardcoded as an optional input above)
        max power is also adjusted by 10/12 to account for the claim in the spec sheet that the tower is designed to meet demands with two cells turned off 
        '''
        
        max_airflow__CFM *= 15 #TODO: this had to be significantly increased to make this system work
        
        return max_airflow__CFM * (power_draw_single_fan__HP)**3/max_power_draw__HP**3 * efficiency # multiply by efficiency
    