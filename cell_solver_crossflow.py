
import numpy as np
import matplotlib.pyplot as plt
from util import Util
from parameters import HalfCell_crossflow_parameters


def m_flow_water_i__kg_s(C_water_i__kg_m3, halfcell_parameters):
    return C_water_i__kg_m3 * halfcell_parameters.v_water_absolute__m_s * halfcell_parameters.discretized_area_waterside__m2

def H_flow_water_i__W(C_water_i__kg_m3, T_water_i__K, halfcell_parameters, util):
    m_i__kg_s = m_flow_water_i__kg_s(C_water_i__kg_m3, halfcell_parameters)
    Cp__J_kg_K = halfcell_parameters.water_properties.liquid_heat_capacity__J_kg_K(T_water_i__K)
    return m_i__kg_s * Cp__J_kg_K * (T_water_i__K -  util.T_ref__K)

def m_flow_vapor_i__kg_s(C_vapor_i__kg_m3, halfcell_parameters):
    return C_vapor_i__kg_m3 * halfcell_parameters.v_air_cell__m_s * halfcell_parameters.discretized_area_vaporside__m2

def H_flow_vapor_i__W(C_vapor_i__kg_m3, T_air_i__K, halfcell_parameters, util):
    v_flow_air__m3_s = halfcell_parameters.v_air_cell__m_s * halfcell_parameters.discretized_area_vaporside__m2
    
    mass_fraction = util.mass_fraction_vapor(C_vapor_i__kg_m3, T_air_i__K, halfcell_parameters.P_amb__Pa)
    air_density__kg_m3 = halfcell_parameters.moist_air_properties.density__kg_m3(halfcell_parameters.P_amb__Pa, T_air_i__K, mass_fraction)
    Cp__J_kg_K = halfcell_parameters.moist_air_properties.heat_capacity(T_air_i__K, mass_fraction)
    return air_density__kg_m3 * v_flow_air__m3_s * Cp__J_kg_K * (T_air_i__K -  util.T_ref__K)

def m_flow_evap_i__kg_s(C_water_i__kg_m3, C_vapor_i__kg_m3, T_water_i__K, T_air_i__K, halfcell_parameters, util):
    '''
    calculates the mass flowrate of the evaporation for a given discretized point
    this is the mass source term for the two different streams (air and water)
    C_water_i is the (kg water)/(m^3 total volume) of water for the discretized point
    C_water_i is the (kg vapor)/(m^3 total volume) of vapor for the discretized point
    T_water_i__K is the temperature of water for the discretized point
    T_air_i__K is the temperature of the moist air for the discretized point
    '''
    
    SA_droplet__m2 = util.get_sphere_SA(halfcell_parameters.droplet_diameter__m/2)
    Vol_droplet__m3 = util.get_sphere_Volume(halfcell_parameters.droplet_diameter__m/2)
    density_water_i__kg_m3 = halfcell_parameters.water_properties.liquid_density__kg_m3(T_water_i__K)
    
    # this basically divides all water into droplets and multiiplies the number of droplets by the surface area of one droplet
    # we will need this after we calculate the mass flux coming off a single water droplet
    num_droplets = C_water_i__kg_m3 * halfcell_parameters.discretized_volume__m3 / (Vol_droplet__m3 * density_water_i__kg_m3)
    SA_water_total__m2 = SA_droplet__m2 * num_droplets
    
    # calculate temperature dependant moist air properties
    D_AB__m2_s = util.D_AB_vapor_air__m2_s(T_air_i__K, halfcell_parameters.P_amb__Pa)
    mass_fraction_vapor_i = util.mass_fraction_vapor(C_vapor_i__kg_m3, T_air_i__K, halfcell_parameters.P_amb__Pa)
    air_kinematic_viscocity__m2_s = halfcell_parameters.moist_air_properties.kinematic_viscocity__m2_s(T_air_i__K, mass_fraction_vapor_i)
    
    # nondimensional numbers
    Re_D = util.Re(halfcell_parameters.v_water_airspeed__m_s, halfcell_parameters.droplet_diameter__m, air_kinematic_viscocity__m2_s)
    Sc_air = util.Sc(air_kinematic_viscocity__m2_s, D_AB__m2_s) 
    Sh_air = util.Nu_correlation(Re_D, Sc_air) 
    
    # find the mass transfer coefficient, h_m, using the heat transfer correlation
    h_m__m_s = util.h_mass_transfer(Sh_air, halfcell_parameters.droplet_diameter__m, D_AB__m2_s) 
    
    # we need to caluculate the driving force, which is mostly just the concentration gradient
    # find the concentration of vapor at the surface of the water drop, using vapor pressures
    # use ideal gas law to solve for n/V, then convert with molecular weight to kg/m^3
    C_vapor_droplet_surface_i__kg_m3 = halfcell_parameters.water_properties.vapor_pressure__Pa(T_water_i__K) / (util.R * T_water_i__K) * util.MW_water__kg_mol
    
    m_flow_evap_i = SA_water_total__m2 * h_m__m_s * (C_vapor_droplet_surface_i__kg_m3 - C_vapor_i__kg_m3)
    
    return m_flow_evap_i

def Enthalpy_convective_i__W(C_water_i__kg_m3, C_vapor_i__kg_m3, T_water_i__K, T_air_i__K, halfcell_parameters, util):
    '''
    Calculates the enthalphy flux due to convective heat transfer coming from the water droplets for a given volume
    '''
    SA_droplet__m2 = util.get_sphere_SA(halfcell_parameters.droplet_diameter__m/2)
    Vol_droplet__m3 = util.get_sphere_Volume(halfcell_parameters.droplet_diameter__m/2)
    density_water_i__kg_m3 = halfcell_parameters.water_properties.liquid_density__kg_m3(T_water_i__K)
    
    # this basically divides all water into droplets and multiiplies the number of droplets by the surface area of one droplet
    # we will need this after we calculate the mass flux coming off a single water droplet
    num_droplets = C_water_i__kg_m3 * halfcell_parameters.discretized_volume__m3 / (Vol_droplet__m3 * density_water_i__kg_m3)
    SA_water_total__m2 = SA_droplet__m2 * num_droplets
    
    # calculate temperature dependant moist air properties
    k_thermal_conductivity = halfcell_parameters.moist_air_properties.k_thermal_cond__W_m_K(T_air_i__K)
    mass_fraction_vapor_i = util.mass_fraction_vapor(C_vapor_i__kg_m3, T_air_i__K, halfcell_parameters.P_amb__Pa)
    air_kinematic_viscocity__m2_s = halfcell_parameters.moist_air_properties.kinematic_viscocity__m2_s(T_air_i__K, mass_fraction_vapor_i)
    air_dynamic_viscocity__kg_m_s = halfcell_parameters.moist_air_properties.dynamic_viscocity__kg_m_s(T_air_i__K, mass_fraction_vapor_i)
    
    # nondimensional numbers
    Re_D = util.Re(halfcell_parameters.v_water_airspeed__m_s, halfcell_parameters.droplet_diameter__m, air_kinematic_viscocity__m2_s)
    Cp_air__J_kg_K = halfcell_parameters.moist_air_properties.heat_capacity(T_air_i__K, mass_fraction_vapor_i)
    Pr_air = util.Pr(Cp_air__J_kg_K, air_dynamic_viscocity__kg_m_s, k_thermal_conductivity) 
    Nu_air = util.Nu_correlation(Re_D, Pr_air) 
    
    # find the heat transfer coefficient, h_m, using the heat transfer correlation
    h__W_m2_K = util.h_heat_transfer(Nu_air, halfcell_parameters.droplet_diameter__m, k_thermal_conductivity)
    
    H_convective_i__W = SA_water_total__m2 * h__W_m2_K * (T_water_i__K - T_air_i__K)
    
    return H_convective_i__W

def solve_halfcell(halfcell_parameters, util):
    
    n_x = halfcell_parameters.n_x
    n_z = halfcell_parameters.n_z
    C_water_0__kg_m3 = halfcell_parameters.C_water_0__kg_m3
    C_vapor_in__kg_m3 = halfcell_parameters.C_vapor_in__kg_m3
    T_water_0__K = halfcell_parameters.T_water_0__K
    T_air_in__K = halfcell_parameters.T_air_in__K
    water_properties = halfcell_parameters.water_properties
    moist_air_properties = halfcell_parameters.moist_air_properties
    
    discretized_area_waterside__m2 = halfcell_parameters.discretized_area_waterside__m2
    discretized_area_vaporside__m2 = halfcell_parameters.discretized_area_vaporside__m2
    v_water_absolute__m_s = halfcell_parameters.v_water_absolute__m_s
    v_air_cell__m_s = halfcell_parameters.v_air_cell__m_s
    
    vol_flowrate_water_discretized__m3_s = v_water_absolute__m_s * discretized_area_waterside__m2
    vol_flowrate_vapor_discretized__m3_s = v_air_cell__m_s * discretized_area_vaporside__m2
    
    # handle the case where the cell is turned off (basically, allow the water to pass through unchanged):
    if halfcell_parameters.turned_off:
        avg_Temperature = T_water_0__K
        m_flow_water_out__kg_s =  m_flow_water_i__kg_s(C_water_0__kg_m3, halfcell_parameters) * n_x
        C_water_matrix__kg_m3 = np.ones([n_z, n_x]) * C_water_0__kg_m3
        C_vapor_matrix__kg_m3 = np.ones([n_z, n_x]) * C_vapor_in__kg_m3
        T_water_i_matrix__K = np.ones([n_z, n_x]) * T_water_0__K
        T_air_i_matrix__K = np.ones([n_z, n_x]) * T_air_in__K
        return avg_Temperature, m_flow_water_out__kg_s, C_water_matrix__kg_m3, C_vapor_matrix__kg_m3, T_water_i_matrix__K, T_air_i_matrix__K
    
    # take average at end to get the finite volume values. Not going to mess with the simultaneous solution.
    C_water_matrix__kg_m3 = np.zeros([n_z + 1, n_x + 1])
    C_vapor_matrix__kg_m3 = np.zeros([n_z + 1, n_x + 1])
    T_water_i_matrix__K = np.zeros([n_z + 1, n_x + 1])
    T_air_i_matrix__K = np.zeros([n_z + 1, n_x + 1])
    
    for z_i in range(n_z + 1):
        for x_i in range(n_x + 1):
            if z_i > 0:
                entering_C_water__kg_m3 = C_water_matrix__kg_m3[z_i - 1, x_i]
                entering_T_water__K = T_water_i_matrix__K[z_i - 1, x_i]
            else:
                entering_C_water__kg_m3 = C_water_0__kg_m3
                entering_T_water__K = T_water_0__K
                    
            if x_i > 0:
                entering_C_vapor__kg_m3 = C_vapor_matrix__kg_m3[z_i, x_i - 1]
                entering_T_air__K = T_air_i_matrix__K[z_i, x_i - 1]
            else:
                entering_C_vapor__kg_m3 = C_vapor_in__kg_m3
                entering_T_air__K = T_air_in__K
            
            water_in_top__kg_s = m_flow_water_i__kg_s(entering_C_water__kg_m3, halfcell_parameters)
            vapor_in_left__kg_s = m_flow_vapor_i__kg_s(entering_C_vapor__kg_m3, halfcell_parameters)
            water_evap__kg_s = m_flow_evap_i__kg_s(entering_C_water__kg_m3, 
                                                   entering_C_vapor__kg_m3, 
                                                   entering_T_water__K, 
                                                   entering_T_air__K, 
                                                   halfcell_parameters, util)
            enthalpy_evap__W = water_evap__kg_s * water_properties.enthalpy_vaporization__J_kg(entering_T_water__K)
            
            enthalpy_evap_convective__W = Enthalpy_convective_i__W(entering_C_water__kg_m3, 
                                                        entering_C_vapor__kg_m3, 
                                                        entering_T_water__K,
                                                        entering_T_air__K,
                                                        halfcell_parameters, util)
            
            water_out_bottom__kg_s = water_in_top__kg_s - water_evap__kg_s
            vapor_out_right__kg_s = vapor_in_left__kg_s + water_evap__kg_s
            
            enthalpy_change_water__W = -enthalpy_evap__W - enthalpy_evap_convective__W
            enthalpy_change_air__W =  enthalpy_evap_convective__W 
            
            delta_T_water__K = enthalpy_change_water__W / ( water_out_bottom__kg_s * 
                                                           water_properties.liquid_heat_capacity__J_kg_K(entering_T_water__K) )
            
            
#            vapor_mass_fraction = get_vapor_mass_fraction(entering_C_vapor__kg_m3, entering_T_air__K,  halfcell_parameters, util)
            vapor_mass_fraction = util.mass_fraction_vapor(entering_C_vapor__kg_m3, entering_T_air__K, halfcell_parameters.P_amb__Pa)
            
            # need to account for the fact that the vapor is not the same temperature as the air
            adjustedT_air__K = vapor_mass_fraction * entering_T_water__K + (1 - vapor_mass_fraction) * entering_T_air__K
            
            air_out_right__kg_s = vapor_out_right__kg_s / vapor_mass_fraction
            delta_T_air__K = enthalpy_change_air__W / ( air_out_right__kg_s * moist_air_properties.heat_capacity(entering_T_air__K,vapor_mass_fraction) )
            
            if z_i > 0:
                C_water_matrix__kg_m3[z_i, x_i] = water_out_bottom__kg_s / vol_flowrate_water_discretized__m3_s
                T_water_i_matrix__K[z_i, x_i] = entering_T_water__K + delta_T_water__K
            else:
                C_water_matrix__kg_m3[z_i, x_i] = entering_C_water__kg_m3
                T_water_i_matrix__K[z_i, x_i] = entering_T_water__K
            
            if x_i > 0:
                C_vapor_matrix__kg_m3[z_i, x_i] = vapor_out_right__kg_s / vol_flowrate_vapor_discretized__m3_s
                T_air_i_matrix__K[z_i, x_i] = adjustedT_air__K + delta_T_air__K 
            else:
                C_vapor_matrix__kg_m3[z_i, x_i] = entering_C_vapor__kg_m3
                T_air_i_matrix__K[z_i, x_i] = entering_T_air__K
                
    
    # find the averages of the adjacent points and store those in an array of size [n_z, n_x]. Need to average along both axis to get the center of the cubes
    C_water_matrix__kg_m3 = np.average([C_water_matrix__kg_m3[:-1],C_water_matrix__kg_m3[1:]],axis=0)
    C_water_matrix__kg_m3 = np.average([C_water_matrix__kg_m3[:,:-1],C_water_matrix__kg_m3[:,1:]],axis=0)
    
    C_vapor_matrix__kg_m3 = np.average([C_vapor_matrix__kg_m3[:-1],C_vapor_matrix__kg_m3[1:]],axis=0)
    C_vapor_matrix__kg_m3 = np.average([C_vapor_matrix__kg_m3[:,:-1],C_vapor_matrix__kg_m3[:,1:]],axis=0)
    
    T_water_i_matrix__K = np.average([T_water_i_matrix__K[:-1],T_water_i_matrix__K[1:]],axis=0)
    T_water_i_matrix__K = np.average([T_water_i_matrix__K[:,:-1],T_water_i_matrix__K[:,1:]],axis=0)
    
    T_air_i_matrix__K = np.average([T_air_i_matrix__K[:-1],T_air_i_matrix__K[1:]],axis=0)
    T_air_i_matrix__K = np.average([T_air_i_matrix__K[:,:-1],T_air_i_matrix__K[:,1:]],axis=0)
    
    
    # get the total mass flow of watercoming out the bottom of the cell
    m_flow_water_out__kg_s = np.sum( m_flow_water_i__kg_s(C_water_matrix__kg_m3[-1,:], halfcell_parameters) )
    # get the average temperature:
    avg_Temperature = np.average(T_water_i_matrix__K[-1,:])
    
    return avg_Temperature, m_flow_water_out__kg_s, C_water_matrix__kg_m3, C_vapor_matrix__kg_m3, T_water_i_matrix__K, T_air_i_matrix__K


# code to test the various methods here
if __name__ == '__main__':  
    plotting_on = True
    
    num_reps = 12
    
    for i in range(num_reps):
        util = Util()
        if i%2 == 0:
            power_draw_single_fan__HP = 183* 0.001
        else:
            power_draw_single_fan__HP = 183* 1.0 # 150 is "design" condition and 183 is supposedly the max, but we need about double that to get the flow rate we need for realistic temperatures
        halfcell_parameters = HalfCell_crossflow_parameters(util, power_draw_single_fan__HP, hotbed_flowrate__gal_min = 180000/24)
        
        avg_Temperature, m_flow_water_out__kg_s, C_water_matrix__kg_m3, C_vapor_matrix__kg_m3, T_water_i_matrix__K, T_air_i_matrix__K = solve_halfcell(halfcell_parameters, util)
        
        print('\nmass flowrate of water coming out of cell (kg/s):',round(m_flow_water_out__kg_s,2))
        print('Average water temperature coming out of cell (K):',round(avg_Temperature,2))
    
    if plotting_on:
        plt.figure()
        z_array = np.linspace(0,1,halfcell_parameters.n_z)
        x_array = np.linspace(0,1,halfcell_parameters.n_x)
    #    for x_i in range(int(halfcell_parameters.n_x/10)):
    #        x_i *= 10
        for x_i in range(halfcell_parameters.n_x):
            if halfcell_parameters.n_x > 20 and x_i%int(halfcell_parameters.n_x/20) != 0:
                continue # this will keep plotting to no more than about 100 lines
            plt.subplot(1, 4, 1)
            plt.title('C water')
            plt.plot(C_water_matrix__kg_m3[:,x_i],z_array, label='x = '+str(x_i))
            plt.legend()
            
            plt.subplot(1, 4, 2)
            plt.title('T water')
            plt.plot(T_water_i_matrix__K[:,x_i],z_array, label='x = '+str(x_i))
            plt.legend()
            
        plt.subplot(1, 4, 1)
        plt.gca().invert_yaxis()
        plt.xlabel('Concentration, kg/m^3')
        plt.ylabel('Fraction of Cell Height')
        plt.subplot(1, 4, 2)
        plt.gca().invert_yaxis()
        plt.xlabel('Temperature, K')
        plt.ylabel('Fraction of Cell Height')
        
        for z_i in range(halfcell_parameters.n_z):
            plt.subplot(1, 4, 3)
            plt.title('C vapor')
    #        plt.gca().invert_yaxis()
            plt.plot(x_array, C_vapor_matrix__kg_m3[z_i,:], label='z = '+str(z_i))
            plt.legend()
            
            plt.subplot(1, 4, 4)
            plt.title('T air')
    #        plt.gca().invert_yaxis()
            plt.plot(x_array, T_air_i_matrix__K[z_i,:], label='z = '+str(z_i))
            plt.legend()
            
            plt.show()
        
        plt.subplot(1, 4, 3)
        plt.ylabel('Concentration, kg/m^3')
        plt.xlabel('Fraction of Cell Width')
        plt.subplot(1, 4, 4)
        plt.ylabel('Temperature, K')
        plt.xlabel('Fraction of Cell Width')
    
    