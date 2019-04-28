
import numpy as np
import matplotlib.pyplot as plt
from util import Util
from parameters import Basin_parameters, HalfCell_crossflow_parameters

from cell_solver_crossflow import solve_halfcell

class Basin:
    
    def __init__(self, fan_array_init__percent, 
                 T_amb__degF = 65, Relative_humidity__percent = 20, T_hotbed__degF = 92, 
                 hotbed_flowrate__gal_min = np.ones(12) * 180000/24,
                 fan_motor_efficiencies = np.ones(12)):
        
        # set unchangeable parameters
        self.fan_motor_efficiencies_array = fan_motor_efficiencies
        self.hotbed_flowrate_array__gal_min = hotbed_flowrate__gal_min
        
        # set changeable parameters
        self.update_distrubances(T_amb__degF, 
                                 Relative_humidity__percent, 
                                 T_hotbed__degF)
        
        # get parameters:
        self.util = Util()
        self.basin_parameters = Basin_parameters(self.util)
        self.n_x = self.basin_parameters.n_x
        self.n_y = self.basin_parameters.n_y
        self.CFL = self.basin_parameters.CFL
        
        self.deltaX = self.basin_parameters.basin_length__m / (self.n_x-1)
        self.deltaY = self.basin_parameters.basin_height__m / (self.n_y-1)
        
        fan_array_init__hp = fan_array_init__percent * self.basin_parameters.fan_max__hp
        
#        print('fan_array_init__hp',fan_array_init__hp)
        
        # update (initialize) fans
        self.update_fans(fan_array_init__hp)
        
        # initalize basin Temperature profile
        basin_T_init = self.basin_parameters.T_basin_init__K * np.ones(self.n_x * self.n_y)
        basin_T_init[self.util.get_horizontal_indexes(self.n_x, self.n_y)] = self.Dirichlet_array
        self.basin_T = basin_T_init
        
        self.average_outlet_T__K = np.average(self.basin_T[self.util.get_vertical_indexes(self.n_x, self.n_y)])
        
        # initalize array for storing solutions
        self.basin_T_array = []
        self.basin_T_array.append(self.basin_T)
        self.normalized_C_water_arrays_array = []
        self.normalized_C_vapor_arrays_array = []
        self.normalized_T_water_arrays_array = []
        self.normalized_T_air_arrays_array = []
        
        self.normalized_C_water_arrays_array.append(self.normalized_C_water_arrays)
        self.normalized_C_vapor_arrays_array.append(self.normalized_C_vapor_arrays)
        self.normalized_T_water_arrays_array.append(self.normalized_T_water_arrays)
        self.normalized_T_air_arrays_array.append(self.normalized_T_air_arrays)
        
        # time counter
        self.t_i__s = 0
        
#        print('\nComplete')
        
    def update_distrubances(self, T_amb__degF, 
                            Relative_humidity__percent, 
                            T_hotbed__degF):
        '''
        these are ambient and plant conditions that are subject to change
        '''
        self.T_amb__degF = T_amb__degF
        self.Relative_humidity__percent = Relative_humidity__percent
        self.T_hotbed__degF = T_hotbed__degF
        
    def update_fans(self, fan_array__hp):
        # get inital conditions:
        ( self.Dirichlet_array, m_flow_cell_total__kg_s, 
          self.normalized_C_water_arrays, self.normalized_C_vapor_arrays, 
          self.normalized_T_water_arrays, self.normalized_T_air_arrays, 
          z_array ) = self.get_cell_info(fan_array__hp, self.n_x)
        T_avg_in = np.average(self.Dirichlet_array)
        V_flow_water_cell__m3_s = m_flow_cell_total__kg_s / self.basin_parameters.water_properties.liquid_density__kg_m3(T_avg_in)
        self.v_x_array__m_s, self.v_y_array__m_s = self.get_timestep_conditions(V_flow_water_cell__m3_s, self.n_x, self.n_y, self.basin_parameters)
        
        #calculate the delta_t for both axes and use the more restrictive one
        delta_t_x__s = self.CFL * self.deltaX / np.max( np.abs(self.v_x_array__m_s) )
        delta_t_y__s = self.CFL * self.deltaY / np.max( np.abs(self.v_y_array__m_s) )
        self.delta_t__s = np.min([delta_t_x__s,delta_t_y__s])
        
    def advance_one_timestep(self, fans_sp_array__percent, recalculate_fans = True):
      
        # determine whether or not to recalculate cell conditions (must happen after we reach the cutoff time)
        if recalculate_fans:
            fans_sp_array__hp = fans_sp_array__percent * self.basin_parameters.fan_max__hp
#            print('Time =', round(self.t_i__s / self.util.s_per_hr,3),'hours. Implementing new fan setpoints.')
            # update fans
            self.update_fans(fans_sp_array__hp)
            
        # advance the basin profile:
        self.basin_T = self.Forward_Euler(self.basin_T, self.v_x_array__m_s, self.v_y_array__m_s,
                                          self.deltaX, self.deltaY, self.delta_t__s, self.n_x, self.n_y, 
                                          self.basin_parameters.water_properties,
                                          self.Dirichlet_array)
        self.average_outlet_T__K = np.average(self.basin_T[self.util.get_vertical_indexes(self.n_x, self.n_y)])
        self.basin_T_array.append(self.basin_T)
        self.normalized_C_water_arrays_array.append(self.normalized_C_water_arrays)
        self.normalized_C_vapor_arrays_array.append(self.normalized_C_vapor_arrays)
        self.normalized_T_water_arrays_array.append(self.normalized_T_water_arrays)
        self.normalized_T_air_arrays_array.append(self.normalized_T_air_arrays)
        self.t_i__s += self.delta_t__s
        
    def get_timestep_conditions(self, V_flow_water_cell__m3_s, n_x, n_y, basin_parameters):
        
        # TODO: need to rework the math to account for cells being shut off. This assumes a linear velocity profile, which only works for uniform velocity coming out of all cells
        # should be a piecewise linear function where the slope is propoertional the to flowrate
        
        v_water_basin_in__m_s  = V_flow_water_cell__m3_s / (basin_parameters.basin_length__m * basin_parameters.basin_width__m)
        v_water_basin_out__m_s = V_flow_water_cell__m3_s / (basin_parameters.basin_height__m * basin_parameters.basin_width__m)
        v_water_basin_out__m_s *= -1 # the water is moving to the left
        
        # normal spacing
        self.v_x_array__m_s = np.linspace(v_water_basin_out__m_s,0,n_x)
        self.v_y_array__m_s = np.linspace(v_water_basin_in__m_s,0,n_y)
        
        return self.v_x_array__m_s, self.v_y_array__m_s 
    
    def dT_dt(self, T, index, deltaX, deltaY, D_function, v_x, v_y, n_x, n_y):
        x_i, y_j    = self.util.get_xi_yi(index, n_x, n_y)
        index_x_i1  = self.util.get_index(x_i + 1, y_j,     n_x, n_y)
        index_y_jm1 = self.util.get_index(x_i,     y_j - 1, n_x, n_y)
        
        # only advection upwind spacing, remember that v_x is negative and v_y is positive
        # impose Nueman boundaries of zero on 1 side
        if x_i == n_x - 1: # right
            index_x_i1 = index
        return  ( -v_x[x_i] * ( T[index_x_i1] - T[index] ) / (deltaX)
                 - v_y[y_j] * ( T[index] - T[index_y_jm1] ) / (deltaY) )
    
    def Forward_Euler(self, solution_current, v_x_array, v_y_array,deltaX,deltaY,delta_t,n_x, n_y, water_properties,Dirichlet_array):
        
        D_function = lambda T : water_properties.liquid_thermal_diffusivity(T)
        new_solution = np.zeros(len(solution_current))
        for index in range(len(solution_current)):
            x_i, y_j = self.util.get_xi_yi(index,n_x, n_y)
            
            # apply Dirichlet condition:
            if y_j == 0:
                new_solution[index] = Dirichlet_array[x_i]
            else:
                new_solution[index] = solution_current[index] + delta_t * self.dT_dt(solution_current,index,deltaX,deltaY,D_function,v_x_array, v_y_array,n_x,n_y)
            
        return new_solution
    
    
    def get_cell_info(self, fan_array__hp, n_x):
        num_cells = len(fan_array__hp)
        temperature_array = np.zeros(num_cells)
        m_flow_cell_total__kg_s = 0
        normalized_C_water_arrays = []
        normalized_C_vapor_arrays = []
        normalized_T_water_arrays = []
        normalized_T_air_arrays = []
        
        for i, fan__hp in enumerate(fan_array__hp):
#            print('Working on cell',i+1,'of',num_cells)
            
            halfcell_parameters = HalfCell_crossflow_parameters(self.util, fan__hp, 
                                                                hotbed_flowrate__gal_min = self.hotbed_flowrate_array__gal_min[i],
                                                                fan_motor_efficiency = self.fan_motor_efficiencies_array[i],
                                                                T_amb__degF = self.T_amb__degF,
                                                                Relative_humidity__percent = self.Relative_humidity__percent,
                                                                T_hotbed__degF = self.T_hotbed__degF)
            
            (avg_Temperature, m_flow_water_out__kg_s, C_water_matrix__kg_m3, 
             C_vapor_matrix__kg_m3, T_water_i_matrix__K, T_air_i_matrix__K) = solve_halfcell(halfcell_parameters, self.util)
            
            #collapse halfcell martices into arrays:
            
            C_water_in__kg_m3 = C_water_matrix__kg_m3[0,0]
            T_water_i_in__K = T_water_i_matrix__K[0,0]
            # average along x axis:
            C_water_array__kg_m3 = np.average(C_water_matrix__kg_m3, axis = 1)
            T_water_i_array__K = np.average(T_water_i_matrix__K, axis = 1)
            
            
            C_vapor_in__kg_m3 = C_vapor_matrix__kg_m3[0,0]
            T_air_i_in__K = T_air_i_matrix__K[0,0]
            # just take the exit conditions:
            C_vapor_array__kg_m3 = C_vapor_matrix__kg_m3[:,-1]
            T_air_i_array__K = T_air_i_matrix__K[:,-1]
            
            m_flow_cell_total__kg_s += m_flow_water_out__kg_s * 2 # need both halves of the cell
            
            normalized_C_water_arrays.append(C_water_array__kg_m3/C_water_in__kg_m3)
            normalized_C_vapor_arrays.append(C_vapor_array__kg_m3/C_vapor_in__kg_m3)
            normalized_T_water_arrays.append(T_water_i_array__K/T_water_i_in__K)
            normalized_T_air_arrays.append(T_air_i_array__K/T_air_i_in__K)
            
            temperature_array[i] = avg_Temperature
        
        Dirichlet_array = np.zeros(n_x)
        for i in range(n_x):
            cell = int(i/n_x * num_cells)
            Dirichlet_array[i] = temperature_array[cell]
        
        z_array = np.linspace(0,1,halfcell_parameters.n_z) # needed for plotting
        
        return Dirichlet_array, m_flow_cell_total__kg_s, normalized_C_water_arrays, normalized_C_vapor_arrays, normalized_T_water_arrays, normalized_T_air_arrays, z_array

    def get_cooled_water_temp__F(self):
        return self.util.degK_to_degF( self.average_outlet_T__K )

    def get_total_power__kW(self, fan_array__percent):
        return np.sum(fan_array__percent * self.basin_parameters.fan_max__hp * self.util.kw_per_hp)

if __name__ == '__main__':
    
    util = Util()
    
    fan_array__percent = np.ones(12) * 0.85
    
    basin = Basin(fan_array__percent)
    
    time_array_hr = []
    basin_temp_array_F = []
    
    # record results
    time_array_hr.append(basin.t_i__s / util.s_per_hr)
    basin_temp_array_F.append( util.degK_to_degF( basin.average_outlet_T__K) )
    
    percent_counter = 0
    percent_increment = 0.1
    time_final__s = 24 * util.s_per_hr
    
    while basin.t_i__s < time_final__s:
        if basin.t_i__s > percent_counter * time_final__s:
            print(str(int(round(100*percent_counter))) + '% complete (time =', round(basin.t_i__s/util.s_per_hr,2),'hours)' )
            
            if round(100*percent_counter) == 50:
                fan_array__percent = np.ones(12) * 0.5
                basin.advance_one_timestep(fan_array__percent, recalculate_fans=True)
                # record results
                time_array_hr.append(basin.t_i__s / util.s_per_hr)
                basin_temp_array_F.append( util.degK_to_degF( basin.average_outlet_T__K) )
                percent_counter += percent_increment
                continue
            
            percent_counter += percent_increment
            
        basin.advance_one_timestep(fan_array__percent, recalculate_fans=False)
        # record results
        time_array_hr.append(basin.t_i__s / util.s_per_hr)
        basin_temp_array_F.append( util.degK_to_degF( basin.average_outlet_T__K) )
    
    plt.figure()
    plt.plot(time_array_hr, basin_temp_array_F, '--k', lw=2, label='average')
    plt.title('Basin Outlet Temperatures vs Time')
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (degF)')
    plt.legend()
    plt.grid()
    plt.show()
        