
from material_properties import Water_properties, Moist_air_properties

class Basin_parameters:
    
    def __init__(self, util):
        
        # create water properties object
        water_properties = Water_properties()
        
        # values determined by system
        basin_length__ft = 481
        basin_width__ft = 51.5
        basin_height__ft = 10
        fan_max__hp = 182.9528 # adjusted the maxflow in the util script to allow for higher flowrates.
        
        # values I pick
        T_basin_init__degF = 70
        n_x = 30
        n_y = 8
        CFL = 1
        
        # unit conversions
        basin_length__m = basin_length__ft * util.m_per_ft
        basin_width__m = basin_width__ft * util.m_per_ft
        basin_height__m = basin_height__ft * util.m_per_ft
        T_basin_init__K = util.degF_to_degK(T_basin_init__degF)
        
        # store all the values
        self.basin_length__m    = basin_length__m
        self.basin_width__m     = basin_width__m
        self.basin_height__m    = basin_height__m
        self.T_basin_init__K    = T_basin_init__K
        self.water_properties   = water_properties
        self.n_x = n_x
        self.n_y = n_y
        self.CFL = CFL
        self.fan_max__hp = fan_max__hp
    
class HalfCell_crossflow_parameters:
    '''
    note: this only works for half of the cell. Do the calculations and multiply by two for the other half
    '''
    
    def __init__(self, util, power_draw_single_fan__HP, # max power draw is 182.9528 HP
                 hotbed_flowrate__gal_min = 180000/24, fan_motor_efficiency = 1, # I have divided by 24 because each cell has 2 sides. Both sides are the same though
                 P_amb__Pa = 86700, T_amb__degF = 65, Relative_humidity__percent = 20, T_hotbed__degF = 100): # note: I'm using the TMY pressure; Huntington is really 81854
        
        self.fan_motor_efficiency = fan_motor_efficiency
        
        if power_draw_single_fan__HP == 0:
            self.turned_off = True
        else:
            self.turned_off = False
        
        power_draw_single_fan_minimum__HP = 36 # breaks with anything less than 18 HP for the default parameters, so use 30 HP to be safe (~16%) anything above zero will be rounded up to this
        power_draw_single_fan__HP = max(power_draw_single_fan__HP, power_draw_single_fan_minimum__HP)
        
        # create water and air properties objects
        water_properties = Water_properties()
        moist_air_properties = Moist_air_properties()
        
        # values defined by system specs:
        tower_length__ft = 481
        evaporation_zone_width__ft = 4 # might need to account for the two sides somewhere else. 
        tower_height__ft = 29.75
        num_fans = 12
        
        # values I can tune:
        
        ## details on raindrop size: https://www.ems.psu.edu/~fraser/Bad/BadRain.html
        ## speed of falling raindrop: https://journals.ametsoc.org/doi/pdf/10.1175/1520-0450%281969%29008<0249%3ATVORA>2.0.CO%3B2
        
        # I orginally throught this should be about half a centimeter, but now I don't know. (about half a centimeter (must be greater than 0.5 mm, http://wxguys.ssec.wisc.edu/2013/09/10/how-fast-do-raindrops-fall/))
        droplet_diameter__m = 0.007 * .6  # pre-tuned = 0.007 m
        
        # maximum speed of the falling water. The specs say it has 1.33 ft of free fall height. Assuming it stops before falling and then falls only due to gravity w/o air resistence, this gives an average free fall speed of sqrt(distance*gravity/2) = 1.41 m/s 
        v_water_absolute__m_s = 1.41 * 0.65 # pre-tuned = 1.41
        # remember that the hotbed_flowrate__gal_min can also be specified
        # consider adding a heat transfer effeciency correction factor
        # consider adding a wind direction effeciency correction factor
        n_x = 100 # previous: 500 = gives temperatures of 306.21 and 290.06. 100 keeps it within 0.5 K
        n_z = 5
        
        # unit conversions:
        hotbed_flowrate__m3_s = hotbed_flowrate__gal_min * util.m3_s_per_gal_min
        tower_length__m = tower_length__ft * util.m_per_ft
        evaporation_zone_width__m = evaporation_zone_width__ft * util.m_per_ft
        tower_height__m = tower_height__ft * util.m_per_ft
        flowrate_out__m3_s = util.flowrate_out__CFM(power_draw_single_fan__HP, efficiency = self.fan_motor_efficiency) * util.m3_s_per_CFM
        T_hotbed__K = util.degF_to_degK(T_hotbed__degF)
        T_amb__K = util.degF_to_degK(T_amb__degF)
        T_water_0__K = T_hotbed__K
        T_air_in__K = T_amb__K
        
        # additional calculations:
        m_flow_hotbed__kg_s = hotbed_flowrate__m3_s * water_properties.liquid_density__kg_m3(T_hotbed__K)
        cell_face_area_total__m2 = tower_height__m * tower_length__m
        
        # need to divide by 2 since there are two sides to each cell, but only one fan
        v_air_cell__m_s = flowrate_out__m3_s / cell_face_area_total__m2 / 2 
        
        v_water_airspeed__m_s = (v_water_absolute__m_s**2 + v_air_cell__m_s**2)**(1/2) # for cross flow, the velocity vectors are perpendicular, so the net velocity is found with the pythagorean theorem
        
        SA_evaporation_zone_horizontal__m2 = evaporation_zone_width__m * tower_length__m / num_fans
        deltaZ__m = tower_height__m / n_z # volumetric descretization
        deltaX__m = evaporation_zone_width__m / n_x # volumetric descretization
        C_water_0__kg_m3 = m_flow_hotbed__kg_s / (v_water_absolute__m_s * SA_evaporation_zone_horizontal__m2)
        
        C_vapor_in__kg_m3 = water_properties.P_vapor_ambient__Pa(Relative_humidity__percent, T_amb__K) / (util.R * T_amb__K) * util.MW_water__kg_mol
        
        discretized_volume__m3          = deltaX__m * deltaZ__m * tower_length__m / num_fans
        discretized_area_waterside__m2  = deltaX__m *             tower_length__m / num_fans
        discretized_area_vaporside__m2  =             deltaZ__m * tower_length__m / num_fans
        
#        internal parameters
        self.deltaX__m              = deltaX__m
        self.n_x                    = n_x
        self.droplet_diameter__m    = droplet_diameter__m
        self.v_water_airspeed__m_s  = v_water_airspeed__m_s
        self.v_water_absolute__m_s  = v_water_absolute__m_s
        self.P_amb__Pa              = P_amb__Pa
        self.water_properties       = water_properties
        self.moist_air_properties   = moist_air_properties
        self.tower_length__m        = tower_length__m
        self.discretized_volume__m3 = discretized_volume__m3
        self.discretized_area_waterside__m2 = discretized_area_waterside__m2
        self.discretized_area_vaporside__m2 = discretized_area_vaporside__m2
        self.deltaZ__m              = deltaZ__m
        self.n_z                    = n_z
        self.v_air_cell__m_s        = v_air_cell__m_s
        
#        boundary parameters
        self.C_water_0__kg_m3       = C_water_0__kg_m3
        self.T_water_0__K           = T_water_0__K
        self.T_air_in__K             = T_air_in__K
        self.C_vapor_in__kg_m3      = C_vapor_in__kg_m3
        