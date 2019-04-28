
import numpy as np

class Water_properties:
    
    def __init__(self):
        
        # unit conversions
        Pa_per_kPa = 1000
        CtoK = 273.15
        J_per_kJ = 1000
        W_per_mW = 1/1000
        
        self.MW_water__kg_mol = 18.02/1000
        
        # data source: https://www.engineeringtoolbox.com/water-properties-d_1573.html
        # format, each row is: Temperature(degC), Vapor Pressure(kPa), Hvap(J/mol), Hvap(kJ/kg), Hvap(Wh/kg), Hvap(BTU(IT)/lb)
        water_data = np.array([ [0.01, 0.61165, 45054, 2500.9, 694.69, 1075.2],
                                [2, 0.70599, 44970, 2496.2, 693.39, 1073.2],
                                [4, 0.81355, 44883, 2491.4, 692.06, 1071.1],
                                [10, 1.2282, 44627, 2477.2, 688.11, 1065.0],
                                [14, 1.5990, 44456, 2467.7, 685.47, 1060.9],
                                [18, 2.0647, 44287, 2458.3, 682.86, 1056.9],
                                [20, 2.3393, 44200, 2453.5, 681.53, 1054.8],
                                [25, 3.1699, 43988, 2441.7, 678.25, 1049.7],
                                [30, 4.2470, 43774, 2429.8, 674.94, 1044.6],
                                [34, 5.3251, 43602, 2420.3, 672.31, 1040.5],
                                [40, 7.3849, 43345, 2406.0, 668.33, 1034.4],
                                [44, 9.1124, 43172, 2396.4, 665.67, 1030.3],
                                [50, 12.352, 42911, 2381.9, 661.64, 1024.0],
                                [54, 15.022, 42738, 2372.3, 658.97, 1019.9],
                                [60, 19.946, 42475, 2357.7, 654.92, 1013.6],
                                [70, 31.201, 42030, 2333.0, 648.06, 1003.0],
                                [80, 47.414, 41579, 2308.0, 641.11, 992.26],
                                [90, 70.182, 41120, 2282.5, 634.03, 981.30],
                                [96, 87.771, 40839, 2266.9, 629.69, 974.59],
                                [100, 101.42, 40650, 2256.4, 626.78, 970.08],
                                [110, 143.38, 40167, 2229.6, 619.33, 958.56],
                                [120, 198.67, 39671, 2202.1, 611.69, 946.73],
                                [140, 361.54, 38630, 2144.3, 595.64, 921.88],
                                [160, 618.23, 37508, 2082.0, 578.33, 895.10],
                                [180, 1002.8, 36286, 2014.2, 559.50, 865.95],
                                [200, 1554.9, 34944, 1939.7, 538.81, 833.92],
                                [220, 2319.6, 33462, 1857.4, 515.94, 798.54],
                                [240, 3346.9, 31804, 1765.4, 490.39, 758.99],
                                [260, 4692.3, 29934, 1661.6, 461.56, 714.36],
                                [280, 6416.6, 27798, 1543.0, 428.61, 663.37],
                                [300, 8587.9, 25304, 1404.6, 390.17, 603.87],
                                [320, 11284, 22310, 1238.4, 344.00, 532.42],
                                [340, 14601, 18507, 1027.3, 285.36, 441.66],
                                [360, 18666, 12967, 719.8, 199.9, 309.5],
                                [373.946, 22064, 0, 0.0, 0.0, 0.0] ])
        
        # data source: https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html
        # format, each row is: Temperature(degC), Density(gm/cm^3), Density(kg/m^3), ... (and a bunch of other stuff)
        water_data2 = np.array([[0.1, 0.9998495, 999.85, 1.9400, 62.4186, 8.3441, 9.8052, 62.419, -0.68],
                                [1, 0.9999017, 999.90, 1.9401, 62.4218, 8.3446, 9.8057, 62.422, -0.50],
                                [4, 0.9999749, 999.97, 1.9403, 62.4264, 8.3452, 9.8064, 62.426, 0.003],
                                [10, 0.9997000, 999.70, 1.9397, 62.4094, 8.3429, 9.8037, 62.409, 0.88],
                                [15, 0.9991026, 999.10, 1.9386, 62.3719, 8.3379, 9.7978, 62.372, 1.51],
                                [20, 0.9982067, 998.21, 1.9368, 62.3160, 8.3304, 9.7891, 62.316, 2.07],
                                [25, 0.9970470, 997.05, 1.9346, 62.2436, 8.3208, 9.7777, 62.244, 2.57],
                                [30, 0.9956488, 995.65, 1.9319, 62.1563, 8.3091, 9.7640, 62.156, 3.03],
                                [35, 0.9940326, 994.03, 1.9287, 62.0554, 8.2956, 9.7481, 62.055, 3.45],
                                [40, 0.9922152, 992.22, 1.9252, 61.9420, 8.2804, 9.7303, 61.942, 3.84],
                                [45, 0.99021, 990.21, 1.9213, 61.8168, 8.2637, 9.7106, 61.817, 4.20],
                                [50, 0.98804, 988.04, 1.9171, 61.6813, 8.2456, 9.6894, 61.681, 4.54],
                                [55, 0.98569, 985.69, 1.9126, 61.5346, 8.2260, 9.6663, 61.535, 4.86],
                                [60, 0.98320, 983.20, 1.9077, 61.3792, 8.2052, 9.6419, 61.379, 5.16],
                                [65, 0.98055, 980.55, 1.9026, 61.2137, 8.1831, 9.6159, 61.214, 5.44],
                                [70, 0.97776, 977.76, 1.8972, 61.0396, 8.1598, 9.5886, 61.040, 5.71],
                                [75, 0.97484, 974.84, 1.8915, 60.8573, 8.1354, 9.5599, 60.857, 5.97],
                                [80, 0.97179, 971.79, 1.8856, 60.6669, 8.1100, 9.5300, 60.667, 6.21],
                                [85, 0.96861, 968.61, 1.8794, 60.4683, 8.0834, 9.4988, 60.468, 6.44],
                                [90, 0.96531, 965.31, 1.8730, 60.2623, 8.0559, 9.4665, 60.262, 6.66],
                                [95, 0.96189, 961.89, 1.8664, 60.0488, 8.0274, 9.4329, 60.049, 6.87],
                                [100, 0.95835, 958.35, 1.8595, 59.8278, 7.9978, 9.3982, 59.828, 7.03],
                                [110, 0.95095, 950.95, 1.8451, 59.3659, 7.9361, 9.3256, 59.366, 8.01],
                                [120, 0.94311, 943.11, 1.8299, 58.8764, 7.8706, 9.2487, 58.876, 8.60],
                                [140, 0.92613, 926.13, 1.7970, 57.8164, 7.7289, 9.0822, 57.816, 9.75],
                                [160, 0.90745, 907.45, 1.7607, 56.6503, 7.5730, 8.8990, 56.650, 11.0],
                                [180, 0.88700, 887.00, 1.7211, 55.3736, 7.4024, 8.6985, 55.374, 12.3],
                                [200, 0.86466, 864.66, 1.6777, 53.9790, 7.2159, 8.4794, 53.979, 13.9] ])
        
        # data source: https://www.engineeringtoolbox.com/water-dynamic-kinematic-viscosity-d_596.html
        # format, each row is: Temperature(degC), Pressure(MPa), Dynamic Viscocity(Pa s) ... (and a bunch of other stuff),...Kinematic viscocity(m^2/s *10^-6)
        water_data3 = np.array([[0.01, 0.000612, 0.0017914, 1.79140, 3.7414, 1.7918],
                                [10, 0.0012, 0.0013060, 1.30600, 2.7276, 1.3065],
                                [20, 0.0023, 0.0010016, 1.00160, 2.0919, 1.0035],
                                [25, 0.0032, 0.0008900, 0.89004, 1.8589, 0.8927],
                                [30, 0.0042, 0.0007972, 0.79722, 1.6650, 0.8007],
                                [40, 0.0074, 0.0006527, 0.65272, 1.3632, 0.6579],
                                [50, 0.0124, 0.0005465, 0.54650, 1.1414, 0.5531],
                                [60, 0.0199, 0.0004660, 0.46602, 0.9733, 0.4740],
                                [70, 0.0312, 0.0004035, 0.40353, 0.8428, 0.4127],
                                [80, 0.0474, 0.0003540, 0.35404, 0.7394, 0.3643],
                                [90, 0.0702, 0.0003142, 0.31417, 0.6562, 0.3255],
                                [100, 0.101, 0.0002816, 0.28158, 0.5881, 0.2938],
                                [110, 0.143, 0.0002546, 0.25461, 0.5318, 0.2677],
                                [120, 0.199, 0.0002320, 0.23203, 0.4846, 0.2460],
                                [140, 0.362, 0.0001966, 0.19664, 0.4107, 0.2123],
                                [160, 0.618, 0.0001704, 0.17043, 0.3559, 0.1878],
                                [180, 1.00, 0.0001504, 0.15038, 0.3141, 0.1695],
                                [200, 1.55, 0.0001346, 0.13458, 0.2811, 0.1556],
                                [220, 2.32, 0.0001218, 0.12177, 0.2543, 0.1449],
                                [240, 3.35, 0.0001111, 0.11106, 0.2320, 0.1365],
                                [260, 4.69, 0.0001018, 0.10181, 0.2126, 0.1299],
                                [280, 6.42, 0.0000936, 0.09355, 0.1954, 0.1247],
                                [300, 8.59, 0.0000859, 0.08586, 0.1793, 0.1206],
                                [320, 11.3, 0.0000783, 0.07831, 0.1636, 0.1174],
                                [340, 14.6, 0.0000703, 0.07033, 0.1469, 0.1152],
                                [360, 18.7, 0.0000603, 0.06031, 0.1260, 0.1143] ])
        
        # data source: https://www.engineeringtoolbox.com/specific-heat-capacity-water-d_660.html
        # format, each row is: Temperature(degC), ... (and a bunch of other stuff),...Isobaric Heat Capacity-Cp (kJ / (kg K)), ... (two other columns)
        water_data4 = np.array([[0.01, 75.981, 4.2174, 0.001172, 1.0073, 76.026, 4.2199, 0.001172, 1.0079],
                                [10, 75.505, 4.1910, 0.001164, 1.0010, 75.586, 4.1955, 0.001165, 1.0021],
                                [20, 74.893, 4.1570, 0.001155, 0.9929, 75.386, 4.1844, 0.001162, 0.9994],
                                [25, 74.548, 4.1379, 0.001149, 0.9883, 75.336, 4.1816, 0.001162, 0.9988],
                                [30, 74.181, 4.1175, 0.001144, 0.9834, 75.309, 4.1801, 0.001161, 0.9984],
                                [40, 73.392, 4.0737, 0.001132, 0.9730, 75.300, 4.1796, 0.001161, 0.9983],
                                [50, 72.540, 4.0264, 0.001118, 0.9617, 75.334, 4.1815, 0.001162, 0.9987],
                                [60, 71.644, 3.9767, 0.001105, 0.9498, 75.399, 4.1851, 0.001163, 0.9996],
                                [70, 70.716, 3.9252, 0.001090, 0.9375, 75.491, 4.1902, 0.001164, 1.0008],
                                [80, 69.774, 3.8729, 0.001076, 0.9250, 75.611, 4.1969, 0.001166, 1.0024],
                                [90, 68.828, 3.8204, 0.001061, 0.9125, 75.763, 4.2053, 0.001168, 1.0044],
                                [100, 67.888, 3.7682, 0.001047, 0.9000, 75.950, 4.2157, 0.001171, 1.0069],
                                [110, 66.960, 3.7167, 0.001032, 0.8877, 76.177, 4.2283, 0.001175, 1.0099],
                                [120, 66.050, 3.6662, 0.001018, 0.8757, 76.451, 4.2435, 0.001179, 1.0135],
                                [140, 64.306, 3.5694, 0.000992, 0.8525, 77.155, 4.2826, 0.001190, 1.0229],
                                [160, 62.674, 3.4788, 0.000966, 0.8309, 78.107, 4.3354, 0.001204, 1.0355],
                                [180, 61.163, 3.3949, 0.000943, 0.8109, 79.360, 4.4050, 0.001224, 1.0521],
                                [200, 59.775, 3.3179, 0.000922, 0.7925, 80.996, 4.4958, 0.001249, 1.0738],
                                [220, 58.514, 3.2479, 0.000902, 0.7757, 83.137, 4.6146, 0.001282, 1.1022],
                                [240, 57.381, 3.1850, 0.000885, 0.7607, 85.971, 4.7719, 0.001326, 1.1397],
                                [260, 56.392, 3.1301, 0.000869, 0.7476, 89.821, 4.9856, 0.001385, 1.1908],
                                [280, 55.578, 3.0849, 0.000857, 0.7368, 95.285, 5.2889, 0.001469, 1.2632],
                                [300, 55.003, 3.0530, 0.000848, 0.7292, 103.60, 5.7504, 0.001597, 1.3735],
                                [320, 54.819, 3.0428, 0.000845, 0.7268, 117.78, 6.5373, 0.001816, 1.5614],
                                [340, 55.455, 3.0781, 0.000855, 0.7352, 147.88, 8.2080, 0.002280, 1.9604],
                                [360, 59.402, 3.2972, 0.000916, 0.7875, 270.31, 15.004, 0.004168, 3.5836] ])
        
        # data source: https://www.engineeringtoolbox.com/water-vapor-d_979.html
        # format, each row is: Temperature(K), Specific Heat-Cp (kJ/kg K)
        water_data5 = np.array([[175, 1.850],
                                [200, 1.851],
                                [225, 1.852],
                                [250, 1.855],
                                [275, 1.859],
                                [300, 1.864],
                                [325, 1.871],
                                [350, 1.880],
                                [375, 1.890],
                                [400, 1.901],
                                [450, 1.926],
                                [500, 1.954],
                                [550, 1.984],
                                [600, 2.015],
                                [650, 2.047],
                                [700, 2.080],
                                [750, 2.113],
                                [800, 2.147],
                                [850, 2.182],
                                [900, 2.217],
                                [950, 2.252],
                                [1000, 2.288] ])
        
        # data source: https://www.engineeringtoolbox.com/saturated-steam-properties-d_101.html
        # format, each row has values for saturated steam: Absolute Pressure(kPa), Evaporation Temperature(degC), ... (a bunch of others)
        water_data6 = np.array([[0.8, 3.8, 160, 0.00626, 15.8, 2493, 2509, 9.058],
                                [2.0, 17.5, 67.0, 0.0149, 73.5, 2460, 2534, 8.725],
                                [5.0, 32.9, 28.2, 0.0354, 137.8, 2424, 2562, 8.396],
                                [10.0, 45.8, 14.7, 0.0682, 191.8, 2393, 2585, 8.151],
                                [20.0, 60.1, 7.65, 0.131, 251.5, 2358, 2610, 7.909],
                                [28, 67.5, 5.58, 0.179, 282.7, 2340, 2623, 7.793],
                                [35, 72.7, 4.53, 0.221, 304.3, 2327, 2632, 7.717],
                                [45, 78.7, 3.58, 0.279, 329.6, 2312, 2642, 7.631],
                                [55, 83.7, 2.96, 0.338, 350.6, 2299, 2650, 7.562],
                                [65, 88.0, 2.53, 0.395, 368.6, 2288, 2657, 7.506],
                                [75, 91.8, 2.22, 0.450, 384.5, 2279, 2663, 7.457],
                                [85, 95.2, 1.97, 0.507, 398.6, 2270, 2668, 7.415],
                                [95, 98.2, 1.78, 0.563, 411.5, 2262, 2673, 7.377],
                                [100, 99.6, 1.69, 0.590, 417.5, 2258, 2675, 7.360] ])
    
        # data source: https://www.engineeringtoolbox.com/water-liquid-gas-thermal-conductivity-temperature-pressure-d_2012.html
        # format, each row is: Temperature(degC), Thermal conductivity (mW/ m K), ....
        water_data7 = np.array([[0.01, 555.75, 0.4779, 0.3211],
                                [10, 578.64, 0.4975, 0.3343],
                                [20, 598.03, 0.5142, 0.3455],
                                [30, 614.50, 0.5284, 0.3551],
                                [40, 628.56, 0.5405, 0.3632],
                                [50, 640.60, 0.5508, 0.3701],
                                [60, 650.91, 0.5597, 0.3761],
                                [70, 659.69, 0.5672, 0.3812],
                                [80, 667.02, 0.5735, 0.3854],
                                [90, 672.88, 0.5786, 0.3888],
                                [99.6, 677.03, 0.5821, 0.3912] ])
        
        # unpack and convert to proper units
        self.temperature_data__K                     = water_data[:,0] + CtoK
        self.vapor_pressure_data__Pa                 = water_data[:,1] * Pa_per_kPa # get the data in terms of Pa
        self.Hvap_data__J_kg                         = water_data[:,3] * J_per_kJ # only need the data in terms of mass (not mols)
        
        self.temperature_data2__K                    = water_data2[:,0] + CtoK
        self.water_density_data2__kg_m3              = water_data2[:,2] # get the data in terms of kg/m^3
        
        self.temperature_data3__K                    = water_data3[:,0] + CtoK
        self.water_dynamic_viscocity_data3__kg_m_s   = water_data3[:,2]  # get the data in terms of Pa_s
        self.water_kinematic_viscocity_data3__m2_s   = water_data3[:,-1] *1e-6 # get the data in terms of m^2/s
        
        self.temperature_data4__K                    = water_data4[:,0] + CtoK
        self.water_Cp_data4__J_kg_K                  = water_data4[:,-3] * J_per_kJ # get the data in terms of J/kg K
        
        self.temperature_data5__K                    = water_data5[:,0]
        self.vapor_Cp_data5__J_kg_K                  = water_data5[:,1] * J_per_kJ # get the data in terms of J/kg K
        
        self.P_sat_steam_data6__Pa                   = water_data6[:,0] * Pa_per_kPa # get the data in terms of Pa
        self.temperature_data6__K                    = water_data6[:,1] + CtoK
        
        self.temperature_data7__K                    = water_data7[:,0] + CtoK
        self.water_k_thermal_cond_data7__W_m_K      = water_data7[:,1] * W_per_mW # get the data in terms of W/m K
        
    ##############################################################
    # interploation methods for temperature dependant properties #
    ##############################################################
    
    def enthalpy_vaporization__J_kg(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data__K, self.Hvap_data__J_kg)
    
    def vapor_pressure__Pa(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data__K, self.vapor_pressure_data__Pa)
    
    def liquid_density__kg_m3(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data2__K, self.water_density_data2__kg_m3)
    
    def liquid_dynamic_viscocity__kg_m_s(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data3__K, self.water_dynamic_viscocity_data3__kg_m_s)
    
    def liquid_kinematic_viscocity__m2_s(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data3__K, self.water_kinematic_viscocity_data3__m2_s)
    
    def liquid_heat_capacity__J_kg_K(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data4__K, self.water_Cp_data4__J_kg_K)
    
    def vapor_heat_capacity__J_kg_K(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data5__K, self.vapor_Cp_data5__J_kg_K)

    def Pressure_saturated__Pa(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data6__K, self.P_sat_steam_data6__Pa)

    def P_vapor_ambient__Pa(self, RH, T_amb):
        P_sat_water = self.Pressure_saturated__Pa(T_amb)
        return RH/100 * P_sat_water

    def liquid_k_thermal_cond__W_m_K(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data7__K, self.water_k_thermal_cond_data7__W_m_K)

    def liquid_thermal_diffusivity(self, Temperature__K):
        '''
        thermal diffusivity is D = k/(rho * cp)
        '''
        return self.liquid_k_thermal_cond__W_m_K(Temperature__K) / ( self.liquid_density__kg_m3(Temperature__K) * self.liquid_heat_capacity__J_kg_K(Temperature__K) )

class Dry_air_properties:
    def __init__(self):
        
        # unit conversions
        CtoK = 273.15
        J_per_kJ = 1000
        W_per_mW = 1/1000
        
        self.MW_air__kg_mol = 28.97/1000
        
        # data source: https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
        # format, each row is: Temperature(degC), Dynamic Viscocity(Pa s *10^-6) ... (and a bunch of other stuff),...Kinematic viscocity(m^2/s *10^-6), Kinematic viscocity(ft^2/s *10^-4)
        air_data = np.array([   [-75, 13.18, 0.01318, 0.2753, 0.03188, 7.40, 0.796],
                                [-50, 14.56, 0.01456, 0.3041, 0.03523, 9.22, 0.992],
                                [-25, 15.88, 0.01588, 0.3317, 0.03842, 11.18, 1.203],
                                [-15, 16.40, 0.01640, 0.3425, 0.03966, 12.01, 1.292],
                                [-10, 16.65, 0.01665, 0.3477, 0.04028, 12.43, 1.338],
                                [-5, 16.90, 0.01690, 0.3530, 0.04089, 12.85, 1.383],
                                [0, 17.15, 0.01715, 0.3582, 0.04149, 13.28, 1.430],
                                [5, 17.40, 0.01740, 0.3633, 0.04209, 13.72, 1.477],
                                [10, 17.64, 0.01764, 0.3685, 0.04268, 14.16, 1.524],
                                [15, 17.89, 0.01789, 0.3735, 0.04327, 14.61, 1.573],
                                [20, 18.13, 0.01813, 0.3786, 0.04385, 15.06, 1.621],
                                [25, 18.37, 0.01837, 0.3836, 0.04443, 15.52, 1.671],
                                [30, 18.60, 0.01860, 0.3885, 0.04500, 15.98, 1.720],
                                [40, 19.07, 0.01907, 0.3983, 0.04614, 16.92, 1.822],
                                [50, 19.53, 0.01953, 0.4080, 0.04725, 17.88, 1.925],
                                [60, 19.99, 0.01999, 0.4175, 0.04835, 18.86, 2.030],
                                [80, 20.88, 0.02088, 0.4361, 0.05051, 20.88, 2.248],
                                [100, 21.74, 0.02174, 0.4541, 0.05260, 22.97, 2.473],
                                [125, 22.79, 0.02279, 0.4760, 0.05513, 25.69, 2.765],
                                [150, 23.80, 0.02380, 0.4971, 0.05758, 28.51, 3.069],
                                [175, 24.78, 0.02478, 0.5176, 0.05995, 31.44, 3.384],
                                [200, 25.73, 0.02573, 0.5374, 0.06225, 34.47, 3.710],
                                [225, 26.66, 0.02666, 0.5567, 0.06448, 37.60, 4.047],
                                [300, 29.28, 0.02928, 0.6115, 0.07083, 47.54, 5.117],
                                [412, 32.87, 0.03287, 0.6865, 0.07952, 63.82, 6.869],
                                [500, 35.47, 0.03547, 0.7409, 0.08581, 77.72, 8.366],
                                [600, 38.25, 0.03825, 0.7988, 0.09252, 94.62, 10.19],
                                [700, 40.85, 0.04085, 0.8532, 0.09883, 112.6, 12.12],
                                [800, 43.32, 0.04332, 0.9047, 0.1048, 131.7, 14.17],
                                [900, 45.66, 0.04566, 0.9535, 0.1104, 151.7, 16.33],
                                [1000, 47.88, 0.04788, 1.000, 0.1158, 172.7, 18.59],
                                [1100, 50.01, 0.05001, 1.045, 0.1210, 194.6, 20.95] ])
        
        # data source: https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
        # format, each row is: Temperature(K), ... (and a bunch of other stuff),...Isobaric heat capacity-Cp (kJ/kg K), ...(4 other columns)
        air_data2 = np.array([  [60, -213, -352, 0.03398, 1.173, 0.0003258, 0.2802, 0.2287, 0.05506, 1.901, 0.000528, 0.45405, 0.37071, 1.621],
                                [78.79, -194, -318, 0.03044, 1.051, 0.0002919, 0.2510, 0.2050, 0.05599, 1.933, 0.000537, 0.46169, 0.37695, 1.839],
                                [81.61, -192, -313, 0.02172, 0.7500, 0.0002083, 0.1791, 0.1463, 0.03154, 1.089, 0.000303, 0.26010, 0.21237, 1.452],
                                [100, -173, -280, 0.02109, 0.7280, 0.0002022, 0.1739, 0.1420, 0.03012, 1.040, 0.000289, 0.24833, 0.20276, 1.428],
                                [120, -153, -244, 1.02109, 0.7211, 0.0002003, 0.1722, 0.1406, 1.03012, 1.022, 0.000284, 0.24410, 0.19930, 1.009],
                                [140, -133, -208, 0.02081, 0.7184, 0.0001996, 0.1716, 0.1401, 0.02937, 1.014, 0.000282, 0.24219, 0.19774, 1.411],
                                [160, -113, -172, 0.02077, 0.7172, 0.0001992, 0.1713, 0.1399, 0.02928, 1.011, 0.000281, 0.24147, 0.19716, 1.410],
                                [180, -93.2, -136, 0.02076, 0.7166, 0.0001991, 0.1712, 0.1397, 0.02920, 1.008, 0.000280, 0.24076, 0.19657, 1.407],
                                [200, -73.2, -99.7, 0.02075, 0.7163, 0.0001990, 0.1711, 0.1397, 0.02917, 1.007, 0.000280, 0.24052, 0.19638, 1.406],
                                [220, -53.2, -63.7, 0.02075, 0.7163, 0.0001990, 0.1711, 0.1397, 0.02914, 1.006, 0.000279, 0.24028, 0.19618, 1.404],
                                [240, -33.2, -27.7, 0.02075, 0.7164, 0.0001990, 0.1711, 0.1397, 0.02914, 1.006, 0.000279, 0.24028, 0.19618, 1.404],
                                [260, -13.2, 8.3, 0.02076, 0.7168, 0.0001991, 0.1712, 0.1398, 0.02914, 1.006, 0.000279, 0.24028, 0.19618, 1.403],
                                [273.2, 0.0, 32.0, 0.02077, 0.7171, 0.0001992, 0.1713, 0.1398, 0.02914, 1.006, 0.000279, 0.24028, 0.19618, 1.403],
                                [280, 6.9, 44.3, 0.02078, 0.7173, 0.0001993, 0.1713, 0.1399, 0.02914, 1.006, 0.000279, 0.24028, 0.19618, 1.402],
                                [288.7, 15.6, 60.0, 0.02078, 0.7175, 0.0001993, 0.1714, 0.1399, 0.02914, 1.006, 0.000279, 0.24030, 0.19620, 1.402],
                                [300, 26.9, 80.3, 0.02080, 0.7180, 0.0001994, 0.1715, 0.1400, 0.02915, 1.006, 0.000280, 0.24036, 0.19625, 1.402],
                                [320, 46.9, 116, 0.02083, 0.7192, 0.0001998, 0.1718, 0.1403, 0.02917, 1.007, 0.000280, 0.24052, 0.19638, 1.400],
                                [340, 66.9, 152, 0.02087, 0.7206, 0.0002002, 0.1721, 0.1405, 0.02923, 1.009, 0.000280, 0.24100, 0.19677, 1.400],
                                [360, 86.9, 188, 0.02092, 0.7223, 0.0002006, 0.1725, 0.1409, 0.02926, 1.010, 0.000281, 0.24123, 0.19696, 1.398],
                                [380, 107, 224, 0.02098, 0.7243, 0.0002012, 0.1730, 0.1412, 0.02931, 1.012, 0.000281, 0.24171, 0.19735, 1.397],
                                [400, 127, 260, 0.02105, 0.7266, 0.0002018, 0.1735, 0.1417, 0.02937, 1.014, 0.000282, 0.24219, 0.19774, 1.396],
                                [500, 227, 440, 0.02150, 0.7424, 0.0002062, 0.1773, 0.1448, 0.02983, 1.030, 0.000286, 0.24597, 0.20083, 1.387],
                                [600, 327, 620, 0.02213, 0.7641, 0.0002123, 0.1825, 0.1490, 0.03044, 1.051, 0.000292, 0.25103, 0.20496, 1.375],
                                [700, 427, 800, 0.02282, 0.7877, 0.0002188, 0.1881, 0.1536, 0.03114, 1.075, 0.000299, 0.25675, 0.20963, 1.365],
                                [800, 527, 980, 0.02351, 0.8117, 0.0002255, 0.1939, 0.1583, 0.03183, 1.099, 0.000305, 0.26249, 0.21432, 1.354],
                                [900, 627, 1160, 0.02415, 0.8338, 0.0002316, 0.1991, 0.1626, 0.03247, 1.121, 0.000311, 0.26772, 0.21858, 1.344],
                                [1100, 827, 1520, 0.02525, 0.8716, 0.0002421, 0.2082, 0.1700, 0.03356, 1.159, 0.000322, 0.27675, 0.22596, 1.329],
                                [1500, 1227, 2240, 0.02673, 0.9230, 0.0002564, 0.2204, 0.1800, 0.03505, 1.210, 0.000336, 0.28901, 0.23597, 1.311],
                                [1900, 1627, 2960, 0.02762, 0.9535, 0.0002649, 0.2277, 0.1859, 0.03593, 1.241, 0.000345, 0.29631, 0.24193, 1.301] ])
        
        # data source: https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
        # format, each row is: Temperature(degC), Thermal conductivity (mW/m K) ... (and a bunch of other stuff)
        air_data3 = np.array([  [-190, 7.82, 0.00672, 0.00452],
                                [-150, 11.69, 0.01005, 0.00675],
                                [-100, 16.20, 0.01393, 0.00936],
                                [-75, 18.34, 0.01577, 0.01060],
                                [-50, 20.41, 0.01755, 0.01179],
                                [-25, 22.41, 0.01927, 0.01295],
                                [-15, 23.20, 0.01995, 0.01340],
                                [-10, 23.59, 0.02028, 0.01363],
                                [-5, 23.97, 0.02061, 0.01385],
                                [0, 24.36, 0.02094, 0.01407],
                                [5, 24.74, 0.02127, 0.01429],
                                [10, 25.12, 0.02160, 0.01451],
                                [15, 25.50, 0.02192, 0.01473],
                                [20, 25.87, 0.02225, 0.01495],
                                [25, 26.24, 0.02257, 0.01516],
                                [30, 26.62, 0.02289, 0.01538],
                                [40, 27.35, 0.02352, 0.01580],
                                [50, 28.08, 0.02415, 0.01623],
                                [60, 28.80, 0.02477, 0.01664],
                                [80, 30.23, 0.02599, 0.01746],
                                [100, 31.62, 0.02719, 0.01827],
                                [125, 33.33, 0.02866, 0.01926],
                                [150, 35.00, 0.03010, 0.02022],
                                [175, 36.64, 0.03151, 0.02117],
                                [200, 38.25, 0.03289, 0.02210],
                                [225, 39.83, 0.03425, 0.02301],
                                [300, 44.41, 0.03819, 0.02566],
                                [412, 50.92, 0.04378, 0.02942],
                                [500, 55.79, 0.04797, 0.03224],
                                [600, 61.14, 0.05257, 0.03533],
                                [700, 66.32, 0.05702, 0.03832],
                                [800, 71.35, 0.06135, 0.04122],
                                [900, 76.26, 0.06557, 0.04406],
                                [1000, 81.08, 0.06971, 0.04685],
                                [1100, 85.83, 0.07380, 0.04959] ])

        self.temperature_data__K                  = air_data[:,0] + CtoK
        self.air_dynamic_viscocity_data__kg_m_s   = air_data[:,1] * 1e-6  # get the data in terms of Pa_s
        self.air_kinematic_viscocity_data__m2_s   = air_data[:,-2] *1e-6  # get the data in terms of m^2/s
        
        self.temperature_data2__K                 = air_data2[:,0]
        self.air_Cp_data2__J_kg_K                 = air_data2[:,-5] * J_per_kJ # get the data in terms of J/kg K
        
        self.temperature_data3__K                 = air_data3[:,0]
        self.air_k_thermal_cond_data3__W_m_K      = air_data3[:,1] * W_per_mW # get the data in terms of W/m K
        
        
    def dynamic_viscocity__kg_m_s(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data__K, self.air_dynamic_viscocity_data__kg_m_s)
    
    def kinematic_viscocity__m2_s(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data__K, self.air_kinematic_viscocity_data__m2_s)
    
    def heat_capacity__J_kg_K(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data2__K, self.air_Cp_data2__J_kg_K)
    
    def k_thermal_cond__W_m_K(self, Temperature__K):
        return np.interp(Temperature__K, self.temperature_data3__K, self.air_k_thermal_cond_data3__W_m_K)
    

class Moist_air_properties:
    def __init__(self):
        self.dry_air_properties = Dry_air_properties()
        self.water_properties = Water_properties()
        
    def dynamic_viscocity__kg_m_s(self, Temperature__K, vapor_mass_fraction):
        # TODO need to get reliable data for viscocities of water vapor and account for this here. Dry air is a great approximation in the meantime
        return self.dry_air_properties.dynamic_viscocity__kg_m_s(Temperature__K)
    
    def kinematic_viscocity__m2_s(self, Temperature__K, vapor_mass_fraction):
        # TODO need to get reliable data for viscocities of water vapor and account for this here.  Dry air is a great approximation in the meantime
        return self.dry_air_properties.kinematic_viscocity__m2_s(Temperature__K)
    
    def k_thermal_cond__W_m_K(self, Temperature__K):
        # TODO need to get reliable data for thermal conductivity of water vapor and account for this here.  Dry air is a great approximation in the meantime
        return self.dry_air_properties.k_thermal_cond__W_m_K(Temperature__K)
    
    def density__kg_m3(self, Pressure__Pa, Temperature__K, vapor_mass_fraction):
        R__J_mol_K = 8.31451
        density__mol_m3 = Pressure__Pa / (R__J_mol_K * Temperature__K)
        MW_air__kg_mol = self.dry_air_properties.MW_air__kg_mol # 28.97/1000
        MW_water__kg_mol = self.water_properties.MW_water__kg_mol # 18.02/1000
        MW_avg__kg_mol = vapor_mass_fraction * MW_water__kg_mol + (1-vapor_mass_fraction) * MW_air__kg_mol
        return density__mol_m3 * MW_avg__kg_mol
    
    def heat_capacity(self, Temperature__K, vapor_mass_fraction):
        air_mass_fraction = 1 - vapor_mass_fraction
        return ( air_mass_fraction * self.dry_air_properties.heat_capacity__J_kg_K(Temperature__K) + 
                 vapor_mass_fraction * self.water_properties.vapor_heat_capacity__J_kg_K(Temperature__K) )
        