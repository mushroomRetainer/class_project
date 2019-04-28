
import numpy as np
from PID_controller import PID
from basin_solver import Basin
import matplotlib.pyplot as plt
from neural_network_optimization import Neural_netork
from util import Util
from data_reader import TMY_data, training_data
import csv


def random_return_temperatures(num_values, random_seed = 12345, lower_bound = 70, upper_bound = 100):
    np.random.seed(random_seed)
    return np.random.rand(num_values) * (upper_bound - lower_bound) + lower_bound
    

# all arrays will be indexed by hour. Use "int(time_hours)" of the basin clock to get the appropriate index (note: basin clock is in seconds)

def main(NN_status = 'learn', TMY_quantity = 'day', training_file = None):
    
    current_time__hr = 0
#    time_f__hr = 8 # calculated based on the length of the TMY data
    PID_frequency = 0.1 # 0.1 seems to be what the PID needs. 0.2 causes oscillations
    NNO_frequency = 1
    num_fans = 12
    
    # setpoint:
    cooled_water_sp__F = 70
    
    # unique cell parameters:
    total_halfcell_flow = 180000/2 # cut in half since we are using half cells and doubling the flow afterward
    percent_flow = np.array([11,10,8,10,7,8,6,9,9,8,9,5])/100 # percent allowcation of flow
#    hotbed_flowrate__gal_min = np.ones(num_fans) * 180000/24 # equal flow
    hotbed_flowrate__gal_min = percent_flow * total_halfcell_flow
    fan_motor_efficiencies = [0.95, 0.99, 0.86, 0.84, 0.94, 0.91, 0.77, 0.76, 0.94, 0.76, 0.95, 0.77] # randomly generated (0.75-1), but rounded and hardcoded for consistency
    
    # distrubances:
    T_amb_array__degF, Relative_humidity_array__percent = TMY_data(Util(), quantity = TMY_quantity)
    time_f__hr = len(T_amb_array__degF)
    T_hotbed_array__degF = random_return_temperatures(time_f__hr)
#    T_amb_array__degF = np.ones(time_f__hr) * 65
#    Relative_humidity_array__percent = np.ones(time_f__hr) * 20
#    T_hotbed_array__degF = np.ones(time_f__hr) * 100 #92
    
    # initialize system:
    fan__percent = 0.90
    fan_array__percent = np.ones(num_fans) * fan__percent
    basin = Basin(fan_array__percent, 
                  T_amb__degF = T_amb_array__degF[0], 
                  Relative_humidity__percent = Relative_humidity_array__percent[0], 
                  T_hotbed__degF = T_hotbed_array__degF[0], 
                  hotbed_flowrate__gal_min = hotbed_flowrate__gal_min,
                  fan_motor_efficiencies = fan_motor_efficiencies)
    
    # initialize neural network object (also used for creating random biases):
    neural_netork = Neural_netork(NN_status, num_fans)
    
    if NN_status == 'on':
        # need to train nueral network
        X, Y = training_data(training_file)
        neural_netork.train(X, Y)

    # set up PID: (0.75 and 0.4 are a good pair, but the approach is a bit slow)
    K_c = -0.03867 *0.75 # basic Tuning = -0.03867
    tau_I = 0.42 * 0.35 # basic Tuning = 0.42
    tau_D = None
    u_lb = 0
    u_ub = 1
    u_bias = fan__percent
    pid = PID(cooled_water_sp__F, K_c, tau_I, tau_D, u_bias=u_bias, u_lb=u_lb, u_ub=u_ub) #, invert_error=True
    
    # set up arrays to store results
    time_array__hr = []
    basin_temp_array__F = []
    PID_output_array__percent = []
    fan_bias_array_array__percent = []
    fan_array_array__percent = []
    total_power_array__kW = []
    T_amb_granular_array__degF = []
    Relative_humidity_granular_array__percent = []
    T_hotbed_granular_array__degF = []
    
    
    while current_time__hr < time_f__hr:
        print('Percent Complete: ' + str(int(round(current_time__hr/time_f__hr * 100))) + '%. Time = ',int(round(current_time__hr)),'hours')
        current_hour_index = int(current_time__hr)
        basin.update_distrubances(T_amb__degF = T_amb_array__degF[current_hour_index], 
                                  Relative_humidity__percent = Relative_humidity_array__percent[current_hour_index], 
                                  T_hotbed__degF = T_hotbed_array__degF[current_hour_index])
        
        # get neural network output or random biases (depending on NN_status)
        fan_bias_array__percent = neural_netork.get_biases(basin_temp__F = cooled_water_sp__F, 
                                                           T_amb__F = T_amb_array__degF[current_hour_index], 
                                                           RH_percent = Relative_humidity_array__percent[current_hour_index], 
                                                           T_hotbed__F = T_hotbed_array__degF[current_hour_index])
        
        # need to update clamping values so the PID can still force the system all the way on or all the way off
#        print('fan_bias_array__percent:',fan_bias_array__percent)
        new_lb = 0 - max(fan_bias_array__percent)
        new_ub = 1 - min(fan_bias_array__percent)
        pid.u_lb = new_lb
        pid.u_ub = new_ub

        NNO_recalc_time__hr = current_time__hr + NNO_frequency
        while current_time__hr < NNO_recalc_time__hr:
            
            # get PID output
            basin_temp__F = basin.get_cooled_water_temp__F()
            PID_output = pid.output(basin_temp__F, PID_frequency)
            
            # advance and update fans
            fan_array__percent = fan_bias_array__percent + PID_output
            
            # clamp the biased outputs:
            fan_array__percent[fan_array__percent>1] = 1
            fan_array__percent[fan_array__percent<0] = 0
            
            basin.advance_one_timestep(fan_array__percent, recalculate_fans = True)
            current_time__hr = basin.t_i__s / basin.util.s_per_hr
            
            # advance without updating fans (the basin takes several advance calls to catch up)
            PID_recalc_time__hr = current_time__hr + PID_frequency
            while current_time__hr < PID_recalc_time__hr:
                basin.advance_one_timestep(fan_array__percent, recalculate_fans = False)
                current_time__hr = basin.t_i__s / basin.util.s_per_hr
            
            # record values:
            time_array__hr.append( current_time__hr )
            basin_temp_array__F.append( basin_temp__F )
            PID_output_array__percent.append( PID_output )
            fan_array_array__percent.append( fan_array__percent )
            fan_bias_array_array__percent.append( fan_bias_array__percent )
            total_power_array__kW.append( basin.get_total_power__kW(fan_array__percent) )
            T_amb_granular_array__degF.append(T_amb_array__degF[current_hour_index])
            Relative_humidity_granular_array__percent.append(Relative_humidity_array__percent[current_hour_index])
            T_hotbed_granular_array__degF.append(T_hotbed_array__degF[current_hour_index])
    
    plt.figure()
    plt.plot(time_array__hr, basin_temp_array__F, '-b', label='basin outlet Temp')
    plt.plot(time_array__hr, T_amb_granular_array__degF, '--b', label='ambient air Temp')
    plt.plot(time_array__hr, T_hotbed_granular_array__degF, '-r', label='tower inlet Temp')
    plt.plot([0, time_array__hr[-1]], [cooled_water_sp__F, cooled_water_sp__F], '-k', label='setpoint')
    plt.title('Basin Temperatures vs Time')
    plt.xlabel('Time (hours)')
    plt.ylabel('Temperature (degF)')
    plt.ylim([50,110])
    plt.legend()
    plt.grid()
    plt.show()
    
    fan_array_array__percent = np.array(fan_array_array__percent)
    plt.figure()
    for i in range(num_fans):
        plt.plot(time_array__hr, fan_array_array__percent[:,i], label='fan #'+str(i+1))
    plt.plot(time_array__hr, PID_output_array__percent, '--k', lw=2, label='PID output')
    plt.title('Fan Percent vs Time')
    plt.xlabel('Time (hours)')
    plt.ylabel('Fan Percent')
    plt.legend()
    plt.grid()
    plt.show()
    
    plt.figure()
    plt.plot(time_array__hr, total_power_array__kW)
    plt.title('Tower Power vs Time')
    plt.xlabel('Time (hours)')
    plt.ylabel('Tower Power (kW)')
    plt.grid()
    plt.show()
    
    time_array__hr = np.array(time_array__hr)
    total_power_array__kW = np.array(total_power_array__kW)
    delta_time_array__hr = time_array__hr[1:] - time_array__hr[:-1]
    total_power_kWh = np.sum(delta_time_array__hr * total_power_array__kW[:-1]) # ignore the last point since the simulation stops there
    
    print('Total energy consumed:',round(total_power_kWh,2),'kWh')
    

    # write csv of results
    print('Writing csv file with results')
    
    fan_bias_array_array__percent = np.array(fan_bias_array_array__percent)
    
    with open('results_'+NN_status+'_'+TMY_quantity+'.csv', 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter=',',
                               quotechar='|', quoting=csv.QUOTE_MINIMAL)
        headers = ['time (hr)',
                   'PID_output',
                   'fan_bias_1',
                   'fan_bias_2',
                   'fan_bias_3',
                   'fan_bias_4',
                   'fan_bias_5',
                   'fan_bias_6',
                   'fan_bias_7',
                   'fan_bias_8',
                   'fan_bias_9',
                   'fan_bias_10',
                   'fan_bias_11',
                   'fan_bias_12',
                   'basin_temp (degF)',
                   'T_amb (degF)',
                   'RH (percent)',
                   'T_hotbed (degF)',
                   'total_power (kW)'
                    ]
        csvwriter.writerow(headers)
        for i in range(len(time_array__hr)):
            row = [time_array__hr[i],
                   PID_output_array__percent[i], # not for NN, just for comparison purposees if needed
                   fan_bias_array_array__percent[i,0],
                   fan_bias_array_array__percent[i,1],
                   fan_bias_array_array__percent[i,2],
                   fan_bias_array_array__percent[i,3],
                   fan_bias_array_array__percent[i,4],
                   fan_bias_array_array__percent[i,5],
                   fan_bias_array_array__percent[i,6],
                   fan_bias_array_array__percent[i,7],
                   fan_bias_array_array__percent[i,8],
                   fan_bias_array_array__percent[i,9],
                   fan_bias_array_array__percent[i,10],
                   fan_bias_array_array__percent[i,11],
                   basin_temp_array__F[i],
                   T_amb_granular_array__degF[i],
                   Relative_humidity_granular_array__percent[i],
                   T_hotbed_granular_array__degF[i],
                   total_power_array__kW[i]
                   ]
            csvwriter.writerow(row)

if __name__ == "__main__":
    main(NN_status = 'learn', TMY_quantity = 'day')
    main(NN_status = 'on', TMY_quantity = 'day', training_file = 'results_learn_month.csv')