
import pandas as pd
import numpy as np

def TMY_data(util, quantity = 'day'):
    '''
    can choose 'day', 'week', 'month', or 'year'
    all start on May 1 (except year obviously)
    '''
    data = pd.read_csv('Hanksville.csv')
    
    temperature_array__degC = data['Dry-bulb (C)'].values.astype(np.float)
    
    temperature_array__degF = util.degC_to_degF(temperature_array__degC)
    relativeHumidity_array__percent = data['RHum (%)'].values.astype(np.float)
    dates = data['Date (MM/DD/YYYY)'].values
    times = data['Time (HH:MM)'].values
    
    if quantity == 'year':
        print('Using one year of TMY data')
        return temperature_array__degF, relativeHumidity_array__percent
    elif quantity == 'day':
        start = 2880
        stop = 2904
    elif quantity == 'week':
        start = 2880
        stop = 3048
    elif quantity == 'month':
        start = 2880
        stop = 3624
    else:
        start = 2880
        stop = 2883
        
    print('Using one',quantity,'of TMY data, from',dates[start],times[start],'through',dates[stop-1],times[stop-1])
    
    return temperature_array__degF[start:stop], relativeHumidity_array__percent[start:stop]

def training_data(filename):
    data = pd.read_csv(filename)
    Y = data['total_power (kW)'].values
    X = data.drop(['time (hr)','PID_output','total_power (kW)'], axis=1).astype(np.float) 
    return X, Y

if __name__ == "__main__":
    from util import Util
    
    temp_1, rh_1 = TMY_data(Util(),quantity = 'day')
    temp_2, rh_2 = TMY_data(Util(),quantity = 'week')
    temp_3, rh_3 = TMY_data(Util(),quantity = 'month')
    temp_4, rh_4 = TMY_data(Util(),quantity = 'year')