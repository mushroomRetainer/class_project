
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.animation as animation

def plot_T_profile(solution,n_x, n_y, basin_length__m, basin_height__m,
                   title,normalized_C_water_arrays, normalized_C_vapor_arrays, 
                   normalized_T_water_arrays, normalized_T_air_arrays, z_array, util):
    
    nrows, ncols = n_y, n_x
    solution = util.degK_to_degF(solution) # convert to degF for plotting
    grid = solution.reshape((nrows, ncols))
    
    x_min = 0
    x_max = basin_length__m
    y_min = 0
    y_max = basin_height__m
    plt.figure()
    
    plt.subplot(241)
    plt.title('C_water')
    plt.ylabel('Nondimensional Cell Height')
    plt.plot(np.transpose(normalized_C_water_arrays), z_array)
    plt.xlabel('Ratio to inlet')
    
    plt.subplot(242)
    plt.title('C_vapor')
    plt.gca().set_yticklabels(['']*10)
    plt.plot(np.transpose(normalized_C_vapor_arrays), z_array)
    plt.xlabel('Ratio to inlet')
    
    plt.subplot(243)
    plt.title('T_water')
    plt.gca().set_yticklabels(['']*10)
    plt.plot(np.transpose(normalized_T_water_arrays), z_array)
    plt.xlabel('Ratio to inlet')
    
    plt.subplot(244)
    plt.title('T_vapor')
    plt.gca().set_yticklabels(['']*10)
    plt.plot(np.transpose(normalized_T_air_arrays), z_array)
    plt.xlabel('Ratio to inlet')
    
    plt.subplot(212)
    
    plt.imshow(grid, extent=(x_min, x_max, y_max*10, y_min), # flip the y-axis
               interpolation='nearest', cmap=cm.viridis)
    plt.title(title)
    plt.colorbar()
    plt.xlabel('x-coordinate (m)')
    plt.ylabel('y-coordinate (m*10, to make graph visible)')
    
    plt.show()
    plt.pause(0.0001) #allows you to see plots as they are plotted, rather than all at the end
    
#def get_animation_info():
#    
#    fig = plt.figure("Animation")
#    
#    ax = fig.add_subplot(111)
#    plt.xlabel('x-coordinate (Length/deltaX = '+str(round(basin_length__m))+' m/'+str(round(deltaX,2))+' m)')
#    plt.ylabel('y-coordinate (Height/deltaY = '+str(round(basin_height__m))+' m/'+str(round(deltaY,2))+' m)')
#    plt.title('Temperature Profile (Dark Blue = 65 degF, Yellow =  100 degF)')
    