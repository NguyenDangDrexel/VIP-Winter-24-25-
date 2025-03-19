"""
1. Run simulation for every 0.2 seconds 
2. Get the altitude, time from simulation 
3. Send altitude,time to Arduino + receive response 
4. Python: if we do not receive control_signal from Arduino: keep the latest signal. To prevent overshooting, if control_signal the same for 2 or 3 control_interval => set it to 0 
"""
"""
                                            READ ME
    THIS SIMULATION IS RUNNING 
    NOW, TRY TO SEND DATA TO ARDUINO 
... 





"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os  
import random 
from tqdm.auto import tqdm 
from scipy.integrate import odeint 
from scipy.signal import medfilt 
import statistics 

import time 
import serial 

#-------------------------------- import neccessary modules ------------------ 

import sys
path = r"C:\Drexel\Drexel\2023\Courses\Summer 2024 - SGN\VIP program - balloon project\Weather-Balloon-Drexel\NEBP_project"
sys.path.append(path) 

#-------------------------------------------------------------------------------------------- 
# IDEA FOR TESTING ARDUINO CODE: WRITE THE CODE THAT WE IMPLEMENT ON ARDUINO IN PYTHON => TEST it 
N = 1
S_xy = S_x = S_y = S_xx=0 
Sum_xx = Sum_xy= sum_int = 0
control_interval = 20 
Kp = 0.4 
Ki = 0 
Kd = 0.1
time_constant = 1142 
target_velocity = 0 
sum_int = 0
error = 0 
error_prev =0
sum_int =0
d_error =0
control_signal =0
altitude_sp_min = 20000
altitude_sp = 28000
altitude_sp_max = 30000
avg_velocity = 7 
pressure_chamber_slope = 1
"""

x = t[index]
    y = altitude 
    S_xy += x*y
    S_x += x 
    S_y +=y 
    S_xx += x**2 

Sum_xx = S_xx - (S_x**2 / N) 
Sum_xy = S_xy - (S_x * S_y)/N 

velocity = Sum_xy / Sum_xx
"""
# the function in can access global varible, but it cannot change this variable if we don't specify that "global x" inside the function  
def Controller_Arduino (altitude_i,time_i,PID): 
    global  S_xy, S_x, S_y, S_xx, Sum_xx, Sum_xy, sum_int,control_signal 
    global control_interval, Kp, Ki, Kd, time_constant, target_velocity
    global N 
    global error,sum_int,error_prev
    time_constant = int ((altitude_sp - altitude_sp_min)/ avg_velocity) 
    if (PID !='Active'): # if we use this way => MUST CHECK THE NUMBER OF DATA POINT EXTREMELY CAREFUL WHEN IMPLEMENT INTO ARDUINO 
    # if (N <=101 ):
        # print (f"Altitude_{N}")

        x = time_i 
        y = altitude_i 
        S_xy += x*y
        S_x += x 
        S_y +=y 
        S_xx += x*x 
        N +=1 
        if ((altitude_sp_min <=altitude_i< altitude_sp)):
                # velocity_sp = avg_velocity - abs((altitude_i - altitude_sp_min) / time_constant) + 0 
                velocity_sp = (altitude_sp - altitude_i)/time_constant
                result = velocity_sp
        else:
            result = 0 
    # elif (PID == 'Active'):
    else :
        print (f"Number of data points received: {N}")
        Sum_xx = S_xx - ((S_x*S_x)/(N-1) ) # S_xx - ((S_x*S_x)/(N) ) THERE IS AN EXTREMELY BIG DIFFERENT BETWEEN N-1 AND N => MUST BE CAREFUL WHEN IMPLEMENT. IN REAL LIFE, IT MUST BE N, NOT N-1 
        Sum_xy = S_xy - ((S_x * S_y)/(N-1))
        velocity = (Sum_xy/Sum_xx)*pressure_chamber_slope 
        estimated_altitude = altitude_i
  
    #   // PID CONTROLLER 
        if (estimated_altitude <altitude_sp_min):
            control_signal =0 
            velocity_sp = 0 
      
        elif ((altitude_sp_min <=estimated_altitude< altitude_sp)):

            velocity_sp = avg_velocity - abs((estimated_altitude - altitude_sp_min) / time_constant) + target_velocity 
            # // CALCULATE CONTROL SIGNAL 
            error = velocity-velocity_sp  
            sum_int += error * control_interval  
            d_error = (error - error_prev) / control_interval  
            control_signal =  Kp * error + Ki * sum_int + Kd * d_error 
            control_signal = max(0.0, min(1.0, control_signal)) 
            error_prev = error 
        
    #   // when balloon go beyond setpoint
        else:
            velocity_sp =0 
            # // CALCULATE CONTROL SIGNAL 
            error = velocity-velocity_sp  
            sum_int += error * control_interval  
            d_error = (error - error_prev) / control_interval  
            control_signal =  Kp * error + Ki * sum_int + Kd * d_error 
            control_signal = max(0.0, min(1.0, control_signal)) 
            error_prev = error 

        S_xy = S_x = S_y = S_xx = Sum_xx = Sum_xy = 0.0
        N = 1
        velocity = round (velocity,2)
        result = [control_signal, velocity] 

    return result 

time_space = np.linspace(0, 1000, 100)  # 100 data points for x from 0 to 10.2

# Generate noise (normal distribution)
noise = np.random.normal(0, 1, 100)  # Mean 0, standard deviation 1

# Compute y = ax + noise
altitude = (5.5 * time_space) + noise 

v_lx = [] 
# Testing linear regression 
# for index, time_i in enumerate (time_space):
#     velocity  = Controller_Arduino (altitude_i= altitude [index],time_i = time_i, PID = '')
#     v_lx.append (velocity)
#     # print (v_lx)

# import matplotlib.pyplot as plt 

# plt.plot (time_space,v_lx)
# plt.show ()


