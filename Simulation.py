"""
Author: Ben Dang - Nguyen Dang 
Email: ppd34@drexel.edu 


This is the simulation of the balloon. Only run on the laptop 

"""
# ---------------------- instruction to use simulation ---------------- 
"""
    y = odeint(venting, y0, [0, delta_t], args=(u,wind_speed  )) # => real value of altitude 
    inputs: 
    delta_t: time step 
    u: control signal 
    wind_speed: CAUTION: DO NOT PUT WIND_SPEED FOR EVERY TIME STEP. USE time_i % 10 == 0 to insert wind_speed every 10 seconds 
        
"""


avg_velocity = 7 
# Initial volume of Helium - 
V0 = 6.6 # m3
# mass of the balloon (helium inside + payload + latex)
m = 5.5
# Discharged rate constant - needed 
discharged_constant = 0.2


import matplotlib.pyplot as plt
import numpy as np
import random 
from tqdm.auto import tqdm 
from scipy.integrate import odeint 

# ---------------------------- READ ME --------------------------------------------


#-------------------------------- import neccessary modules ------------------ 

from Models.ThermodynamicModels import Balloon
from Models.Atmospheric_models.AtmosphericModel import Pressure, Temperature, Density 
from Controller.PIDPython import Controller_Arduino # this the python script 
#-------------------------------------------------------------------------------------------- 

# --------------------------------physics parameters -----------------------------

# Molar mass of Helium (kg/mol) <=> rho (kg/m^3)
muy_he = 0.004 #(kg/mol)
# Density of Helium 
rho_he = 0.1784 # kg/m^3 
# Drag coefficient  
Cd = 0.44
# Gas constant 
R = 8.31 
# Gravitational constant 
g = 9.8 

# -------------------------------- balloon parameters -----------------------------
# Surface tension pressure; from 100 to 200 Pa 
difference_pressure = 200 # (Pa)

# Initial mol of helium 
initial_mol_helium = Pressure (0)* V0 / (R *Temperature (0))

# -------------------------------- air_drag parameters -----------------------------
C = 0.5 # m_air = C * rho * V  
include_m_air  = True # True if we want to 

# -------------------------------- venting parameters -----------------------------

# Diameter of venting system (m): 
diameter = 0.095 
# Area of the venting 
area_of_orifice = 1/4 * np.pi * diameter**2 * 1 

  # 
# -------------------------------- simulation parameters -----------------------------
# case 1: constant volume balloon 
# case 2: constant temperature 
# case 3: adiabatic 
# case 4: temperature inside = temperature ouside 

initial_params = [101325,initial_mol_helium,293, 0] # [P0,n0,T0,delta_P]

def venting(y:list, t,u:float, wind_speed) -> list :
    """
    parameters 

    Returns: list 
        return the velocity, acceleration, dn/dt at time t 
        [ dydt = [xdot, xddot, ndot] ]
    
    CAUTION: DO NOT USE CASE 3. CASE 3 STILL WORKS, BUT ... NOV 5 2024 
    
    """
    # vector state : y = [x,xdot,n]
    x = y[0] 
    xdot = y[1]
    n = y[2] # mol of Helium 
    
    # Choose case to do simulation 
    case = 4
    if case == 2: 
        initial_params =[101325,initial_mol_helium,293, 0] 
        T = Balloon.iosthermal_balloon (x,initial_params= initial_params)['Temperature']
        P = Balloon.constant_volume_balloon (x,initial_params= initial_params)['Pressure']
        V = R * T *n /P # ONLY for this case 
        A = Balloon.iosthermal_balloon (x,initial_params= initial_params)['Area']
    elif case == 3: # this case is not accurate => do not use 
        gamma = 1.4 
        initial_params =[101325,n,293, 0] 
        # P = Balloon.iosthermal_balloon (x,initial_params= initial_params)['Pressure']
        P = Pressure(x) 
        V = ((n/initial_mol_helium)* (initial_params[0]/P))** (1/gamma) * (initial_mol_helium*R*293/101325)
        A = np.pi * ((3/4) * (V/np.pi))**(2/3) 
    elif case  == 4: 
        initial_params = [101325,initial_mol_helium,293, 0] # [P0,n0,T0,delta_P]
        P = Pressure (x)
        T = Balloon.simultaneous_temperature_balloon (x,initial_params= initial_params)['Temperature'] # T = t(h)
        V = R * T *n /P 
        A = np.pi * ((3/4) * (V/np.pi))**(2/3) 


    # Use ideal gas law to define volume 
    # V = n*R*T/ P  # NOTE: in this line, n changes => cannot use the V from Balloon, which is the closed balloon 

    # air_drag 
    m_air = C * Density(x) * V 
    
    # F = ma - dynamic equation 
    rel_velocity = xdot + wind_speed # add it from 5 to 10 seconds 
    if xdot > 0: 
        xddot = (-m * g + Density(x) * V * g) / (m + m_air) - 0.5 * Density(x) * rel_velocity**2 * A * Cd  / (m + m_air)  
    else: 
        xddot = (-m * g + Density(x) * V * g) / (m + m_air) + 0.5 * Density(x) * xdot**2 * A * Cd / (m + m_air) 
    # adding noise to Cd to simulate vertical movement of air - 1 option (xdot-noise) **2. This noise can be random, or the function of altitude/time 
    # Define the function to calculate volume discharged rate 
    Q = u*-discharged_constant * area_of_orifice * ((2 * R) / (muy_he) * difference_pressure) ** (1/2) * ((Temperature(x)) / (Pressure(x))) ** (1/2)
    rho_he_var = (P*muy_he)/(R*T)
    # Define ndot = Q * rho / muy 
    ndot = (Q * rho_he_var)/ muy_he 
    # P, T, n => V => diameter (drag equation?) ??? 
    # dydt = [x,xdot,n]
    if x <=35000: 
        dydt = [xdot, xddot, ndot] 
    else: 
        dydt = [0,0,ndot]

    return dydt 

# Resolution of sensor 
def resolution_sensor (h, resolution = 0.2): 
    """
    This function simulates the resolution of the altitude sensor. 
    The default resolution is 0.2 meter 
    Params: 
        h: noised altitude 
        resolution: resolution of the sensor 

    Returns: 
        h_resolution 
    """
    h_resolution = round (h/resolution) * resolution 
    return h_resolution
# ------------------------------------------------------------------------------------------------------
# Random noise of sensor 
def random_noise_sensor (h:float,amplitude:int):
    """
    This function simulates the random noise of sensor by using random.randint (begin,end) method. 
    This method generates random number in the range begin -> end.
    Params: 
        h: altitude from odeint (), which is the "true" altitude of the balloon 
        amplitude
    Return: 
        h_noise: altitude + random noise 
    """
    h += random.randint (-1,1) * amplitude  
    return h 
# ---------------------------------------------------------------------------------------------------------
# random vertical velocity 
def random_vertical_velocity (amplitude = 1):
    """
    This function generates the vertical velocity
    Params: 
        amplitude: default value = 1 
    """

    vertical_velocity = random.randint (-1,1) * amplitude
    return vertical_velocity


# ---------------------------------- RUN SIMULATION ---------------------------- 
# ---------- PARAMS 
ODEINT_initial_params = [0,0,initial_mol_helium]# [x,xdot,n]
time_params = [5000,1] #[time_simulation,time_step] 
PID_params = [0.4,0,0.1] # [kp,ki,kd] 
control_params = [20,0,3] # [control_interval, velocity_final_setpoint, min_time_of_valve] 
altitude_setpoints = [20000,28000,30000] # [al_min, al_target,al_max ] 
noise = [0.2,3] # sensor noise, wind noise 
time_window = 20 #seconds 

altitude_sp_min = altitude_setpoints[0]

altitude_sp = altitude_setpoints[1]

# ------------------------------------------------------------- 
# response time of the sensor
sensor_pulses = 5 # (Hz)
# PID time parameters

simulation_time = 6000                                 # final time for simulation;  
delta_t = 0.2                                              # time distep for odeint - do not change 
nsteps = int (simulation_time/delta_t) + 1             # number of time steps 
ts = np.arange(0, simulation_time,delta_t)  
u = 0 
# ---------------------- list for plotting 
raw_velocity = [] # store the real value 
raw_altitude = [] 
# lists for storing the results
venting_status = [] # u = valve % open
v_list = []
sps = []
sps2 =[]
v_list_measure =[]
altitude_lx = []
venting_time = np.zeros (nsteps) 
# ----------------------------------------- TARGET SETPOINT ------------------ 

# -------------------- USE FOR LOOPS TO RUN SIMULATION 
y0 = ODEINT_initial_params 
velocity_Arduino = 0 
# time_values = np.arange(0, 5000, 0.2)
for i in tqdm (ts):# 0.1 seconds i:0.2;0.4;0.6=> 5hz . 
    if i==0: 
        estimated_altitude =0 
    if estimated_altitude < altitude_sp_min:
            u =0 
            velocity_sp = 10
    elif altitude_sp_min <= y[-1][0] < altitude_sp:
        velocity_sp = avg_velocity- abs ((estimated_altitude  - altitude_sp_min))/1142   + 0  # this setpoint change with time 

    # ----------------------------- read the value of u from Arduino 
    altitude_t = round(y0[0], 2)
    
    if int (i *5) % 100 ==0 and int (i*5) !=0: # Since i is 0.2;0.4;0.6 => o = 22.2 => (int i) = 20 => wrong 
        # Event -trigger PID 
        # simulate what happend when pressure chamber reach minimum pressure 
        PID = 'Active' 
        
        u,velocity_Arduino = Controller_Arduino (altitude_i = altitude_t,time_i = i,PID = PID )
        print (f'Control Signal: {u}')
        Controller_Arduino (altitude_i = altitude_t,time_i = i,PID  = ' ') # call this function to acculmulate altitude, time 

    else:
        Controller_Arduino (altitude_i = altitude_t,time_i = i,PID  = ' ')
    if i ==0:
        altitude_lx.append (altitude_t)
        # sps.append (velocity_sp2)
        sps2.append (velocity_sp)
        v_list.append (y0[1]) # store the velocity for plottin 
        venting_status.append (u)
        v_list_measure.append (velocity_Arduino)
    else: 
        altitude_lx.append (altitude_t)
        # sps.append (velocity_sp2)
        sps2.append (velocity_sp)
        v_list.append (y[-1][1]) # store the velocity for plottin 
        venting_status.append (u)
        v_list_measure.append (velocity_Arduino) 

    # ------------------------------- add some noise to the wind speed ------------------ 
    if i % 5 ==0: # 2 => too much noise
        period = 100
        omega = 2*np.pi/period
        wind_speed = 0
        for n in range (1,4,1): # try to add multiple sine functions with different period together 
            wind_speed += (np.sin(omega*n *i)) * random_vertical_velocity (0.5) # now, only one 
    # wind_speed =0
# ------------------- run the odeint () 
    
    # wind_speed = 0 
    y = odeint(venting, y0, [0, delta_t], args=(u,wind_speed  )) # => real value of altitude
# --------------------------------- update the current states 
    y0 =[y[-1][0],y[-1][1],y[-1][2]] 
    # real value of the altitude
    altitude = y[-1][0] 
# # --------------------- noise from the sensor -----------------------
    # random noise of the sensor
    altitude = random_noise_sensor (altitude,amplitude= 0.2) 
    # resolution of the sensor (rounding) 
    altitude = resolution_sensor (altitude)
    # altitude = resolution_sensor (altitude)
    raw_altitude.append (altitude)
# --------------------------------------------------------- 
    estimated_altitude = altitude 



import matplotlib.patches as patches
start_point1 = (0,0)
start_point2 = (0,20000)
start_point3 = (0,28000)

width = 6000 
height1 = 20000
height2 = 8000
height3 = 2000
# ------------------------------------------- PLOTTING -------------------------------- 
fig, (ax3, ax1,ax2) = plt.subplots(3, 1, figsize=(10, 8))
start_time = 0
end_time = 5000 * 10
simulation_time = time_params[0]
if end_time <= simulation_time: 
    end_time = end_time 
else: 
    end_time = simulation_time

start_time = 17500
# Plot the first set of data on the first subplot
ax1.plot(ts[start_time:], v_list[start_time:], label='velocity')
ax1.plot(ts[start_time:], v_list_measure[start_time:], label='average velocity',color = 'red')


ax1.plot (ts[start_time:],sps2[start_time:],label = 'velocity setpoint')  
# ax1.plot (ts[start_time:],sps2[start_time:],label = 'setpoint from simulation') # PROBLEM??? 

# ax1.plot (ts[start_time:],estimated_velocity[start_time:],label = 'estimated velocity')

ax1.set_xlabel('Time (s)',fontsize = 20)
ax1.set_ylabel('Velocity (m/s)',fontsize = 20 )
ax1.legend(fontsize= 20,loc='upper right')
ax1.grid(True)
plt.setp(ax1.get_xticklabels(), fontsize=20)
plt.setp(ax1.get_yticklabels(), fontsize=20)

# Plot the second set of data on the second subplot
 
ax2.plot(ts[start_time:], venting_status[start_time:], color='orange')
ax2.set_xlabel('Time (s)',fontsize = 20)
ax2.set_ylabel('Control Signal ',fontsize = 20)
# ax2.legend(fontsize = 20 )
ax2.grid(True)
# Plotting the third graph 
ax3.plot(ts, altitude_lx, label='Altitude (km)', color='red')
# ax3.plot(ts[start_time:], estimated_altitude[start_time:], label=' estimated altitude', color='blue')
# -------------------------- plotting venting region - start from here ---------------------
# add the control regions 
# ---------------- Region 1 - green region 
rect1 = patches.Rectangle(start_point1, width, height1, 
                          facecolor='green', edgecolor=None, alpha=0.7)
ax3.add_patch(rect1)
# ---------------- Region 2 - Yellow region 
rect2 = patches.Rectangle(start_point2, width, height2, 
                          facecolor='yellow', edgecolor=None, alpha=0.7)
ax3.add_patch(rect2)
# ---------------- Region 3 - Red region 
rect3 = patches.Rectangle(start_point3, width, height3, 
                          facecolor='red', edgecolor=None, alpha=0.7)
ax3.add_patch(rect3)
# Add text 
ax3.text(1000,20000,  r'Free Rise ($u_{sp} = u_{max}$)', color='black', 
        fontsize=20, ha='left', va='top')  # Center the text
ax3.text(1000,28000, 'Active Control', color='black', 
        fontsize=20, ha='left', va='top')
ax3.text(1000,30000, r'Above Setpoint ($u_{sp} = 0$)', color='black', 
        fontsize=20, ha='left', va='bottom')

# Set custom y-tick positions and labels
yticks = [0, 10000, 20000, 30000]
ytick_labels = ['0','10','20','30']

# Apply the custom y-ticks and labels
ax3.set_yticks(yticks)
ax3.set_yticklabels(ytick_labels, fontsize=20)

# Create the dash line 
x= np.linspace (0,6000,50)
y1 = np.linspace (20000,20000,50)
y2 = np.linspace (28000,28000,50)
y3 = np.linspace (30000,30000,50)
# -------------------------- plotting venting region - end here ---------------------
ax3.plot (x,y1, color = 'blue', linestyle = '-')
ax3.plot (x,y2, color = 'blue', linestyle = '-')
ax3.plot (x,y3, color = 'blue', linestyle = '-')

ax3.set_xlabel('Time (s)',fontsize = 20)
ax3.set_ylabel('Altitude (km)',fontsize = 20)
ax3.set_ylim(0, 32000)
# ax3.legend(fontsize = 20,loc='upper left')
ax3.grid(True)
plt.tight_layout()

plt.setp(ax2.get_xticklabels(), fontsize=20)
plt.setp(ax2.get_yticklabels(), fontsize=20)
plt.setp(ax3.get_xticklabels(), fontsize=20)
plt.setp(ax3.get_yticklabels(), fontsize=20)





plt.show() 

# Region 1: , color green - can go freely  
start_point1 = (0,0)
width = 6000 
height1 = 20000
# Region 2: yellow light: hit the break slightly 
start_point2 = (0,20000)
width = 6000 
height2 = 8000
# Region 3: yellow red light: hit the break as hard as you can 
start_point2 = (0,20000)
start_point3 = (0,28000)

width = 6000 
height1 = 20000
height2 = 8000
height3 = 2000

