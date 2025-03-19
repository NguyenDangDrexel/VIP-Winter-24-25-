# Plot the Temperature, Pressure vs altitude 
from AtmosphericModel import Pressure, Temperature , Density
import numpy as np 
import plotly.express as px 
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import matplotlib.pyplot as plt 
# 1hPa = 100 Pa 
start = 0 
stop = 30
step = 0.1
Pressure_lx = []
Temperature_lx = [] 
Density_lx = [] 
altitude = np.arange (start = start, stop =stop ,step =step)
# Pressure_lx = 
for h in altitude: 
    h = h * 1000
    Pressure_lx.append (Pressure (h)/100)
    Temperature_lx.append (Temperature(h)-273)
    Density_lx.append (Density(h))
""" Create subplots using plotly 
    # fig = make_subplots(rows=1, cols=3, subplot_titles=('Pressure vs Altitude', 'Temperature(K) vs Altitude', 'Density vs Altitude'))

    # # Add scatter plot traces for pressure
    # fig.add_trace(go.Scatter(x= altiude, y=Pressure_lx, mode='markers', name='Pressure', marker=dict(color='blue')), row=1, col=1)

    # # Add scatter plot traces for temperature
    # fig.add_trace(go.Scatter(x= altiude, y=Temperature_lx, mode='markers', name='Temperature', marker=dict(color='red')), row=1, col=2)

    # fig.add_trace(go.Scatter(x= altiude, y=Density_lx, mode='markers', name='Temperature', marker=dict(color='red')), row=1, col=3)

    # # Update layout
    # fig.update_xaxes(title_text='Altitude (m)', row=1, col=1)
    # fig.update_yaxes(title_text='Pressure (Pa)', row=1, col=1)
    # fig.update_xaxes(title_text='Altitude (m)', row=1, col=2)
    # fig.update_yaxes(title_text='Temperature (K)', row=1, col=2)
    # fig.update_xaxes(title_text='Altitude (m)', row=1, col=3)
    # fig.update_yaxes(title_text='Density (kg/m3)', row=1, col=3)

    # fig.update_layout(title='Pressure and Temperature vs Altitude')
    # fig.show() 
""" 


# Create a figure using matplotlib 
fig, axs = plt.subplots(2,1, figsize=(8, 14))  # 3 rows, 1 column
# fig.suptitle ("US Standard Atmosphere Model")
# fig.subplots_adjust(hspace=0.5)  # Increase vertical space

# First subplot:
axs[0].plot(altitude, Pressure_lx, label="Pressure", color="blue")
# axs[0].set_title("Pressure versus Altitude")
axs[0].set_xlabel("Altitude (km)", fontsize = 20 )
axs[0].set_ylabel("Pressure (hPa)", fontsize = 20)
axs[0].legend()
axs[0].grid(True)
plt.setp(axs[0].get_xticklabels(), fontsize=20)
plt.setp(axs[0].get_yticklabels(), fontsize=20)
plt.setp(axs[1].get_xticklabels(), fontsize=20)
plt.setp(axs[1].get_yticklabels(), fontsize=20)




# Second subplot: 
# axs[1].plot(altitude, Density_lx, label="Density", color="green")
# # axs[1].set_title("Temperature vs Altitude")
# axs[1].set_xlabel("Altitude (km)")
# axs[1].set_ylabel("Density (Kg/m$^3$)")
# axs[1].legend()
# axs[1].grid(True)

# Third subplot: 
axs[1].plot(altitude, Temperature_lx, label="Temperature", color="red")
# axs[1].set_title("Temperature versus Altitude")
axs[1].set_xlabel("Altitude (km)", fontsize = 20)
axs[1].legend()
axs[1].grid(True)
axs[1].set_ylabel("Temperature (C)", fontsize = 20)



# Adjust layout
plt.tight_layout()

# Show the figure
plt.show()