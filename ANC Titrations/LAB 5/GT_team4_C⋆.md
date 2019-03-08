````python

from aguaclara.core.units import unit_registry as u
import aguaclara.core.materials as mat
import aguaclara.core.constants as con
import aguaclara.core.utility as ut
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.core.physchem as pc


import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate, integrate


#need to know average water temperature, barometric pressure(partial pressure of oxygen in atmospheres)
C = epa.O2_sat(1*u.atm , 295*u.kelvin)

C





from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import collections
import os
from pathlib import Path

# This code is included because there is a bug in the version of this code that is in epa.
def aeration_data(DO_column, dirpath):

    #return the list of files in the directory
    filenames = os.listdir(dirpath)
    #extract the flowrates from the filenames and apply units
    airflows = ((np.array([i.split('.', 1)[0] for i in filenames])).astype(np.float32))
    #sort airflows and filenames so that they are in ascending order of flow rates
    idx = np.argsort(airflows)
    airflows = (np.array(airflows)[idx])*u.umole/u.s
    filenames = np.array(filenames)[idx]

    filepaths = [os.path.join(dirpath, i) for i in filenames]
    #DO_data is a list of numpy arrays. Thus each of the numpy data arrays can have different lengths to accommodate short and long experiments
    # cycle through all of the files and extract the column of data with oxygen concentrations and the times
    DO_data=[epa.column_of_data(i,0,DO_column,-1,'mg/L') for i in filepaths]
    time_data=[(epa.column_of_time(i,0,-1)).to(u.s) for i in filepaths]
    aeration_collection = collections.namedtuple('aeration_results','filepaths airflows DO_data time_data')
    aeration_results = aeration_collection(filepaths, airflows, DO_data, time_data)
    return aeration_results
#delete 550
# The column of data containing the dissolved oxygen concentrations

DO_column = 2
dirpath = "/Users/Emily Liu/Desktop/CEE 4530/Aeration"
filepaths, airflows, DO_data, time_data = aeration_data(DO_column, dirpath)


# Plot the raw data

for i in range(airflows.size):
  plt.plot(time_data[i], DO_data[i],'-')
plt.xlabel(r'$time (s)$')
plt.ylabel(r'Oxygen concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(airflows.magnitude)
plt.show()

#delete data that is less than 2 or greater than 6 mg/L

DO_min = 2 * u.mg/u.L
DO_max = 6 * u.mg/u.L

DO_data[12]
idx_start = (np.abs(DO_data[5]-DO_min)).argmin()
for i in range(airflows.size):
  if i == 12:
    print('Omit these data')
  else:
    idx_start = (np.abs(DO_data[i]-DO_min)).argmin()
    idx_end = (np.abs(DO_data[i]-DO_max)).argmin()
    time_data[i] = time_data[i][idx_start:idx_end] - time_data[i][idx_start]
    DO_data[i] = DO_data[i][idx_start:idx_end]


slopes = np.zeros(24)
for i in range(airflows.size):
  if i == 12:
    print('Omit these data')
  else:
    x= time_data[i][0]-time_data[i]
    y= np.log((C-DO_data[i])/(C-DO_data[i][0]))
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    slopes[i] =slope

slopes

###We pick a representative data set###
"We pick the data set for flow rate = 800"
flowrate = airflows[19]
DO = DO_data[19]
t = np.linspace(0,70) *u.second
kvl = slopes[18] * 1/u.second
c0 = DO_data[19][0]
#C is C*
t0 = time_data[19][0]
t19 = time_data[19]

conc = C -(C-c0)*np.exp(-kvl*(t-t0))

plt.plot(t, conc, 'r-', t19, DO, 'b-' )
plt.savefig("/Users/Emily Liu/Desktop/CEE 4530/lab5q5")
###Plotting kvl as a function of air flow rate#####

plt.plot(airflows, slopes, 'go')
plt.savefig("/Users/Emily Liu/Desktop/CEE 4530/lab5q6")

####plotting OTE as a function of flow rate###
OD = 6 * u.mg/u.L #the oxygen deficit
fo2 =0.21 #molar fraction of
V = 750 *u.mL #volume of the reactor
MWO2 = 32 *u.g/u.mole
OTE = V*slopes*OD/(fo2*airflows*MWO2)
plt.plot(airflows, OTE, 'ko')
plt.savefig("/Users/Emily Liu/Desktop/CEE 4530/lab5q7")
