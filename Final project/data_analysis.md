```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

DO_column = 3 #index of DO column

###Trial 1
DO_first_row_t1 = 4
dirpath_trial1 = "/Users/Emily/Desktop/4530/Final project/datalog 5-5-2019.xls"
DO_data_t1=epa.column_of_data(dirpath_trial1, DO_first_row, 3, -1, 'mg/L')
DO_time_data_t1 = (epa.column_of_time(dirpath_trial1,DO_first_row,-1)).to(u.s)
# Plot the raw data
plt.plot(DO_time_data_t1, DO_data_t1, '-')
plt.title('DO concentration over time, 5/5')
plt.xlabel('time (s)')
plt.ylabel('DO concentration (mg/L)')
plt.savefig('/Users/Emily/Desktop/4530/Final project/5-5_data.png')
plt.show()

###Trial 2
DO_first_row_t2 = 309
dirpath_trial2 = "/Users/Emily/Desktop/4530/Final project/datalog 5-6-2019.xls"
DO_data_t2=epa.column_of_data(dirpath_trial2, DO_first_row_t2, 3, -1, 'mg/L')
DO_time_data_t2 = (epa.column_of_time(dirpath_trial2,DO_first_row_t2,-1)).to(u.s)
# Plot the raw data
plt.plot(DO_time_data_t2, DO_data_t2, '-')
plt.title('DO concentration over time, 5/6')
plt.xlabel('time (s)')
plt.ylabel('DO concentration (mg/L)')
plt.savefig('/Users/Emily/Desktop/4530/Final project/5-6_data.png')
plt.show()
