
```python
from aguaclara.core.units import unit_registry as u
u.define('equivalent = mole = eq')
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.core.utility as ut
from scipy import optimize
import numpy as np
import math
from scipy import special
import aguaclara.core.physchem as pc
from aguaclara.research.procoda_parser import *
from scipy.optimize import curve_fit
import collections
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

Gran_datat0 = '/Users/Emily/Desktop/4530/ANC Titrations/lab 3 titration time 0.tsv'


# The epa.Gran function imports data from your Gran data file as saved by ProCoDA.
# The epa.Gran function assigns all of the outputs in one statement
V_titrant, pH, V_Sample, Normality_Titrant, V_equivalent, ANC = epa.Gran(Gran_data)
ANC0 = ANC
#plotting the pH vs. equivalent volume of Titrant
#Question: do we have to call Gran to get the pH and V Titrant

plt.plot(V_titrant, pH, 'r')
plt.axvline(x = 0.6, color = 'y')
plt.axvline(x = 0.8, color = 'y')
plt.xlabel("Titrant Volume (mL)")
plt.ylabel("pH")
plt.savefig('/Users/Emily/github/el545-4530/ANC Titrations/pH vs V_titrant')
plt.tight_layout()
plt.show()


#Define the gran function.
def F1(V_sample,V_titrant,pH):
  return (V_sample + V_titrant)/V_sample * epa.invpH(pH)
#Create an array of the F1 values.
V_titrant
F1_data = F1(V_Sample,V_titrant,pH)
#By inspection I guess that there are 4 good data points in the linear region.
N_good_points = 3
#use scipy linear regression. Note that the we can extract the last n points from an array using the notation [-N:]
slope, intercept, r_value, p_value, std_err = stats.linregress(V_titrant[-N_good_points:],F1_data[-N_good_points:])
#reattach the correct units to the slope and intercept.
intercept = intercept*u.mole/u.L
slope = slope*(u.mole/u.L)/u.mL
V_eq = -intercept/slope
ANC_sample = V_eq*Normality_Titrant/V_Sample
print('The r value for this curve fit is', r_value)
print('The equivalent volume was', V_eq)
print('The acid neutralizing capacity was', ANC_sample.to(u.meq/u.L))

#The equivalent volume agrees well with the value calculated by ProCoDA.
#create an array of points to draw the linear regression line
x=[V_eq.magnitude,V_titrant[-1].magnitude ]
y=[0,(V_titrant[-1]*slope+intercept).magnitude]
#Now plot the data and the linear regression
plt.plot(V_titrant, F1_data,'o')
V_titrant
F1_data
plt.plot(x, y,'r')
plt.xlabel('Titrant Volume (mL)')
plt.ylabel('Gran function (mole/L)')
plt.legend(['data'])
plt.savefig('/Users/Emily/github/el545-4530/ANC Titrations/time 0 ')
plt.tight_layout()
plt.show()


#Plotting the measured ANC
Kw = 10**(-14) * (u.mole/u.L)**2
K1_carbonate = 10**(-6.37)*u.mol/u.L
K2_carbonate = 10**(-10.25)*u.mol/u.L
K_Henry_CO2 = 10**(-1.5) * u.mole/(u.L*u.atm)
P_CO2 = 10**(-3.5) * u.atm


def invpH(pH):
    #This function calculates inverse pH
    #Parameters

    return 10**(-pH)*u.mol/u.L


def alpha0_carbonate(pH):
    #This function calculates the fraction of total carbonates of the form
    #H2CO3
    #Parameters


    alpha0_carbonate = 1/(1+(K1_carbonate/invpH(pH)) *
                            (1+(K2_carbonate/invpH(pH))))
    return alpha0_carbonate


def alpha1_carbonate(pH):
    #This function calculates the fraction of total carbonates of the form
    #HCO3-
    #Parameters

    alpha1_carbonate = 1/((invpH(pH)/K1_carbonate) + 1 +
                          (K2_carbonate/invpH(pH)))
    return alpha1_carbonate


def alpha2_carbonate(pH):
    #This function calculates the fraction of total carbonates of the form
    #CO3-2
    #Parameters

    alpha2_carbonate = 1/(1+(invpH(pH)/K2_carbonate) *
                            (1+(invpH(pH)/K1_carbonate)))
    return alpha2_carbonate


def ANC_closed(pH, Total_Carbonates):
    #Acid neutralizing capacity (ANC) calculated under a closed system where
    #there are no carbonates exchanged with the atmosphere during the
    #experiment. Based on pH and total carbonates in the system.

    return (Total_Carbonates * (u.eq/u.mol * alpha1_carbonate(pH) +
            2 * u.eq/u.mol * alpha2_carbonate(pH)) +
            1 * u.eq/u.mol * Kw/invpH(pH) - 1 * u.eq/u.mol * invpH(pH))

def ANC_open(pH):
    #Acid neutralizing capacity (ANC) calculated under an open system, based
    #on pH.
    return ANC_closed(pH, P_CO2*K_Henry_CO2/alpha0_carbonate(pH))


data_file_path = "/Users/Emily/Desktop/4530/Acid Rain/lab2 acid rain.tsv"
#Now we create a pandas dataframe with the data in the file


time = epa.column_of_data(data_file_path,1,1,-1,"")
pH = epa.column_of_data(data_file_path,1,0,-1,"")


conservativeANC = [0]
ANC_input = -0.001 *u.eq/u.L
ANC_o = 50 * (u.eq/u.L)*10**(-6)
conservativeANC = (ANC_input)*(1-(2.7182**(-time)))+ (ANC_o)*(2.7182**(-time))


closedANC = [0]
ANC_input = -0.001 *u.eq/u.L
ANC_out = 50 * (u.eq/u.L)*10**(-6)
TC= 1.854*(u.eq/u.L)*10**(-3)
closedANC = ANC_closed(pH, TC)


openANC = np.array([0])
ANC_input = -0.001 *u.eq/u.L
ANC_out = 50 * (u.eq/u.L)*10**(-6)
openANC= ANC_open(pH)


residencetime = 609.137

t_over_theta = np.array([0*60, 5*60, 10*60, 15*60, 20*60])/residencetime
ANC0 = ANC0/(u.mol/u.L)*(u.eq/u.L)
data5 = "/Users/Emily/github/el545-4530/ANC Titrations/lab 3 titration time 5.tsv"
V_titrant, pH, V_Sample, Normality_Titrant, V_equivalent, ANC = epa.Gran(data5)
ANC5 = ANC/(u.mol/u.L)*(u.eq/u.L)
data10 = '/Users/Emily/github/el545-4530/ANC Titrations/lab 3 time 10.tsv'
V_titrant, pH, V_Sample, Normality_Titrant, V_equivalent, ANC = epa.Gran(data10)
ANC10 = ANC/(u.mol/u.L)*(u.eq/u.L)
data15 = '/Users/Emily/github/el545-4530/ANC Titrations/lab3 t 15 first try.tsv'
V_titrant, pH, V_Sample, Normality_Titrant, V_equivalent, ANC = epa.Gran(data15)
ANC15 = ANC/(u.mol/u.L)*(u.eq/u.L)
data20 = '/Users/Emily/github/el545-4530/ANC Titrations/lab 3 t 20.tsv'
V_titrant, pH, V_Sample, Normality_Titrant, V_equivalent, ANC = epa.Gran(data20)
ANC20 = ANC/(u.mol/u.L)*(u.eq/u.L)
anc_data = [ANC0, ANC5, ANC10, ANC15, ANC20]
anc_data = (np.stack(anc_data))*(u.eq/u.L)
anc_data


ancCOMP = [anc_data, t_over_theta]
plt.plot(time, conservativeANC, 'r')
plt.scatter(t_over_theta, anc_data, c='g', marker = 'o')
plt.plot(time, closedANC, 'b')
plt.plot(time, openANC, 'y')
plt.xlabel('t/theta')
plt.ylabel('ANC')
plt.legend(['conservative ANC','closed ANC','open ANC', 'Data'])
plt.grid(True)
plt.savefig('/Users/Emily/github/el545-4530/ANC Titrations/part3')
plt.show()
