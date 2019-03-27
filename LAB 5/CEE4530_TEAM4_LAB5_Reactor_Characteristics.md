##Lab 5: Reactor Characteristics
######Emily Liu el545
######Haonan Kuang hk876
###Introduction
Reactors are defined by boundaries --real or imaginary-- that physically confine the processes within. Chlorine contactor tanks are designed to maximize the contact time between chlorine and pathogens before water is delivered to consumers. Hence, the goal is to maximize that time. It is important to know a reactor’s level of mixing and residence time. Mixing levels characterize reactor categories: completely mixed flow reactors (CMFR), plug flow reactors (PFR), and flow with dispersion reactors (FDR). PFRs are not attainable in real life, and so reactors are either CMFR or FDR; we use these models in our experiment.

Reactors can be studied by examining the response curve of the effluent concentration after adding a spike or pulse of a tracer in the influent or after a step function change in input concentration.  As a function of dimensionless time $t^{\star}$,$\frac{t}{\theta}$, these response curves are known as E curves and F curves for pulse and step inputs respectively. The E curve is the exit age distribution and is the actual output of a tracer from a reactor, while the F curve is the cumulative age distribution and is the cumulative fraction of tracer that has exited a reactor at time $t^{\star}$.  The curves are related as such:
$$F_{(t^{\star})} = \int^{t^{\star}}_{0}E_{(t^{\star})}dt^{\star}$$
In a CMFR, a mass balance on a conservative tracer gives the following differential equation:
$$V_{r}\frac{dC}{dt} = (C_{in} -C)Q$$
where Q is the volumetric flow rate and $V_{r}$ is the volume of the reactor. If a pulse of tracer is discharged directly into a reactor such that the initial concentration is  $C_{0} = \frac{C_{tr}V{r}}{V_{r}}$ and the input concentration is zero --$C_{in = 0}$-- the solution to the differential equation in terms of dimensionless time is:
$$E_{(t^{\star})} = \frac{C_{(t^{\star})}V_{r}}{C_{tr}V_{tr}} = e^{(-t^{\star})}$$

Analytical solutions that describe FDRs are more difficult to come across than CMFRs, and oftentimes a parameter describing the dispersion is fit to the data rather than predicted prior to the experiment. For arbitrary mixing levels, we can use a one-dimensional advection-dispersion equation under open boundary conditions or a CMFR in series under closed boundary conditions. Under open boundary conditions, dispersion is possible across a boundary, and some of the tracer introduced at the reactor inlet can be carried upstream. Under closed boundary conditions, the reactor has a diffusion or dispersion coefficient that differs from those of the entrance and exit.

The governing differential equation for a one-dimensional advection-dispersion equation under open boundary conditions is as follows:
$$\frac{\partial C}{\partial t} = -U\frac{\partial C}{\partial x} + D_{d}\frac{\partial^{2}C}{\partial x^{2}}$$
where C is the concentration of a conservative substance, U is the average fluid velocity in the x direction, $D_{d}$ is the longitudinal dispersion coefficient, and t is time.
The solution of the governing differential equation under complete mixing in the y-z plane and advective and dispersive transport only in the x direction for any x and t (after t=0) is:
$$C(x,t) = \frac{M}{A\sqrt{4\pi{}D_{d}t}}exp[\frac{{-x^{'}}^{2}}{4D_{d}t}]$$
where M is the mass of conservative material in the spike, $D_{d}$ is the axial dispersion coefficient ([$\frac{L^{2}}{T}$]), $x^{'} = x - Ut$ (U is the longitudinal advective velocity in the reactor), and A is the cross-sectional area of the reactor. In the above equation, we can measure dispersion, and obtain a maximum value of C at t=x/U, and at this time $C(x,t) = \frac{M}{A\sqrt{4\pi{}D_{d}t}}$.
We can normalize the above equation and obtain its dimensionless form:
$$E_{(t^{\star})}=\sqrt{\frac{Pe}{4\pi{}t^{\star}}}exp[\frac{-(1-t^{\star})^{2}Pe}{4t^{\star}}]$$
The Peclet number Pe is the ratio of advective to dispersive transport. If M and A are known, $D_{d}$ can be estimated. $D_{d}$ can be normalized by dividing into a velocity and a length, such that:
$$Pe = \frac{UL}{D_{d}}$$
where L is the length of the reactor and U is the mean velocity. We can estimate the Peclet number by utilizing the variance of tracer concentration, $\sigma{}_{t}^{2}$. For open systems, this relationship is given by
$$\sigma{}_{t}^{2} = (\frac{2}{Pe}+\frac{8}{Pe^{2}})\theta^{2}$$
For closed systems, the relationship is:
$$\sigma{}_{t}^{2} = (\frac{2}{Pe}-\frac{2}{Pe^{2}}(1-e^{-Pe}))\theta^{2}$$

With CMFRs in series under closed boundary conditions, the concentration of the $N^{th}$ reactor can be described non-dimensionally as follows:
$$E_{N(t^{\star})} = \frac{N^{N}}{(N-1)!}(t^{\star})^{N-1}e^{(t^{\star})}$$

###Objective
In this experiment, our goal is to create a reactor design that will maximize the time it takes for water to travel from the tank influent to the effluent. We use tracer studies --in this case the conservative tracer is red dye-- to accurately estimate the effective contact time.
The data obtained from each reactor design will be analyzed and fit to a one dimensional advection-dispersion model under open boundary conditions and a model of CMFR reactors in series under closed boundary conditions.

###Procedures
After setup the reactor, make sure the connections of the tubes are tight and the flow goes in the right direction. Configure ProCoDA and set it to monitor both the photometer voltage and the photometer concentration. Calibrate one of them using the photometer calibration and leave the other uncalibrated as volts, which could help to do a post experiment recalibration of the photometer. Accurately measure the reactor volume and calculate the volume of red dye # 40 to get a maximum concentration of approximately 30 mg/L near the influent of the reactor. Start to log data exactly when adding the red dye and record the red dye concentration near the effluent as a function of time. Don’t stop log data until the majority of the red dye has exited. For the baffle design, install various configurations of perforated or staggered baffles. These baffles should be installed so that each compartment has the same volume. In order to get more accurate results, seal the gap between baffles and reactor walls with duct tape. This setup would look as follows:
![linear](https://raw.githubusercontent.com/el545/el545-4530/master/image1.jpeg)
 Afterwards, repeat the same procedures of the CMFR to collect data.

###Results and Discussion
n our experiment, we collected data for CMFR, two baffles, three baffles, and four baffles. For the first CMFR trial, when we compare the measured dye data and the CMFR model on the same graph, they have very different trends, as shown in figure 1. However, this bias could be the result that the residence time in the CMFR model does not alter with the measurement.  On the other hand, after analyzing the data through multivariable nonlinear regression for each baffle trial, we obtain the best fit between the experimental data and the two models (AD and CMFR) where the sum of the squared errors has been minimized. From these graphs (figure 2, 3, and 5), the measured dye concentration is very close to both models for each one. Peclet numbers and tracer residence times are calculated respectively for each baffle design. Our Peclect numbers are 7.877, 11.854, and 9.535, respectively. In figure 3, there is a biased point. This point occurs because when we were doing the trial, we accidentally taped the photometer and let air in. After we deleted the biased point, we plot the graph again in figure 4 and the Peclet number is 12.139. In order to get a maximum $t^{\star}$ that close to 1, which means a perfect PFR, we extract the time when F equals or closest to 0.1 and find $t^{\star}$ by applying the residence time, which is 290 seconds, 330 seconds, and 260 seconds respectively. Therefore, we got the $t^{\star}$ to be 0.9138, 0.9394, and 0.9808 for trials 1, 2, and 3 respectively. They are close to 1, namely our baffling conditions are near perfect.



```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.core.utility as ut
import numpy as np
import matplotlib.pyplot as plt

##Finding the concentration for which F = 0.1
F = 0.1
cstar = F*30 #for some odd reason, the find_nearest function below won't work with units
#Note that the cstar value has units of mg/L

def find_nearest(array, value):
  ##This function was found on this page --https://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
  #It's the first answer
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def tstar(time_data, concentration_data, concentration,residence_time):
  ##Given a concentration, residence time, arrays of concentration and time,this function finds
  #the time that corresponds to that concentration
  for i in range(concentration_data.size):
    if concentration_data[i] == concentration:
      print('The concentration closest to 3.0 mg/L is', concentration)
      print("The value for tstar is ", time_data[i]/residence_time)
#The following file is from a CMFR
print('This is our CMFR data')
CMFR_path = 'https://raw.githubusercontent.com/hk876/CEE-4530/master/reactor/lab5_team4_cmfr.xls'

# find the row after the last note in the file. This assumes that the last note marks the beginning of the test.
epa.notes(CMFR_path)
CMFR_firstrow = epa.notes(CMFR_path).last_valid_index() + 1
CMFR_firstrow

#I eliminate the beginning of the data file because this is a CMFR and the
#first data was taken before the dye reached the sensor. Note that eliminating
#some data from the beginning of the data file will not change the analysis
#except in the estimate of the initial tracer mass.#
CMFR_firstrow = CMFR_firstrow + 10
CMFR_time_data = (epa.column_of_time(CMFR_path,CMFR_firstrow,-1)).to(u.s)

CMFR_concentration_data = epa.column_of_data(CMFR_path,CMFR_firstrow,1,-1,'mg/L')

#You should use real measured values!#
CMFR_V = 1.984*u.L
CMFR_Q = 266.9994 * u.mL/u.min

#here we set estimates that we will use as starting values for the curve fitting
CMFR_theta_hydraulic = (CMFR_V/CMFR_Q).to(u.s)
CMFR_C_bar_guess = np.max(CMFR_concentration_data)

#The Solver_CMFR_N will return the initial tracer concentration,
#residence time, and number of reactors in series.
#This experiment was for a single reactor and so we expect N to be 1!
CMFR_CMFR = epa.Solver_CMFR_N(CMFR_time_data, CMFR_concentration_data, CMFR_theta_hydraulic, CMFR_C_bar_guess)
#use dot notation to get the 3 elements of the tuple that are in CMFR.
print('The model estimated mass of tracer injected was',ut.round_sf(CMFR_CMFR.C_bar*CMFR_V ,2) )
print('The model estimate of the number of reactors in series was', CMFR_CMFR.N)
print('The tracer residence time was',ut.round_sf(CMFR_CMFR.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(CMFR_CMFR.theta/CMFR_theta_hydraulic).magnitude)

##Finding the concentration for which F = 0.1
concCMFR = find_nearest(CMFR_concentration_data, cstar)*(u.mg/u.L)
tstar(CMFR_time_data,CMFR_concentration_data, concCMFR, ut.round_sf(CMFR_CMFR.theta ,2))
#create a model curve given the curve fit parameters.

CMFR_CMFR_model = CMFR_CMFR.C_bar * epa.E_CMFR_N(CMFR_time_data/CMFR_CMFR.theta,CMFR_CMFR.N)
plt.plot(CMFR_time_data.to(u.min), CMFR_concentration_data.to(u.mg/u.L),'r')
plt.plot(CMFR_time_data.to(u.min), CMFR_CMFR_model,'b')

plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model'])
plt.savefig('/Users/Emily Liu/Desktop/CEE 4530/LAB 5/cmfr', bbox_inches = 'tight')
plt.show()

#Load a data file for a reactor with  two baffles.

two_baffle_path = 'https://raw.githubusercontent.com/hk876/CEE-4530/master/reactor/lab%205%20trial%201%20redo.xls'
two_baffle_firstrow = epa.notes(two_baffle_path).last_valid_index() + 1
two_baffle_time_data = (epa.column_of_time(two_baffle_path,two_baffle_firstrow,-1)).to(u.s)
two_baffle_concentration_data = epa.column_of_data(two_baffle_path,two_baffle_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomoly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#
print('Trial 1')
two_baffle_concentration_data = two_baffle_concentration_data - two_baffle_concentration_data[0]
two_baffle_V = 1.984*u.L
two_baffle_Q = 380 * u.mL/u.min
two_baffle_theta_hydraulic = (two_baffle_V/two_baffle_Q).to(u.s)
two_baffle_C_bar_guess = np.max(two_baffle_concentration_data)/2
#use solver to get the CMFR parameters
two_baffle_CMFR = epa.Solver_CMFR_N(two_baffle_time_data, two_baffle_concentration_data, two_baffle_theta_hydraulic, two_baffle_C_bar_guess)
two_baffle_CMFR.C_bar
two_baffle_CMFR.N
two_baffle_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
two_baffle_CMFR_model = (two_baffle_CMFR.C_bar*epa.E_CMFR_N(two_baffle_time_data/two_baffle_CMFR.theta, two_baffle_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
two_baffle_AD = epa.Solver_AD_Pe(two_baffle_time_data, two_baffle_concentration_data, two_baffle_theta_hydraulic, two_baffle_C_bar_guess)
two_baffle_AD.C_bar
two_baffle_AD.Pe
two_baffle_AD.theta

print('The model estimated mass of tracer injected was',ut.round_sf(two_baffle_AD.C_bar*two_baffle_V ,2) )
print('The model estimate of the Peclet number was', two_baffle_AD.Pe)
print('The tracer residence time was',ut.round_sf(two_baffle_AD.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(two_baffle_AD.theta/two_baffle_theta_hydraulic).magnitude)

#Create the advection dispersion model curve based on the solver parameters
two_baffle_AD_model = (two_baffle_AD.C_bar*epa.E_Advective_Dispersion((two_baffle_time_data/two_baffle_AD.theta).to_base_units(), two_baffle_AD.Pe)).to(u.mg/u.L)

##Finding tstar such that F = 0.1
concTRIAL1 = find_nearest(two_baffle_concentration_data, cstar)*(u.mg/u.L)
tstar(two_baffle_time_data,two_baffle_concentration_data, concTRIAL1,ut.round_sf(two_baffle_AD.theta ,2) )

#Plot the data and the two model curves.
##NOT SURE ABOUT THE POINTS -- THEY LOOK A LITTLE ODD SO WE REPLACED THEM WITH DASHED LINES##
plt.plot(two_baffle_time_data.to(u.s), two_baffle_concentration_data.to(u.mg/u.L),'r--')
plt.plot(two_baffle_time_data.to(u.s), two_baffle_CMFR_model,'b')
plt.plot(two_baffle_time_data.to(u.s), two_baffle_AD_model,'g')
plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])

plt.savefig('/Users/Emily Liu/Desktop/CEE 4530/LAB 5/trial1', bbox_inches = 'tight')
plt.show()


#Trial 2: 4 reactors of equal volume and 3 baffles#
print('Trial 2')
three_baffle_path = 'https://raw.githubusercontent.com/hk876/CEE-4530/master/reactor/lab%205%20%20trial%202.xls'
three_baffle_firstrow = epa.notes(three_baffle_path).last_valid_index() + 1
three_baffle_time_data = (epa.column_of_time(three_baffle_path,three_baffle_firstrow,-1)).to(u.s)
three_baffle_concentration_data = epa.column_of_data(three_baffle_path,three_baffle_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomoly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

three_baffle_concentration_data = three_baffle_concentration_data - three_baffle_concentration_data[0]
three_baffle_V = 2.03625*u.L
three_baffle_Q = 380 * u.mL/u.min
three_baffle_theta_hydraulic = (three_baffle_V/three_baffle_Q).to(u.s)
three_baffle_C_bar_guess = np.max(three_baffle_concentration_data)/2
#use solver to get the CMFR parameters
three_baffle_CMFR = epa.Solver_CMFR_N(three_baffle_time_data, three_baffle_concentration_data, three_baffle_theta_hydraulic, three_baffle_C_bar_guess)
three_baffle_CMFR.C_bar
three_baffle_CMFR.N
three_baffle_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
three_baffle_CMFR_model = (three_baffle_CMFR.C_bar*epa.E_CMFR_N(three_baffle_time_data/three_baffle_CMFR.theta, three_baffle_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
three_baffle_AD = epa.Solver_AD_Pe(three_baffle_time_data, three_baffle_concentration_data, three_baffle_theta_hydraulic, three_baffle_C_bar_guess)
three_baffle_AD.C_bar
three_baffle_AD.Pe
three_baffle_AD.theta

print('The model estimated mass of tracer injected was',ut.round_sf(three_baffle_AD.C_bar*three_baffle_V ,2) )
print('The model estimate of the Peclet number was', three_baffle_AD.Pe)
print('The tracer residence time was',ut.round_sf(three_baffle_AD.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(three_baffle_AD.theta/three_baffle_theta_hydraulic).magnitude)

#Create the advection dispersion model curve based on the solver parameters
three_baffle_AD_model = (three_baffle_AD.C_bar*epa.E_Advective_Dispersion((three_baffle_time_data/three_baffle_AD.theta).to_base_units(), three_baffle_AD.Pe)).to(u.mg/u.L)

#Finding tstar such that F = 0.1
conctrial2 = find_nearest(three_baffle_concentration_data, cstar)*(u.mg/u.L)
tstar(three_baffle_time_data,three_baffle_concentration_data, conctrial2,ut.round_sf(three_baffle_AD.theta ,2))



#Plot the data and the two model curves.
##NOT SURE ABOUT THE POINTS -- THEY LOOK A LITTLE ODD SO WE REPLACED THEM WITH DASHED LINES##
plt.plot(three_baffle_time_data.to(u.s), three_baffle_concentration_data.to(u.mg/u.L),'r--')
plt.plot(three_baffle_time_data.to(u.s), three_baffle_CMFR_model,'b')
plt.plot(three_baffle_time_data.to(u.s), three_baffle_AD_model,'g')
plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])

plt.savefig('/Users/Emily Liu/Desktop/CEE 4530/LAB 5/trial2', bbox_inches = 'tight')
plt.show()


##Trial 2 with the outlying data removed##
print('Trial 2 sans outliers')
threer_baffle_path = 'https://raw.githubusercontent.com/hk876/CEE-4530/master/reactor/lab%205%20%20trial%202%20delete.xls'
threer_baffle_firstrow = epa.notes(threer_baffle_path).last_valid_index() + 1
threer_baffle_time_data = (epa.column_of_time(threer_baffle_path,threer_baffle_firstrow,-1)).to(u.s)
threer_baffle_concentration_data = epa.column_of_data(threer_baffle_path,threer_baffle_firstrow,1,-1,'mg/L')


threer_baffle_concentration_data = threer_baffle_concentration_data - threer_baffle_concentration_data[0]
threer_baffle_V = 2.03625*u.L
threer_baffle_Q = 380 * u.mL/u.min
threer_baffle_theta_hydraulic = (threer_baffle_V/threer_baffle_Q).to(u.s)
threer_baffle_C_bar_guess = np.max(threer_baffle_concentration_data)/2
#use solver to get the CMFR parameters
threer_baffle_CMFR = epa.Solver_CMFR_N(threer_baffle_time_data, threer_baffle_concentration_data, threer_baffle_theta_hydraulic, threer_baffle_C_bar_guess)
threer_baffle_CMFR.C_bar
threer_baffle_CMFR.N
threer_baffle_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
threer_baffle_CMFR_model = (threer_baffle_CMFR.C_bar*epa.E_CMFR_N(threer_baffle_time_data/threer_baffle_CMFR.theta, threer_baffle_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
threer_baffle_AD = epa.Solver_AD_Pe(threer_baffle_time_data, threer_baffle_concentration_data, threer_baffle_theta_hydraulic, threer_baffle_C_bar_guess)
threer_baffle_AD.C_bar
threer_baffle_AD.Pe
threer_baffle_AD.theta

print('The model estimated mass of tracer injected was',ut.round_sf(threer_baffle_AD.C_bar*threer_baffle_V ,2) )
print('The model estimate of the Peclet number was', threer_baffle_AD.Pe)
print('The tracer residence time was',ut.round_sf(threer_baffle_AD.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(threer_baffle_AD.theta/threer_baffle_theta_hydraulic).magnitude)

#Create the advection dispersion model curve based on the solver parameters
threer_baffle_AD_model = (threer_baffle_AD.C_bar*epa.E_Advective_Dispersion((threer_baffle_time_data/threer_baffle_AD.theta).to_base_units(), threer_baffle_AD.Pe)).to(u.mg/u.L)


#Finding tstar such that F = 0.1
conctrial2r = find_nearest(threer_baffle_concentration_data, cstar)*(u.mg/u.L)
tstar(threer_baffle_time_data,threer_baffle_concentration_data, conctrial2r,ut.round_sf(threer_baffle_AD.theta ,2))



#Plot the data and the two model curves.
##NOT SURE ABOUT THE POINTS -- THEY LOOK A LITTLE ODD SO WE REPLACED THEM WITH DASHED LINES##
plt.plot(threer_baffle_time_data.to(u.s), threer_baffle_concentration_data.to(u.mg/u.L),'r--')
plt.plot(threer_baffle_time_data.to(u.s), threer_baffle_CMFR_model,'b')
plt.plot(threer_baffle_time_data.to(u.s), threer_baffle_AD_model,'g')
plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])

plt.savefig('/Users/Emily Liu/Desktop/CEE 4530/LAB 5/trial2w_no_outlier', bbox_inches = 'tight')
plt.show()



###Trial 3: Four baffles
print('Trial 3')
four_baffle_path = 'https://raw.githubusercontent.com/hk876/CEE-4530/master/reactor/lab%205%20trial%203.xls'
four_baffle_firstrow = epa.notes(four_baffle_path).last_valid_index() + 1
four_baffle_time_data = (epa.column_of_time(four_baffle_path,four_baffle_firstrow,-1)).to(u.s)
four_baffle_concentration_data = epa.column_of_data(four_baffle_path,four_baffle_firstrow,1,-1,'mg/L')



four_baffle_concentration_data = four_baffle_concentration_data - four_baffle_concentration_data[0]
four_baffle_V = 2.025*u.L
four_baffle_Q = 380 * u.mL/u.min
four_baffle_theta_hydraulic = (four_baffle_V/four_baffle_Q).to(u.s)
four_baffle_C_bar_guess = np.max(four_baffle_concentration_data)/2
#use solver to get the CMFR parameters
four_baffle_CMFR = epa.Solver_CMFR_N(four_baffle_time_data, four_baffle_concentration_data, four_baffle_theta_hydraulic, four_baffle_C_bar_guess)
four_baffle_CMFR.C_bar
four_baffle_CMFR.N
four_baffle_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
four_baffle_CMFR_model = (four_baffle_CMFR.C_bar*epa.E_CMFR_N(four_baffle_time_data/four_baffle_CMFR.theta, four_baffle_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
four_baffle_AD = epa.Solver_AD_Pe(four_baffle_time_data, four_baffle_concentration_data, four_baffle_theta_hydraulic, four_baffle_C_bar_guess)
four_baffle_AD.C_bar
four_baffle_AD.Pe
four_baffle_AD.theta

print('The model estimated mass of tracer injected was',ut.round_sf(four_baffle_AD.C_bar*four_baffle_V ,2) )
print('The model estimate of the Peclet number was', four_baffle_AD.Pe)
print('The tracer residence time was',ut.round_sf(four_baffle_AD.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(four_baffle_AD.theta/four_baffle_theta_hydraulic).magnitude)

#Create the advection dispersion model curve based on the solver parameters
four_baffle_AD_model = (four_baffle_AD.C_bar*epa.E_Advective_Dispersion((four_baffle_time_data/four_baffle_AD.theta).to_base_units(), four_baffle_AD.Pe)).to(u.mg/u.L)


##finding tstar such that F=0.1
conctrial3 = find_nearest(four_baffle_concentration_data, cstar)*(u.mg/u.L)
tstar(four_baffle_time_data,four_baffle_concentration_data, conctrial3,ut.round_sf(four_baffle_AD.theta ,2))

#Plot the data and the two model curves.
##NOT SURE ABOUT THE POINTS -- THEY LOOK A LITTLE ODD SO WE REPLACED THEM WITH DASHED LINES##
plt.plot(four_baffle_time_data.to(u.s), four_baffle_concentration_data.to(u.mg/u.L),'r--')
plt.plot(four_baffle_time_data.to(u.s), four_baffle_CMFR_model,'b')
plt.plot(four_baffle_time_data.to(u.s), four_baffle_AD_model,'g')
plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])

plt.savefig('/Users/Emily Liu/Desktop/CEE 4530/LAB 5/trial3', bbox_inches = 'tight')
plt.show()
