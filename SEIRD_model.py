#!/usr/bin/env python
# coding: utf-8



from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# # Let’s run the basic SIR model

# describe the model
def deriv(y, t, N, beta, k, delta, alpha, rho):
    S, E, I, R, D = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - delta * E
    dIdt = delta * E - (1 - alpha ) * k * I - alpha * rho* I
    dRdt = (1 - alpha) * k * I
    dDdt = alpha * rho * I #death
    return dSdt, dEdt, dIdt, dRdt, dDdt

#death % per age group in sweden (dead/confirmed sick)
alpha_by_agegroup = {
    "0-29": 0.00001, "30-59": 0.0001, "60-89": 0.08, "89+": 0.305
}

#swedish age distribution
proportion_of_agegroup = {
    "0-29": 0.36, "30-59": 0.385, "60-89": 0.245, "89+": 0.01
}


# describe the parameters
N =  1000000        # population 2470
beta = 2    #expected amount of people an infected person infects per day
k = 1/7
delta = 1/10 #incubation period 
rho = 1/9 #days until death
alpha = sum(alpha_by_agegroup[i] * proportion_of_agegroup[i] 
            for i in list(alpha_by_agegroup.keys()))
S0, E0, I0, R0, D0 = N-1, 1, 0, 0, 0  # initial conditions: one infected, rest susceptible



t = np.linspace(0, 99, 100) # Grid of time points (in days)
y0 = S0, E0, I0, R0, D0 # Initial conditions vector

# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, k, delta, alpha, rho))
S, E, I, R, D = ret.T




def plotsir(t, S, E, I, R, D):
  f, ax = plt.subplots(1,1,figsize=(10,4))
  #ax.plot(t, S, 'b', alpha=0.7, linewidth=2, label='Susceptible')
  #ax.plot(t, E, 'y', alpha=0.7, linewidth=2, label='Exposed')
  #ax.plot(t, I, 'r', alpha=0.7, linewidth=2, label='Infected')
  #ax.plot(t, R, 'g', alpha=0.7, linewidth=2, label='Recovered')
  ax.plot(t, D, 'k', alpha=0.7, linewidth=2, label='Dead')
  #ax.plot(t, S+E+I+R+D, 'k--', alpha=0.7, linewidth=2, label='Total')
  
  ax.set_xlabel('Time (days)')

  ax.yaxis.set_tick_params(length=0)
  ax.xaxis.set_tick_params(length=0)
  ax.grid(b=True, which='major', c='w', lw=2, ls='-')
  legend = ax.legend()
  legend.get_frame().set_alpha(0.5)
  for spine in ('top', 'right', 'bottom', 'left'):
      ax.spines[spine].set_visible(False)
  #plt.savefig("Plot.png")
  plt.show();


# plot the graph



plotsir(t, S, E, I, R, D)





