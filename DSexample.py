from HBMBEMlib import DirichletSignorini as DS
from HBMBEMlib import HBM_method as HBM
import matplotlib.pyplot as plt
import scipy.optimize as sop
import numpy as np

## Define system of interest
bar = DS.DSBarSystem()
bar.g0 = 0.5;

## Define DSinuation parameters and initialize result list
Ta = np.arange(3.8,3.01,-0.001)
result = [];
m = 40
last_admissible_soln = np.zeros(m+1); 
last_admissible_soln[1] = 1.5;

## Sequential continuation loop

for i in range(len(Ta)):
    
    newpoint = DS.HbmBemContinuationPoint(Ta[i],barsys = bar,rho =40,harmonics = m) 
    cosmatrix = newpoint.calcCosmatrix()
    nonlinearfunc = lambda disp: HBM.HBM_u(newpoint,disp,cosmatrix)
    
    result_temp = sop.root(nonlinearfunc,last_admissible_soln,method='hybr')
    newpoint.displacementL = result_temp.x
    newpoint.forceL = newpoint.calcForce(newpoint.displacementL)

    HBM.calcEnergyResidualDS(newpoint)
    result.append(newpoint)
    if newpoint.energy>1e-5 and newpoint.res < 1:
        last_admissible_soln = newpoint.displacementL

## Post processing
energy = []
res = []
for item in result:
    energy.append(item.energy)
    res.append(item.res)
    
energy = np.array(energy)
res = np.array(res)

adm = res<1e-2 #admissible solution set
result_adm = np.array(result)[adm]
energy_adm = energy[adm]
Ta_adm = np.array(Ta)[adm]
res_adm = res[adm]

## Plot
# energy plot
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.scatter(2*np.pi/Ta_adm,energy_adm)
ax.set_yscale('log')
plt.ylim([1e-1,1e3])

