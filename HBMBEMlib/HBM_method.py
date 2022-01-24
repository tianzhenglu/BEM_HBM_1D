import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  
  
def calcResidual(point,cosmatrix,displacement,force):
    # calculate residual
    #   Warning:The input displacement and force are arrays. Output are SQUEEZED to array.
    #           This awkward implementation is introduced by scipy.root, which automatically
    #           squeeze dimension.
    displacementInTime = cosmatrix@(displacement.reshape(-1,1)) 
    forceInTime = cosmatrix@(force.reshape(-1,1))
    maxpart = point.rho*(displacementInTime - point.barsys.g0)- forceInTime
    maxresult = np.maximum(maxpart,np.zeros_like(maxpart))
    signorini = forceInTime + maxresult
    return (displacementInTime.squeeze(),forceInTime.squeeze(), signorini.squeeze())
    

def HBM_u(point,displacement,cosmatrix):
    # Calculate HBM of signorini residual
    # Use displacement harmonis as input
    force = point.calcForce(displacement)
    res_in_time = calcResidual(point,cosmatrix,displacement,force)[2]
    res_HBM = cosmatrix.transpose() @ res_in_time / point.nt
    return res_HBM

def HBM_p(point,force,cosmatrix):
    # Calculate HBM of signorini residual
    # Use force harmonis as input
    displacement = point.calcDisplacement(force)
    res_in_time = calcResidual(point,cosmatrix,displacement,force)[2]
    res_HBM = cosmatrix.transpose() @ res_in_time / point.nt
    return res_HBM

def calcDisplacementField(point):
    # Calculate displacemetn field during the priod
    # Calc force/displacement at 0 if not calculated
    if (not point.force0[0]) and (not point.displacement0[0]):
        point.calc0()
    xa = np.arange(point.nx + 1)*point.dx #collocation points in space
    ka = np.arange(point.harmonics+1)*point.Omega/point.barsys.c
    ux = np.zeros((point.nx+1,point.nt))
    L = point.barsys.L
    cosmatrix = point.calcCosmatrix()
    for i in range(point.nx+1):
        x = xa[i]
        ux_amp = np.zeros(point.harmonics+1)
        ux_amp[1:] = -np.sin(ka[1:]*x)*point.force0[1:]/(2*ka[1:])\
            -np.sin(ka[1:]*(L-x))*point.forceL[1:]/(2*ka[1:])\
            + np.cos(ka[1:]*(L-x))*point.displacementL[1:]/2\
            + np.cos(ka[1:]*x)*point.displacement0[1:]/2
        ux_amp[0] = x/L*(point.displacementL[0]-point.displacement0[0])+point.displacement0[0]
        ux[[i],:] =  ux_amp@cosmatrix.transpose()
    return ux
    
def calcEnergyResidualDS(point):
    # calculate energy and reisudal for DS systems
    # ONLY USE AFTER calculated all displacement/force attributes
    ux = calcDisplacementField(point)
    point.energy = np.asscalar((ux[1:,0]-ux[:-1,0])@(ux[1:,[0]]-ux[:-1,[0]])/2./point.dx)
    cosmatrix = point.calcCosmatrix()
    signorini =  calcResidual(point,cosmatrix,point.displacementL,point.forceL)[2]
    point.res = np.asscalar(sum(signorini*signorini)/point.nt)

def plot2D(point):
    # Plot displacmentL, forceL, signoriniL versus time
    # Only use after calc0()
    plt.figure(1)
    ta = np.arange(point.nt)*point.dt
    cosmatrix = point.calcCosmatrix()
    displacementInTime,forceInTime, signorini =  calcResidual(point,cosmatrix,point.displacementL,point.forceL)
    plt.plot(ta,displacementInTime)
    plt.plot(ta,forceInTime)
    plt.plot(ta,signorini)
    
def plot3D(point):
    # Plot displacment field versus space and time
    # Only use after calc0()
    ux=calcDisplacementField(point)
    fig = plt.figure(2)
    ax = fig.add_subplot(111, projection='3d')
    X,Y = np.meshgrid(np.arange(point.nt)*point.dt,np.arange(point.nx+1)*point.dx)
    ax.plot_surface(X,Y,ux)
    
    
if __name__ == "__main__":
    pass
    