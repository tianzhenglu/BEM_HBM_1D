import numpy as np
import matplotlib.pyplot as plt

class DSBarSystem():
    # Class of system of interest: Dirichlet-Signorini bar
    def __init__(self,L=1.,density=1.,E=1.,g0 = 1):
        self.L = L
        self.density = density
        self.E = E
        self.c = np.sqrt(self.E/self.density)
        self.g0 = g0

class HbmBemContinuationPoint():
    # This is a conitnuation point for HBMBEM, cos term only
    # This is continuation for DS system
    def __init__(self,T,barsys= DSBarSystem(),harmonics = 20,cd = 20,rho=4):
        # Continuation Properties
        self.T = T #Period
        self.Omega = 2*np.pi/T #Targeted frequency
        self.cd = cd #Number of Reimann sum, n_h in the paper
        self.rho = rho #Signorini nonsmooth condition coefficient, rho >0
        self.harmonics = harmonics #Number of truncation harmonics
        self.barsys = barsys #properties of the bar
        self.nt = self.cd*self.harmonics
        self.dt = self.T/self.nt
        self.nx = 100;
        self.dx = self.barsys.L/self.nx;
        # Continuation results
        self.res = None
        self.energy = None
        self.displacementL = np.zeros(1+harmonics) #cos coefficient at x = L
        self.displacement0 = np.zeros(1+harmonics)
        self.forceL = np.zeros(1+harmonics) 
        self.force0 = np.zeros(1+harmonics) 
        self.upConversion = np.zeros(1+harmonics)
    
    def calcUpConversion(self):
        # This method calculates relationship between u and p
        # Method is based on Nonlinear Dynamics paper, equation 25
        # displacementL = upConversion * forceL
        kn = np.arange(self.harmonics+1)*self.Omega/self.barsys.c
        
        self.upConversion[1:] =  np.tan(kn[1:]*self.barsys.L)/kn[1:]
        self.upConversion[0] = self.barsys.L            
    
    def calcDisplacement(self,force):
        # Calculate displacementL amplitudes from forceL
        if not self.upConversion[0]:
            self.calcUpConversion()
        return self.upConversion * force
        
    def calcForce(self,displacement):
        # Calculate forceL amplitudes from displacementL
        if not self.upConversion[0]:
            self.calcUpConversion()
        return displacement / self.upConversion
        
    def calc0(self):
        # Calculate values at x = 0
        kn = np.arange(self.harmonics+1)*self.Omega/self.barsys.c
        self.force0[1:] = -self.displacementL[1:]/(np.sin(kn[1:]*self.barsys.L)/kn[1:])
        self.force0[0] = -self.forceL[0]
        
    def calcCosmatrix(self):
        # Calculate inverse fourier matrix. It converges amplitudes of harmonics
        # to time domain signal
        ta = np.arange(self.nt)*self.dt
        cosmatrix = np.zeros((self.harmonics*self.cd,self.harmonics+1))
    
        for i in range(self.harmonics+1):
            cosmatrix[:,i] = np.cos(ta*i*self.Omega)
        return cosmatrix
    

           
        
if __name__ == "__main__":
    pass
#    a = HbmBemContinuationPoint(T,barsys = bar)
#    a.displacementL = displacement
#    a.forceL = a.calcForce(a.displacementL)
#    cosmatrix = a.calcCosmatrix()
#    plt.plot(a.calcResidual(cosmatrix,a.displacementL,a.forceL))
