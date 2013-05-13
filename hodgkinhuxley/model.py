import matplotlib.pylab as plt
import numpy as np
from scipy.integrate import odeint

from gates import MGate, HGate, NGate

class Model:

    def __init__(self):

        #Gate functions
        self.mgate = MGate()
        self.hgate = HGate()
        self.ngate = NGate()
        
        #Channel conductances
        self.g_n = 120
        self.g_k = 36
        self.g_l = 0.3
        
        #Nernst Potentials
        self.v_n = 115
        self.v_k = -12
        self.v_l = 10.6 
        
        #Membrane Capacitance 
        self.Cm = 1
    
    
    
      
    def plot_steady_state(self):
       	#Plot steady state gate functions over a range of membrane potentials
    	m_inf = np.vectorize(self.mgate.inf)
    	h_inf = np.vectorize(self.hgate.inf)
    	n_inf = np.vectorize(self.ngate.inf)
    	
    	v = np.linspace(-60, 120, 181)
    	
    	#Plot Graph
    	plt.figure()
    	plt.plot(v, m_inf(v), label=r'$m_{\infty} (v)$')
    	plt.plot(v, h_inf(v), label=r'$h_{\infty} (v)$')
    	plt.plot(v, n_inf(v), label=r'$n_{\infty} (v)$')
    	plt.xlabel('Membrane Potential (mV)')
    	plt.title('Steady State Gate Functions')
    	plt.legend(loc=0)
    	plt.show()
    	
    def plot_time_constants(self):
        #Plot time constants over a range of membrane potentials	
    	m_tau = np.vectorize(self.mgate.tau)
    	h_tau = np.vectorize(self.hgate.tau)
    	n_tau = np.vectorize(self.ngate.tau)
    	
    	v = np.linspace(-60, 120, 181)
    	
    	plt.figure()
    	plt.plot(v, m_tau(v), label=r'$\tau_m (v)$')
    	plt.plot(v, h_tau(v), label=r'$\tau_h (v)$')
    	plt.plot(v, n_tau(v), label=r'$\tau_v (v)$')
    	plt.xlabel('Membrane Potential (mV)')
    	plt.title('Gate Time Constants')
    	plt.legend(loc=0)
    	plt.show()
    	   	
    def f(self, y, t):
    
        v = y[0]
        m = y[1]
        h = y[2]
        n = y[3]
        
        # Channel Currents
        i_k = -self.g_k*(n**4)*(v-self.v_k)
        i_n = -self.g_n*(m**3)*h*(v-self.v_n)
        i_l = -self.g_l*(v-self.v_l)
        
        # Applied Current        
        i_app = self.I[min((t/self.dt),((self.T[-1]/self.dt)-1))]        
        
        # System of Differential Equations
        dv =  (i_k+i_n+i_l+i_app)/self.Cm
        dm = self.mgate.dif_eq(m,v)
        dh = self.hgate.dif_eq(h,v)
        dn = self.ngate.dif_eq(n,v)
        
        return [dv,dm,dh,dn]
            	
    def simulate(self, I, dt):
    
        self.dt = dt
        self.I = I
        
        # Compute time grid
        self.T = np.arange(0,(len(self.I)+1)*self.dt,self.dt)
         
        # Initial Conditions
        v0 = 0.0
        m0 = self.mgate.inf(v0)
        h0 = self.hgate.inf(v0)
        n0 = self.ngate.inf(v0)
        y0 = [v0, m0, h0, n0]
        
        # Solve Equations
        soln = odeint(self.f, y0, self.T)
        v = soln[:,0]
        m = soln[:,1]
        h = soln[:,2]
        n = soln[:,3]
        
        # Plot graph of action potential
        plt.figure()
        plt.plot(self.T,v)
        plt.xlabel('Time (mS)')
        plt.ylabel('Membrane Potential (mV)')
        plt.title('Action Potential Plot')
        
        # Plot graph of gate variables
        plt.figure()
        plt.plot(self.T,m,label='m')
        plt.plot(self.T,h,label='h')
        plt.plot(self.T,m,label='n')
        plt.xlabel('Time (mS)')
        plt.title('Gate variables')
        plt.legend(loc=0)
        
        plt.show()
    
    
    
	