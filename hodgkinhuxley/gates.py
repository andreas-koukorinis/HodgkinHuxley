from numpy import exp

class Gate:
    # Parent class for channel gates, provides routines for calculating
    # steady state values and time constants, and evaluating differential
    # equations
    def inf(self,v):
        a = self.alpha(v)
        b = self.beta(v)
        return a/(a+b)

    def tau(self,v):
        a = self.alpha(v)
        b= self.beta(v)
        return 1/(a+b)

    def dif_eq(self,d,v):
        a = self.alpha(v)
        b = self.beta(v)
        return a * (1-d) - b*d

class MGate(Gate):
    # Defines alpha and beta functions for the activating gates on the sodium channel
    def alpha(self,v):
    	if v==25:
    		return 1.
    	else:
        	return  0.1 * (25 -v)/(exp((25-v)/10)-1)
    
    def beta(self,v):
        return 4 * exp(-v/18)

class HGate(Gate):
    # Defines alpha and beta functions for the inactivating gate on the sodium channel
    def alpha(self,v):
        return 0.07*exp(-v/20)

    def beta(self,v):
        return 1/(exp((30-v)/10)+1)

class NGate(Gate):
    # Defines alpha and beta functions for the activating gate on the potassium channel
    def alpha(self,v):
    	if v==10:
    		return 0.1
    	else:
        	return 0.01 * (10-v) /(exp((10-v)/10)-1)

    def beta(self,v):
        return 0.125*exp(-v/80)    
