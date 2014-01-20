# -*- coding: utf-8 -*-
"""
Centrifugal Pump design tool with EA/LOX engine
Only use for Sigle-stage certrifugal pump.
It is limited the use of approximate culcuration.
This script is New BSD License.
cf. Basic Design on Pumps by Hirohisa Takeda
"""
__author__ = "Takahiro Inagawa <ina111meister_at_gmail.com>"
__status__ = "beta"
__version__ = "0.0.1"
__date__    = "20 January 2014"

import numpy as np
import scipy as sp
from scipy.optimize import fsolve
from scipy.interpolate import UnivariateSpline
#from scipy.integrate import odeint
import matplotlib.pyplot as plt
#from matplotlib.font_manager import FontProperties
#fp = FontProperties(fname=r'C:\WINDOWS\Fonts\TakaoPGothic.ttf') # HG ゴシック系

class Rocket(object):
    """
    Hydrocarbon Rocket Engine Class
    It output mass flow rate when "thrust, Isp and O/F" are determined. 
    """
    def __init__(self, thrust, Isp, OFratio):
        self.g = 9.8066 #Gravitational acceleration
        self.thrust = thrust
        self.Isp = Isp
        self.OFratio = OFratio
        
        # It must change when propellant is changed.
        self.density_LOX = 1140.0 # kg/m3
        self.density_EA = 789.0 # kg/m3
        self.density_water = 1000.0 # kg/m3
        
    def culc(self):        
        self.mdot_proppelant = self.thrust / self.Isp / self.g
        self.mdot_oxidant = self.mdot_proppelant * (self.OFratio / (self.OFratio + 1))
        self.flow_oxidant = self.mdot_oxidant / self.density_LOX
        self.mdot_fuel = self.mdot_proppelant * (1 / (self.OFratio + 1))
        self.flow_fuel = self.mdot_fuel / self.density_EA
        
    def show(self):
        print "==== Rocket Engine ===="
        print "Thrust:\t\t%d (N)" % (self.thrust)
        print "Isp:\t\t%d (sec)" % (self.Isp)
        print "O/F:\t\t%.2f " % (self.OFratio)
        print "Proppelant mass flow:\t%.2f (kg/s)" % (self.mdot_proppelant)
        print "Oxidant mass flow:\t%.2f (kg/s)" % (self.mdot_oxidant)
        print "Oxidant flow Q:\t%.4f (m3/s)" % (self.flow_oxidant)
        print "Fuel mass flow:\t%.2f (kg/s)" % (self.mdot_fuel)
        print "Fuel flow Q:\t%.4f (m3/s)" % (self.flow_fuel)
        
class Pump(object):
    """
    Centrifugal Pump Class
    """
    def __init__(self, rocket, working_fluid = "fuel"):
        """
        @arg rocket: rocket class
        @arg working_fluid = "oxidant", "fuel" or "water"
        """
        self.g = rocket.g
        if working_fluid == "oxidant":
            self.density = rocket.density_LOX
        elif working_fluid == "fuel":
            self.density = rocket.density_EA
        elif working_fluid == "water":
            self.density = rocket.density_water
        self.rpm = 16000
        self.press_out = 3.0 #MPa
        self.flow_q = 0.04 #m3/s
        
        self.head_total = self.press_out / self.density / rocket.g * (10 ** 6)
        
        self.press_in = 0.1 #MPa
        self.head_loss = 3 #m
        self.head_v = 1.046862 #m
        self.height_pump = 0 #m
        """ Impeller size """
        self.D2 = 0.1
        """ efficiency, design constant """
        self.eta_h = 1.0
        self.eta_v = 0.9
        self.eta_m = 0.9
        """ Impeller inlet """
        self.rambda1 = 0.25
        self.rambda2 = 1.1
        self.hubratio = 0.5
        """ Impeler outlet """        
        self.Z = 5
        self.beta2 = 25
        self.R1R2ratio = 0.3
        
    def culc_Ns(self):
        """
        Culcurate Specific speed Ns and type number
        """
        self.omega = self.rpm / 60 * 2 * np.pi
        self.flow_Q = self.flow_q * 60
        self.Ns = self.rpm * (self.flow_Q ** 0.5) / (self.head_total ** 0.75)
        self.type_number = self.omega * self.flow_q**0.5 / (self.g * self.head_total)**0.75
        
    def design_size(self):
        """
        Culcurate inner diameter and radius
        """
        self.D1 = self.D2 * self.R1R2ratio
        self.R2 = self.D2 / 2
        self.R1 = self.D1 / 2
        
    def design_inlet(self):
        """
        design impeller inlet
        """
        self.head_th = self.head_total / self.eta_h
        self.Ns_th = self.rpm * (self.flow_Q ** 0.5) / (self.head_th ** 0.75)
        
        self.head_a = self.press_in / self.g / self.density * 10**6
        self.NPSH = self.head_a - self.height_pump - self.head_loss - self.head_v
        
        self.tanbeta1 = (self.rambda1 / 2 / (self.rambda1 + self.rambda2)) **(0.5)
        self.beta1 = np.rad2deg(np.arctan(self.tanbeta1))
        self.S = 757.7 * (1 - self.hubratio**2)**0.5 \
                 / ((self.rambda1**0.5) * (self.rambda1+self.rambda2)**0.25)
        self.head_sv = (1.0/60)**2 * (4*np.pi)**(2./3) / 2 / self.g \
                       * (self.rpm * np.sqrt(self.flow_Q) / (np.sqrt(1 - self.hubratio**2)) * self.tanbeta1) \
                       * (self.rambda2 + self.rambda1 + self.rambda1 / self.tanbeta1 **2) 
        self.coef_head_sv = (self.rambda1 / np.sin(np.deg2rad(self.beta1))**2 \
                            + self.rambda2) / 2 / self.g
        self.v1 = self.rpm **(2./3) * self.flow_Q **(1./3) / 57.35
        self.Ds = 1.1 / np.sqrt(1 - self.hubratio**2) * (self.flow_Q / self.rpm)**(1./3)
        self.Dh = self.Ds * self.hubratio
    
    def _func(self, x):
        return self.sigma * x**3 - 0.5 * x - self.chi**2 / np.tan(np.deg2rad(self.beta2))
        
    def design_outlet(self):
        """
        design impeller outlet
        """
        x1 = [110, 200, 300, 400, 500, 600, 800, 1000, 1200]
        y1 = [5.2, 4.0, 3.2, 2.8, 2.5, 2.4, 2.0, 1.8, 1.75]
        y2 = UnivariateSpline(x1,y1, k=5,s=1)
        self.sqrtD2B2 = y2(self.Ns_th)
        self.B2 = self.D2 / self.sqrtD2B2 **2
        """ Slip factor """
        self.epsilon_limit = 1 / np.exp(8.16 * np.sin(np.deg2rad(self.beta2)) / self.Z)
        if self.R1R2ratio < self.epsilon_limit:
            self.sigma = 1 - np.sqrt(np.sin(np.deg2rad(self.beta2))) / self.Z **0.7
        else:
            self.sigma = (1 - np.sqrt(np.sin(np.deg2rad(self.beta2))) / self.Z **0.7) \
                         * (1 - ((self.R1R2ratio - self.epsilon_limit) / (1 - self.epsilon_limit))**3)
        self.k = 1 - self.sigma
        """ Outlet Non-dimensional velocity """
        self.chi = self.Ns_th * self.sqrtD2B2 / 2320
        self.Ku = fsolve(self._func, 1)
        self.Kc = 0.5 / self.Ku
        self.Km = (self.chi / self.Ku)**2
        self.Ka = (self.Kc**2 + self.Km**2)**0.5
        self.Kw = self.Km / np.sin(np.deg2rad(self.beta2))
        self.alpha2 = np.rad2deg(np.arctan(self.Km / self.Kc))
        """ Outlet velocity """
        self.sqrt2gHth = np.sqrt(2 * self.g * self.head_th)
        self.u2 = self.Ku * self.sqrt2gHth
        self.vm2 = self.Km * self.sqrt2gHth
        self.vu2 = self.Kc * self.sqrt2gHth
        self.v2 = self.Ka * self.sqrt2gHth
        self.w2 = self.Kw * self.sqrt2gHth
            
    def _func2(self, x):
        return self.vu2 * self.R2 * x**2 - self.flow_Q/120 * x - self.R3 * self.flow_Q/60
    
    def design_volute(self):
        """
        design volute casing
        """
        self.R3 = self.R2 * 1.05
        self.B3 = fsolve(self._func2, 1)
        self.T = self.B3
        self.Hth_impeller = self.u2\
                            * (self.sigma * self.u2 - self.vm2 / np.tan(np.deg2rad(self.beta2))) / self.g
        self.Hth_volute = np.pi * self.rpm\
                          / (1800*self.g*self.B3*np.log(1+self.T/self.R3))* self.flow_Q        
        self.flow_Q_max = 60 * self.sigma * self.u2 * (1/(np.pi*self.D2*self.B2*np.tan(np.deg2rad(self.beta2)))\
                          + 1/(self.R2*self.B3*np.log(1+self.T/self.R3)))**-1
        
    def show(self):
        print "==== Centrifugal Pump ====="
        print "Working fluid densigy:\t%d (kg/m3)" % (self.density)
        print "Rotation speed N:\t%d (rpm)" % (self.rpm)
        print "Potation speed N:\t%.1f (rad/s)" % (self.omega)
        print "Pump out pressure:\t%.1f (MPa)" % (self.press_out)
        print "Mass flow Q:\t%.2f (m3/min)" % (self.flow_Q)
        print "Mass flow q:\t%.2f (m3/sec)" % (self.flow_q)
        print "Total Head H:\t%.2f (m)" % (self.head_total)
        print "Specific Speed Ns:\t%.1f (m min)" % (self.Ns)
        print "Type Number K:\t%.2f" % (self.type_number)
        print "---- Impeller inlet design ----"
        print "NPSH:\t\t%.2f (m)" % (self.NPSH)
        print "tan(beta1):\t\t%.3f" % (self.tanbeta1)
        print "Inlet beta1:\t%.2f (deg)" % (self.beta1)
        print "suction sp. speed S:\t%.2f" % (self.S)
        print "Required NPSH Hsv:\t%.2f (m)" % (self.head_sv)
        print "Coef of Hsv: \t%.2f" % (self.coef_head_sv)
        print "Inlet v1:\t\t%.3f (m/s)" % (self.v1)
        print "Diameter shroud Ds:\t%.2f (mm)" % (self.Ds * 1000)
        print "Diameter hub Dh:\t%.2f (mm)" % (self.Dh * 1000)
        print "---- Impeller outlet design ----"
        print "sqrt(D2/B2):\t%.2f" % (self.sqrtD2B2)
        print "Theoretical Ns_th:\t%.2f" % (self.Ns_th)
        print "Slip factor sigma:\t%.2f" % (self.sigma)
        print "Parameter chi:\t%.2f" % (self.chi)
        print "Outer diameter D2:\t%.2f (mm)" % (self.D2 * 1000)
        print "Outer width B2:\t%.2f (mm)" % (self.B2 * 1000)
        print "Ku:\t\t%.3f" % (self.Ku)
        print "Kc:\t\t%.3f" % (self.Kc)
        print "Km:\t\t%.3f" % (self.Km)
        print "Ka:\t\t%.3f" % (self.Ka)
        print "Kw:\t\t%.3f" % (self.Kw)
        print "Outer alpha2:\t%.2f (deg)" % (self.alpha2)
        print "---- Volute casing design ----"
        print "Volute Radius R3:\t%.2f (mm)" % (self.R3 * 1000)
        print "Volute Width B3:\t%.2f (mm)" % (self.B3 * 1000)
        print "Volute Hight T:\t%.2f (mm)" % (self.T * 1000)
        print "Hth(Impeller):\t%.2f (m)" % (self.Hth_impeller)
        print "Hth(Volute):\t%.2f (m)" % (self.Hth_volute)
        print "Q when max efficiency:\t%.2f (m3/min)" % (self.flow_Q_max)
        
    def plot_Non_dimensional_velocity_dia(self):
        """
        plot Non-dimensional velocity diagram
        """
        a = self.sqrt2gHth
        plt.figure()
        plt.plot([0,a*self.Ku],[0,0],'k-')
        plt.plot([0,a*self.Kc],[0,a*self.Km],'k-')
        plt.plot([a*self.Kc,a*self.Kc + a*self.k * self.Ku],[a*self.Km,a*self.Km],'k-')
        plt.plot([a*self.Kc + a*self.k*self.Ku, a*self.Ku],[a*self.Km,0],'k-')
        plt.axes().set_aspect('equal','datalim')
        plt.xlabel("circumferential direction velocity (m/s)")
        plt.ylabel("diametrical direction velocity (m/s)")
        plt.title("velocity diagram")
        plt.show()
        
    def plot_slip_factor_dia(self):
        """
        plot slip factor diagram
        """
        plt.figure()
        beta = 25
        beta = np.deg2rad(beta)
        array_R1R2 = np.linspace(0,1)
        array_Z = np.array([2,4,8,16,32])
        array_epsilon = np.zeros(len(array_Z))
        array_sigma = np.zeros((len(array_Z),len(array_R1R2)))
        for z in range(len(array_Z)):
            array_epsilon[z] = 1 / np.exp(8.16 * np.sin(beta) / array_Z[z])
            for i in range(len(array_R1R2)):
                if array_R1R2[i] < array_epsilon[z]:
                    array_sigma[z][i] = 1 - np.sqrt(np.sin(beta)) / array_Z[z] **0.7
                else:
                    array_sigma[z][i] = (1 - np.sqrt(np.sin(beta)) / array_Z[z] **0.7) \
                         * (1 - ((array_R1R2[i] - array_epsilon[z]) / (1 - array_epsilon[z]))**3)
            plt.plot(array_R1R2, array_sigma[z],label="Z=%d"%(array_Z[z]))
        plt.xlabel("R1/R2")
        plt.ylabel("sigama")
        plt.title(r"Slip factor diagram ($\beta_2$=25deg)")
        plt.axis([0.0,1.0,0.0,1.0])
        plt.legend(loc=3)
        plt.grid()
#        plt.axes().set_aspect('equal')
        plt.show()
        
    def plot_D2B2forNs_dia(self):
        """
        plot Approximate sqrt{D2B2} value for Ns
        """
        x = [110, 200, 300, 400, 500, 600, 800, 1000, 1200]
        y = [5.2, 4.0, 3.2, 2.8, 2.5, 2.4, 2.0, 1.8, 1.75]
        spline = UnivariateSpline(x,y, k=5,s=1)
        Ns = np.linspace(110,1200)
        D2B2 = spline(Ns)
        plt.figure()
        plt.loglog(Ns,D2B2)
        plt.xlabel("Ns")
        plt.ylabel(r" $\sqrt{D_2/B_2} $" )
        plt.title(r"Approximate $\sqrt{D_2/B_2}$ Value for Ns (Z=6)")
        plt.xlim(100,1200)
        plt.ylim(1,7)
        plt.grid(True,which="both")
        plt.xticks(np.arange(100,1200,100),
                   [100,200,300,400,500,600,"",800,"",1000,"",1200])
        plt.yticks(np.arange(1,7,1),
                   [1,2,3,4,5,6,7])
        plt.show()
        
def main():
    rocket = Rocket(200000, 230, 1.6)
    rocket.culc()
    rocket.show()
    
    pump = Pump(rocket,"fuel")
    pump.culc_Ns()
    pump.design_size()
    pump.design_inlet()
    pump.design_outlet()
    pump.design_volute()
    pump.show()
    pump.plot_Non_dimensional_velocity_dia()
    pump.plot_slip_factor_dia()
    pump.plot_D2B2forNs_dia()


if __name__ == '__main__':
    main()