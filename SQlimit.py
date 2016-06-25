# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 17:54:16 2016

@author: marcus
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from  scipy.integrate import cumtrapz
plt.style.use('ggplot')
#plt.style.use(['ggplot','dark_background'])

try: # use seaborn to choose better colors 
    import seaborn as sns
    with_sns = True
except: # if you don't have seaborn
    with_sns = False



k = 1.38064852e-23   # m^2 kg s^-2 K^-1, Boltzmann constant
h = 6.62607004e-34   # m^2 kg s^-1     , planck constant
c = 2.99792458e8     # m s^-1          , speed of light
eV= 1.6021766208e-19 # joule           , eV to joule
q = 1.6021766208e-19 # C               , elemental charge


# http://rredc.nrel.gov/solar/spectra/am1.5/
ref_solar = pd.read_csv("ASTMG173.csv", header=1) # nm vs W m^-2 nm^-1
# data range: 280nm to 4000nm, 0.31eV to 4.42857 eV
WL, solar_per_nm = ref_solar.iloc[:,0], ref_solar.iloc[:,2] # WL (nm), W*m-2*nm-1

E = 1240.0/WL # eV
# jacobian transformation, W m^-2 eV^-1
solar_per_E= solar_per_nm*(eV/1e-9)*h*c/(eV*E)**2 

Es = np.arange(0.32, 4.401, 0.001)

# linear interpolation to get an equally spaced spectrum
AM15 = np.interp(Es, E[::-1], solar_per_E[::-1]) # W m^-2 eV^-1
AM15flux = AM15/(Es*eV) # number of photon m^-2 eV^-1


class SQlim(object):
    def __init__(self, T=300, EQE_EL=1.0, intensity=1.0):
        """
        T: temperature in K
        EQE_EL: radiative efficiency (EL quantum yield)
        intensity: light concentration, 1.0 = one Sun, 100 mW/cm^2
        """
        self.T = T
        self.EQE_EL=EQE_EL
        self.intensity = intensity
        self.Es = np.arange(0.32, 4.401, 0.001)
        self.Calculate()
    
    def Calculate(self):
        self.J0  = self.cal_E_J0()
        self.Jsc = self.cal_E_Jsc()
        self.Voc = self.cal_E_Voc()
        self.PCE = self.cal_E_PCE()
        self.FF  = self.PCE/(self.Voc*self.Jsc)*100.0*self.intensity  # unit:%
        return None

    def cal_E_Jsc(self):
        fluxcumm = cumtrapz(AM15flux[::-1], self.Es[::-1], initial=0)
        fluxaboveE = fluxcumm[::-1]*-1*self.intensity
        Jsc = fluxaboveE*q*0.1 # mA/cm^2  (0.1: from A/m2 to mA/cm2)
        return Jsc

    def cal_E_J0(self):
        '''
        Calculate and return E vs J0, the dark saturation current
        J0 = q * (integrate(phi dE) from E to infinity)  / EQE_EL
        phi is the black body radiation at T (flux vs energy)

        '''
        phi=2*np.pi*(self.Es*eV)**2*eV/(h**3*c**2)/(np.exp(self.Es*eV/(k*self.T))-1)
        #fluxcumm = sp.integrate.cumtrapz(phi[::-1], self.Es[::-1], initial=0)
        fluxcumm = cumtrapz(phi[::-1], self.Es[::-1], initial=0)
        fluxaboveE = fluxcumm[::-1]*-1
        J0 = fluxaboveE*q*0.1/self.EQE_EL # mA/cm^2  (0.1: from A/m2 to mA/cm2)
        return J0

    def cal_E_Voc(self):
        '''
        Calculate and return E vs Voc
        Voc = (kT/q)*[ ln(Jsc/J0 + 1 )]
        '''
        return (k*self.T/q)*np.log((self.Jsc/self.J0)+1)

    def cal_E_PCE(self):
        PCE = []
        for i, E in enumerate(self.Es):
            V = np.arange(0, E, 0.001)
            J =-1*self.Jsc[i] +self.J0[i]*(np.exp(q*V/(k*self.T))-1)
            PCE.append(-1*np.min(J*V)/self.intensity)
        return PCE

    def Sim_JV(self, Eg, Vstep=0.001, plot=False, Vmin=-0.5):
        if not self.Es[0] <= Eg <= self.Es[-1]:
            print "invalid bandgap \nvalid range: 0.32 to 4.4"
            return
        V = np.arange(Vmin, Eg, Vstep)
        paras = self.get_paras(Eg, toPrint=False)
        J0, Jsc = paras["J0"], paras["Jsc"]
        J =-1*Jsc + J0*(np.exp(q*V/(k*self.T))-1)
        mask = (J<=200)
        if plot:
            title = "Theoretical JV for Eg={0} eV".format(round(Eg,2))
            plt.figure()
            plt.plot(V[mask], J[mask],'r')
            plt.plot([-1,Eg],[0,0], 'k')
            plt.plot([0,0], [-2*Jsc,200],'k')
            plt.ylim(-1.5*Jsc, min(100, 1.5*Jsc))
            plt.xlim(-0.5,Eg)
            plt.xlabel("Voltage (V)", fontsize=16)
            plt.ylabel("Current density (mA/$\mathregular{cm^2}$)", fontsize=16)
            plt.tick_params(labelsize=16)
            plt.title(title)
            plt.tight_layout()
        return np.vstack([V,J]) # col1: V, col2:J
        
        
    def get_paras(self, Eg, toPrint=True):
        ''' input Eg, return or print the corresponding parameters
        '''
        if not self.Es[0] <= Eg <= self.Es[-1]:
            print "invalid bandgap \nvalid range: 0.32 to 4.4"
            return
        Voc = np.interp([Eg], self.Es, self.Voc)[0]
        Jsc = np.interp([Eg], self.Es, self.Jsc)[0]
        FF  = np.interp([Eg], self.Es, self.FF)[0]
        PCE = np.interp([Eg], self.Es, self.PCE)[0]
        J0  = np.interp([Eg], self.Es, self.J0)[0]
        if toPrint:
            print "Bandgap: %.3f eV ;"  % Eg, "J0=%.3g mA/cm^2" % J0, "\n"
            print "Voc = %.4g \t V"       % Voc
            print "Jsc = %.4g \t mA/cm^2" % Jsc
            print "FF  = %.2f \t %%"      % FF
            print "PCE = %.3f \t %%"      % PCE
            return
        para={}
        para["Voc"]=Voc
        para["Jsc"]=Jsc
        para["FF"] =FF
        para["PCE"]=PCE
        para["J0"]=J0
        return para

    def saveall(self, savename="SQ limit"):
        result=pd.DataFrame()
        result["Bandgap (eV)"]=self.Es
        result["Voc (V)"] = self.Voc
        result["Jsc (mA/cm^2)"] = self.Jsc
        result["FF (%)"] = self.FF
        result["PCE (%)"] = self.PCE
        result["J0 (mA/cm^2)"] = self.J0
        result.to_csv(savename+".csv", index=False)
        return result

    def plotall(self, xlims=(0.32,3.0)):

        fig, ax = plt.subplots(2,2, sharex=True)
        axs=[(0,0), (0,1), (1,0), (1,1)]
        ys = [self.Voc, self.Jsc, self.FF, self.PCE]
        ylabel=["Voc (V)","Jsc (mA/$\mathregular{cm^2}$)","FF (%)","PCE (%)"]
        for i, a in enumerate(axs):
            ax[axs[i]].plot(self.Es, ys[i])
            ax[axs[i]].set_ylabel(ylabel[i])
            plt.setp(ax[axs[i]].get_xticklabels(), visible=True)
            ax[axs[i]].set_xlabel("Bandgap (eV)")       
        ax[(0, 0)].set_xlim(xlims)
        plt.tight_layout()
        return

############################################################################
###################### Useful functions for plotting results #######


def E_loss(Eg, SQ=SQlim(), xmin=300, xmax=2500, savefig=False):
    """
    input bandgap Eg, plot the energy loss and the available energy
    return None
    
    """
    if Eg>4.4 or Eg<0.32:
        print "invalid bandgap \nvalid range: 0.32 to 4.4"
        return None
        
    xmax = max(xmax, 1240.0/Eg)
    WLs=np.arange(280,4001,1.0)
    AM15nm=np.interp(WLs, WL, solar_per_nm)
    plt.figure(figsize=(8,4.5))
    ax = plt.gca()
    # color options: darkgoldenrod, darkorange, yellow, black
    colors = {'therm':'lightcoral', 'extract':'gold', 
              'avail':'LightSkyBlue', 'trans':'grey'}

    para=SQ.get_paras(Eg, toPrint=False)
    factor = para["Voc"]*para["FF"]/100.0/Eg
    extract = AM15nm/(1240.0/WLs)*Eg*((1240.0/WLs)>=Eg)
    Eavail = extract*factor
    therm = AM15nm*((1240.0/WLs)>=Eg)

    extractloss = extract - Eavail
    thermloss = therm-extract
    transloss = AM15nm*((1240.0/WLs)<Eg)
    
    ax.fill_between(WLs, 0, transloss, linewidth=0, facecolor=colors['trans'])
    ax.fill_between(WLs, 0, therm, linewidth=0, facecolor=colors['therm'])
    ax.fill_between(WLs, 0, extract, linewidth=0, facecolor=colors['extract'])
    ax.fill_between(WLs, 0, Eavail, linewidth=0, facecolor=colors['avail'])
    
    
    E_tot=np.sum(AM15nm)
    E_pct= {'trans':np.sum(transloss)/E_tot, 'therm':np.sum(thermloss)/E_tot,
            'extract': np.sum(extractloss)/E_tot, 'avail':np.sum(Eavail)/E_tot}


    legends = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i], edgecolor=None)
               for i in ['trans', 'therm', 'extract', 'avail'] ]  
    
    labels = [
    "{:.1f}% Not Absorbed".format(100.0*E_pct['trans']),
    "{:.1f}% Thermalization Loss".format(100.0*E_pct['therm']),
    "{:.1f}% Extraction Loss".format(100.0*E_pct['extract']),
    "{:.1f}% Available Energy".format(100.0*E_pct['avail'])]

    ax.plot([Eg],[0])
    ax.legend(legends, labels, frameon=False, 
              fontsize=14, loc="upper right").draggable()

    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("Irradiance  (W $\mathregular{m^{-2}\ nm^{-1}}$)", size=18)
    ax.set_xlabel("Wavelength (nm)",size=18)
    ax.tick_params(labelsize=16)

    plt.tight_layout()
    if savefig:
        plt.savefig("available_E.pdf", transparent=True)

    losses = pd.DataFrame()
    losses["Wavelength"] = WLs
    losses["Thermalization Loss"] = thermloss
    losses["Extraction Loss"] = extractloss
    losses["Not Absorbed"] = transloss
    losses["Available"] = Eavail
    
    return losses



def available_E(Egs, SQ=SQlim(), 
                E_MPP=True, xmin=300, xmax=2500, savefig=False):
    """
    plot the theoretical maximum available energies from a series of 
    solar cells with different Egs. 
    Note: this is NOT the theoretical efficiencies for two-terminal 
    tandem cells. It's more like mechanically stacked tandem cells.
    
    Egs: an array-like object of bandgaps; float/int are accepted too (one Eg)
    
    E_MPP: whether to scale to MPP or not 
           False: Eavail = Eg * Jsc
           True : Eavail = Voc * Jsc *FF
    """
    # 1-J : 1.337
    # 2-J : (1.63,0.96) or (1.8, 1.1)
    # 3-J : (1.82, 1.16, 0.71)
    
    try: # if input Eg is a float or integer
        l = len(Egs)
    except:
        l, Egs = 1, [Egs]
    EgMax, Egmin = max(Egs), min(Egs)
    if EgMax>4.4 or Egmin<0.32:
        print "invalid bandgap \nvalid range: 0.32 to 4.4"
        return None
    xmax = max(xmax, 1240.0/Egmin)
    xmin = min(xmin, 1240.0/EgMax)
    Egs  = sorted(list(Egs) + [4.5]) # add a dummy 4.5 eV, a dummy Jsc 0.0


    WLs=np.arange(280,4001,1.0)
    AM15nm=np.interp(WLs, WL, solar_per_nm)
    plt.figure(figsize=(8,4.5))
    ax = plt.gca()
    # color options: darkgoldenrod, darkorange, yellow, black
    solarcolor = "gold"
    ax.fill_between(WLs, 0, AM15nm, linewidth=0.0, facecolor=solarcolor)
               
    factor = 1.0 
    E_subcell = pd.DataFrame()
    E_subcell["WL"] = WLs
    E_subcell["Solar"] = AM15nm
    PCEsubcell, tot_E = [], 1000.0

    Jscs = [SQ.get_paras(E, toPrint=False)["Jsc"] for E in Egs[:-1]] + [0.0]

    if with_sns:
        colors = sns.color_palette("husl", len(Egs))
    else:
        cm = plt.get_cmap('gist_rainbow') # gist_rainbow
        cNorm     = matplotlib.colors.Normalize(0, 1.2*(l-1))
        scalarMap = matplotlib.cm.ScalarMappable(norm=cNorm, cmap=cm)
        colors = [scalarMap.to_rgba(i) for i in xrange(l)]
        
    for n, Eg in enumerate(Egs[:-1]):
        color = colors[n]
        if E_MPP:
            SQ_E = SQlim(intensity=(Jscs[n]-Jscs[n+1])/Jscs[n])
            para = SQ_E.get_paras(Eg, toPrint=False)
            factor = para["Voc"]*para["FF"]/100.0/Eg
        mask = ((1240.0/WLs)>=Eg) * (Egs[n+1]>=(1240.0/WLs))
        Eavail = factor*AM15nm/(1240.0/WLs)*Eg*mask
        E_subcell["Eg="+str(round(Eg,3))] = Eavail
        PCEsubcell.append(100*np.sum(Eavail)/tot_E)
        ax.fill_between(WLs, 0, Eavail, facecolor=color, linewidth=0.0)
        #ax.fill_between(WLs, 0, Eavail, facecolor=color, linewidth=0.2)
                         
    ax.set_xlim(xmin, xmax)
    ax.set_ylabel("Irradiance  (W $\mathregular{m^{-2}\ nm^{-1}}$)", size=18)
    ax.set_xlabel("Wavelength (nm)",size=18)
    ax.tick_params(labelsize=16)
    
    legends = [plt.Rectangle((0, 0), 1, 1, facecolor=colors[i], 
                             edgecolor=None) for i in range(l)[::-1] ]
    if l==1:
        labels = ["Eg={0} eV, PCE={1:.1f}%".format(Egs[0], PCEsubcell[0])]
    else:
        labels = ["Cell {0}, Eg={1:.2f} eV, PCE={2:.1f}%".format(l-i, Egs[i], 
                   PCEsubcell[i])  for i in range(l)[::-1]]
    ax.legend(legends, labels, frameon=False, loc="upper right", 
              fontsize=14).draggable() 

    plt.tight_layout()
    if savefig:
        fname = "available_E_tandem"
        for E in Egs[:-1]:
            fname += "_{:.2f}".format(E) #str(round(Es,2))
        plt.savefig(fname+".pdf", transparent=True)
    plt.show()
    return E_subcell, PCEsubcell

   

def Concentration(Suns=[1,10,100,1000]):
    for sun in sorted(Suns,reverse=True):
        plt.plot(Es, SQlim(intensity=sun).PCE, label = str(sun)+" sun")
    plt.xlabel("Bandgap (eV)", size=28)
    plt.ylabel("Efficiency (%)", size=28)
    plt.xlim(0.3,3.0)
    plt.tick_params(labelsize=24)
    ax=plt.gca()
    ax.legend(loc='upper right').draggable()
    plt.tight_layout()
    

def Js_tandem(Es):
    P=[0]+[SQ.get_paras(E, toPrint=False)["Jsc"] 
           for E in sorted(Es,reverse=True)]
    J=[p-P[n]for n,p in enumerate(P[1:])][::-1]
    return np.round(np.array(J),2)


if __name__=="__main__":
    SQ = SQlim()
    

### Examples  ####

# plot Voc, Jsc, FF, PCE
SQ.plotall() # plot Voc, Jsc, FF, PCE
    
# get Voc, Jsc, FF, PCE for a specific bandgap
#SQ.get_paras(1.337)  # print out the values
#SQ.get_paras(1.337, toPrint=False) # get a dictionary of parameters
    

# plot available energy and loss for a specific bandgap  
#E_loss(1.3)  
    
# plot available energy for cells with many different bandgap materials
#available_E([0.9, .6,1.5,1.2]) 

### calculate and plot JV curve
#SQ.Sim_JV(1.3, plot=True)

    


