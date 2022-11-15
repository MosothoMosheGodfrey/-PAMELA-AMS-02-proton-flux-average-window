

from __future__ import division
import warnings,datetime, os.path,csv,  matplotlib
warnings.simplefilter(action = "ignore", category = RuntimeWarning)
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

# Setting up standard style....
from matplotlib import rcParams
plt.rcParams['agg.path.chunksize'] = 10000
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Computer Modern Roman'
plt.rcParams['font.sans-serif'] = 'Computer Modern Roman'
plt.rcParams['font.monospace'] = 'Computer Modern Roman'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams["axes.labelweight"] = "bold"#"medium"
plt.rcParams['axes.labelsize'] =23
plt.rcParams["axes.linewidth"] = 1.5
plt.rcParams['axes.labelcolor'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'  
plt.rcParams['axes.facecolor'] = 'white'  
#plt.rcParams['axes.facecolor'] = 'black'
plt.rcParams['axes.edgecolor'] = 'black'
plt.rcParams['font.size'] =23
plt.rcParams['legend.fontsize'] =23
plt.rcParams['legend.fancybox'] = True
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams['lines.antialiased'] = True
matplotlib.rc('text',usetex=True)
plt.rcParams['text.color'] = 'black'
plt.rcParams['text.hinting'] = 'auto'
rcParams['mathtext.default'] = 'regular'
plt.rcParams['mathtext.fontset'] = 'cm'

plt.rcParams['mathtext.default'] = 'regular'
plt.rcParams['mathtext.default'] = 'rm' # sf
rcParams['pdf.compression'] = 0
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.color'] = 'gray'
plt.rcParams['grid.linestyle'] = '--'
plt.rcParams['grid.linewidth'] = '1.2'
plt.rcParams['grid.alpha'] = 0.2
# visible ticks
plt.rcParams['xtick.labelsize'] =20
plt.rcParams['ytick.labelsize'] =20
plt.rcParams["xtick.major.width"] = 1.5
plt.rcParams["ytick.major.width"] = 1.5
plt.rcParams["xtick.minor.width"] = 1.5
plt.rcParams["ytick.minor.width"] = 1.5
rcParams["xtick.minor.visible"] = True
rcParams["ytick.minor.visible"] = True
rcParams["xtick.top"] = False
rcParams["ytick.right"] = True
rcParams["xtick.major.size"] = 12
rcParams["xtick.major.pad"] = 10
rcParams["ytick.major.size"] = 12
rcParams["ytick.major.pad"] = 10
rcParams["xtick.minor.size"] = 6
rcParams["ytick.minor.size"] = 6
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams['xtick.color'] = 'black'
plt.rcParams['ytick.color'] = 'black'
 
# Save figures as PDFs
from matplotlib.backends.backend_pdf import PdfPages
FigureZ  =  r"./Figures/"
pdfsav_AverDays= PdfPages(FigureZ+"AverDays_.pdf")
pdfsav_FluxComp= PdfPages(FigureZ+"FluxComp2014_.pdf")
  
  
# ================================================================================================================ #
# ============  Note ! Here we :                                                                             ======# 
# ============                 1. ...calculate the data averaging window of the AMS-02 & PAMELA observation  ======# 
# ============                 2. ...compare the flux data energy ranges & trim data only where both         ======#
# ============                       experiments overlap                                                     ======# 
# ================================================================================================================ #

# Declare variable/constant_values,etc.....
A,Z,E0 = np.array(1,float),np.array(1,float),np.array(0.938,float)   #  E0 in  GeV/n.


# Function to extract B-OR-C rotation dates & type experiment....
def Date_Extraction(Rlines):
    Date_start =str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-1]) +"/"+str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-2]) +"/"+str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-3]) 
    Date_end =str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-1]) +"/"+str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-2]) +"/"+str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-3]) 
    Date =Date_start+"-"+Date_end
    Exper =Rlines[6].split()[3]
    return Date,Exper
# Plot the Duration between two Dates....
def date_difference(d1, d2):
    rD = (d2-d1).days
    rD +=1    
    return rD
# Converts energy (in GeV) to rigidity (in GV)
def Prot_energy2rigidity(epn):
    epn  =  np.array(epn)
    Pz =  np.array(A/Z)*np.sqrt(np.array(epn,float)*np.array(epn,float) +2*np.array(epn,float)*np.array(E0,float))
    return Pz
#                     <<<<<<<<<<<< Mosotho & Strauss 2020 LIS >>>>>>>>>>>>>>>         #
def LIS20P(En_): # Units   (m2.sr.s.GeV)^-1 .
    En_  =  np.array(En_)
    #jLISP = (20.98)*(((0.729 + En_**(0.795))**-3.77641509434))*(En_**0.179)
    jLISP = (20.98)*(((0.777*E0  + En_**(0.795))**-3.77641509434))*(En_**0.179)
    return jLISP*1000 
# Convert Energy spectra to rigidity spectra.... [see  Boschinia et al. 2022,https://doi.org/10.1016/j.asr.2022.03.026,Eq. (2)]

def dJdP(T,_Flux):
    Rn  =  np.array(Prot_energy2rigidity(np.array(T)))
    #Intens  =  np.array(np.array((A/Z)*(T + E0)*np.power( np.power(np.power(T,2) + 2*T*E0,0.5),-1.0))* np.array(_Flux))
    Xz  =   np.array((A/Z)*(A/Z)*Rn/np.sqrt(Rn*Rn*(Z/A)*(Z/A) + E0))
    Xz  =   np.array((A/Z)*(A/Z)*Rn/( T+ E0))
    Intens  =   Xz* np.array(_Flux)
    return Rn,Intens

def dTdP(T):
    Rtn  =  np.array((Z/A)*(np.power(np.power(T,2) + 2*T*E0,0.5))/(T + E0))
    return Rtn
def dTdP_bins(T_low,T_high):
    P_high,P_low  =  np.array(Prot_energy2rigidity(np.array(T_high)),Prot_energy2rigidity(np.array(T_low)))
    Rtn_  =  np.array((T_high - T_low)/np.array(P_high - P_low))
    return Rtn_


# ======================================== Find : Energy/rigidity range of interest !  ======================================= #

# Use only overlaping dataset within the same rigidity range.............. 

# PAMELA Energy range.... Extracted from actual data. Visit   https://tools.ssdc.asi.it/CosmicRays/
KE_ams=[0.49246985,0.62080371,0.76383548,0.92531532,1.10505384,1.3027742,1.52278626,1.7649741,2.03390915,2.32952873,2.65171712,3.00518957,3.38987119,3.81056494,4.27213138,4.77453613,5.3177229,5.90657771,6.54601423,7.23599804,7.98145873,8.7873403,9.65361668,10.60019067,11.59714435,12.64442777,13.79189627,15.03955751,16.38741036,17.83544785,19.38365948,21.03203286,22.83051581,24.77911189,26.82784852,29.02668495,31.37561701,33.87463914,36.57372942,39.47288769,42.57211197,45.8713991,49.42073652] 
# AMS-02 Energy range.... Extracted from actual data. Visit   https://tools.ssdc.asi.it/CosmicRays/
KE_pam=[0.0885,0.1,0.1075,0.115,0.125,0.14,0.155,0.165,0.18,0.19995,0.215,0.23,0.25,0.275,0.3,0.325,0.35,0.375,0.41,0.445,0.48,0.52,0.56,0.6,0.645,0.695,0.745,0.8,0.86,0.925,0.99,1.055,1.13,1.21,1.295,1.38,1.47,1.57,1.67,1.775,1.89,2.01,2.135,2.265,2.405,2.55,2.7,2.86,3.03,3.21,3.395,3.59,3.905,4.355,4.85,5.395,5.99,6.645,7.365,8.155,9.025,9.98,11.025,12.17,13.43,14.815,16.33,17.99,19.81,21.805,24.,26.405,29.04,31.93,35.095,38.57] 
KEmin  = max([min(KE_pam),min(KE_ams)])
KEmax  = min([max(KE_pam),max(KE_ams)])
KE  = []
for YE in sorted(KE_pam +KE_ams):
    if KEmin  <=  YE <=  KEmax :
        KE.append(np.array(YE,float))
    else:
        pass
##Same_Rigidity_range=np.array(Same_Rigidity_rangeX)
Same_Rigidity_range=np.array(Prot_energy2rigidity(np.array(KE,float)))
# Rigidity boundaries.........
Pmin_  =  min(Same_Rigidity_range)
Pmax_  =  max(Same_Rigidity_range)
# Convert LIS spectra from energy to rigidity.......
T=np.logspace(np.log10(0.09),np.log10(55),200)
LIS_p  =  LIS20P(T)
Prigid_,LIS_p  =  dJdP(T,LIS_p) 
Flux_boundry  =  [min(LIS_p),max(LIS_p)]
# ================ Figure [1]........
Fig_Flux_0 =plt.figure(figsize =(14,8),facecolor="w",edgecolor="black")
plt.subplots_adjust(left =None,bottom =None,right =None,top =None,wspace =0,hspace =0.01)
axSubs2 =Fig_Flux_0.add_subplot(111)
Fig_Flux_0.text(0.5,0.04,"Rigidity, P [GV]",ha ="center",color ="k",fontsize=25)
axSubs2.plot([Pmin_,Pmin_],[min(Flux_boundry),max(Flux_boundry)],linewidth  =  2,color = "gray",linestyle = "--")
axSubs2.plot([Pmax_,Pmax_],[min(Flux_boundry),max(Flux_boundry)],linewidth  =  2,color = "gray",linestyle = "--")
axSubs2.fill([Pmin_,Pmax_,Pmax_,Pmin_],[min(Flux_boundry),min(Flux_boundry),max(Flux_boundry),max(Flux_boundry)],fill = False,hatch = "\/\/\/",alpha = 0.3,color = "silver")
axSubs2.annotate(r"\textbf{Rigidity range of interest : "+ str(round(Pmin_,2))+"-"+str(round(Pmax_,2))+ " GV}",xy =(6,min(Flux_boundry)+0.1),xycoords ="data",textcoords ="data",ha ="center",color ="black" )
axSubs2.annotate("",xy =(Pmin_,min(Flux_boundry)+0.05),xytext =(Pmax_,min(Flux_boundry)+0.05),arrowprops =dict(arrowstyle ="<->"),color ="black")
axSubs2.set_ylabel(r"${\mathcal{J}}\left({\rm{P}},t\right)$ $\left[ {\rm{GV}.\rm{m^2.s.sr}}\right]^{-1} $",fontsize=25)
# LIS....
Cb4,=axSubs2.loglog(Prigid_,LIS_p,label = "LIS",linewidth  =  2,color = "black",linestyle = "--")

 


 
##  Directory : PAMELA/AMS Fluxes.....
Dir_ = r"./Data_/PAMELA_AMS02/"
f  =  sorted(os.listdir(Dir_))# Get all text files
 
Date_TimeFull,ExpRiment  =[],[] 
for Q0,_clr2 in zip(f,plt.cm.rainbow(np.linspace(0,1,len(f)))):
    Q  =  Dir_ +Q0
    with open(str(Q),"r") as f:
        line  = f.readlines()
        lines  = line[11:-2]
    FullDate = Date_Extraction(line)[0]
    Date_TimeFull.append(FullDate)
    ExpR_ = Date_Extraction(line)[-1]
    ExpRiment.append(ExpR_)
    YearSart,MonthStart =FullDate.split("-")[0].split("/")[-1],FullDate.split("-")[0].split("/")[-2] 

    EnergY,Energl,Energh,J,Flux_Exp_l,Flux_Exp_h  = [],[],[],[],[],[]
    if ExpR_  == "PAMELA":
        for hh in lines:
            J.append((float(hh.split()[3])))
            Flux_Exp_l.append(float(hh.split()[4]))
            Flux_Exp_h.append(float(hh.split()[5]))
            EnergY.append(float(hh.split()[0]))
            Energl.append(float(hh.split()[1]))
            Energh.append(float(hh.split()[2]))
    elif ExpR_  == "AMS-02":
        for hh in lines:
            J.append((float(hh.split()[3])))
            Flux_Exp_l.append(float(hh.split()[4]))
            Flux_Exp_h.append(float(hh.split()[5]))
            EnergY.append(float(hh.split()[0]))
            Energl.append(float(hh.split()[1]))
            Energh.append(float(hh.split()[2]))
    else:
        pass
    # Convert to numpy array...
    EnergY  =  np.array(np.array(EnergY,str),float)
    Flux_Exp_  =  np.array(J,float)
    Flux_Exp_h  =  np.array(np.array(Flux_Exp_h,str),float)
    Flux_Exp_l  =  np.array(np.array(Flux_Exp_l,str),float)
    # Convert Flux: From Energy to Rigidity...
    if ExpR_  == "PAMELA": 
        Full_Rigidity,Full_djdp_Rigidity  =   dJdP(EnergY,Flux_Exp_) 
        _,Full_djdp_h_Rigidity  =   dJdP(EnergY,Flux_Exp_h) 
        _,Full_djdp_l_Rigidity  =   dJdP(EnergY,Flux_Exp_l) 
        axSubs2.set_ylim(min(Flux_boundry),max(Full_djdp_Rigidity)+300)
    elif ExpR_  == "AMS-02":
        Full_Rigidity,Full_djdp_Rigidity,Full_djdp_h_Rigidity,Full_djdp_l_Rigidity = EnergY,Flux_Exp_,Flux_Exp_h,Flux_Exp_l
    else:        
        print(" --- Re-check your data directory if empty ! ---")
        sys.exit()
    Full_Rigidity,Full_djdp_Rigidity =(list(Full_Rigidity)),(list(Full_djdp_Rigidity))
    Full_djdp_l_Rigidity,Full_djdp_h_Rigidity =(list(Full_djdp_l_Rigidity)),(list(Full_djdp_h_Rigidity))
    
    # ========================== Data-trimming ! Only where both PAMELA and AMS-02 overlap.....  ! 
    # Interpolate overlaping values from the AMS/PAMELA energy ranges.....
    IndxFlux  =  interpolate.splrep(Full_Rigidity,Full_djdp_Rigidity)
    FLUX_Data_SameRigidity_  =  np.array(interpolate.splev(np.array(Same_Rigidity_range),IndxFlux))
    IndxFlux_l  =  interpolate.splrep(Full_Rigidity,Full_djdp_l_Rigidity)
    FLUX_Data_SameRigidity_l=  np.array(interpolate.splev(np.array(Same_Rigidity_range),IndxFlux_l))
    IndxFlux_h  =  interpolate.splrep(Full_Rigidity,Full_djdp_h_Rigidity)
    FLUX_Data_SameRigidity_h=  np.array(interpolate.splev(np.array(Same_Rigidity_range),IndxFlux_h))
    
    ### =================== Plot 2014 Rigidity Flux  ================== !
    if YearSart ==  "2014"  :
        if   MonthStart == "01":
            if ExpR_  == "PAMELA": 
                axSubs2.errorbar(Full_Rigidity,Full_djdp_Rigidity,xerr=None,yerr = np.array([np.array(Full_djdp_l_Rigidity),np.array(Full_djdp_h_Rigidity)]) ,fmt = "o",mfc = "white",**dict(ecolor = "gray",color = "gray",capsize = 2.5,elinewidth = 0.8,linestyle = "none"))
                Cb0,=[axSubs2.errorbar(Same_Rigidity_range,FLUX_Data_SameRigidity_,label = str(ExpR_)+" : "+"[Interpolated]",xerr=None,yerr = np.array([np.array(FLUX_Data_SameRigidity_l),np.array(FLUX_Data_SameRigidity_h)]) ,fmt = "o",mfc = "blue",**dict(ecolor = "blue",color = "blue",capsize = 2.5,elinewidth = 0.8,linestyle = "none"))]
            else:
                axSubs2.errorbar(Full_Rigidity,Full_djdp_Rigidity,xerr=None,yerr = np.array([np.array(Full_djdp_l_Rigidity),np.array(Full_djdp_h_Rigidity)]) ,fmt = "*",markersize=10,mfc = "white",**dict(ecolor = "gray",color = "gray",capsize = 2.5,elinewidth = 0.8,linestyle = "none"))
                Cb1,=[axSubs2.errorbar(Same_Rigidity_range,FLUX_Data_SameRigidity_,label = str(ExpR_)+"  "+"[Interpolated]",xerr=None,yerr = np.array([np.array(FLUX_Data_SameRigidity_l),np.array(FLUX_Data_SameRigidity_h)]),fmt = "*",markersize=10,mfc = "red",**dict(ecolor = "red",color = "red",capsize = 2.5,elinewidth = 0.8,linestyle = "none"))]    
linez_2=[Cb4,Cb0,Cb1]
axSubs2.legend(linez_2,[l.get_label() for l in linez_2],title=' Jan. 2014',shadow=False,frameon=False,ncol=1,loc=9,bbox_to_anchor=(0.54,0.5),labelspacing=0.1,handlelength=1,handletextpad=0.5,columnspacing=0.7,markerscale=1.2) 
pdfsav_FluxComp.savefig(Fig_Flux_0,dpi =250)
pdfsav_FluxComp.close()
 
 
# ================ Figure [2]........
Fig_Dates =plt.figure(figsize =(14,8),facecolor="w",edgecolor="black")
plt.subplots_adjust(left  = 0.06,bottom  = 0.11,right  = 0.98,top  = 0.98,wspace  = None,hspace  = None)
ax_ =Fig_Dates.add_subplot(111)
Ftime, AverWind,EXP =[],[],[]
for hh, ExpRi  in zip(Date_TimeFull, ExpRiment ):
    d1 = datetime.datetime.strptime( str(hh.split(",")[0].split("-")[0])  ,"%d/%m/%Y") 
    d2 = datetime.datetime.strptime( str(hh.split(",")[0].split("-")[-1])  ,"%d/%m/%Y") 
    Tim = [ str(hh.split(",")[0].split("-")[-1])]
    Diff = [ date_difference(d1, d2)]
    CombinedTIME = [datetime.datetime.strptime(d,"%d/%m/%Y") for d in Tim]
    if str(hh.split(",")[0].split("-")[-1]) =='26/07/2006' and str(ExpRi) !=str("AMS-02"):       
        plt.plot(CombinedTIME, Diff , color ='red',ls ="none" ,marker ="s",markeredgecolor='red',markersize=5,label = "PAMELA")
        ax_.plot([CombinedTIME,CombinedTIME], [[-1],Diff] , color ='red',ls ="-" ,label = None)
    elif str(hh.split(",")[0].split("-")[-1]) =='09/05/2017' and str(ExpRi) ==str("AMS-02"):
        ax_.plot(CombinedTIME, Diff , color ='blue',ls ="none" ,marker ="s",markeredgecolor='blue',markersize=5,label = "AMS-02")
        ax_.plot([CombinedTIME,CombinedTIME], [[-1],Diff] , color ='blue',ls ="-" ,label = None)
    else:
        if str(ExpRi) ==str("AMS-02"): CLR ='blue'
        else: CLR ='red'
        ax_.plot( CombinedTIME,Diff , color =CLR,ls ="none" ,marker ="s",markeredgecolor=CLR,markersize=5,label = None)
        ax_.plot([CombinedTIME,CombinedTIME], [[-1],Diff] , color =CLR,ls ="-" ,label = None)
    AverWind.append(Diff[0])
    Ftime.append(CombinedTIME[0])
    EXP.append(ExpRi)
ax_.legend(shadow=False,frameon=False,ncol=2,loc=1,fontsize=20,labelspacing=0.1,handlelength=1,handletextpad=0.5,columnspacing=0.7,markerscale=1.2)
ax_.set_xlim(datetime.datetime.strptime("2006-01-01","%Y-%m-%d"),datetime.datetime.strptime("2018-01-01","%Y-%m-%d"))
ax_.set_ylim( 0,31) 
ax_.set_xlabel( r"$ t_f$  [day/month/year]", color="k" ) 
ax_.set_ylabel( r"$\Delta t = t_f - t_s$  [total number of days] ", color="k" ) 
ax_.fill_between([datetime.datetime.strptime("2011-06-10","%Y-%m-%d"),datetime.datetime.strptime("2014-02-11","%Y-%m-%d")],[40,40], hatch="*", color='gray',linewidth=1, alpha=0.2,ls='-')
ax_.annotate(r"Average window : $t_s-$Initial time window, $t_f-$ Final time window",xy =(datetime.datetime.strptime("2008-12-05","%Y-%m-%d"),30),xycoords ="data",textcoords ="data",ha ="center",color ="k",fontsize=15)

pdfsav_AverDays.savefig(Fig_Dates,dpi =250)
pdfsav_AverDays.close()
 

import datetime
#==================== * Create output file to save data  =====================#
Heading  = "# Generated data from TOA Flux Measurements \n"
Heading += "# AMS-02 & PAMELA proton Flux data obtained from SSDC cosmic rays database www.ssdc.asi.it \n"
Heading += "# File generated date ====>> "+str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))+"\n" 
Heading += "# Format : Final date-time [d.m.y], Time difference [days], Experiment  " + "\n"
# Create a text file
Data_name = str("./Data_/Flux_TimeOfMeasurements.txt")
with open(Data_name,"w") as S: 
    S.write(Heading)
    S.close()
with open(Data_name,"a") as fx:
    writerx = csv.writer(fx,delimiter=",")
    writerx.writerows(zip(list(Ftime),list(AverWind),list(EXP)))
    fx.close()

plt.show()
 


