

from __future__ import division
import warnings,datetime, os.path,csv,  matplotlib
warnings.simplefilter(action = "ignore", category = RuntimeWarning)
import matplotlib.pyplot as plt
import numpy as np

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

# ================================================================================================================ #
# ============  Note ! Here we calculate the data averaging window of the AMS-02 & PAMELA observation =============# 
# ================================================================================================================ #

# Function to extract B-OR-C rotation dates & type experiment....
def Date_Extraction(Rlines):
    Date_start =str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-1]) +"/"+str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-2]) +"/"+str(Rlines[6].split(';')[1].split(',')[0].split()[0].split("-")[-3]) 
    Date_end =str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-1]) +"/"+str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-2]) +"/"+str(Rlines[6].split(';')[1].split(',')[-1].split()[0].split("-")[-3]) 
    Date =Date_start+"-"+Date_end
    Exper =Rlines[6].split()[3]
    return Date,Exper


# Save figures as PDFs 
FigureZ  =  r"./Figures/"
pdfsav_AverDays= PdfPages(FigureZ+"AverDays_.pdf")
  
 
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
 
# Plot the Duration between two Dates....
def date_difference(d1, d2):
    rD = (d2-d1).days
    rD +=1    
    return rD

 
# ================ Figure [1]........
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
 


