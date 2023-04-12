import uproot
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt

events1 = uproot.open("~/tafc/Tafc_Project/ProjectResultsSipmSebastiao.root:DetectorEfficiency1;1")

Data = events1.arrays(["NSIMPs", "Nphotons_detected", "Nphotons_total"])
print(Data['NSIMPs'])
NSIMPs = [1,2,3,4,5,6,7,8]
sum_total = [0] * len(NSIMPs)
sum_accepted = [0] * len(NSIMPs)
for y in Data:
    x=y['NSIMPs']
    if x==8:
        sum_total[7]+=y["Nphotons_total"]
        sum_accepted[7]+=y["Nphotons_detected"]
    elif x==7:
        sum_total[6]+=y["Nphotons_total"]
        sum_accepted[6]+=y["Nphotons_detected"]
    elif x==6:
        sum_total[5]+=y["Nphotons_total"]
        sum_accepted[5]+=y["Nphotons_detected"]
    elif x==5:
        sum_total[4]+=y["Nphotons_total"]
        sum_accepted[4]+=y["Nphotons_detected"]
    elif x==4:
        sum_total[3]+=y["Nphotons_total"]
        sum_accepted[3]+=y["Nphotons_detected"]
    elif x==3:
        sum_total[2]+=y["Nphotons_total"]
        sum_accepted[2]+=y["Nphotons_detected"]
    elif x==2:
        sum_total[1]+=y["Nphotons_total"]
        sum_accepted[1]+=y["Nphotons_detected"]
    elif x==1:
        sum_total[0]+=y["Nphotons_total"]
        sum_accepted[0]+=y["Nphotons_detected"]

ratio= [0] * len(NSIMPs)
print(sum_accepted)
print(sum_total)
for i in range(len(sum_total)):
    ratio[i]=100*sum_accepted[i]/sum_total[i]
print(ratio)
plt.plot(NSIMPs, ratio)

#naming the x axis
plt.xlabel('Number SIPMs')
#naming the y axis
plt.ylabel('Detector Efficiency Scintillator 1 (%)')
#
# #giving a title to my graph
plt.title('Detector Efficiency for Scintillator 1 as a function of number os SIPMs')

plt.show()
