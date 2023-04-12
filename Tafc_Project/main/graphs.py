import uproot
import pandas as pd
from IPython.display import display
import matplotlib.pyplot as plt

events1 = uproot.open("~/ProjectResultsMariana.root:GeomEfficiency;1")
Data = events1.arrays(['distance', 'Nmuons_total', 'Nmuons_accepted'])
distances = [5,10,15,20,25]
print(len(distances))
sum_total = [0] * len(distances)
sum_accepted = [0] * len(distances)
for y in Data:
    x=y['distance']
    if x==25:
        sum_total[4]+=y['Nmuons_total']
        sum_accepted[4]+=y['Nmuons_accepted']
    elif x==20:
        sum_total[3]+=y['Nmuons_total']
        sum_accepted[3]+=y['Nmuons_accepted']
    elif x==15:
        sum_total[2]+=y['Nmuons_total']
        sum_accepted[2]+=y['Nmuons_accepted']
    elif x==10:
        sum_total[1]+=y['Nmuons_total']
        sum_accepted[1]+=y['Nmuons_accepted']
    elif x==5:
        sum_total[0]+=y['Nmuons_total']
        sum_accepted[0]+=y['Nmuons_accepted']
ratio= [0] * len(distances)
for i in range(5):
    ratio[i]=100*sum_accepted[i]/sum_total[i]
print(ratio)
plt.plot(distances, ratio)

# naming the x axis
plt.xlabel('Distances (cm)')
# naming the y axis
plt.ylabel('Geometrical Efficiencies (%)')

# giving a title to my graph
plt.title('Geometrical Efficiency as a function of Distance')

# function to show the plot
plt.show()
