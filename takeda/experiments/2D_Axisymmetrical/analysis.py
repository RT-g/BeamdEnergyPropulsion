import pandas as pd
import numpy as np
from math import pi
import matplotlib.pyplot as plt

columns = pd.read_csv('21021605/r_scale.csv').columns.tolist()
index = pd.read_csv('21021605/time_scale.csv', names=['time'])
df = pd.read_csv('21021605/plateau_pressure_data.csv', names=columns)
df = pd.concat([index, df], axis=1)
l = len(df)
F_list = []

for i in range(l):
    s = df.loc[i]
    F = 0
    for j in range(1,len(s)):
        r = float(s.index[j])
        p = s[j]
        # p = complex(s[j].replace('i','j')).real
        F += p*2*pi*r*0.005
    F_ave = F/pi/2.5**2
    F_list.append(F_ave)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
ax.plot(index.time.tolist(), F_list, color='blue', label='pressure')
ax.set_xlabel('t')
ax.set_ylabel('p_plateau')
plt.savefig('t_plateau.png')
plt.show()
