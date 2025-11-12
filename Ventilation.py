
DEBUGGING = True  # set true to enable some debug graphs and console outputs/dumps 


#user defined variables

H = [3,3,10]                     #room height
A_eff = [0.5,0.5,0.5]        #effective vent area
S = [100,100,30]                  #floor area
n = [5,5,0]                     #people
w = 100                         #power input per person
T_day = 16                      #max daytime temp
T_night = 10                    #min nighttime temp
people_leave_work = True

'''
Room numbering system:
                 _______
                 |
                 |     |
  _______________|     |
 |                     |
 |               |     |
 |       0       |     |
                 |  2  |
 |_______________|     |
 |                     |
 |               |     |
 |       1       |     |
                 |      
 |_______________|_____|

'''

#fundamental variables
c = 0.15                        #plume constant
B = 0.0281*(w/1000)             #buoyancy input per person

def get_ext_temp(t):
    return ((T_night-T_day)/2)*cos(2*pi*t/(3600*24))+(T_day+T_night)/2
def temp_from_g(g):
    return T_night+(g/g_real)*(273+T_night)
def g_from_temp(T):
    return g_real*(T-T_night)/(273+T_night)

#setup variables
h = [0.5*H[0],0.5*H[1],0.9*H[2]]    #interface height - needs initial value < H
g_real = 9.81                       #gravity on earth
dt = 1                              #timestep in seconds
steps = 24*3600*7                   #number of steps
t_start = 0                         #start time, in seconds past midnight
g_h = [g_from_temp(T_night+1),g_from_temp(T_night+1),g_from_temp(T_night+1)]             #effective gravity of the hot layer - needs initial value so g isn't 0
g_c = [0,0,0]                       #effective gravity of the cold layer

Q_out = [0,0,0]
Q_in = [0,0,0]
B_out = [0,0,0]
B_in = [0,0,0]
old_h = [0,0,0]
#plotting variables
ts = [i/3600 for i in range(t_start,steps*dt,dt)] # /3600 to convert to hrs for nice plot
Ths = [[],[],[]]
Tcs = [[],[],[]]
hs = [[],[],[]]

import numpy as np
from math import sqrt, cos, pi
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

def dump():
    print("h:   ", h)
    print("g_h: ", g_h)
    print("g_c: ", g_c)


for i in range(steps):
    t=i*dt
    T_ext = get_ext_temp(t)
    old_h = h.copy()
    for rm in range(2):
        #rooms
        #calculate volume flux, buoyancy in and out of hot layer, if people are in work and interface h is below the ceiling

        #########if h[rm]/H[rm]< 0.999  and (people_leave_work == False or (people_leave_work and (9*3600 <  t%(24*3600) < 17*3600)) ): - is this equiv?
        if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600) and h[j]/H[j]>0.999):
            Q_out[rm] = A_eff[rm]*sqrt((g_h[rm])*(H[rm]-h[rm]))
            B_out[rm] = g_h[rm]*Q_out[rm]

            if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
                Q_in[rm] = 0
                B_in[rm] = 0
            else:
                Q_in[rm] = c*n[rm]*B**(1/3)*h[rm]**(5/3)
                B_in[rm] = B*n[rm] + g_c[rm]*Q_in[rm] 
            
            #apply
            
            h[rm] -= dt*(Q_in[rm] - Q_out[rm])/S[rm]
            g_h[rm] = dt*((g_h[rm]*(H[rm]-old_h[rm])*S[rm] + B_in[rm] - B_out[rm])/((H[rm]-h[rm])*S[rm]))
            g_c[rm] = dt*(g_c[rm]*old_h[rm]*S[rm] + g_from_temp(get_ext_temp(t))*Q_out[rm]-g_c[rm]*Q_in[rm])/(h[rm]*S[rm])
        else:
            Q_out[rm] = 0

    #chimney
    if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600) and h[2]/H[2]>0.99):
        rm = 2 #beneficial for debug, don't touch
        Q_out[rm] = A_eff[rm]*sqrt((g_h[rm]-g_c[rm])*(H[rm]-h[rm]))
        B_out[rm] = g_h[rm]*Q_out[rm]

        if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
            Q_in[rm] = Q_out[0] + Q_out[1]
            B_in[rm] = B_out[0] + B_out[1]
        else:
            Q_in[rm] = Q_out[0] + Q_out[1] + c*n[rm]*B**(1/3)*h[rm]**(5/3)
            B_in[rm] = B_out[0] + B_out[1] + B*n[rm] + g_c[rm]*c*n[rm]*B**(1/3)*h[rm]**(5/3)
        
        #apply
        
        h[rm] -= dt*(Q_in[rm] - Q_out[rm])/S[rm]
        g_h[rm] = dt*((g_h[rm]*(H[rm]-old_h[rm])*S[rm] + B_in[rm] - B_out[rm])/((H[rm]-h[rm])*S[rm]))
        if Q_out[rm] > Q_out[0] + Q_out[1]:
            if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
                g_c[rm] = dt*(g_c[rm]*old_h[rm]*S[rm] + g_from_temp(get_ext_temp(t))*(Q_out[rm]-(Q_out[0] + Q_out[1])))/(h[rm]*S[rm])
            else:
                g_c[rm] = dt*(g_c[rm]*old_h[rm]*S[rm] + g_from_temp(get_ext_temp(t))*(Q_out[rm]-(Q_out[0] + Q_out[1]))-g_c[rm]*c*n[rm]*B**(1/3)*h[rm]**(5/3))/(h[rm]*S[rm])
        if h[rm]<0.01:
            print("hot layer in chimney got too low")
            dump()
            break

    #gather data for graphing 
    Ths[0].append(temp_from_g(g_h[0]))
    Tcs[0].append(temp_from_g(g_c[0]))
    hs[0].append(h[0])
    Ths[1].append(temp_from_g(g_h[1]))
    Tcs[1].append(temp_from_g(g_c[1]))
    hs[1].append(h[1])
    Ths[2].append(temp_from_g(g_h[2]))
    Tcs[2].append(temp_from_g(g_c[2]))
    hs[2].append(h[2])

#plot graphs

# debug graphs
if DEBUGGING:

    ext_temp_over_day = []
    for i in ts:
        ext_temp_over_day.append(get_ext_temp(i*3600))

    plt.title("temp over the day")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Temp (deg C)")
    plt.plot(ts,ext_temp_over_day)
    plt.yticks(range(0, 30, 2))
    plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # makes x-axis tickers every 24 hrs
    plt.legend()
    plt.show()


plt.title("Bottom floor layer temps.")
plt.xlabel("Time (hrs)")
plt.ylabel("Temp (deg C)")
plt.plot(ts,Ths[1], label="Hot layer")
plt.plot(ts, Tcs[1], label="Cold layer")
plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # makes x-axis tickers every 24 hrs
plt.legend()
plt.show()#temps of layers in floor 1

plt.title("Bottom floor interface height")
plt.xlabel("Time (hrs)")
plt.ylabel("Interface height (m)")
plt.plot(ts,hs[1])
plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # as above
plt.legend()
plt.show()#interface height, rm 1

plt.title("Chimney layer temps.")
plt.xlabel("Time (hrs)")
plt.ylabel("Temp (deg C)")
plt.plot(ts,Ths[2], label="Hot layer")
plt.plot(ts, Tcs[2], label="Cold layer")
plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # as above
plt.legend()
plt.show()#temps of layers in chimney

plt.title("Chimney interface height")
plt.xlabel("Time (hrs)")
plt.ylabel("Interface height (m)")
plt.plot(ts,hs[2])
plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # as above
plt.legend()
plt.show()#interface height, chimney
