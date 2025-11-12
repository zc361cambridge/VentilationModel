import numpy as np
from math import sqrt, cos, pi
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

DEBUGGING = True  # set true to enable some debug graphs and console outputs/dumps 


#user defined variables

H = [3,3,10]                    #room height
A_eff = [0.5,0.5,0.5]           #effective vent area
S = [20,20,30]                  #floor area
num_people = [1,1,0]            #people
w = 100                         #power input per person
T_day = 16                      #max daytime temp (deg C)
T_night = 10                    #min nighttime temp (deg C)
people_leave_work = True        #if true it will simulate people only being in the office 9-5
varying_ext_temp = True         #simulates day/night cycle with cosine between day/night temps., else avg. is used constantly

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
B = 0.0281*(w/1000)             #buoyancy input per person = 0.00281 for 100W per person

def get_ext_temp(t):
    if varying_ext_temp:
        return ((T_night-T_day)/2)*cos(2*pi*t/(3600*24))+(T_day+T_night)/2
    else:
        return (T_night+T_day)/2
        
def temp_from_g(g, t):
    return get_ext_temp(t)+(g/g_real)*(273+get_ext_temp(t))
def g_from_temp(T, t):
    return g_real*(T-get_ext_temp(t))/(273+get_ext_temp(t))

#setup variables
h = [0.999999*H[0],0.999999*H[1],0.999999*H[2]]           #interface height, set at roof but needs a minimum amount of hot layer
g_real = 9.81                                       #gravity on earth
dt = 1                                             #timestep in seconds
steps = 24*3600*7                                   #number of steps
t_start = 0                                         #start time, in seconds past midnight
g_h = [g_from_temp(get_ext_temp(t_start), t_start), 
       g_from_temp(get_ext_temp(t_start), t_start),
       g_from_temp(get_ext_temp(t_start), t_start)] #effective gravity of the hot layer (at t_0 it is 0, equiv to ext temp.)
g_c = [0,0,0]                                       #effective gravity of the cold layer (always 0 in displacement flow)

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

def dump():
    print("h:   ", h)
    print("g_h: ", g_h)
    print("g_c: ", g_c)


for i in range(steps):
    t=i*dt
    old_h = h.copy()
    for j in range(2):
        #rooms
        #calculate volume flux, buoyancy in and out of hot layer, if 
        #########if h[rm]/H[rm]< 0.999  and (people_leave_work == False or (people_leave_work and (9*3600 <  t%(24*3600) < 17*3600)) ): - is this equiv?
        
        if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600) and h[j]/H[j]>0.999):
            Q_out[j] = A_eff[j]*sqrt((g_h[j]-g_c[j])*(H[j]-h[j]))
            B_out[j] = g_h[j]*Q_out[j]

        Q_out[j] = A_eff[j]*sqrt((g_h[j])*(H[j]-h[j]))

        if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
            n = [0,0,0] # make sure there are no plumes if there are no people
        else:
            n = num_people
        
        Q_in[j] = c*n[j]*B**(1/3)*h[j]**(5/3)
        
        #apply
        
        h[j] -= dt*(Q_in[j] - Q_out[j])/S[j]

        # software needs a minimum thickness of hot layer
        # if it is at a minimum assume the temp of it (g_h) is 0 (i.e. assume it's not there)
        if h[j] > 0.999999*H[j]:
            h[j] = 0.999999*H[j]
            g_h[j] = 0

        if n != [0,0,0]:                            # only applies when there is a plume, otherwise g_h won't change
            g_h[j] = B**(2/3) / (c * h[j]**(5/3))   # i think technically a steady-state equation but works somewhat in this situation 
        #g_h[j] = dt*((g_h[j]*(H[j]-old_h[j])*S[j] + B_in[j] - B_out[j])/((H[j]-h[j])*S[j]))


    #gather data for graphing 
    Ths[0].append(temp_from_g(g_h[0],t))
    Tcs[0].append(temp_from_g(g_c[0],t))
    hs[0].append(h[0])
    Ths[1].append(temp_from_g(g_h[1],t))
    Tcs[1].append(temp_from_g(g_c[1],t))
    hs[1].append(h[1])
    Ths[2].append(temp_from_g(g_h[2],t))
    Tcs[2].append(temp_from_g(g_c[2],t))
    hs[2].append(h[2])
'''
    #chimney
    if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600) and h[2]/H[2]>0.99):
        j = 2 #beneficial for debug, don't touch
        Q_out[j] = A_eff[j]*sqrt((g_h[j]-g_c[j])*(H[j]-h[j]))
        B_out[j] = g_h[j]*Q_out[j]

        if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
            Q_in[j] = Q_out[0] + Q_out[1]
            B_in[j] = B_out[0] + B_out[1]
        else:
            Q_in[j] = Q_out[0] + Q_out[1] + c*n[j]*B**(1/3)*h[j]**(5/3)
            B_in[j] = B_out[0] + B_out[1] + B*n[j] + g_c[j]*c*n[j]*B**(1/3)*h[j]**(5/3)
        
        #apply
        
        h[j] -= dt*(Q_in[j] - Q_out[j])/S[j]
        g_h[j] = dt*((g_h[j]*(H[j]-old_h[j])*S[j] + B_in[j] - B_out[j])/((H[j]-h[j])*S[j]))
        if Q_out[j] > Q_out[0] + Q_out[1]:
            if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
                g_c[j] = dt*(g_c[j]*old_h[j]*S[j] + g_from_temp(get_ext_temp(t))*(Q_out[j]-(Q_out[0] + Q_out[1])))/(h[j]*S[j])
            else:
                g_c[j] = dt*(g_c[j]*old_h[j]*S[j] + g_from_temp(get_ext_temp(t))*(Q_out[j]-(Q_out[0] + Q_out[1]))-g_c[j]*c*n[j]*B**(1/3)*h[j]**(5/3))/(h[j]*S[j])
        if h[j]<0.01:
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
'''

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
plt.show()#interface height, flr 1

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
'''
