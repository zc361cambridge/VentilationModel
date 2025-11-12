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
h = [0.5*H[0],0.5*H[1],0.9*H[2]]               #interface height
g_real = 9.81                   #gravity on earth
dt = 1                         #timestep in seconds
steps = 24*3600*7                 #number of steps
t_start = 0                     #start time, in seconds past midnight
g_h = [g_from_temp(T_night+1),g_from_temp(T_night+1),g_from_temp(T_night+1)]             #effective gravity of the hot layer
g_c = [0,0,0]                   #effective gravity of the cold layer

Q_out = [0,0,0]
Q_in = [0,0,0]
B_out = [0,0,0]
B_in = [0,0,0]
old_h = [0,0,0]
#plotting variables
ts = [i/3600 for i in range(t_start,steps*dt,dt)]
Ths = [[],[],[]]
Tcs = [[],[],[]]
hs = [[],[],[]]

from math import sqrt, cos, pi
from matplotlib import pyplot as plt

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
        if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600) and h[j]/H[j]>0.999):
            Q_out[j] = A_eff[j]*sqrt((g_h[j]-g_c[j])*(H[j]-h[j]))
            B_out[j] = g_h[j]*Q_out[j]

            if people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600):
                Q_in[j] = 0
                B_in[j] = 0
            else:
                Q_in[j] = c*n[j]*B**(1/3)*h[j]**(5/3)
                B_in[j] = B*n[j] + g_c[j]*Q_in[j] 
            
            #apply
            
            h[j] -= dt*(Q_in[j] - Q_out[j])/S[j]
            g_h[j] = dt*((g_h[j]*(H[j]-old_h[j])*S[j] + B_in[j] - B_out[j])/((H[j]-h[j])*S[j]))
            g_c[j] = dt*(g_c[j]*old_h[j]*S[j] + g_from_temp(get_ext_temp(t))*Q_out[j]-g_c[j]*Q_in[j])/(h[j]*S[j])
        else:
            Q_out[j] = 0

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

#plot graphs
plt.plot(ts,Ths[1])
plt.plot(ts, Tcs[1])
plt.show()#temps of layers in floor 1

plt.plot(ts,hs[1])
plt.show()#interface height, flr 1

plt.plot(ts,Ths[2])
plt.plot(ts, Tcs[2])
plt.show()#temps of layers in chimney

plt.plot(ts,hs[2])
plt.show()#interface height, chimney
