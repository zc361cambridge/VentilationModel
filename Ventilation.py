'''Comprehensive ventialtion simulation for 2 storey atrium building'''
#imports
from math import sqrt, cos, pi, sin
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator

#user defined variables
H = [3,3]                       #room height
H_a = 9                         #height of atrium
vents = [[1,1],[1,1]]           #vent areas, as [low,high], or equivalently: [in,out]
vents_a = [4,3]
h_v = [3,6]                     #vent height above ext. vent, in atrium
S = [60,60]                     #floor area
S_a = 36
n = [30,30]                     #num of people
w = 100                         #power input per person
T_day = 20                      #max daytime temp
T_night = 5                     #min nighttime temp
people_leave_work = True        #Do people go home outside of 9-5
peak_solar = [20000,20000]      #peak solar heating, in watts. Try 1000*S.
sunrise = 6                     #sunrise, time in 24h
sunset = 18

'''
Room numbering system:
                 _______
                 |
                 |     |
  _______________|     |
 |                     |
 |               |     |
 |       1       |     |
                 |  a  |
 |_______________|     |
 |                     |
 |               |     |
 |       0       |     |
                 |      
 |_______________|_____|
'''

#fundamental variables
c = 0.15                        #plume entrainment constant
alpha = c                       #conversion to half line plume
g_real = 9.81
rho = 1.225

#function definitions
def get_ext_temp(t):
    '''t is time in seconds since start of the test. return is Celcius'''
    return ((T_night-T_day)/2)*cos(2*pi*t/(3600*24))+(T_day+T_night)/2
def temp_from_g(g):
    '''returns the temp of air from its effective gravity. mostly for readability'''
    return T_night+(g/g_real)*(273+T_night)
def g_from_temp(T):
    '''returns effective gravity of air by temp, relative to ext air at midnight'''
    return g_real*(T-T_night)/(273+T_night)
def get_A_eff(a,b):
    '''used in setup, calculates effective vent area from real vent areas'''
    return sqrt(2)*a*b/sqrt(a**2+b**2)
def get_z(M,Q,B):
    '''takes momentum flux, volume flux, buoyancy flux, 
    returns plume vertical origin correct
    method taken from Hunt and Kaye 2001
    '''
    gamma = 5*Q**2*B/(4*alpha*M**(5/2))
    if gamma > 88:
        return 0.853*gamma**(-0.2)*5*Q/(6*alpha*M**0.5)
    L = 2**(-1.5)*alpha**(-0.5)*M**(0.75)/B**0.5
    return 1.057*L
def get_dp(d, g_c, g_h, t, h, H):
    '''get pressure drop across room due to chimney effect
    ***WARNING***: can be negative
    '''
    net = g_c*h + g_h*(H-h) - g_from_temp(get_ext_temp(t))*H #net delta p in chimney
    if d>h:
        return g_c*h + g_h*(d-h) - g_from_temp(get_ext_temp(t))*d - net/2
    return g_c*d - g_from_temp(get_ext_temp(t))*d - net/2
def sign(n):
    if n<0:
        return -1
    return 1
def dump(h, g_h, g_c):
    '''returns some basic info in event of handled exception. For debugging'''
    print("h:   ", h)
    print("g_h: ", g_h)
    print("g_c: ", g_c)
def w_to_B(watts):
    '''converts watts to buoyancy input'''
    return 0.0281*(watts/1000)
def get_sun(t, rise, fall):
    '''returns the fraction of the peak solar radiation for a set time'''
    t = (t%(3600*24))/3600
    if t<rise or t>fall:
        return 0
    return sin(pi*(t-rise)/(fall-rise))

#setup variables
B = w_to_B(w)                       #buoyancy input per person
dt = 1                              #timestep in seconds
steps = 24*3600*2
A_eff = [get_A_eff(vents[i][0],vents[i][1]) for i in range(2)]      #effective vent area
A_eff_a = get_A_eff(vents_a[0],vents_a[1])


#main function
def main():
    #simulation variables
    h = 0.99*H_a                    #chimney interface height
    old_h = 0
    g = [0,0]                       #effective gravity, rooms
    flow = [1,1,1]                  #for debugging
    #0: into hot layer, 1: h plume, 2: c plume, 3: inhale from h, 4 inhale from c
    Q = [0,0]                       #volumetric flow rate, rooms
    B_out = [0,0]                   #buoyancy flux out of room
    B_in = [0,0]
    g_h = 0                         #effective gravity of hot layer, chimney
    g_c = 0
    z = 0
    dp = [0,0]
    throughflow = 0

    #plotting variables
    ts = [] # /3600 is to convert to hrs for nice plot
    Ts = [[],[]]
    Ths = []
    Tcs = []
    flows = [[],[],[]]
    hs = []
    dps = [[],[]]
    qs = [[],[]]
    deltas = [[],[]]
    throughflows = []
    for i in range(steps):
        #variables that must be reset
        mixing = True
        half_mixing = False
        force_mixing = False
        B_in_c = 0
        B_out_c = 0
        B_in_h = 0
        B_out_h = 0
        Q_in_c = 0
        Q_out_c = 0
        Q_in_h = 0
        Q_out_h = 0

        t = i*dt
        old_h = h
        g_ext = g_from_temp(get_ext_temp(t))

        for j in range(2):
            #rooms
            dp[j] = get_dp(h_v[j], g_c, g_h, t, h, H_a)
            dps[j].append(dp[j])
            #dp = 0
            Q[j] = A_eff[j]*sqrt(abs((g[j]-g_ext)*H[j]-dp[j]/rho)) #absolute value
            B_out[j] = Q[j]*g[j]
            if (g[j]-g_ext)*H[j]-dp[j]/rho>0:
                #correct flow direction
                B_in[j] = Q[j]*g_ext

                #effect on chimney:
                Q_in_c -= Q[j]
                if h_v[j] > h:
                    #direct into hot layer
                    flow[j] = 0
                    Q_in_h += Q[j]
                    B_in_h += B_out[j]

                else:
                    #plume physics
                    if g[j] > g_c:
                        #hot plume
                        flow[j] = 1
                        if B_out[j] > 0:
                            z = get_z(Q[j]/vents[j][1],Q[j],B_out[j])
                            Q_in_h += alpha*B_out[j]**(1/3)*((h-h_v[j]+z)**(5/3)-(z)**(5/3))+Q[j]
                            B_in_h += B_in[j]+g_c*alpha*B_out[j]**(1/3)*((h-h_v[j]+z)**(5/3)-(z)**(5/3))
                            Q_out_c += alpha*B_out[j]**(1/3)*((h-h_v[j]+z)**(5/3)-(z)**(5/3))
                            B_out_c += g_c*(alpha*B_out[j]**(1/3)*((h-h_v[j]+z)**(5/3)-(z)**(5/3)))
                            mixing = False
                        else:
                            dump(h, g_h, g_c)
                            raise ValueError(f"negative B_out: {B_out[j]}")
                    else:
                        #cold plume
                        flow[j] = 2
                        B_in_c += B_out[j]
            else:
                #reverse flow
                Q_in_c += Q[j]
                if h_v[j] > h:
                    #inhaling from hot layer
                    flow[j] = 3
                    B_in[j] = Q[j]*g_h

                    #effect on chimney:
                    B_out_h += B_in[j]
                    Q_out_h += Q[j]
                else:
                    #inhaling from cold layer
                    flow[j] = 4
                    B_in[j] = Q[j]*g_c

                    #effect on chimney
                    B_out_c += B_in[j]
                    Q_out_c += Q[j]

            if not (people_leave_work and (t%(24*3600) < 9*3600 or  t%(24*3600) > 17*3600)):
                #if people are in the building, add their heat
                B_in[j] += B*n[j]

            B_in[j] += get_sun(t,sunrise,sunset)*w_to_B(peak_solar[j])
            g[j] += dt*(B_in[j]-B_out[j])/(H[j]*S[j])
            flows[j].append(flow[j]+j/10)

        #chimney
        if Q_in_h < throughflow and (not mixing) and h/H_a >=0.99:
            half_mixing = True

        if mixing and g_c >= g_h:
            force_mixing = True
            h = 0.99*H_a

        if (h/H_a <0.99 or (not mixing)) and not force_mixing:
            #displacement case
            old_g_c = g_c
            old_h = h
            throughflow = sign(g_h*(H_a-h)+g_c*h-g_ext*H_a)*A_eff_a*sqrt(abs(g_h*(H_a-h)+g_c*h-g_ext*H_a))
            if not half_mixing:
                if throughflow > 0:
                    #case 1: normal displacement flow
                    flow[2] = 1
                    Q_out_h += throughflow
                    B_out_h += Q_out_h*g_h
                    Q_in_c += abs(throughflow)
                else:
                    #case 2: top vent is backing up
                    flow[2] = 2
                    Q_in_c += abs(throughflow)
                    Q_out_c += abs(throughflow)
            else:
                if throughflow > 0:
                    #case 3: half mixing flow: top layer is very small
                    #so vent outputs hot plume input, and some cold layer
                    #distinct from mixing, where plume would mix into room
                    flow[2] = 3
                    Q_out_c += throughflow - Q_in_h
                    Q_in_c += throughflow - Q_in_h
                    Q_in_h = 0
                    Q_out_h = 0
                    B_in_h = 0
                    B_out_h = 0
                else:
                    #edge case 4
                    #very unlikely to occur at all, equivalent to case 2
                    flow[2] = 4
                    Q_in_c += abs(throughflow)
                    Q_out_c += abs(throughflow)
                    print("edge case")

            h -= dt*(Q_in_h - Q_out_h)/S_a
            g_h = (g_h*(H_a - old_h)*S_a + dt*(B_in_h - B_out_h))/((H_a - h)*S_a)
            if Q_in_c > 0:
                g_c = (g_c*old_h*S_a + dt*(g_from_temp(get_ext_temp(t))*Q_in_c + B_in_c - B_out_c))/(h*S_a)
            else:
                g_c = (g_c*old_h*S_a + dt*(B_in_c - B_out_c - abs(Q_in_c)*g_c))/(h*S_a)

            if half_mixing:
                g_h = g_c
            deltas[0].append(t)
            deltas[1].append(g_c-old_g_c)

            mixing = False

            #error handling
            if old_g_c-g_c >1e-13 and g_ext>old_g_c and flow[0] != 2 and flow[1] != 2:
                print("Second law violation: cold layer is cooling on its own!")
                dump(h, g_h, g_c)
                break
            if h<0:
                print("hot layer flowing out through lower vent, assumptions no longer hold!")
                dump(h, g_h, g_c)
                break
            if h>H_a:
                print(f"somehow, your hot layer has negative size at time {t}")
                dump(h, g_h, g_c)
                break

        if mixing:
            #mixing case
            g_h = (g_c*h + g_c*(H_a-h))/H_a #normalise temp
            g_c = g_h
            old_g_c = g_c
            flow[2] = 0
            h = 0.99*H_a
            # there has to be a virtual hot/cold layer, so that
            # when mixing stops, we can go back to displacement

            throughflow = sign((g_c-g_ext)*H_a)*A_eff_a*sqrt(abs(g_c-g_ext)*H_a)
            Q_out_c = abs(throughflow) - Q_in_h
            B_out_c = Q_out_c*g_c
            B_in_c += Q_out_c*g_ext + B_in_h

            g_c += dt*(B_in_c-B_out_c)/(H_a*S_a)
            g_h = g_c
            deltas[0].append(t)
            deltas[1].append(g_c - old_g_c)

        #gather data for graphing
        ts.append(t/3600)
        Ts[0].append(temp_from_g(g[0]))
        Ts[1].append(temp_from_g(g[1]))
        Ths.append(temp_from_g(g_h))
        Tcs.append(temp_from_g(g_c))
        hs.append(h)
        qs[0].append(Q[0])
        qs[1].append(Q[1])
        flows[2].append(flow[2]+0.2)
        throughflows.append(throughflow)

    #plot graphs

    plt.title("Room temps.")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Temp (deg C)")
    plt.plot(ts,Ts[0], label="Room 0")
    plt.plot(ts,Ts[1], label="Room 1")
    plt.plot(ts, [get_ext_temp(i*3600) for i in ts], label="External")
    plt.gca().xaxis.set_major_locator(MultipleLocator(24)) # makes x-axis tickers every 24 hrs
    plt.legend()
    plt.show()

    plt.title("Chimney layer temps.")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Temp (deg C)")
    plt.plot(ts, Ths, label="Hot layer")
    plt.plot(ts, Tcs, label="Cold layer")
    plt.plot(ts, [get_ext_temp(i*3600) for i in ts], label="External")
    plt.gca().xaxis.set_major_locator(MultipleLocator(24))
    plt.legend()
    plt.show()

    plt.title("Throughflow")
    plt.xlabel("Time (hrs)")
    plt.ylabel("Throughflow (m^3/s)")
    plt.plot(ts,throughflows, label="chimney")
    plt.plot(ts,qs[1], label="room 1")
    plt.plot(ts,qs[0], label="room 0")
    plt.gca().xaxis.set_major_locator(MultipleLocator(24))
    plt.legend()
    plt.show()

    plt.title("Delta p")
    plt.xlabel("Time (hrs)")
    plt.ylabel("dp (pa)")
    plt.plot(ts,dps[1], label="room 1")
    plt.plot(ts,dps[0], label="room 0")
    plt.gca().xaxis.set_major_locator(MultipleLocator(24))
    plt.legend()
    plt.show()

#using best practices: because I know someone is going to read this code
if __name__ == "__main__":
    main()
