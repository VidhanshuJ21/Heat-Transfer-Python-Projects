import numpy as np
import plotly.graph_objects as go

# Geometry
L = 0.018
N = 181
dx = L/(N-1)
x = np.linspace(0,L,N)

# Materials: Glass | Air | Glass
glass = (x <= 0.004) | (x >= 0.014)

k   = np.where(glass,0.80,0.025)
rho = np.where(glass,2500,1.2)
cp  = np.where(glass,800,1005)

# Boundary conditions
Ti, To = 20., -10.
hi, ho = 10., 30.

#Initial temperature
T = np.ones(N)*Ti

#Interface conductivity
kf = 2*k[:-1]*k[1:]/(k[:-1]+k[1:])

# Implicit time step
dt = 10.0

# Transient solution
for n in range(3000):

    A = np.zeros((N,N))
    b = T.copy()

    for i in range(1,N-1):

        aW = dt*kf[i-1]/(rho[i]*cp[i]*dx**2)
        aE = dt*kf[i]/(rho[i]*cp[i]*dx**2)

        A[i,i-1] = -aW
        A[i,i]   = 1+aW+aE
        A[i,i+1] = -aE

    a = dt*hi/(rho[0]*cp[0]*dx)
    A[0,0] = 1+a+dt*kf[0]/(rho[0]*cp[0]*dx**2)
    A[0,1] = -dt*kf[0]/(rho[0]*cp[0]*dx**2)
    b[0] = T[0]+a*Ti

    a = dt*ho/(rho[-1]*cp[-1]*dx)
    A[-1,-1] = 1+a+dt*kf[-1]/(rho[-1]*cp[-1]*dx**2)
    A[-1,-2] = -dt*kf[-1]/(rho[-1]*cp[-1]*dx**2)
    b[-1] = T[-1]+a*To

    Tnew = np.linalg.solve(A,b)

    if np.max(abs(Tnew-T)) < 1e-7:
        T = Tnew
        print("Converged at", n, "steps")
        break

    T = Tnew

# 2D temperature visualization
Y = np.linspace(0,30,40)
Z = np.tile(T,(len(Y),1))

fig = go.Figure(go.Heatmap(
    x=x*1000,
    y=Y,
    z=Z,
    colorscale="Turbo",
    zmin=-10,
    zmax=20,
    colorbar=dict(title="Temperature [°C]")
))

fig.add_vline(x=4,line_dash="dash",line_color="black")
fig.add_vline(x=14,line_dash="dash",line_color="black")

fig.add_annotation(x=2,y=28,text="GLASS",showarrow=False)
fig.add_annotation(x=9,y=28,text="AIR",showarrow=False)
fig.add_annotation(x=16,y=28,text="GLASS",showarrow=False)

fig.update_layout(
    title="2D Transient Conduction - Double Pane Window",
    xaxis_title="Window Thickness [mm]",
    yaxis_title="Window Height [mm]"
)

fig.show()

#Numerical results
print("\nTemperature Results")
print("-------------------")
print("T1 =",round(T[np.argmin(abs(x-0.004))],2),"°C")
print("T2 =",round(T[np.argmin(abs(x-0.004))],2),"°C")
print("T3 =",round(T[np.argmin(abs(x-0.014))],2),"°C")
print("T4 =",round(T[-1],2),"°C")