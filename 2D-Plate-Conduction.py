
import numpy as np
import plotly.graph_objects as go

# Plate and mesh
L = 10.0
N = 51
dx = L/(N-1)

# Temperature field
T = np.ones((N,N))*273.0

# Boundary conditions
T[0,:]  = 373.0       # Bottom
T[-1,:] = 273.0       # Top
T[:,0]  = 273.0       # Left
T[:,-1] = 273.0       # Right

# Iterative solution of Laplace equation
for iteration in range(10000):
    Told = T.copy()

    T[1:-1,1:-1] = 0.25*(
        Told[2:,1:-1] + Told[:-2,1:-1] +
        Told[1:-1,2:] + Told[1:-1,:-2]
    )

    if np.max(abs(T-Told)) < 1e-5:
        break

# Numerical results
print("Iterations =", iteration)
print("Maximum Temperature =", T.max(), "K")
print("Minimum Temperature =", T.min(), "K")
print("Center Temperature =", T[N//2,N//2], "K")

# Temperature contour
x = np.linspace(0,L,N)
y = np.linspace(0,L,N)

fig = go.Figure(go.Contour(
    x=x, y=y, z=T,
    colorscale="Turbo",
    contours=dict(showlabels=True),
    colorbar=dict(title="Temperature [K]")
))

fig.update_layout(
    title="2D Steady-State Heat Conduction - Steel Plate",
    xaxis_title="x [m]",
    yaxis_title="y [m]"
)

fig.show()

