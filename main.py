import numpy as np
from matplotlib import pyplot as plt

left = 0.0
right = 1.0
t0 = 0.0
max_time = 3.0
dx = 0.01
dt = 0.0001
u0 = 0.0
rho0 = 1.0
p0 = 50.0
x0 = 0.5
r0 = 0.01
gamma = 1.7
V0 = 2.0
alpha = 1
a = 2.3

cells = (int)((right - left) / dx)

rho = np.full(cells, rho0)
u = np.full(cells, u0)
p = np.zeros(cells)
E = np.zeros(cells)
R = np.zeros(cells)  # Euler
r = np.zeros(cells)  # Lagrange
tmp = np.zeros(cells)
V = np.zeros(cells)
U = np.zeros(cells)

# Gaussian profile.
for i in range(cells):
    p[i] = p0 * np.exp(-1.0 * np.power((i * dx - x0) / r0, 2))

for i in range(cells):
    E[i] = p[i] / (rho[i] * (gamma - 1))
    R[i] = dx * i
    r[i] = dx * i
    V[i] = 1.0 / rho[i]
    U[i] = 1.0 / rho[i]

def k(j):
    du = u[j + 1] - u[j]
    if du >= 0:
        return 0.0
    else:
        return 2.0 * a * a * du * du / (V[j] + U[j])

u[0] = 0.
u[cells - 1] = 0.

p[0] = 0
p[cells - 1] = 0

plt.ion()
fig, (U1, P1, Ro1, R1) = plt.subplots(4)
fig.set_size_inches(11, 12)
fig.suptitle(f"alpha = {alpha}, gamma = {gamma}, Start Vel = {V0}, Start pressure = {p0}, Time = %.2f"%(t0))
#ax.set_title("FPS: %.2f, particle_num: %d, cur_time = %.3f"% ((1.0 / (time.time() - t1 + 1e-6)), sim.get_particle_num(), sim.get_cur_time()))

U1.set_xlim(0, 1)
U1.set_title("Скорость")
lineU, = U1.plot([], [], lw=1)

R1.set_xlim(0, 1)
R1.set_title("Эйлер от Лагранжа")
lineR, = R1.plot([], [], lw=1)

Ro1.set_xlim(0, 0.9)
Ro1.set_title("Плотность")
lineRo, = Ro1.plot([], [], lw=1)

P1.set_xlim(0, 1)
P1.set_title("Давление")
lineP, = P1.plot([], [], lw=1)

while (t0 < max_time):

    for i in range(1, cells - 2):
        u[i] = (u[i] - dt / dx * V0 *
                ((p[i] - p[i - 1]) + (k(i) - k(i - 1))))


    for i in range(cells - 1):
        R[i] = R[i] + u[i] * dt

    U = np.copy(V)
    for i in range(cells - 2):
        V[i] = V0 * (R[i + 1] - R[i]) / (r[i + 1] - r[i])
        rho[i] = 1.0 / V[i]

    V[0] = V[1]
    V[cells - 1] = V[cells - 2]

    for i in range(cells - 2):
        E[i] = E[i] - (p[i] + k(i)) * (V[i] - U[i])

    E[0] = E[1]
    E[cells - 1] = E[cells - 2]

    for i in range(cells - 2):
        p[i] = E[i] * rho[i] * (gamma - 1)

    U1.set_ylim(u.min(), u.max())
    lineU.set_data(r, u)

    R1.set_ylim(R[1:-1].min(), R[1:-1].max())
    lineR.set_data(r, R)

    P1.set_ylim(p.min(), p.max())
    lineP.set_data(r, p)

    Ro1.set_ylim(rho.min(), rho.max())
    lineRo.set_data(r, rho)

    fig.canvas.draw()
    fig.canvas.flush_events()

    t0 += dt
