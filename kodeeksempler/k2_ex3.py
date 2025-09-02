import numpy as np
import matplotlib.pyplot as plt
import math

# Konstanter
m = 0.145               # massen av gjenstanden, kg
k = 1.31*10**(-3)       # luftmotstandstallet, kg/m
g = 9.81                # tyngdeakselerasjon, m/s^2
v0 = 42.7               # startfart, m/s
y0 = 1.20               # starthøyde, m
theta_deg = 40          # vinkelen, i grader

theta = theta_deg * math.pi / 180  # konverterer vinkelen til radianer

# Konstante krefter
G = np.array([0, -m*g])     # tyngden i N


# Variable krefter, utregning av kraftsum og akselerasjon
def a(vel):                   # akselerasjonsfunksjon
    e_v = np.linalg.norm(vel)          # enhetsvektor for farten
    L = -k*np.linalg.norm(vel)**2 * e_v  # luftmotstandsvektor, N
    sum_F = G + L            # vektorsummen av kreftene, N
    aks = sum_F/m            # akselerasjonsvektoren, m/s^2
    return aks


# Startverdier for bevegelsen
r = np.array([0, y0])                         # startposisjon, m
v = np.array([v0*math.cos(theta), v0*math.sin(theta)])  # startfart, m/s
t = 0                                      # starttid, s

# Lister for lagring av verdier
r_liste = [r]
v_liste = [v]

# Simulering av bevegelsen
dt = 0.001         # tidssteg i simuleringen, s

while r[1] >= 0:   # stopper når y = 0
    v = v + a(v)*dt  # regner ut neste fartsvektor
    r = r + v*dt     # regner ut neste posisjonsvektor
    t = t + dt       # går til neste tidspunkt

    # Lagring av 2D-verdier i lister
    r_liste = np.concatenate([r_liste, [r]])
    v_liste = np.concatenate([v_liste, [v]])

# Tegning av graf
plt.plot(r_liste[:, 0], r_liste[:, 1])       # lager grafen
plt.title("Skrått kast med luftmotstand")  # tittel på grafen
plt.xlabel("$x$ / m")                      # x-akse-tittel
plt.ylabel("$y$ / m")                      # y-akse-tittel
plt.grid()                                 # legger til rutenett
plt.show()                                 # viser grafen