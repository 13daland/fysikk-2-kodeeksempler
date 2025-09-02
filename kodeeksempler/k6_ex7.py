import matplotlib.pyplot as plt
import numpy as np

# Konstanter
q = 1.0                     # elektrisk ladning, C
m = 0.010                   # masse, kg
E = np.array([0, -0.0040, 0])  # elektrisk felt, N/C
B = np.array([0, 0, -0.010])   # magnetisk felt, T

# Konstante krefter
F_e = q*E  # elektrisk kraft, N


# Variable krefter, utregning av kraftsum og akselerasjon
def a(v):                       # akselerasjonsfunksjon
    F_m = np.cross(q*v, B)      # magnetisk kraft, N
    sum_F = F_e + F_m           # kraftsum, N
    aks = sum_F / m             # akselerasjon, m/s^2
    return aks                  # returnerer akselerasjonen


# Startverdier for bevegelsen
r = np.array([0, 0, 0])         # startposisjon, m
v = np.array([0.24, 0, 0])      # startfart, m/s
t = 0                           # starttid, s

# Lister for lagring av verdier
x_verdier = [r[0]]
y_verdier = [r[1]]

# Simulering av bevegelsen
dt = 0.001  # tidssteg i simuleringen, s

while t < 14:  # stopper når t = 14
    v = v + a(v)*dt         # regner ut ny fart
    r = r + v*dt            # regner ut ny posisjon
    t = t + dt              # går til neste tidspunkt
    x_verdier.append(r[0])  # legger x inn i listen
    y_verdier.append(r[1])  # legger y inn i listen

# Tegning av graf
plt.plot(x_verdier, y_verdier)  # lager grafen
plt.title("Kryssede felt")      # tittel på grafen
plt.xlabel("$x$ / m")           # navn på x-aksen
plt.ylabel("$y$ / m")           # navn på y-aksen
plt.grid()                      # lager rutenett
plt.show()                      # viser grafen
