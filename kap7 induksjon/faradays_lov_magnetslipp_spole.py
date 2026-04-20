import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --------------------------
# Parametere du kan endre
# --------------------------
N = 10                 # Antall vindinger
B_surface = 0.5        # T (magnetens feltstyrke ved "overflaten" langs aksen)
m = 0.05               # kg (magnetens masse)
h0 = 0.25              # m (slipphøyde over spolens senter, z(0) = h0)
R_coil = 0.03          # m (spolens radius)
R_internal = 2.0       # ohm (indre resistans i spolen)

g = 9.81               # m/s^2
mu0 = 4 * np.pi * 1e-7 # vakuumpermeabilitet
z0 = 0.02              # m, effektiv "halv-lengde" for magneten langs aksen

# --------------------------
# Graf-selektor
# --------------------------
# Skriv ett av: "emf", "v", "a", "h"
plot_type = "emf"   # <-- endre denne til "v", "a" eller "h" ved behov

# --------------------------
# Feltmodell for magneten
# --------------------------
# Juster magnetisk dipolmoment slik at B(z0) = B_surface langs aksen.
# B(z) ≈ (mu0 / (2π)) * (m_mag / z^3) (ideelt), men vi bruker (z^2+z0^2)^(3/2)
# for å unngå singularitet.
# m_mag bestemmes slik at B(z0) = B_surface i den glatte modellen.

# B(z) = C / (z^2 + z0^2)^(3/2), med C = (mu0/(2π))*m_mag
# B(z0) = C / (2*z0^2)^(3/2) = B_surface
# => C = B_surface * (2*z0**2)**(3/2)
# => m_mag = C * (2π)/mu0
C_B = B_surface * (2 * z0**2)**1.5
m_mag = C_B * (2 * np.pi) / mu0

def B_axial(z):
    """
    Magnetfelt langs spolens akse.
    z: avstand fra magnetens sentrum til spolens senter (positiv oppover).
    """
    C = (mu0 / (2 * np.pi)) * m_mag
    return C / (z**2 + z0**2)**1.5

def dB_dz(z):
    """
    Derivert av B_axial(z) mht z.
    B(z) = C * (z^2 + z0^2)^(-3/2)
    dB/dz = -3 C z (z^2 + z0^2)^(-5/2)
    """
    C = (mu0 / (2 * np.pi)) * m_mag
    return -3 * C * z * (z**2 + z0**2)**(-2.5)


# Spolens tverrsnittsareal
A_coil = np.pi * R_coil**2

# --------------------------
# ODE-system: tilstand [z, v]
# --------------------------
def dynamics(t, y):
    z, v = y

    # Tyngdekraft (nedover = negativ)
    F_g = -m * g

    # Induktiv bremsing
    # EMF = -N * A * dB/dz * v
    # F_ind = - (N^2 * A^2 * (dB/dz)^2 / R_internal) * v
    dBdz = dB_dz(z)
    F_ind = - (N**2 * A_coil**2 * dBdz**2 / R_internal) * v

    # Total kraft
    F_tot = F_g + F_ind

    dzdt = v
    dvdt = F_tot / m
    return [dzdt, dvdt]


# --------------------------
# Integrasjon
# --------------------------
# Starttilstand: z(0) = h0 (magnet over spolen), v(0) = 0 (slippes fra ro)
y0 = [h0, 0.0]

# Estimér simuleringstid: fritt fall til ca. spolesenter + litt ekstra
t_free_fall = np.sqrt(2 * h0 / g)
t_span = (0, 1.5 * t_free_fall)

# For bedre oppløsning: spesifiser tidspunkter vi vil ha ut
t_eval = np.linspace(t_span[0], t_span[1], 5000)

sol = solve_ivp(dynamics, t_span, y0, t_eval=t_eval, rtol=1e-7, atol=1e-9)

t = sol.t
z = sol.y[0]  # posisjon
v = sol.y[1]  # hastighet

# --------------------------
# Beregn EMF og akselerasjon
# --------------------------
dBdz_arr = dB_dz(z)
B_arr = B_axial(z)
Phi = B_arr * A_coil

# EMF = -N * dPhi/dt = -N * A * dB/dz * v
emf = -N * A_coil * dBdz_arr * v

# Akselerasjon = dv/dt (kan også fås direkte fra kraftuttrykk)
# Her bruker vi kraftformelen for konsistens:
F_g_arr = -m * g
F_ind_arr = - (N**2 * A_coil**2 * dBdz_arr**2 / R_internal) * v
a = (F_g_arr + F_ind_arr) / m

# --------------------------
# Plot etter valg
# --------------------------
plt.figure(figsize=(8, 4))

if plot_type == "emf":
    plt.plot(t, emf)
    plt.ylabel("emf (V)")
    plt.title("Indusert elektromotorisk spenning")
elif plot_type == "v":
    plt.plot(t, v)
    plt.ylabel("hastighet v (m/s)")
    plt.title("Hastighet som funksjon av tid")
elif plot_type == "a":
    plt.plot(t, a)
    plt.ylabel("akselerasjon a (m/s²)")
    plt.title("Akselerasjon som funksjon av tid")
elif plot_type == "h":
    plt.plot(t, z)
    plt.ylabel("posisjon z (m)")
    plt.title("Posisjon som funksjon av tid (z over spolesenter)")
else:
    # Fallback: vis emf hvis brukeren har skrevet noe annet
    plt.plot(t, emf)
    plt.ylabel("emf (V)")
    plt.title("Indusert elektromotorisk spenning (default)")

plt.xlabel("t (s)")
plt.grid(True)
plt.tight_layout()
plt.show()
