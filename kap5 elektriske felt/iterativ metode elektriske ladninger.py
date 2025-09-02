import matplotlib.pyplot as plt
from e_felt_av_n_ladninger import e_felt


charges = (         # Ladning til faste ladninger, C
    -2e-8,
    -2e-8,
    2e-8,
    #2e-8
)
charges_positions = (   # Posisjoner til faste ladninger, (x, y) m
    (-0.05, 0),
    ( 0.05, 0),
    (0, -0.05),
    #(0,  0.05)
)

p_mass = .1       # Massen til ladningen som skal beveges, kg
p_charge = 4e-7     # Ladningen til fri ladning, C
p_initial_pos = [-.01, .1]   # Startposisjon til fri ladning, (x, y) m
p_initial_vel = [0.15, -0.13]   # startfast til fri ladning (x, y) m/s

dt = 1e-4       # Tidssteg, s
t_max = 40      # Når simuleringen skal termineres, s

t_vals = [0]

t = 0
x, y = p_initial_pos
dx, dy = p_initial_vel

x_vals = [x]        # Lagre x-posisjoner for fri partikkel, m
y_vals = [y]        # Lagre y-posisjoner for fri partikkel, m
dx_vals = [dx]      # Lagre x-hastigheter for fri partikkel, m/s
dy_vals = [dy]      # Lagre y-hastigheter for fri partikkel, m/s
v_vals = [(dx**2+dy**2)**.5]

def vec_size(k1, k2):
    return (k1**2 + k2**2)**.5


while t < t_max:
    ex, ey = e_felt([x, y], charges, charges_positions)

    dx += ex * dt * p_charge / p_mass
    dy += ey * dt * p_charge / p_mass

    dx_vals.append(dx)
    dy_vals.append(dy)
    v_vals.append(vec_size(dx, dy))

    x += dx*dt
    y += dy*dt

    x_vals.append(x)
    y_vals.append(y)

    t_vals.append(t)

    t += dt


plt.plot(x_vals, y_vals, 'k-', label="fri ladning")
plt.plot([x_vals[0], x_vals[-1]], [y_vals[0], y_vals[-1]], 'ko')

for i in range(len(charges)):
    if charges[i] > 0:
        plt.plot(charges_positions[i][0], charges_positions[i][1], 'rP')
    else:
        # print(charges_positions[i], "negativ")
        plt.plot(charges_positions[i][0], charges_positions[i][1], 'b*')

plt.show()
plt.plot(t_vals, v_vals)
plt.show()



print(x, y)
