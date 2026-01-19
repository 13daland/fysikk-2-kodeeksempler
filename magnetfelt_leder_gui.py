import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors
import matplotlib.cm as cm

mu0 = 4*np.pi*1e-7

# ============================
# Fysikk
# ============================

class Wire:
    def __init__(self, pos, direction, I):
        self.pos = np.array(pos, dtype=float)
        d = np.array(direction, dtype=float)
        self.dir = d / np.linalg.norm(d)
        self.I = I

def B_from_wire(r, wire):
    R = r - wire.pos
    cross = np.cross(wire.dir, R)
    d2 = np.dot(cross, cross)
    if d2 < 1e-12:
        return np.zeros(3)
    return mu0 * wire.I / (2*np.pi) * cross / d2

def B_total(r, wires):
    B = np.zeros(3)
    for w in wires:
        B += B_from_wire(r, w)
    return B

# ============================
# App
# ============================

class App:
    def __init__(self, root):
        self.root = root
        root.title("3D magnetfelt fra rette strømledere")

        # --------------------
        # Modell
        # --------------------

        self.wires = [
            Wire([0,0,0], [0,0,1], 10.0)
        ]

        self.q = 1e-6
        self.m = 1e-6
        self.r = np.array([0.05, 0.0, 0.0])
        self.v = np.array([0.0, 2.0, 0.0])

        self.running = False
        self.dt = 0.0002
        self.speed = 1.0

        # Visning
        self.lim = 0.20
        self.min_dist = 0.02
        self.B_scale = 5e-5

        # Tettere rutenett
        self.n = 8   # 8×8×8 → merkbart tettere, men fortsatt håndterbart

        self.build_gui()
        self.init_grid()
        self.update_plot()

    # ============================
    # GUI
    # ============================

    def build_gui(self):
        main = ttk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(main, width=340)
        left.pack(side=tk.LEFT, fill=tk.Y)
        left.pack_propagate(False)

        right = ttk.Frame(main)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Plot
        self.fig = Figure(figsize=(8,8))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # ---------- Ledere ----------
        ttk.Label(left, text="Ny leder", font=("Arial", 11, "bold")).pack(pady=5)

        self.I_entry  = self.make_entry(left, "I [A]", "10")
        self.x_entry  = self.make_entry(left, "x [m]", "0")
        self.y_entry  = self.make_entry(left, "y [m]", "0")
        self.z_entry  = self.make_entry(left, "z [m]", "0")
        self.dx_entry = self.make_entry(left, "dx", "0")
        self.dy_entry = self.make_entry(left, "dy", "0")
        self.dz_entry = self.make_entry(left, "dz", "1")

        ttk.Button(left, text="Legg til leder", command=self.add_wire).pack(pady=4)

        ttk.Separator(left).pack(fill=tk.X, pady=8)

        # ---------- Testladning ----------
        ttk.Label(left, text="Testladning", font=("Arial", 11, "bold")).pack()

        self.q_entry  = self.make_entry(left, "q [C]", "1e-6")
        self.m_entry  = self.make_entry(left, "m [kg]", "1e-6")

        self.rx_entry = self.make_entry(left, "r0 x [m]", "0.05")
        self.ry_entry = self.make_entry(left, "r0 y [m]", "0.0")
        self.rz_entry = self.make_entry(left, "r0 z [m]", "0.0")

        self.vx_entry = self.make_entry(left, "v0 x [m/s]", "0.0")
        self.vy_entry = self.make_entry(left, "v0 y [m/s]", "2.0")
        self.vz_entry = self.make_entry(left, "v0 z [m/s]", "0.0")

        ttk.Button(left, text="Reset partikkel", command=self.reset_particle).pack(pady=5)

        ttk.Separator(left).pack(fill=tk.X, pady=8)

        ttk.Button(left, text="Play / Pause", command=self.toggle).pack()
        ttk.Label(left, text="Hastighet").pack()
        self.speed_scale = ttk.Scale(left, from_=0.1, to=5,
                                     orient=tk.HORIZONTAL,
                                     command=self.set_speed)
        self.speed_scale.set(1.0)
        self.speed_scale.pack(fill=tk.X, padx=10)

    def make_entry(self, parent, label, default):
        f = ttk.Frame(parent)
        f.pack(anchor="w", padx=10)
        ttk.Label(f, text=label, width=12).pack(side=tk.LEFT)
        e = ttk.Entry(f, width=12)
        e.insert(0, default)
        e.pack(side=tk.LEFT)
        return e

    # ============================
    # Handling
    # ============================

    def add_wire(self):
        I  = float(self.I_entry.get())
        x  = float(self.x_entry.get())
        y  = float(self.y_entry.get())
        z  = float(self.z_entry.get())
        dx = float(self.dx_entry.get())
        dy = float(self.dy_entry.get())
        dz = float(self.dz_entry.get())

        self.wires.append(Wire([x,y,z], [dx,dy,dz], I))
        self.update_plot()

    def reset_particle(self):
        self.q = float(self.q_entry.get())
        self.m = float(self.m_entry.get())
        self.r = np.array([
            float(self.rx_entry.get()),
            float(self.ry_entry.get()),
            float(self.rz_entry.get())
        ])
        self.v = np.array([
            float(self.vx_entry.get()),
            float(self.vy_entry.get()),
            float(self.vz_entry.get())
        ])
        self.update_plot()

    def toggle(self):
        self.running = not self.running
        if self.running:
            self.step()

    def set_speed(self, val):
        self.speed = float(val)

    def step(self):
        if not self.running:
            return

        B = B_total(self.r, self.wires)
        a = self.q * np.cross(self.v, B) / self.m

        dt = self.dt * self.speed
        self.v += a * dt
        self.r += self.v * dt

        self.update_plot()
        self.root.after(30, self.step)

    # ============================
    # Plot
    # ============================

    def init_grid(self):
        xs = np.linspace(-self.lim, self.lim, self.n)
        ys = np.linspace(-self.lim, self.lim, self.n)
        zs = np.linspace(-self.lim, self.lim, self.n)
        self.X, self.Y, self.Z = np.meshgrid(xs, ys, zs)
        self.colorbar = None

    def update_plot(self):
        self.ax.clear()

        U = np.zeros_like(self.X)
        V = np.zeros_like(self.Y)
        W = np.zeros_like(self.Z)
        M = np.zeros_like(self.X)

        for idx, _ in np.ndenumerate(self.X):
            r = np.array([self.X[idx], self.Y[idx], self.Z[idx]])

            # Ikke tegn nær leder
            if any(np.linalg.norm(np.cross(w.dir, r - w.pos)) < self.min_dist
                   for w in self.wires):
                continue

            B = B_total(r, self.wires)
            U[idx], V[idx], W[idx] = B
            M[idx] = np.linalg.norm(B)

        # Fast skalering
        U2 = U / self.B_scale
        V2 = V / self.B_scale
        W2 = W / self.B_scale

        L = np.sqrt(U2**2 + V2**2 + W2**2)
        L = np.clip(L, 0, 2.0)
        mask = L > 0
        U2[mask] *= L[mask] / np.sqrt(U2[mask]**2 + V2[mask]**2 + W2[mask]**2)
        V2[mask] *= L[mask] / np.sqrt(U2[mask]**2 + V2[mask]**2 + W2[mask]**2)
        W2[mask] *= L[mask] / np.sqrt(U2[mask]**2 + V2[mask]**2 + W2[mask]**2)

        norm = mcolors.Normalize(vmin=0, vmax=max(np.max(M), 1e-12))
        colors = cm.viridis(norm(M.flatten()))

        self.ax.quiver(self.X, self.Y, self.Z,
                       U2, V2, W2,
                       colors=colors, length=0.05, normalize=False)

        # Colorbar kun én gang
        if self.colorbar is None:
            sm = cm.ScalarMappable(norm=norm, cmap="viridis")
            sm.set_array([])
            self.colorbar = self.fig.colorbar(sm, ax=self.ax, shrink=0.6, pad=0.1)
            self.colorbar.set_label("|B| [T]")
        else:
            self.colorbar.mappable.set_norm(norm)

        # Ledere
        t = np.linspace(-self.lim, self.lim, 2)
        for w in self.wires:
            pts = w.pos[None,:] + t[:,None]*w.dir[None,:]
            self.ax.plot(pts[:,0], pts[:,1], pts[:,2], "r", linewidth=3)

        # Partikkel
        self.ax.scatter(self.r[0], self.r[1], self.r[2], color="blue", s=40)

        self.ax.set_xlim(-self.lim, self.lim)
        self.ax.set_ylim(-self.lim, self.lim)
        self.ax.set_zlim(-self.lim, self.lim)

        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")
        self.ax.set_zlabel("z [m]")
        self.ax.set_title("3D magnetfelt (farge = |B|)")

        self.canvas.draw()

# ============================

root = tk.Tk()
app = App(root)
root.mainloop()
