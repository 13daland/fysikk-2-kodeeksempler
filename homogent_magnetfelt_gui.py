import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D

# ============================
# App
# ============================

class App:
    def __init__(self, root):
        self.root = root
        root.title("Bevegelse av ladning i homogent magnetfelt")

        # --------------------
        # Fysiske parametre
        # --------------------

        # Magnetfelt (T)
        self.B = np.array([0.0, 0.0, 0.01])

        # Partikkel
        self.q = 1e-6      # C
        self.m = 1e-6      # kg
        self.r = np.array([0.0, 0.0, 0.0])
        self.v = np.array([1.0, 0.5, 0.0])

        # --------------------
        # Simulasjon
        # --------------------

        self.running = False

        self.dt = 1e-4     # FAST fysisk tidssteg
        self.steps_per_frame = 1  # styres av slider

        # Visning
        self.lim = 0.2

        # Spor
        self.traj = []

        self.build_gui()
        self.init_plot()
        self.update_plot()

    # ============================
    # GUI
    # ============================

    def build_gui(self):
        main = ttk.Frame(self.root)
        main.pack(fill=tk.BOTH, expand=True)

        left = ttk.Frame(main, width=320)
        left.pack(side=tk.LEFT, fill=tk.Y)
        left.pack_propagate(False)

        right = ttk.Frame(main)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # Plot
        self.fig = Figure(figsize=(8,8))
        self.ax = self.fig.add_subplot(111, projection="3d")
        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # -------- Magnetfelt --------
        ttk.Label(left, text="Magnetfelt B [T]", font=("Arial", 11, "bold")).pack(pady=5)

        self.Bx_entry = self.make_entry(left, "Bx [T]", "0.0")
        self.By_entry = self.make_entry(left, "By [T]", "0.0")
        self.Bz_entry = self.make_entry(left, "Bz [T]", "0.01")

        ttk.Button(left, text="Oppdater B", command=self.update_B).pack(pady=4)

        ttk.Separator(left).pack(fill=tk.X, pady=8)

        # -------- Partikkel --------
        ttk.Label(left, text="Partikkel", font=("Arial", 11, "bold")).pack()

        self.q_entry = self.make_entry(left, "q [C]", "1e-6")
        self.m_entry = self.make_entry(left, "m [kg]", "1e-6")

        self.rx_entry = self.make_entry(left, "r0 x [m]", "0.0")
        self.ry_entry = self.make_entry(left, "r0 y [m]", "0.0")
        self.rz_entry = self.make_entry(left, "r0 z [m]", "0.0")

        self.vx_entry = self.make_entry(left, "v0 x [m/s]", "1.0")
        self.vy_entry = self.make_entry(left, "v0 y [m/s]", "0.5")
        self.vz_entry = self.make_entry(left, "v0 z [m/s]", "0.0")

        ttk.Button(left, text="Reset partikkel", command=self.reset_particle).pack(pady=5)

        ttk.Separator(left).pack(fill=tk.X, pady=8)

        ttk.Button(left, text="Play / Pause", command=self.toggle).pack()

        ttk.Label(left, text="Avspillingshastighet").pack()

        # Slider styrer hvor mange steg vi tar per frame
        self.speed_scale = ttk.Scale(left, from_=1, to=50,
                                     orient=tk.HORIZONTAL,
                                     command=self.set_speed)
        self.speed_scale.set(1)
        self.speed_scale.pack(fill=tk.X, padx=10)

        ttk.Label(left, text="(= hvor mange tidssteg per oppdatering)").pack(pady=4)

    def make_entry(self, parent, label, default):
        f = ttk.Frame(parent)
        f.pack(anchor="w", padx=10)
        ttk.Label(f, text=label, width=14).pack(side=tk.LEFT)
        e = ttk.Entry(f, width=12)
        e.insert(0, default)
        e.pack(side=tk.LEFT)
        return e

    # ============================
    # Handling
    # ============================

    def update_B(self):
        self.B = np.array([
            float(self.Bx_entry.get()),
            float(self.By_entry.get()),
            float(self.Bz_entry.get())
        ])
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
        self.traj = [self.r.copy()]
        self.update_plot()

    def toggle(self):
        self.running = not self.running
        if self.running:
            self.step()

    def set_speed(self, val):
        self.steps_per_frame = int(float(val))

    # ============================
    # Fysikksteg
    # ============================

    def physics_step(self):
        # Lorentz
        a = self.q * np.cross(self.v, self.B) / self.m

        self.v = self.v + a * self.dt
        self.r = self.r + self.v * self.dt

    def step(self):
        if not self.running:
            return

        # Ta mange små fysiske steg før vi tegner
        for _ in range(self.steps_per_frame):
            self.physics_step()

            self.traj.append(self.r.copy())
            if len(self.traj) > 3000:
                self.traj.pop(0)

        self.update_plot()
        self.root.after(30, self.step)

    # ============================
    # Plot
    # ============================

    def init_plot(self):
        self.traj = [self.r.copy()]

    def update_plot(self):
        self.ax.clear()

        # B-vektor (tegn én pil fra origo)
        Bnorm = np.linalg.norm(self.B)
        if Bnorm > 0:
            Bdir = self.B / Bnorm
            self.ax.quiver(0,0,0, Bdir[0], Bdir[1], Bdir[2],
                           length=0.15, color="green", linewidth=3)
            self.ax.text(Bdir[0]*0.16, Bdir[1]*0.16, Bdir[2]*0.16, "B", color="green")

        # Hastighetsvektor
        vnorm = np.linalg.norm(self.v)
        if vnorm > 0:
            vdir = self.v / vnorm
            self.ax.quiver(self.r[0], self.r[1], self.r[2],
                           vdir[0], vdir[1], vdir[2],
                           length=0.1, color="red", linewidth=2)

        # Partikkel
        self.ax.scatter(self.r[0], self.r[1], self.r[2], color="blue", s=40)

        # Spor
        T = np.array(self.traj)
        if len(T) > 1:
            self.ax.plot(T[:,0], T[:,1], T[:,2], color="blue")

        # Akser
        lim = self.lim
        self.ax.set_xlim(-lim, lim)
        self.ax.set_ylim(-lim, lim)
        self.ax.set_zlim(-lim, lim)

        self.ax.set_xlabel("x [m]")
        self.ax.set_ylabel("y [m]")
        self.ax.set_zlabel("z [m]")
        self.ax.set_title("Ladning i homogent magnetfelt")

        self.canvas.draw()

# ============================

root = tk.Tk()
app = App(root)
root.mainloop()
