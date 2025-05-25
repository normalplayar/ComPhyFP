import tkinter as tk
from tkinter import ttk, messagebox
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation

class CombinedOscillatorApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Oscillator & Pendulum Simulator")

        self.params = {
            "Mass": 4.0,
            "Amplitude": 2.0,
            "Length": 2.0,
            "Damping": 0.9,
            "Driven Force": 2.0
        }

        self.entries = {}

        notebook = ttk.Notebook(root)
        notebook.pack(fill="both", expand=True)

        input_frame = ttk.Frame(root)
        input_frame.pack(pady=10)

        for i, (label, val) in enumerate(self.params.items()):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky="e", padx=5, pady=2)
            entry = ttk.Entry(input_frame, width=10)
            entry.insert(0, str(val))
            entry.grid(row=i, column=1, pady=2)
            self.entries[label] = entry

        # KE tab
        self.kinetic_tab = ttk.Frame(notebook)
        notebook.add(self.kinetic_tab, text="Kinetic Energy")

        self.show_normal = tk.BooleanVar(value=True)
        self.show_damped = tk.BooleanVar(value=True)
        self.show_driven = tk.BooleanVar(value=True)

        chk_frame = ttk.Frame(self.kinetic_tab)
        chk_frame.pack(pady=5)
        ttk.Checkbutton(chk_frame, text="Normal", variable=self.show_normal).grid(row=0, column=0, padx=5)
        ttk.Checkbutton(chk_frame, text="Damped", variable=self.show_damped).grid(row=0, column=1, padx=5)
        ttk.Checkbutton(chk_frame, text="Driven", variable=self.show_driven).grid(row=0, column=2, padx=5)

        btn_frame_ke = ttk.Frame(self.kinetic_tab)
        btn_frame_ke.pack(pady=5)
        ttk.Button(btn_frame_ke, text="Run", command=self.run_kinetic_simulation).grid(row=0, column=0, padx=5)
        self.pause_ke_btn = ttk.Button(btn_frame_ke, text="Pause", command=self.toggle_pause_ke, state="disabled")
        self.pause_ke_btn.grid(row=0, column=1, padx=5)
        self.stop_ke_btn = ttk.Button(btn_frame_ke, text="Stop", command=self.stop_kinetic, state="disabled")
        self.stop_ke_btn.grid(row=0, column=2, padx=5)

        self.fig_ke, self.ax_ke = plt.subplots(figsize=(7, 4))
        self.canvas_ke = FigureCanvasTkAgg(self.fig_ke, master=self.kinetic_tab)
        self.canvas_ke.get_tk_widget().pack()
        self.anim_ke = None
        self.paused_ke = False

        # pendulum tab
        self.pendulum_tab = ttk.Frame(notebook)
        notebook.add(self.pendulum_tab, text="Pendulum")

        self.osc_type = tk.StringVar(value="normal")
        radio_frame = ttk.Frame(self.pendulum_tab)
        radio_frame.pack(pady=5)
        for i, mode in enumerate(["normal", "damped", "driven"]):
            ttk.Radiobutton(radio_frame, text=mode.title(), variable=self.osc_type, value=mode).grid(row=0, column=i, padx=5)

        btn_frame_pend = ttk.Frame(self.pendulum_tab)
        btn_frame_pend.pack(pady=5)
        ttk.Button(btn_frame_pend, text="Run", command=self.run_pendulum_simulation).grid(row=0, column=0, padx=5)
        self.pause_pend_btn = ttk.Button(btn_frame_pend, text="Pause", command=self.toggle_pause_pend, state="disabled")
        self.pause_pend_btn.grid(row=0, column=1, padx=5)
        self.stop_pend_btn = ttk.Button(btn_frame_pend, text="Stop", command=self.stop_pendulum, state="disabled")
        self.stop_pend_btn.grid(row=0, column=2, padx=5)

        self.fig_pend, self.ax_pend = plt.subplots(figsize=(6, 6))
        self.ax_pend.set_aspect('equal')
        self.ax_pend.axis('off')
        self.canvas_pend = FigureCanvasTkAgg(self.fig_pend, master=self.pendulum_tab)
        self.canvas_pend.get_tk_widget().pack()

        self.pivot, = self.ax_pend.plot(0, 0, 'ko', markersize=10)
        self.rod, = self.ax_pend.plot([], [], lw=3, color='blue')
        self.bob = plt.Circle((0, 0), 0.15, fc='red')
        self.ax_pend.add_patch(self.bob)

        self.anim_pend = None
        self.paused_pend = False
        self.t_pend = 0

    def get_params(self):
        try:
            mass = float(self.entries["Mass"].get())
            amp = float(self.entries["Amplitude"].get())
            length = float(self.entries["Length"].get())
            damping = float(self.entries["Damping"].get())
            force = float(self.entries["Driven Force"].get())
            if min(mass, amp, length, damping, force) <= 0:
                raise ValueError
            return mass, amp, length, damping, force
        except ValueError:
            messagebox.showerror("Invalid Input", "All values must be positive numbers.")
            return None

    # KE
    def run_kinetic_simulation(self):
        values = self.get_params()
        if not values:
            return
        m, A, L, b, F0 = values
        g = 9.81
        w0 = math.sqrt(g / L)
        wd = w0 * math.sqrt(1 - (b / (2 * m * w0)) ** 2)
        tau = m / b
        dt, frames = 1/60, 600

        x, y_n, y_d, y_dd = [], [], [], []
        self.ax_ke.clear()
        self.ax_ke.set_xlim(0, 10)
        self.ax_ke.set_ylim(0, 50)
        self.ax_ke.set_xlabel("Time (s)")
        self.ax_ke.set_ylabel("Kinetic Energy (J)")
        self.ax_ke.grid(True)

        lines = []
        if self.show_normal.get():
            ln_n, = self.ax_ke.plot([], [], label="Normal")
            lines.append(("normal", ln_n))
        if self.show_damped.get():
            ln_d, = self.ax_ke.plot([], [], label="Damped")
            lines.append(("damped", ln_d))
        if self.show_driven.get():
            ln_dd, = self.ax_ke.plot([], [], label="Driven")
            lines.append(("driven", ln_dd))

        self.ax_ke.legend()

        t = 0
        def update(_):
            nonlocal t
            if self.paused_ke:
                return [line for _, line in lines]
            t += dt
            x.append(t)
            if self.show_normal.get():
                v = -A * w0 * math.sin(w0 * t)
                y_n.append(0.5 * m * v ** 2)
            if self.show_damped.get():
                Ad = A * math.exp(-t / (2 * tau))
                v = -Ad * wd * math.sin(wd * t)
                y_d.append(0.5 * m * v ** 2)
            if self.show_driven.get():
                Adrive = F0 / math.sqrt((m ** 2) * ((w0 ** 2 - wd ** 2) ** 2) + (b ** 2) * wd ** 2)
                v = -Adrive * wd * math.sin(wd * t)
                y_dd.append(0.5 * m * v ** 2)

            for label, line in lines:
                if label == "normal":
                    line.set_data(x, y_n)
                elif label == "damped":
                    line.set_data(x, y_d)
                elif label == "driven":
                    line.set_data(x, y_dd)
            return [line for _, line in lines]

        if self.anim_ke:
            self.anim_ke.event_source.stop()
        self.paused_ke = False
        self.anim_ke = FuncAnimation(self.fig_ke, update, frames=frames, interval=1000/60, blit=False)
        self.canvas_ke.draw()
        self.pause_ke_btn.config(state="normal", text="Pause")
        self.stop_ke_btn.config(state="normal")

    def toggle_pause_ke(self):
        self.paused_ke = not self.paused_ke
        self.pause_ke_btn.config(text="Resume" if self.paused_ke else "Pause")

    def stop_kinetic(self):
        if self.anim_ke:
            self.anim_ke.event_source.stop()
            self.anim_ke = None
        self.paused_ke = False
        self.ax_ke.clear()
        self.canvas_ke.draw()
        self.pause_ke_btn.config(state="disabled", text="Pause")
        self.stop_ke_btn.config(state="disabled")

    # === Pendulum ===
    def run_pendulum_simulation(self):
        values = self.get_params()
        if not values:
            return
        m, A, L, b, F0 = values
        g = 9.81
        w0 = math.sqrt(g / L)
        wd = w0 * math.sqrt(1 - (b / (2 * m * w0)) ** 2)
        dt = 1 / 60
        self.t_pend = 0
        self.paused_pend = False

        self.ax_pend.set_xlim(-3, 3)
        self.ax_pend.set_ylim(-3, 3)

        osc_mode = self.osc_type.get()
        def theta(t):
            if osc_mode == "normal":
                return A * math.cos(w0 * t)
            elif osc_mode == "damped":
                decay = math.exp(-b * t / (2 * m))
                return A * decay * math.cos(wd * t)
            elif osc_mode == "driven":
                Adr = F0 / math.sqrt((m**2) * (w0**2 - w0**2)**2 + (b**2) * (w0**2))
                return Adr * math.cos(w0 * t)
            return 0

        def update(_):
            if self.paused_pend:
                return
            self.t_pend += dt
            angle = theta(self.t_pend)
            x = L * math.sin(angle)
            y = -L * math.cos(angle)
            self.rod.set_data([0, x], [0, y])
            self.bob.center = (x, y)
            return self.rod,

        if self.anim_pend:
            self.anim_pend.event_source.stop()
        self.anim_pend = FuncAnimation(self.fig_pend, update, interval=1000/60)
        self.canvas_pend.draw()
        self.pause_pend_btn.config(state="normal", text="Pause")
        self.stop_pend_btn.config(state="normal")

    def toggle_pause_pend(self):
        self.paused_pend = not self.paused_pend
        self.pause_pend_btn.config(text="Resume" if self.paused_pend else "Pause")

    def stop_pendulum(self):
        if self.anim_pend:
            self.anim_pend.event_source.stop()
            self.anim_pend = None
        self.paused_pend = False
        self.rod.set_data([], [])
        self.canvas_pend.draw()
        self.pause_pend_btn.config(state="disabled", text="Pause")
        self.stop_pend_btn.config(state="disabled")


if __name__ == "__main__":
    root = tk.Tk()
    app = CombinedOscillatorApp(root)
    root.mainloop()