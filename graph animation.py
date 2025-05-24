import math
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

class OscillatorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Oscillator Kinetic Energy Animation")

        self.entries = {}
        self.options = {
            "Mass": 4,
            "Amplitude": 2,
            "Damping": 0.9,
            "Driven Force": 2
        }

        input_frame = ttk.Frame(root)
        input_frame.pack(pady=5)

        for i, (label, default) in enumerate(self.options.items()):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky="e")
            entry = ttk.Entry(input_frame)
            entry.insert(0, str(default))
            entry.grid(row=i, column=1)
            self.entries[label] = entry

        self.show_normal = tk.BooleanVar(value=True)
        self.show_damped = tk.BooleanVar(value=True)
        self.show_driven = tk.BooleanVar(value=True)

        checkbox_frame = ttk.Frame(root)
        checkbox_frame.pack(pady=5)

        ttk.Checkbutton(checkbox_frame, text="Normal Oscillation", variable=self.show_normal).grid(row=0, column=0)
        ttk.Checkbutton(checkbox_frame, text="Damped Oscillation", variable=self.show_damped).grid(row=0, column=1)
        ttk.Checkbutton(checkbox_frame, text="Damped + Driven Oscillation", variable=self.show_driven).grid(row=0, column=2)

        ttk.Button(root, text="Run Simulation", command=self.run_simulation).pack(pady=10)

        self.fig, self.ax = plt.subplots(figsize=(7, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack()

        self.anim = None

    def run_simulation(self):
        try:
            mass = float(self.entries["Mass"].get())
            A = float(self.entries["Amplitude"].get())
            b = float(self.entries["Damping"].get())
            F0 = float(self.entries["Driven Force"].get())
        except ValueError:
            print("Invalid input")
            return

        k = math.pi**2
        w0 = (k / mass)**0.5
        damping_ratio = (b / (2 * mass * w0))**2
        if damping_ratio >= 1:
            damping_ratio = 0.999
        wd = w0 * math.sqrt(1 - damping_ratio)
        tau = mass / b
        frames = 480
        fps = 60

        def velocity(A, w, t): return -A * w * math.sin(w * t)
        def KE(m, v): return 0.5 * m * v**2
        def damped_amp(A0, t, tau): return A0 * math.exp(-t / (2 * tau))
        def driven_amp(F0, m, w0, w, b):
            return F0 / math.sqrt((m**2) * ((w0**2 - w**2)**2) + (b**2) * w**2)

        xdata, y_normal, y_damped, y_driven = [], [], [], []

        for i in range(frames + 1):
            t = i / fps
            xdata.append(t)

            if self.show_normal.get():
                v = velocity(A, w0, t)
                y_normal.append(KE(mass, v))

            if self.show_damped.get():
                Ad = damped_amp(A, t, tau)
                vd = velocity(Ad, wd, t)
                y_damped.append(KE(mass, vd))

            if self.show_driven.get():
                Adrive = driven_amp(F0, mass, w0, wd, b)
                vdrive = velocity(Adrive, wd, t)
                y_driven.append(KE(mass, vdrive))

        self.ax.clear()
        self.ax.set_xlim(0, 8)
        self.ax.set_ylim(0, 25)
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Kinetic Energy (J)")
        self.ax.grid(True)

        lines = []

        if self.show_normal.get():
            line_n, = self.ax.plot([], [], label="Normal")
            lines.append(("normal", line_n))

        if self.show_damped.get():
            line_d, = self.ax.plot([], [], label="Damped")
            lines.append(("damped", line_d))

        if self.show_driven.get():
            line_dd, = self.ax.plot([], [], label="Driven")
            lines.append(("driven", line_dd))

        self.ax.legend()

        def update(frame):
            for label, line in lines:
                if label == "normal":
                    line.set_data(xdata[:frame], y_normal[:frame])
                elif label == "damped":
                    line.set_data(xdata[:frame], y_damped[:frame])
                elif label == "driven":
                    line.set_data(xdata[:frame], y_driven[:frame])
            return [line for _, line in lines]

        if self.anim:
            self.anim.event_source.stop()

        self.anim = FuncAnimation(self.fig, update, frames=frames, interval=25, blit=True)
        self.canvas.draw()

root = tk.Tk()
app = OscillatorGUI(root)
root.mainloop()