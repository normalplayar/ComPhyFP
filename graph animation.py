import math
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from tkinter import messagebox

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

        button_frame = ttk.Frame(root)
        button_frame.pack(pady=10)

        self.run_button = ttk.Button(button_frame, text="Run Simulation", command=self.run_simulation)
        self.run_button.grid(row=0, column=0, padx=5)

        self.pause_button = ttk.Button(button_frame, text="Pause", command=self.toggle_pause, state="disabled")
        self.pause_button.grid(row=0, column=1, padx=5)

        self.stop_button = ttk.Button(button_frame, text="Stop", command=self.stop_animation, state="disabled")
        self.stop_button.grid(row=0, column=2, padx=5)

        self.fig, self.ax = plt.subplots(figsize=(7, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack()

        self.anim = None
        self.paused = False

    def run_simulation(self):
        try:
            mass = float(self.entries["Mass"].get())
            A = float(self.entries["Amplitude"].get())
            b = float(self.entries["Damping"].get())
            F0 = float(self.entries["Driven Force"].get())

            if mass <= 0 or A <= 0 or b <= 0 or F0 <= 0:
                messagebox.showerror("Invalid Input", "All values must be greater than zero.")
                return

        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numeric values.")
            return

        k = math.pi**2
        w0 = (k / mass)**0.5
        damping_ratio = (b / (2 * mass * w0))**2
        if damping_ratio >= 1:
            messagebox.showerror("Input Error", "Damping is too high — leads to overdamping.\nPlease use a smaller value.")
            return
        wd = w0 * math.sqrt(1 - damping_ratio)
        tau = mass / b
        frames = 480
        fps = 60

        def velocity(A, w, t): return -A * w * math.sin(w * t)
        def KE(m, v): return 0.5 * m * v**2
        def damped_amp(A0, t, tau): return A0 * math.exp(-t / (2 * tau))
        def driven_amp(F0, m, w0, w, b):
            return F0 / math.sqrt((m**2) * ((w0**2 - w**2)**2) + (b**2) * w**2)

        self.mass = mass
        self.A = A
        self.b = b
        self.F0 = F0
        self.k = k
        self.w0 = w0
        self.wd = wd
        self.tau = tau
        self.t = 0
        self.dt = 1 / fps
        self.frames = frames

        self.xdata, self.y_normal, self.y_damped, self.y_driven = [], [], [], []

        self.ax.clear()
        self.ax.set_xlim(0, 8)
        self.ax.set_ylim(0, 25)
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Kinetic Energy (J)")
        self.ax.grid(True)

        self.lines = []

        if self.show_normal.get():
            line_n, = self.ax.plot([], [], label="Normal")
            self.lines.append(("normal", line_n))

        if self.show_damped.get():
            line_d, = self.ax.plot([], [], label="Damped")
            self.lines.append(("damped", line_d))

        if self.show_driven.get():
            line_dd, = self.ax.plot([], [], label="Driven")
            self.lines.append(("driven", line_dd))

        self.ax.legend()

        def update(_):
            if self.paused:
                return [line for _, line in self.lines]

            self.t += self.dt
            self.xdata.append(self.t)

            if self.show_normal.get():
                v = -self.A * self.w0 * math.sin(self.w0 * self.t)
                self.y_normal.append(0.5 * self.mass * v**2)

            if self.show_damped.get():
                Ad = self.A * math.exp(-self.t / (2 * self.tau))
                v = -Ad * self.wd * math.sin(self.wd * self.t)
                self.y_damped.append(0.5 * self.mass * v**2)

            if self.show_driven.get():
                Adrive = self.F0 / math.sqrt(
                    (self.mass**2) * ((self.w0**2 - self.wd**2)**2) + (self.b**2) * self.wd**2)
                v = -Adrive * self.wd * math.sin(self.wd * self.t)
                self.y_driven.append(0.5 * self.mass * v**2)

            for label, line in self.lines:
                if label == "normal":
                    line.set_data(self.xdata, self.y_normal)
                elif label == "damped":
                    line.set_data(self.xdata, self.y_damped)
                elif label == "driven":
                    line.set_data(self.xdata, self.y_driven)

            # self.ax.set_xlim(max(0, self.t - 8), self.t)
            return [line for _, line in self.lines]

        if self.anim:
            self.anim.event_source.stop()
        self.paused = False
        self.anim = FuncAnimation(self.fig, update, frames=frames, interval=1000/fps, blit=False)
        self.canvas.draw()

        self.pause_button.config(state="normal", text="Pause")
        self.stop_button.config(state="normal")

    def toggle_pause(self):
        if self.anim is None:
            return

        if self.paused:
            self.paused = False
            self.pause_button.config(text="Pause")
        else:
            self.paused = True
            self.pause_button.config(text="Resume")

    def stop_animation(self):
        if self.anim is None:
            return
        self.anim.event_source.stop()
        self.anim = None
        self.paused = False

        self.ax.clear()
        self.ax.set_xlim(0, 8)
        self.ax.set_ylim(0, 25)
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Kinetic Energy (J)")
        self.ax.grid(True)
        self.canvas.draw()

        self.pause_button.config(state="disabled", text="Pause")
        self.stop_button.config(state="disabled")

root = tk.Tk()
app = OscillatorGUI(root)
root.mainloop()