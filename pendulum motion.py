import tkinter as tk
from tkinter import ttk
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class PendulumGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Pendulum Oscillation Simulator")

        self.params = {
            "Mass": 4.0,
            "Length": 2.0,
            "Damping": 0.9,
            "Driven Force": 2.0
        }

        input_frame = ttk.Frame(root)
        input_frame.pack(pady=10)

        self.entries = {}
        for i, (label, val) in enumerate(self.params.items()):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky="e", padx=5, pady=2)
            entry = ttk.Entry(input_frame, width=10)
            entry.insert(0, str(val))
            entry.grid(row=i, column=1, pady=2)
            self.entries[label] = entry

        self.oscillation_type = tk.StringVar(value="normal")
        radio_frame = ttk.Frame(root)
        radio_frame.pack(pady=5)

        ttk.Radiobutton(radio_frame, text="Normal", variable=self.oscillation_type, value="normal").grid(row=0, column=0, padx=5)
        ttk.Radiobutton(radio_frame, text="Damped", variable=self.oscillation_type, value="damped").grid(row=0, column=1, padx=5)
        ttk.Radiobutton(radio_frame, text="Driven", variable=self.oscillation_type, value="driven").grid(row=0, column=2, padx=5)

        control_frame = ttk.Frame(root)
        control_frame.pack(pady=10)

        ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation).grid(row=0, column=0, padx=5)
        self.pause_button = ttk.Button(control_frame, text="Pause", command=self.pause_animation, state='disabled')
        self.pause_button.grid(row=0, column=1, padx=5)
        self.stop_button = ttk.Button(control_frame, text="Stop", command=self.stop_animation, state='disabled')
        self.stop_button.grid(row=0, column=2, padx=5)

        self.fig, self.ax = plt.subplots(figsize=(6,6))
        self.ax.set_aspect('equal')
        self.ax.axis('off')
        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack()

        self.pivot, = self.ax.plot(0, 0, 'ko', markersize=14)
        self.rod, = self.ax.plot([], [], lw=3, color='blue')
        self.bob_radius = 0.15
        self.bob = plt.Circle((0,0), self.bob_radius, fc='red')
        self.ax.add_patch(self.bob)

        self.anim = None
        self.paused = False
        self.current_time = 0

    def driven_amplitude(self, F0, m, w0, w, b):
        return F0 / math.sqrt((m**2) * (w0**2 - w**2)**2 + (b**2) * (w**2))

    def normal_theta(self, t, A=0.5, w=None, w0=None):
        if w is None:
            w = w0
        return A * math.cos(w * t)

    def damped_theta(self, t, A=0.5, b=None, m=None, w0=None):
        wd = w0 * math.sqrt(1 - (b / (2 * m * w0))**2)
        decay = math.exp(-b * t / (2 * m))
        return A * decay * math.cos(wd * t)

    def driven_theta(self, t, F0=None, m=None, w0=None, w=None, b=None):
        A = self.driven_amplitude(F0, m, w0, w, b)
        return A * math.cos(w * t)

    def run_simulation(self):
        try:
            self.mass = float(self.entries["Mass"].get())
            self.length = float(self.entries["Length"].get())
            self.damping = float(self.entries["Damping"].get())
            self.driven_force = float(self.entries["Driven Force"].get())
        except ValueError:
            print("Invalid input values")
            return

        self.spring_constant = math.pi**2
        self.natural_freq = math.sqrt(self.spring_constant / self.mass)

        self.ax.set_xlim(-self.length - 0.5, self.length + 0.5)
        self.ax.set_ylim(-self.length - 0.5, self.length + 0.5)

        self.bob.radius = min(self.bob_radius, self.length * 0.15)

        self.FPS = 60
        self.dt = 1 / self.FPS
        self.current_time = 0

        self.osc_type = self.oscillation_type.get()
        self.paused = False
        self.pause_button.config(text="Pause", state='normal')
        self.stop_button.config(state='normal')

        def update(_):
            if not self.paused:
                self.current_time += self.dt

                t = self.current_time

                if self.osc_type == 'normal':
                    angle = self.normal_theta(t, w0=self.natural_freq)
                elif self.osc_type == 'damped':
                    angle = self.damped_theta(t, b=self.damping, m=self.mass, w0=self.natural_freq)
                elif self.osc_type == 'driven':
                    angle = self.driven_theta(t, F0=self.driven_force, m=self.mass, w0=self.natural_freq,
                                              w=self.natural_freq, b=self.damping)
                else:
                    angle = 0

                x = self.length * math.sin(angle)
                y = -self.length * math.cos(angle)

                self.rod.set_data([0, x], [0, y])
                self.bob.center = (x, y)

            return self.rod, self.bob

        if self.anim:
            self.anim.event_source.stop()

        self.anim = FuncAnimation(self.fig, update, interval=1000/self.FPS, blit=True)
        self.canvas.draw()

    def pause_animation(self):
        if not self.anim:
            return
        self.paused = not self.paused
        self.pause_button.config(text="Resume" if self.paused else "Pause")

    def stop_animation(self):
        if not self.anim:
            return
        self.anim.event_source.stop()
        self.anim = None
        self.paused = False
        self.current_time = 0
        self.pause_button.config(text="Pause", state='disabled')
        self.stop_button.config(state='disabled')
        self.rod.set_data([0, 0], [0, 0])
        self.bob.center = (0, 0)
        self.canvas.draw()

if __name__ == "__main__":
    root = tk.Tk()
    app = PendulumGUI(root)
    root.mainloop()