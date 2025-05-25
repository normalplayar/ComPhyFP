import tkinter as tk
from tkinter import ttk
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox

class OscillatorKEGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Pendulum and Oscillator Simulator")

        self.params = {
            "Mass": 4.0,
            "Length": 2.0,
            "Damping": 0.9,
            "Driving Force": 2.0
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
        ttk.Radiobutton(radio_frame, text="Damped + Driven", variable=self.oscillation_type, value="dampeddriven").grid(row=0, column=2, padx=5)

        control_frame = ttk.Frame(root)
        control_frame.pack(pady=10)

        ttk.Button(control_frame, text="Run Simulation", command=self.run_simulation).grid(row=0, column=0, padx=5)
        self.pause_button = ttk.Button(control_frame, text="Pause", command=self.pause_animation, state='disabled')
        self.pause_button.grid(row=0, column=1, padx=5)
        self.stop_button = ttk.Button(control_frame, text="Stop", command=self.stop_animation, state='disabled')
        self.stop_button.grid(row=0, column=2, padx=5)

        self.fig, (self.ax_ke, self.ax_pend, self.ax_amp) = plt.subplots(1, 3, figsize=(18, 6))

        self.ax_ke.set_title("Kinetic Energy")
        self.ax_ke.set_xlabel("Time (s)")
        self.ax_ke.set_ylabel("Kinetic Energy (J)")
        self.ax_ke.set_xlim(0, 8)
        self.ax_ke.set_ylim(0, 50)
        self.ax_ke.grid(True)

        self.ax_pend.set_title("Pendulum Oscillation")
        fixed_view_size = 3
        self.ax_pend.set_xlim(-fixed_view_size, fixed_view_size)
        self.ax_pend.set_ylim(-fixed_view_size, fixed_view_size)
        self.ax_pend.set_aspect('equal')
        self.ax_pend.axis('off')

        self.ax_amp.set_title("Amplitude")
        self.ax_amp.set_xlabel("Time (s)")
        self.ax_amp.set_ylabel("Amplitude (rad)")
        self.ax_amp.set_xlim(0, 8)
        self.ax_amp.set_ylim(0, 1.2)
        self.ax_amp.grid(True)

        self.canvas = FigureCanvasTkAgg(self.fig, master=root)
        self.canvas.get_tk_widget().pack()

        self.pivot, = self.ax_pend.plot(0, 0, 'ko', markersize=14)
        self.rod, = self.ax_pend.plot([], [], lw=3, color='blue')
        self.bob_radius = 0.15
        self.bob = plt.Circle((0, 0), 0, fc='red')
        self.ax_pend.add_patch(self.bob)

        self.ke_line, = self.ax_ke.plot([], [], label="Kinetic Energy")
        self.ax_ke.legend()

        self.amp_line, = self.ax_amp.plot([], [], label="Amplitude", color='green')
        self.ax_amp.legend()

        self.anim = None
        self.paused = False
        self.current_time = 0

        self.xdata = []
        self.ydata = []
        self.ampdata = []

    def velocity_normal(self, A, w, t):
        return -A * w * math.sin(w * t)

    def velocity_damped(self, A, wd, t, tau):
        Ad = A * math.exp(-t / (2 * tau))
        return -Ad * wd * math.sin(wd * t)

    def kinetic_energy(self, m, v):
        return 0.5 * m * v ** 2

    def normal_theta(self, t, A=0.5, w=None):
        return A * math.cos(w * t)

    def damped_theta(self, t, A=0.5, b=None, m=None, w0=None):
        wd = w0 * math.sqrt(1 - (b / (2 * m * w0))**2)
        decay = math.exp(-b * t / (2 * m))
        return A * decay * math.cos(wd * t)

    def dampeddriven_theta(self, t, F0=None, m=None, w0=None, w=None, b=None, A0=0.5):
        gamma = b / (2 * m)
        wd = math.sqrt(w0**2 - gamma**2)
        transient = A0 * math.exp(-gamma * t) * math.cos(wd * t)
        A_ss = (F0 / m) / math.sqrt((w0**2 - w**2)**2 + (2 * gamma * w)**2)
        steady_state = A_ss * math.cos(w * t)
        return transient + steady_state

    def amplitude_dampeddriven(self, t, F0=None, m=None, w0=None, w=None, b=None, A0=0.5):
        gamma = b / (2 * m)
        A_ss = (F0 / m) / math.sqrt((w0**2 - w**2)**2 + (2 * gamma * w)**2)
        transient_amp = A0 * math.exp(-gamma * t)
        return transient_amp + A_ss

    def run_simulation(self):
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        try:
            mass = float(self.entries["Mass"].get())
            length = float(self.entries["Length"].get())
            damping = float(self.entries["Damping"].get())
            driving_force = float(self.entries["Driving Force"].get())
            if mass <= 0 or length <= 0 or damping <= 0 or driving_force <= 0:
                messagebox.showerror("Invalid Input", "All values must be greater than zero.")
                return
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter valid numeric values.")
            return

        g = 9.81
        w0 = math.sqrt(g / length)
        damping_ratio = (damping / (2 * mass * w0))**2
        if damping_ratio >= 1:
            damping_ratio = 0.999
        wd = w0 * math.sqrt(1 - damping_ratio)
        tau = mass / damping

        self.mass = mass
        self.length = length
        self.damping = damping
        self.driving_force = driving_force
        self.w0 = w0
        self.wd = wd
        self.tau = tau

        self.A_ke = 2.0

        self.xdata = []
        self.ydata = []
        self.ampdata = []
        self.ke_line.set_data([], [])
        self.amp_line.set_data([], [])
        self.ax_ke.set_xlim(0, 8)
        self.ax_ke.set_ylim(0, 50)
        self.ax_amp.set_xlim(0, 8)
        self.ax_amp.set_ylim(0, 1.2)

        self.rod.set_data([0, 0], [0, 0])
        self.bob.radius = self.bob_radius
        self.bob.center = (0, 0)

        self.current_time = 0
        self.paused = False
        self.pause_button.config(state='normal', text="Pause")
        self.stop_button.config(state='normal')

        osc_type = self.oscillation_type.get()

        def update(frame):
            if not self.paused:
                self.current_time += 1 / 60
                t = self.current_time

                if osc_type == "normal":
                    v = self.velocity_normal(self.A_ke, self.w0, t)
                    amp = 0.5
                elif osc_type == "damped":
                    v = self.velocity_damped(self.A_ke, self.wd, t, self.tau)
                    amp = 0.5 * math.exp(-self.damping * t / (2 * self.mass))
                elif osc_type == "dampeddriven":
                    gamma = self.damping / (2 * self.mass)
                    wd = self.wd
                    angle = self.dampeddriven_theta(t, F0=self.driving_force, m=self.mass, w0=self.w0, w=wd, b=self.damping, A0=0.5)
                    amp = self.amplitude_dampeddriven(t, F0=self.driving_force, m=self.mass, w0=self.w0, w=wd, b=self.damping, A0=0.5)

                    transient_v = (-0.5 * math.exp(-gamma * t) * ( -gamma * math.cos(wd * t) - wd * math.sin(wd * t)))
                    A_ss = (self.driving_force / self.mass) / math.sqrt((self.w0**2 - wd**2)**2 + (2 * gamma * wd)**2)
                    steady_v = -A_ss * wd * math.sin(wd * t)

                    v = transient_v + steady_v

                ke = self.kinetic_energy(self.mass, v)
                self.xdata.append(t)
                self.ydata.append(ke)
                self.ampdata.append(amp)

                self.ke_line.set_data(self.xdata, self.ydata)
                self.ax_ke.set_xlim(0, max(8, t + 1))
                max_ke = max(self.ydata) if self.ydata else 50
                self.ax_ke.set_ylim(0, max(50, max_ke + 10))

                self.amp_line.set_data(self.xdata, self.ampdata)
                self.ax_amp.set_xlim(0, max(8, t + 1))
                max_amp = max(self.ampdata) if self.ampdata else 1.2
                self.ax_amp.set_ylim(0, max(1.2, max_amp + 0.2))

                if osc_type == "normal":
                    theta = self.normal_theta(t, A=0.5, w=self.w0)
                elif osc_type == "damped":
                    theta = self.damped_theta(t, A=0.5, b=self.damping, m=self.mass, w0=self.w0)
                elif osc_type == "dampeddriven":
                    theta = angle

                x = self.length * math.sin(theta)
                y = -self.length * math.cos(theta)

                self.rod.set_data([0, x], [0, y])
                self.bob.center = (x, y)

                self.canvas.draw_idle()

                if t >= 15:
                    self.anim.event_source.stop()
                    self.pause_button.config(state='disabled')
                    self.stop_button.config(state='disabled')

            return self.ke_line, self.rod, self.bob, self.amp_line

        self.anim = FuncAnimation(self.fig, update, blit=False, interval=16, repeat=False)
        self.canvas.draw()

    def pause_animation(self):
        if self.anim is None:
            return
        if self.paused:
            self.anim.event_source.start()
            self.paused = False
            self.pause_button.config(text="Pause")
        else:
            self.anim.event_source.stop()
            self.paused = True
            self.pause_button.config(text="Resume")

    def stop_animation(self):
        if self.anim is not None:
            try:
                self.anim.event_source.stop()
            except Exception:
                pass
            self.anim = None
        self.pause_button.config(state='disabled', text="Pause")
        self.stop_button.config(state='disabled')

        self.xdata.clear()
        self.ydata.clear()
        self.ampdata.clear()
        self.ke_line.set_data([], [])
        self.amp_line.set_data([], [])
        self.rod.set_data([], [])
        self.bob.center = (0, 0)
        self.canvas.draw_idle()

if __name__ == "__main__":
    root = tk.Tk()
    app = OscillatorKEGUI(root)
    root.mainloop()