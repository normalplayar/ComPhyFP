import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.animation import FuncAnimation
from math import sin, cos, radians, degrees, sqrt, atan2

class OscillatorGUI:
    def __init__(self, master):
        self.paused = False
        self.running = False

        self.master = master
        self.master.title("Oscillation Simulator")

        self.create_input_fields()

        # Create 3 subplots
        self.fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(1, 3, figsize=(15, 5))
        self.fig.tight_layout(pad=4.0)

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.master)
        self.canvas.get_tk_widget().grid(row=7, column=0, columnspan=4)

        self.ani = None

    def create_input_fields(self):
        labels = ['Length (m):', 'Mass (kg):', 'Initial Angle (degrees):',
                  'Damping Coefficient:', 'Driving Force Amplitude:', 'Driving Frequency (rad/s):']
        defaults = ['1', '1', '30', '0.1', '0.5', '1']

        self.entries = []
        for i, (label, default) in enumerate(zip(labels, defaults)):
            ttk.Label(self.master, text=label).grid(row=i, column=0, sticky='w')
            entry = ttk.Entry(self.master)
            entry.insert(0, default)
            entry.grid(row=i, column=1)
            self.entries.append(entry)

        self.motion_type = tk.StringVar(value='normal')
        for i, option in enumerate(['normal', 'damped', 'damped + driven']):
            ttk.Radiobutton(self.master, text=option, variable=self.motion_type, value=option).grid(row=i, column=2, sticky='w')

        ttk.Button(self.master, text="Start", command=self.start_simulation).grid(row=6, column=1)

    def start_simulation(self):
        if self.ani:
            self.ani.event_source.stop()

        self.length = float(self.entries[0].get())
        self.mass = float(self.entries[1].get())
        self.angle = radians(float(self.entries[2].get()))
        self.damping = float(self.entries[3].get())
        self.F0 = float(self.entries[4].get())
        self.w_drive = float(self.entries[5].get())
        self.motion = self.motion_type.get()

        self.t = 0
        self.dt = 0.05
        self.theta0 = self.angle
        self.omega0 = sqrt(9.81 / self.length)
        self.gamma = self.damping / (2 * self.mass)
        self.w_d = sqrt(max(self.omega0**2 - self.gamma**2, 0))
        self.A = self.theta0

        self.initial_me = 0.5 * self.mass * (self.A * self.omega0)**2 + self.mass * 9.81 * self.length * (1 - cos(self.A))

        self.x_data = []
        self.ke_data = []
        self.pe_data = []
        self.me_data = []
        self.te_data = []
        self.theta_data = []

        self.ani = FuncAnimation(self.fig, self.update_plot, interval=50)
        self.canvas.draw()

    def normal_theta(self, t, A, w):
        return A * cos(w * t)

    def velocity_normal(self, A, w, t):
        return -A * w * sin(w * t)

    def damped_theta(self, t, A, gamma, wd):
        return A * np.exp(-gamma * t) * np.cos(wd * t)

    def velocity_damped(self, A, gamma, wd, t):
        return A * np.exp(-gamma * t) * (-gamma * np.cos(wd * t) - wd * np.sin(wd * t))

    def dampeddriven_theta(self, t, A, gamma, w0, w_drive, F0, m):
        transient = A * np.exp(-gamma * t) * np.cos(w0 * t)
        denom = sqrt((w0**2 - w_drive**2)**2 + (2 * gamma * w_drive)**2)
        amplitude = (F0 / m) / denom
        phase = atan2(2 * gamma * w_drive, w0**2 - w_drive**2)
        steady_state = amplitude * cos(w_drive * t - phase)
        return transient + steady_state

    def velocity_numerical(self, t):
        delta = 0.001
        theta_now = self.dampeddriven_theta(t, self.A, self.gamma, self.omega0, self.w_drive, self.F0, self.mass)
        theta_prev = self.dampeddriven_theta(t - delta, self.A, self.gamma, self.omega0, self.w_drive, self.F0, self.mass)
        return (theta_now - theta_prev) / delta

    def update_plot(self, frame):
        t = self.t
        A = self.A
        gamma = self.gamma
        w0 = self.omega0
        w_drive = self.w_drive
        F0 = self.F0
        m = self.mass
        L = self.length

        if self.motion == 'normal':
            theta = self.normal_theta(t, A, w0)
            v = self.velocity_normal(A, w0, t)
        elif self.motion == 'damped':
            theta = self.damped_theta(t, A, gamma, self.w_d)
            v = self.velocity_damped(A, gamma, self.w_d, t)
        else:
            theta = self.dampeddriven_theta(t, A, gamma, w0, w_drive, F0, m)
            v = self.velocity_numerical(t)

        x = L * sin(theta)
        y = -L * cos(theta)

        KE = 0.5 * m * (v * L)**2
        PE = m * 9.81 * L * (1 - cos(theta))
        ME = KE + PE
        TE = self.initial_me - ME

        self.x_data.append(t)
        self.ke_data.append(KE)
        self.pe_data.append(PE)
        self.me_data.append(ME)
        self.te_data.append(max(TE, 0))
        self.theta_data.append(degrees(theta))  # In degrees

        while self.x_data and self.x_data[-1] - self.x_data[0] > 8:
            self.x_data.pop(0)
            self.ke_data.pop(0)
            self.pe_data.pop(0)
            self.me_data.pop(0)
            self.te_data.pop(0)
            self.theta_data.pop(0)

        # === Pendulum animation ===
        self.ax1.clear()
        self.ax1.set_xlim(-2.5, 2.5)
        self.ax1.set_ylim(-2.5, 1)
        self.ax1.axis('off')
        self.ax1.plot([0, x], [0, y], lw=2)
        self.ax1.plot(x, y, 'o', markersize=20)
        self.ax1.set_title("Pendulum Motion")

        # === Angular displacement plot ===
        self.ax2.clear()
        self.ax2.plot(self.x_data, self.theta_data, label='Angular Displacement (°)', color='tab:orange')
        self.ax2.set_title("Angular Displacement vs Time")
        self.ax2.set_ylabel("Degrees")
        self.ax2.set_xlabel("Time (s)")
        self.ax2.legend()

        # === Energy plot ===
        self.ax3.clear()
        self.ax3.plot(self.x_data, self.ke_data, label='Kinetic Energy')
        self.ax3.plot(self.x_data, self.pe_data, label='Potential Energy')
        self.ax3.plot(self.x_data, self.me_data, label='Mechanical Energy')
        self.ax3.plot(self.x_data, self.te_data, label='Thermal Energy')
        self.ax3.set_title("Energy vs Time")
        self.ax3.set_xlabel("Time (s)")
        self.ax3.set_ylabel("Energy (J)")
        self.ax3.legend(loc='upper right')

        self.t += self.dt
        self.canvas.draw()

    def pause_animation(self):
        if self.anim is None:
            return
        self.paused = not self.paused
        if self.paused:
            self.pause_button.config(text="Resume")
        else:
            self.pause_button.config(text="Pause")

    def stop_animation(self):
        if self.anim is not None:
            self.anim.event_source.stop()
            self.anim = None
        self.pause_button.config(state='disabled', text="Pause")
        self.stop_button.config(state='disabled')

if __name__ == "__main__":
    root = tk.Tk()
    gui = OscillatorGUI(root)
    root.mainloop()