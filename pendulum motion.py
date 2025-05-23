import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

# normal, damped, driven
oscillation_type = 'driven'

mass = 4
spring_constant = math.pi**2

# for damped, driven
damping_constant = 0.9

# for driven
driven_force = 2.0

length = 2
bob_radius = 0.1

FPS = 60
duration = 8

natural_freq = math.sqrt(spring_constant / mass)
time_constant = mass / damping_constant

def driven_amplitude(F0, m, w0, w, b):
    return F0 / math.sqrt((m**2) * (w0**2 - w**2)**2 + (b**2) * (w**2))

def normal_theta(t, A=0.5, w=None):
    if w is None:
        w = natural_freq
    return A * math.cos(w * t)

def damped_theta(t, A=0.5, b=None, m=None, w0=None):
    if b is None: b = damping_constant
    if m is None: m = mass
    if w0 is None: w0 = natural_freq
    # Damped angular frequency
    wd = w0 * math.sqrt(1 - (b / (2 * m * w0))**2)
    decay = math.exp(-b * t / (2 * m))
    return A * decay * math.cos(wd * t)

def driven_theta(t, F0=None, m=None, w0=None, w=None, b=None):
    if F0 is None: F0 = driven_force
    if m is None: m = mass
    if w0 is None: w0 = natural_freq
    if w is None: w = w0
    if b is None: b = damping_constant
    A = driven_amplitude(F0, m, w0, w, b)
    return A * math.cos(w * t)

fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim(-length - 0.5, length + 0.5)
ax.set_ylim(-length - 0.5, length + 0.5)
ax.set_aspect('equal')
ax.axis('off')

pivot, = ax.plot(0, 0, 'ko', markersize=8)

rod, = ax.plot([], [], lw=3, color='blue')
bob = plt.Circle((0,0), bob_radius, fc='red')
ax.add_patch(bob)

frames = int(FPS * duration)
dt = 1 / FPS

def update(frame):
    t = frame * dt

    if oscillation_type == 'normal':
        angle = normal_theta(t)
    elif oscillation_type == 'damped':
        angle = damped_theta(t)
    elif oscillation_type == 'driven':
        angle = driven_theta(t, F0=driven_force, m=mass, w0=natural_freq, w=natural_freq, b=damping_constant)
    else:
        angle = 0

    x = length * math.sin(angle)
    y = -length * math.cos(angle)

    rod.set_data([0, x], [0, y])
    bob.center = (x, y)
    return rod, bob

anim = FuncAnimation(fig, update, frames=frames, interval=1000/FPS, blit=True)
plt.show()