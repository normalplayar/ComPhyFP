import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

frames = 60
pi = math.pi
e = math.e

def displacement(A, w, t, phaseconstant=0):
    return A * math.cos(w * t + phaseconstant)

def velocity(A, w, t, phaseconstant=0):
    return -A * w * math.sin(w * t + phaseconstant)

def acceleration(A, w, t, phaseconstant=0):
    return -A * (w**2) * math.cos(w * t + phaseconstant)

def kinetic_energy(m, v):
    return 0.5 * m * (v**2)

def potential_energy(k, x):
    return 0.5 * k * (x**2)

def mechanical_energy(k, A):
    return 0.5 * k * (A**2)

def underdamping_displacement(A0, b, m, t, w, phaseconstant=0):
    return A0 * e**(-t * b / (2 * m)) * math.cos(w * t + phaseconstant)

def damping(b, m, w0):
    return (b / (2 * m * w0))**2

def damped_angular_frequency(w0, b, m):
    return w0 * (math.sqrt(1 - damping(b, m, w0)))

def damped_energy(time, time_constant, E0):
    return E0 * e**(-time / time_constant)

def time_constant_calculation(m, b):
    return m / b

def damped_amplitude(A0, time, time_constant):
    return math.sqrt((A0**2) * e**(-time / time_constant))

def driven_amplitude(F0, m, w0, w, b):
    return F0 / (math.sqrt((m**2) * ((w0**2 - w**2)**2) + (b**2) * (w**2)))

xaxis = []
yaxis_kinetic = []
yaxis_damped_kinetic = []
yaxis_damped_driven_kinetic = []

spring_constant = pi**2
mass = 4
amplitude = 2
damping_constant = 0.9
driven_force = 2

natural_frequency = (spring_constant / mass)**0.5
damped_frequency = damped_angular_frequency(natural_frequency, damping_constant, mass)
time_constant = time_constant_calculation(mass, damping_constant)
initial_total_energy = mechanical_energy(spring_constant, amplitude)

for i in range(480 + 1):
    time = i / frames
    xaxis.append(time)

    # Normal kinetic energy
    v = velocity(amplitude, natural_frequency, time)
    yaxis_kinetic.append(kinetic_energy(mass, v))

    # Damped kinetic energy
    damped_amp = damped_amplitude(amplitude, time, time_constant)
    v_damped = velocity(damped_amp, damped_frequency, time)
    yaxis_damped_kinetic.append(kinetic_energy(mass, v_damped))

    # Damped driven kinetic energy
    driven_amp = driven_amplitude(driven_force, mass, natural_frequency, damped_frequency, damping_constant)
    v_driven = velocity(driven_amp, damped_frequency, time)
    yaxis_damped_driven_kinetic.append(kinetic_energy(mass, v_driven))

fig, axis = plt.subplots()
animated_nkinetic, = axis.plot([], [], label="Normal")
animated_dampedkinetic, = axis.plot([], [], label="Damped")
animated_dampeddrivenkinetic, = axis.plot([], [], label="Damped + Driven")

axis.set_xlim([0, 8])
axis.set_ylim([0, 21])
plt.xlabel("Time / s")
plt.ylabel("Kinetic Energy / J")
plt.legend()
plt.grid(True, linewidth=1)

def update(f):
    animated_nkinetic.set_data(xaxis[:f], yaxis_kinetic[:f])
    animated_dampedkinetic.set_data(xaxis[:f], yaxis_damped_kinetic[:f])
    animated_dampeddrivenkinetic.set_data(xaxis[:f], yaxis_damped_driven_kinetic[:f])
    return animated_nkinetic, animated_dampedkinetic, animated_dampeddrivenkinetic

animation = FuncAnimation(fig, update, frames=480, interval=25)
plt.show()