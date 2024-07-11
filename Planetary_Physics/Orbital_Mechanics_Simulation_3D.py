import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Function for calculating the orbit using the Euler-Cromer method in 3D
def orbit(orbital_radius, mass):
    n = 365
    dt = 60 * 60 * 24  # time step (1 day in seconds)
    sun_mass = 2e30  # mass of the Sun in kg
    G = 6.67e-11  # gravitational constant
    orbital_radius = orbital_radius * 1.5e11  # converting AU to meters
    gravitational_param = G * (sun_mass + mass)
    orbital_velocity = np.sqrt(gravitational_param / orbital_radius)

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    vx0 = 0
    vy0 = orbital_velocity
    vz0 = 0
    x0 = orbital_radius
    y0 = 0
    z0 = 0
    sun_pos = [0, 0, 0]

    for i in range(n):
        x[i] = x0
        y[i] = y0
        z[i] = z0
        Rx = x0 - sun_pos[0]
        Ry = y0 - sun_pos[1]
        Rz = z0 - sun_pos[2]
        Rmag = np.sqrt(Rx**2 + Ry**2 + Rz**2)
        Ax = -G * (sun_mass / Rmag**3) * Rx
        Ay = -G * (sun_mass / Rmag**3) * Ry
        Az = -G * (sun_mass / Rmag**3) * Rz

        new_vel_x = vx0 + (Ax * dt)
        new_vel_y = vy0 + (Ay * dt)
        new_vel_z = vz0 + (Az * dt)

        new_xpos = x0 + (new_vel_x * dt)
        new_ypos = y0 + (new_vel_y * dt)
        new_zpos = z0 + (new_vel_z * dt)

        x0 = new_xpos
        y0 = new_ypos
        z0 = new_zpos
        vx0 = new_vel_x
        vy0 = new_vel_y
        vz0 = new_vel_z

    return x / 1.5e11, y / 1.5e11, z / 1.5e11

# Function for calculating the path of a meteoroid in 3D
def meteoroid(Ejection_speed, ejection_angle):
    sun_pos = [0, 0, 0]
    Ejection_speed = Ejection_speed * 1000  # convert km/s to m/s
    ejection_angle = ejection_angle * np.pi / 180  # convert degrees to radians
    n = 365
    dt = 60 * 60 * 24
    sun_mass = 2e30
    G = 6.67e-11
    orbital_radius = 1.5 * 1.5e11  # assuming a starting distance from the Sun in meters
    gravitational_param = G * (sun_mass + 6.4e23)
    orbital_velocity = np.sqrt(gravitational_param / orbital_radius)

    # Calculate initial velocity components
    vx0 = orbital_velocity * np.cos(ejection_angle) + Ejection_speed
    vy0 = orbital_velocity * np.sin(ejection_angle)
    vz0 = 0  # assuming initially no velocity along the z-axis

    x = np.zeros(n)
    y = np.zeros(n)
    z = np.zeros(n)
    x0 = orbital_radius
    y0 = 0
    z0 = 0

    Jupiter_x0 = 5.2 * 1.5e11
    Jupiter_y0 = 0
    Jupiter_z0 = 0
    Jupiter_vx0 = 0
    Jupiter_vy0 = 13e3
    Jupiter_vz0 = 0

    for i in range(n):
        x[i] = x0
        y[i] = y0
        z[i] = z0
        Rx = x0 - sun_pos[0]
        Ry = y0 - sun_pos[1]
        Rz = z0 - sun_pos[2]
        Rmag = np.sqrt(Rx**2 + Ry**2 + Rz**2)

        RxJ = Jupiter_x0 - sun_pos[0]
        RyJ = Jupiter_y0 - sun_pos[1]
        RzJ = Jupiter_z0 - sun_pos[2]
        RmagJ = np.sqrt(RxJ**2 + RyJ**2 + RzJ**2)

        rjx = RxJ - Rx
        rjy = RyJ - Ry
        rjz = RzJ - Rz
        rj = np.sqrt(rjx**2 + rjy**2 + rjz**2)

        AxJ = -G * (1.9e27 / RmagJ**3) * RxJ
        AyJ = -G * (1.9e27 / RmagJ**3) * RyJ
        AzJ = -G * (1.9e27 / RmagJ**3) * RzJ

        Ax = (-G * (sun_mass / Rmag**3) * Rx) - (-G * (1.9e27 / rj**3) * rjx)
        Ay = (-G * (sun_mass / Rmag**3) * Ry) - (-G * (1.9e27 / rj**3) * rjy)
        Az = (-G * (sun_mass / Rmag**3) * Rz) - (-G * (1.9e27 / rj**3) * rjz)

        new_vel_x = vx0 + (Ax * dt)
        new_vel_y = vy0 + (Ay * dt)
        new_vel_z = vz0 + (Az * dt)

        new_vel_xJ = Jupiter_vx0 + (AxJ * dt)
        new_vel_yJ = Jupiter_vy0 + (AyJ * dt)
        new_vel_zJ = Jupiter_vz0 + (AzJ * dt)

        new_xpos = x0 + (new_vel_x * dt)
        new_ypos = y0 + (new_vel_y * dt)
        new_zpos = z0 + (new_vel_z * dt)

        new_xposJ = Jupiter_x0 + (new_vel_xJ * dt)
        new_yposJ = Jupiter_y0 + (new_vel_yJ * dt)
        new_zposJ = Jupiter_z0 + (new_vel_zJ * dt)

        x0 = new_xpos
        y0 = new_ypos
        z0 = new_zpos
        vx0 = new_vel_x
        vy0 = new_vel_y
        vz0 = new_vel_z

        Jupiter_x0 = new_xposJ
        Jupiter_y0 = new_yposJ
        Jupiter_z0 = new_zposJ
        Jupiter_vx0 = new_vel_xJ
        Jupiter_vy0 = new_vel_yJ
        Jupiter_vz0 = new_vel_zJ

    return x / 1.5e11, y / 1.5e11, z / 1.5e11

# Planetary data: orbital radius in AU and mass in kg
planets = {
    'Mercury': {'orbital_radius': 0.39, 'mass': 0.33e24},
    'Venus': {'orbital_radius': 0.723, 'mass': 4.867e24},
    'Earth': {'orbital_radius': 1, 'mass': 5.97e24},
    'Mars': {'orbital_radius': 1.524, 'mass': 0.64171e24},
    'Jupiter': {'orbital_radius': 5.2, 'mass': 1898.19e24},
    'Saturn': {'orbital_radius': 9.58, 'mass': 568.34e24},
    'Uranus': {'orbital_radius': 19.22, 'mass': 86.813e24},
    'Neptune': {'orbital_radius': 30.05, 'mass': 102.413e24},
}

# Calculate orbits for each planet
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')

for planet, data in planets.items():
    x, y, z = orbit(data['orbital_radius'], data['mass'])
    ax.plot(x, y, z, label=planet)

# Calculate and plot the path of the meteoroid
met_x, met_y, met_z = meteoroid(1, 270)  # Example values: speed in km/s, angle in degrees
ax.plot(met_x, met_y, met_z, label='Meteoroid', linestyle='dashed')

ax.scatter(0, 0, 0, label='Sun', color='yellow', s=100)

ax.set_xlim([-5, 5])
ax.set_ylim([-5, 5])
ax.set_zlim([-5, 5])
ax.grid(linestyle='dotted')
ax.set_title('Orbits of Planets and Path of Meteoroid in 3D')
ax.set_xlabel('X-position (AU)')
ax.set_ylabel('Y-position (AU)')
ax.set_zlabel('Z-position (AU)')
ax.legend()

plt.show()
