import numpy as np
import matplotlib.pyplot as plt
#for planet orbits
def orbit(orbital_radius, mass):
  n = 365
  dt = 60 * 60 * 24
  sun_mass = 2e30
  G = 6.67e-11
  orbital_radius = orbital_radius * 1.5e11
  gravitational_param = G * (sun_mass + mass)
  orbital_velocity = np.sqrt(gravitational_param/orbital_radius)

  x = np.zeros(n)
  y = np.zeros(n)
  vx0 = 0
  vy0 = orbital_velocity
  x0 = orbital_radius
  y0 = 0
  sun_pos = [0,0]
  for i in range(n):
    x[i] = x0
    y[i] = y0
    theta = np.arctan2(x0, y0)
    Rx = x0 - (sun_pos[0])
    Ry = y0 - (sun_pos[1])

    Rmag = np.sqrt(Rx**2 + Ry**2)

    Ax = -G * (sun_mass/Rmag**3) * Rx
    Ay = -G * (sun_mass/Rmag**3) * Ry

    new_vel_x = vx0 + (Ax * dt)
    new_vel_y = vy0 + (Ay * dt)

    new_xpos = x0 + (new_vel_x * dt)
    new_ypos = y0 + (new_vel_y * dt)

    x0 = new_xpos
    y0 = new_ypos
    vx0 = new_vel_x
    vy0 = new_vel_y

  return x/1.5e11, y/1.5e11
def meteoroid(Ejection_speed, ejection_angle):
  sun_pos = [0,0]
  Ejection_speed = Ejection_speed * 1000
  ejection_angle = ejection_angle * np.pi/180
  n = 365
  dt = 60 * 60 * 24
  sun_mass = 2e30
  G = 6.67e-11
  orbital_radius = 1.5 * 1.5e11
  gravitational_param = G * (sun_mass + 6.4e23)
  orbital_velocity = np.sqrt(gravitational_param/orbital_radius) + Ejection_speed

  x0 = orbital_radius
  y0 = 0
  vx0 = orbital_velocity * np.cos(ejection_angle)
  vy0 = orbital_velocity * np.sin(ejection_angle)

  Jupiter_x0 = 5.2 * 1.5e11
  Jupiter_y0 = 0
  Jupiter_vx0 = 0
  Jupiter_vy0 = 13e3

  x = np.zeros(n)
  y = np.zeros(n)

  for i in range(n):
    x[i] = x0
    y[i] = y0
    theta = np.arctan2(x0, y0)
    Rx = x0 - (sun_pos[0])
    Ry = y0 - (sun_pos[1])

    RxJ = Jupiter_x0 - (sun_pos[0])
    RyJ = Jupiter_y0 - (sun_pos[1])
    rjx = RxJ - Rx
    rjy = RyJ - Ry
    rj = np.sqrt(rjx**2 + rjy**2)

    Rmag = np.sqrt(Rx**2 + Ry**2)
    RmagJ = np.sqrt(RxJ**2 + RyJ**2)

    AxJ =  (-G * (1.9e27/RmagJ**3)*RxJ)
    AyJ =  (-G * (1.9e27/RmagJ**3)*RyJ)

    Ax = (-G * (sun_mass/Rmag**3) * Rx) - (-G * (1.9e27//rj**3)*rjx)
    Ay = (-G * (sun_mass/Rmag**3) * Ry) - (-G * (1.9e27//rj**3)*rjy)

    new_vel_x = vx0 + (Ax * dt)
    new_vel_y = vy0 + (Ay * dt)

    new_vel_xJ = Jupiter_vx0 + (AxJ * dt)
    new_vel_yJ = Jupiter_vy0 + (AyJ * dt)

    new_xpos = x0 + (new_vel_x * dt)
    new_ypos = y0 + (new_vel_y * dt)

    new_xposJ = Jupiter_x0 + (new_vel_xJ * dt)
    new_yposJ = Jupiter_y0 + (new_vel_yJ * dt)

    x0 = new_xpos
    y0 = new_ypos
    vx0 = new_vel_x
    vy0 = new_vel_y

    Jupiter_x0 = new_xposJ
    Jupiter_y0 = new_yposJ
    Jupiter_vx0 = new_vel_xJ
    Jupiter_vy0 = new_vel_yJ
  return x/1.5e11,y/1.5e11
Jupiter = orbit(5.2, 1.9e27)
Earth = orbit(1, 5.97e24)
Mars = orbit(1.5, 6.4e23)

Met = meteoroid(1, 240) 


fig = plt.figure()
plt.scatter(0,0, label = 'Sun')
plt.plot(Earth[0],Earth[1], label = 'Earth')
plt.plot(Mars[0], Mars[1], label = 'Mars')
plt.plot(Jupiter[0],Jupiter[1], label = 'Jupiter')
plt.plot(Met[0],Met[1],label = 'Meteoroid', linestyle = 'dashed')
plt.ylim(-5.5,5.5)
plt.xlim(-5.5,5.5)
plt.grid(linestyle = 'dotted')
plt.title('Meteoroid Simulation')
plt.xlabel('au')
plt.ylabel('au')
plt.legend()
plt.show()
