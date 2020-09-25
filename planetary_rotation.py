import vpython as vp
import numpy as np

# Simulation parameters
#    Sun
sun_radius = 100
sun_color = vp.color.orange
sun_initial_position = vp.vector(0., 0., 0.)
#    Mars
mars_radius = 25
mars_color = vp.color.red
mars_initial_position = vp.vector(0., 8.*sun_radius, 0.)
mars_initial_velocity = vp.vector(-0.5, 0., 0.)

# Universal constant
gravitational_constant = 1.  # 6.67e-11 N m^2 / kg^2

# Define animation parameters
animation_time_step = 0.01  # seconds
rate_of_animation = 1 / animation_time_step
time_step = 1.0
stop_time = 5000.

# Set initial conditions
time = 0.

# Create spheres
sun = vp.sphere(radius=sun_radius, pos=sun_initial_position, color=sun_color)
mars = vp.sphere(radius=mars_radius, pos=mars_initial_position, color=mars_color, make_trail=True)

# Define properties of objects
mars.mass = 1.
sun.mass = (sun_radius**3 / mars_radius**3) * mars.mass
mars.velocity = mars_initial_velocity

# Calculate the separation of the objects, its magnitude, and unit vector
separation = mars.pos - sun.pos
separation_squared = separation.x**2 + separation.y**2 + separation.z**2
separation_unit_vector = separation / np.sqrt(separation_squared)

# Calculate the gravitational force of one object on the other
# F_g = G * m1 * m2 / d^2
gravitational_force = -(gravitational_constant * mars.mass * sun.mass / separation_squared) * separation_unit_vector
mars.acceleration = gravitational_force / mars.mass

# Move time forward
while time < stop_time:
    vp.rate(rate_of_animation)

    # Kinematic equation for one time step (constant acceleration)
    mars.pos += 0.5 * mars.acceleration * time_step**2 + mars.velocity * time_step
    mars.velocity += mars.acceleration * time_step

    # Update separation, force, and acceleration
    separation = mars.pos - sun.pos
    separation_squared = separation.x**2 + separation.y**2 + separation.z**2
    separation_unit_vector = separation / np.sqrt(separation_squared)
    gravitational_force = -(gravitational_constant * mars.mass * sun.mass / separation_squared) * separation_unit_vector
    mars.acceleration = gravitational_force / mars.mass

    # Move time forward
    time += time_step
