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
mars_initial_position1 = vp.vector(0., 5.*sun_radius, 0.)
mars_initial_position2 = vp.vector(1., 3*sun_radius, 1.)

mars_initial_velocity = vp.vector(-0.275, 0., 0.)
mars_initial_velocity1 = vp.vector(0.25, 0., 0.)
mars_initial_velocity2 = vp.vector(-0.26, 0., 0.)

# must change further arguments to implement different initial conditions

# Universal constant
gravitational_constant = 1.  # 6.67e-11 N m^2 / kg^2

# Define animation parameters
animation_time_step = 0.001  # seconds
rate_of_animation = 1 / animation_time_step
time_step = 1.
stop_time = 25000.


# Plots of analytic solution
time_rate = 0.1

# mars_initial
mars_position = (((mars_initial_velocity.x**2) + (mars_initial_velocity.y**2)) / 2) * time_step

position_plot_field = vp.graph(title='X_Position vs. Time', xtitle=r't [s]', ytitle='x [km]', fast=False)
position_points = vp.gdots(graph=position_plot_field, color=vp.color.green, size=0.1)
mars_position_curve = vp.gcurve(graph=position_plot_field,
                               data=[[0, mars_position], [stop_time, mars_position]],
                               color=vp.color.magenta)
velocity_plot_field = vp.graph(title='X_Velocity vs. Time', xtitle=r't [s]', ytitle='v [km/s]', fast=False)
mars_velocity_points = vp.gdots(graph=velocity_plot_field, color=vp.color.black, size=0.1)

# mars_initial1
mars1_position = (((mars_initial_velocity1.x**2) + (mars_initial_velocity1.y**2)) / 2) * time_step

position1_plot_field = vp.graph(title='X_Position1 vs. Time', xtitle=r't [s]', ytitle='x [km]', fast=False)
position1_points = vp.gdots(graph=position1_plot_field, color=vp.color.green, size=0.1)
mars1_position_curve = vp.gcurve(graph=position1_plot_field,
                               data=[[0, mars1_position], [stop_time, mars1_position]],
                               color=vp.color.cyan)
velocity1_plot_field = vp.graph(title='X_Velocity1 vs. Time', xtitle=r't [s]', ytitle='v [km/s]', fast=False)
mars1_velocity_points = vp.gdots(graph=velocity1_plot_field, color=vp.color.blue, size=0.1)

# mars_initial2
mars2_position = (((mars_initial_velocity2.x**2) + (mars_initial_velocity2.y**2)) / 2) * time_step

position2_plot_field = vp.graph(title='X_Position2 vs. Time', xtitle=r't [s]', ytitle='x [km]', fast=False)
position2_points = vp.gdots(graph=position2_plot_field, color=vp.color.green, size=0.1)
mars2_position_curve = vp.gcurve(graph=position2_plot_field,
                               data=[[0, mars2_position], [stop_time, mars2_position]],
                               color=vp.color.red)
velocity2_plot_field = vp.graph(title='X_Velocity2 vs. Time', xtitle=r't [s]', ytitle='v [km/s]', fast=False)
mars2_velocity_points = vp.gdots(graph=velocity2_plot_field, color=vp.color.blue, size=0.1)

# Set initial conditions
time = 0.

# Create spheres
sun = vp.sphere(radius=sun_radius, pos=sun_initial_position, color=sun_color)
mars = vp.sphere(radius=mars_radius, pos=mars_initial_position, color=mars_color, make_trail=True)
mars1 = vp.sphere(radius=mars_radius, pos=mars_initial_position1, color=mars_color)
mars2 = vp.sphere(radius=mars_radius, pos=mars_initial_position2, color=mars_color)

# Define properties of objects
mars.mass = 1.
sun.mass = (sun_radius**3 / mars_radius**3) * mars.mass

mars.velocity = mars_initial_velocity
mars1.velocity = mars_initial_velocity1
mars2.velocity = mars_initial_velocity2

# Calculate the separation of the objects, its magnitude, and unit vector
separation = mars.pos - sun.pos
separation1 = mars1.pos - sun.pos
separation2 = mars2.pos - sun.pos

separation_squared = separation.x**2 + separation.y**2 + separation.z**2
separation1_squared = separation1.x**2 + separation1.y**2 + separation1.z**2
separation2_squared = separation2.x**2 + separation2.y**2 + separation2.z**2

separation_unit_vector = separation / np.sqrt(separation_squared)
separation_unit_vector1 = separation1 / np.sqrt(separation1_squared)
separation_unit_vector2 = separation2 / np.sqrt(separation2_squared)

# Calculate the gravitational force of one object on the other
# F_g = G * m1 * m2 / d^2
gravitational_force = -(gravitational_constant * mars.mass * sun.mass / separation_squared) * separation_unit_vector
gravitational_force1 = -(gravitational_constant * mars.mass * sun.mass / separation1_squared) * separation_unit_vector1
gravitational_force2 = -(gravitational_constant * mars.mass * sun.mass / separation2_squared) * separation_unit_vector2

mars.acceleration = gravitational_force / mars.mass
mars1.acceleration = gravitational_force1 / mars.mass
mars2.acceleration = gravitational_force2 / mars.mass

# Move time forward
while time < stop_time:
    vp.rate(rate_of_animation)

    # Kinematic equation for one time step (constant acceleration)
    mars.pos += 0.5 * mars.acceleration * time_step**2 + mars.velocity * time_step
    mars1.pos += 0.5 * mars1.acceleration * time_step**2 + mars1.velocity * time_step
    mars2.pos += 0.5 * mars2.acceleration * time_step**2 + mars2.velocity * time_step
    mars.velocity += mars.acceleration * time_step
    mars1.velocity += mars1.acceleration * time_step
    mars2.velocity += mars2.acceleration * time_step

    # Update separation, force, and acceleration
    separation = mars.pos - sun.pos
    separation1 = mars1.pos - sun.pos
    separation2 = mars2.pos - sun.pos

    separation_squared = separation.x**2 + separation.y**2 + separation.z**2
    separation1_squared = separation1.x**2 + separation1.y**2 + separation1.z**2
    separation2_squared = separation2.x**2 + separation2.y**2 + separation2.z**2

    separation_unit_vector = separation / np.sqrt(separation_squared)
    separation_unit_vector1 = separation1 / np.sqrt(separation1_squared)
    separation_unit_vector2 = separation2 / np.sqrt(separation2_squared)

    gravitational_force = -(gravitational_constant * mars.mass * sun.mass / separation_squared) * separation_unit_vector
    gravitational_force1 = -(gravitational_constant * mars.mass * sun.mass / separation1_squared) * separation_unit_vector1
    gravitational_force2 = -(gravitational_constant * mars.mass * sun.mass / separation2_squared) * separation_unit_vector2

    mars.acceleration = gravitational_force / mars.mass
    mars1.acceleration = gravitational_force1 / mars.mass
    mars2.acceleration = gravitational_force2 / mars.mass

    mars_position_curve.plot(pos=(time, mars.pos.x))
    mars_velocity_points.plot(pos=(time, (mars.velocity.x)**2 + (mars.velocity.y)**2))

    mars1_position_curve.plot(pos=(time, mars1.pos.x))
    mars1_velocity_points.plot(pos=(time, (mars1.velocity.x)**2 + (mars1.velocity.y)**2))

    mars2_position_curve.plot(pos=(time, mars2.pos.x))
    mars2_velocity_points.plot(pos=(time, (mars2.velocity.x)**2 + (mars2.velocity.y)**2))

    # Move time forward
    time += time_step
