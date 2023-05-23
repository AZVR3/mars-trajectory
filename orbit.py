# Import modules
import turtle
import math

# Define constants
G = 6.674e-11  # Gravitational constant
AU = 1.496e11  # Astronomical unit
SCALE = 250 / AU  # Scale factor for turtle graphics
TIME_STEP = 100000  # Time step in seconds (increased by a factor of 100)


# Define classes
class Body:
    """A class to represent a celestial body"""

    def __init__(self, name, mass, position, velocity, color):
        self.name = name  # Name of the body
        self.mass = mass  # Mass of the body in kg
        self.position = position  # Position of the body in m
        self.velocity = velocity  # Velocity of the body in m/s
        self.color = color  # Color of the body for turtle graphics
        self.turtle = turtle.Turtle()  # Turtle object to draw the body
        self.turtle.shape("circle")  # Set the shape to circle
        self.turtle.color(color)  # Set the color
        self.turtle.penup()  # Lift the pen up
        self.turtle.goto(position[0] * SCALE, position[1] * SCALE)  # Go to the initial position
        self.turtle.pendown()  # Put the pen down

    def attraction(self, other):
        """Calculate the gravitational force exerted by another body"""
        if self is other:
            raise ValueError("Attraction of object %r to itself requested" % self.name)
        dx = other.position[0] - self.position[0]  # x distance between the bodies
        dy = other.position[1] - self.position[1]  # y distance between the bodies
        dz = other.position[2] - self.position[2]  # z distance between the bodies

        d = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)  # Distance between the bodies

        if d == 0:
            raise ValueError("Collision between objects %r and %r" % (self.name, other.name))

        f = G * self.mass * other.mass / (d ** 2)  # Magnitude of the force

        theta = math.atan2(dy, dx)  # Angle of the force
        phi = math.atan2(dz, math.sqrt(dx ** 2 + dy ** 2))  # Angle of the force
        fx = math.cos(theta) * math.cos(phi) * f  # x component of the force
        fy = math.sin(theta) * math.cos(phi) * f  # y component of the force
        fz = math.sin(phi) * f  # z component of the force

        return fx, fy, fz

    def update(self, bodies):
        """Update the position and velocity of the body"""
        fx_total = 0.0  # Total x force acting on the body
        fy_total = 0.0  # Total y force acting on the body
        fz_total = 0.0  # Total z force acting on the body

        for body in bodies:  # For each other body in the system
            if body is self:
                continue  # Skip itself
            fx, fy, fz = self.attraction(body)  # Calculate the force exerted by the other body
            fx_total += fx  # Add to the total x force
            fy_total += fy  # Add to the total y force
            fz_total += fz  # Add to the total z force

        ax = fx_total / self.mass  # x acceleration of the body
        ay = fy_total / self.mass  # y acceleration of the body
        az = fz_total / self.mass  # z acceleration of the body

        self.velocity[0] += ax * TIME_STEP  # Update x velocity based on acceleration and time step
        self.velocity[1] += ay * TIME_STEP  # Update y velocity based on acceleration and time step
        self.velocity[2] += az * TIME_STEP  # Update z velocity based on acceleration and time step

        self.position[0] += self.velocity[0] * TIME_STEP  # Update x position based on velocity and time step
        self.position[1] += self.velocity[1] * TIME_STEP  # Update y position based on velocity and time step
        self.position[2] += self.velocity[2] * TIME_STEP  # Update z position based on velocity and time step

    def draw(self):
        """Draw the body using turtle graphics"""
        self.turtle.goto(self.position[0] * SCALE, self.position[1] * SCALE)  # Go to the new position


class Spaceship(Body):
    """A class to represent a spaceship"""

    def __init__(self, name, mass, position, velocity, color):
        super().__init__(name, mass, position, velocity, color)
        self.turtle.shape("triangle")

    def launch(self):

        if abs(self.position[0]) < AU / 4 and abs(self.position[1]) < AU / 4:
            print("Launching spaceship from Earth towards Mars using Hohmann transfer")
            r1 = math.sqrt(earth.position[0] ** 2 + earth.position[1] ** 2 + earth.position[2] ** 2)
            r2 = math.sqrt(mars.position[0] ** 2 + mars.position[1] ** 2 + mars.position[2] ** 2)
            a = (r1 + r2) / 2
            v1_circ = math.sqrt(G * sun.mass / r1)
            v1_ellip = math.sqrt(G * sun.mass * (2 / r1 - 1 / a))
            dv1 = v1_ellip - v1_circ
            theta_earth_mars_optimal = math.acos(1 - r1 / a * (dv1 / v1_circ) ** 2)
            theta_earth_mars_current = math.atan2(earth.position[1] - mars.position[1],
                                                  earth.position[0] - mars.position[0])
            theta_earth_mars_diff = theta_earth_mars_optimal - theta_earth_mars_current
            if theta_earth_mars_diff > 0:
                print(f"Wait for {theta_earth_mars_diff * 180 / math.pi:.2f} degrees before launching")
                return
            else:
                print(f"Launch now with a velocity boost of {dv1:.2f} m/s")
                vx_boost = dv1 * math.cos(theta_earth_mars_current)
                vy_boost = dv1 * math.sin(theta_earth_mars_current)
                vz_boost = 0
                vx_new = earth.velocity[0] + vx_boost
                vy_new = earth.velocity[1] + vy_boost
                vz_new = earth.velocity[2] + vz_boost
                self.velocity = [vx_new, vy_new, vz_new]


# Create objects for the sun, earth and Mars using data from NASA (https://nssdc.gsfc.nasa.gov/planetary/factsheet/)
sun = Body("Sun", 1.98892e30, [0.0, 0.0, 0.0], [0.000000000000e+00, 0.000000000000e+00, -7.902486461457e-12], "yellow")
earth = Body("Earth", 5.9742e24, [-9.906191795391e+10, 1.071492238645e+11, -6.724107636888e+06],
             [-1.989124794559e+04, -1.848356583065e+04, -7.509779157144e-01], "blue")
mars = Body("Mars", 6.4171e23, [-2.279366267344e+11, -3.394709621470e+10, 5.026069609079e+09],
            [3.449470346970e+03, -2.111870894296e+04, -5.051854701265e+02], "red")

# Create a spaceship object with initial position and velocity slightly different from earth's
spaceship = Spaceship("Spaceship", 10000,
                      [-9.906191795391e+10 + 10000, 1.071492238645e+11, -6.724107636888e+06 + AU / 10],
                      [-1.989124794559e+04 - 10000 + 1000, -1.848356583065e+04 + 5000, -7.509779157144e-01], "white")

# Create a list of all bodies in the system
bodies = [sun, earth, mars, spaceship]

# Set up turtle graphics
turtle.bgcolor("black")
turtle.screensize(750, 750)
turtle.title("Orbit Simulation")

# Create a turtle object to write text on screen
text_turtle = turtle.Turtle()
text_turtle.hideturtle()
text_turtle.penup()
text_turtle.goto(300, -300)
text_turtle.color("white")

# Initialize elapsed time and speed variables
elapsed_time = 0
speed = 0

# Run simulation loop
while True:
    for body in bodies:
        body.update(bodies)
        body.draw()

    spaceship.launch()

    turtle.update()

    dx = mars.position[0] - spaceship.position[0]
    dy = mars.position[1] - spaceship.position[1]
    dz = mars.position[2] - spaceship.position[2]

    d = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    if d < AU / 10:
        print("Spaceship reached Mars")
        print(f"Final elapsed time: {elapsed_time / 3600:.2f} hours")
        print(f"Final speed: {speed:.2f} m/s")
        break

    elapsed_time += TIME_STEP
    speed = math.sqrt(spaceship.velocity[0] ** 2 + spaceship.velocity[1] ** 2 + spaceship.velocity[2] ** 2)

    text_turtle.clear()
    text_turtle.write(f"Time: {elapsed_time / 3600:.2f} hours\nSpeed: {speed:.2f} m/s", font=("Arial", 16, "normal"))
