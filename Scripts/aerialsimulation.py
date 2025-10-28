import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.animation as animation

# Simulation: Aerial Camera Coverage Simulation

# Simulation parameters
dt = 0.1  # Time step (seconds)
total_time = 10  # Total simulation time (seconds)
S = 5  # Speed of the aircraft (units/sec)
FOV_angle = 30  # Field of View (degrees)
gimbal_max_angle = 45  # Max gimbal rotation (degrees)
gimbal_speed = 20  # Gimbal rotation speed (degrees/sec)

# Initialize variables
aircraft_x = 0
aircraft_y = 100
aircraft_trail = []
camera_angles = []

fig, ax = plt.subplots(figsize=(8, 6))
ax.set_xlim(-100, 500)
ax.set_ylim(0, 200)
plt.ion()  # Turn on interactive mode

# Function to update the animation
def update(frame):
    global aircraft_x
    ax.clear()
    ax.set_xlim(-100, 500)
    ax.set_ylim(0, 200)
    
    # Move aircraft
    aircraft_x += S * dt
    aircraft_trail.append((aircraft_x, aircraft_y))
    
    # Rotate gimbal
    angle = gimbal_max_angle * np.sin(frame * dt * np.radians(gimbal_speed))
    camera_angles.append(angle)
    
    # Plot aircraft
    ax.scatter(aircraft_x, aircraft_y, c='blue', label="Aircraft")
    
    # Plot trail
    if len(aircraft_trail) > 1:
        trail_x, trail_y = zip(*aircraft_trail)
        ax.plot(trail_x, trail_y, 'b--')
    
    # Plot FoV as a sector
    fov_radians = np.radians(FOV_angle)
    gimbal_radians = np.radians(angle)
    arc = patches.Wedge((aircraft_x, aircraft_y), 50, np.degrees(gimbal_radians - fov_radians/2),
                         np.degrees(gimbal_radians + fov_radians/2), color='red', alpha=0.3)
    ax.add_patch(arc)
    
    ax.set_title(f"Time: {frame * dt:.1f}s, Gimbal Angle: {angle:.1f}°")
    ax.legend()
    
    plt.draw()
    plt.pause(0.01)

ani = animation.FuncAnimation(fig, update, frames=int(total_time/dt), interval=dt*1000, blit=False)
plt.show()
