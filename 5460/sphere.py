import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Parameters
r0 = np.pi/2               # starting colatitude (radians)
num_balls = 5
T_values = np.linspace(0, 4, num_balls)  # each ball begins rolling at t=T

# Precompute constant
tan_const = np.tan(r0/4)

# Convert r -> (x,y) on the unit circle
def x_of_r(r):
    return np.sin(r)

def y_of_r(r):
    return np.cos(r)

# Closed‐form solution for delta = t - T
def r_formula(delta):
    return 4 * np.arctan(tan_const * np.exp(-delta))

# Piecewise r(t;T)
def r_t(t, T):
    return np.where(t < T, r0, r_formula(t - T))

# Set up figure
fig, ax = plt.subplots(figsize=(6,6))
ax.set_title("Balls Rolling Up a Hemisphere")
ax.set_xlabel("x")
ax.set_ylabel("y")

# Plot the hemisphere arc
r_vals = np.linspace(0, r0, 400)
ax.plot(x_of_r(r_vals), y_of_r(r_vals), 'k-', lw=2, label='hemisphere: x=sin(r), y=cos(r)')

# Create balls
colors = plt.cm.viridis(np.linspace(0,1,num_balls))
balls = []
for c in colors:
    b = plt.Circle((x_of_r(r0), y_of_r(r0)), 0.02, color=c)
    balls.append(b)
    ax.add_patch(b)

# Time annotation
time_text = ax.text(0.05, 0.95, '', transform=ax.transAxes)

# **New**: square plotting window and equal aspect
ax.set_aspect('equal')
ax.set_xlim(0, 1.1)
ax.set_ylim(0.0, 1.1)

# Init
def init():
    for b in balls:
        b.center = (x_of_r(r0), y_of_r(r0))
    time_text.set_text('')
    return balls + [time_text]

# Animate
def animate(frame):
    t = frame / 30.0
    for i, T in enumerate(T_values):
        rpos = r_t(t, T)
        balls[i].center = (x_of_r(rpos), y_of_r(rpos))
    time_text.set_text(f"t = {t:.1f} s")
    return balls + [time_text]

# Build and save
ani = FuncAnimation(fig, animate, frames=30*20, init_func=init,
                    interval=33, blit=True)

writer = PillowWriter(fps=30)
ani.save('hemisphere_rollup_cartesian.gif', writer=writer)
print("Saved hemisphere_rollup_cartesian.gif")

plt.legend()
plt.show()
