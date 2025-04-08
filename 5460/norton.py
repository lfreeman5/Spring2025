import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

# Dome profile
def h_of_r(r):
    return -(2/3) * r**1.5

# Time it takes to roll from r=1 to r=0:
t0 = 144**0.25  # = (144)^(1/4) ≈ 3.464

# r(t; T): sits at 1 until t=T, then follows reversed profile,
# reaching 0 at t = T + t0, then stays at 0.
def r_t(t, T):
    return np.where(
        t < T,
        1.0,
        np.where(
            t > T + t0,
            0.0,
            (1/144) * (T + t0 - t)**4
        )
    )

# Animation setup
fig, ax = plt.subplots(figsize=(6,4))
ax.set_xlabel('r')
ax.set_ylabel('h')
ax.set_title("Balls Rolling Up Norton's Dome")

# Plot dome outline
r_vals = np.linspace(0, 1.0, 400)
ax.plot(r_vals, h_of_r(r_vals), 'k-', lw=2, label="Norton's Dome")

# Ball delays
num_balls = 7
T_values = np.linspace(0, 4, num_balls)  # each ball waits a bit longer before rolling

# Create ball patches
balls = []
colors = plt.cm.viridis(np.linspace(0,1,num_balls))
for c in colors:
    b = plt.Circle((1.0, h_of_r(1.0)), 0.01, color=c)
    balls.append(b)
    ax.add_patch(b)

# Axes limits
ax.set_xlim(-0.05, 1.05)
ax.set_ylim(h_of_r(1.0)-0.05, 0.05)

# Time label
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    for b in balls:
        b.center = (1.0, h_of_r(1.0))
    time_text.set_text('')
    return balls + [time_text]

def animate(frame):
    t = frame / 30.0  # seconds
    for i, T in enumerate(T_values):
        rpos = r_t(t, T)
        balls[i].center = (rpos, h_of_r(rpos))
    time_text.set_text(f"t = {t:.2f} s")
    return balls + [time_text]

ani = FuncAnimation(fig, animate, frames=30*20, init_func=init,
                    interval=33, blit=True)

# Save as GIF
writer = PillowWriter(fps=30)
ani.save('nortons_dome_rollup.gif', writer=writer)
print("Saved nortons_dome_rollup.gif")

plt.legend()
plt.show()
