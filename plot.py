import pandas as pd
import matplotlib.pyplot as plt

# Load CSV data
data = pd.read_csv("output.csv")

# Extract columns
radius = data["r"]
y1 = data["alpha_norm(r)"]
y2 = data["psi(r)"]
y3 = data["u(r)"]

# Create subplots
# plt.figure(figsize=(10, 6))

# Y1 and Y2 vs Time
plt.plot(radius, y1, label=r"Lapse $\alpha$", color='blue')
plt.plot(radius, y2, label=r"Psi $\Psi$", color='red')
plt.plot(radius, y3, label=r"Scalar $U$", color='purple')
plt.title("Scalar Bubble")
plt.xlabel("Radius")
plt.ylabel("Fields")
plt.legend()
plt.grid(True)
plt.ylim([-0.1,1.5])



# Show the plots
plt.show()
