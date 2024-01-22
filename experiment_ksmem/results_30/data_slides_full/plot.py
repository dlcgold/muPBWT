import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["xtick.labelsize"] = 5
# Assuming your data is in a CSV file named 'exe-time.csv'
# If your data is in a DataFrame, you can skip the read_csv step
# Create a list of color maps to get distinct colors
color_maps = ["tab10", "tab20", "tab20b", "tab20c"]

# Initialize an empty list to store colors
colors = []

# Iterate through color maps and get colors
for cmap in color_maps:
    colors.extend(plt.get_cmap(cmap).colors)

# Select the first 22 distinct colors
colors = colors[:22]
w21 = [3.3, 6.4, 24, 97, 196, 477]
w22 = [3.4, 6.6, 23, 93, 194, 518]
df = pd.read_csv("exe-time.csv")
# df = df.sort_values(by=["panel"])
# Convert max_mem to gigabytes
df["max_mem_gb"] = df["max_mem"] / 1024 / 1024

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# Plot for wall_clock
for i, (panel, group) in enumerate(df.groupby("panel")):
    axes[0].plot(group["tool"], group["wall_clock"], label=f"{panel}", color=colors[i])

axes[0].set_xticks(df["tool"].unique())
axes[0].set_xlabel("K")
axes[0].set_ylabel("Wall Clock Time (s)")
axes[0].set_title("ksmem computing  Wall Clock Time vs Tool")
# axes[0].legend()

# Plot for max_mem
for i, (panel, group) in enumerate(df.groupby("panel")):
    axes[1].plot(group["tool"], group["max_mem_gb"], label=f"{panel}", color=colors[i])

axes[1].set_xticks(df["tool"].unique())
axes[1].set_xlabel("K")
axes[1].set_ylabel("Max Memory (GB)")
axes[1].set_title("ksmem computing Max Memory")
# axes[1].legend()

# fig.title("ksmem computing")
# Adjust layout for better spacing
labels = [i for i in range(1, 23)]
fig.legend(
    loc="lower center", ncol=22, fontsize=8, bbox_to_anchor=(0.5, -0.06), labels=labels
)  # Adjust position as needed # Adjust position as needed
plt.tight_layout()

# Save the combined plot to a PDF file
plt.savefig("exe.pdf", dpi=500, bbox_inches="tight")

# Show the combined plot
# plt.show()

df = pd.read_csv("make-time.csv")

# Convert max_mem to gigabytes
df["max_mem_gb"] = df["max_mem"] / 1024 / 1024

# Create subplots
fig, axes = plt.subplots(1, 2, figsize=(15, 5))

# Plot for wall_clock
for i, (panel, group) in enumerate(df.groupby("panel")):
    axes[0].plot(group["tool"], group["wall_clock"], label=f"{panel}", color=colors[i])

axes[0].set_xticks(df["tool"].unique())
axes[0].set_xlabel("K")
axes[0].set_ylabel("Wall Clock Time (s)")
axes[0].set_title("building Wall Clock Time")
# axes[0].legend()

# Plot for max_mem
for i, (panel, group) in enumerate(df.groupby("panel")):
    axes[1].plot(group["tool"], group["max_mem_gb"], label=f"{panel}", color=colors[i])

axes[1].set_xticks(df["tool"].unique())
axes[1].set_xlabel("K")
axes[1].set_ylabel("Max Memory (GB)")
axes[1].set_title("building Max Memory")
# axes[1].legend()


# fig.title("building")
# Adjust layout for better spacing
labels = [i for i in range(1, 23)]
fig.legend(
    loc="lower center", ncol=22, fontsize=8, bbox_to_anchor=(0.5, -0.06), labels=labels
)  # Adjust position as needed
plt.tight_layout()
# Save the combined plot to a PDF file
plt.savefig("make.pdf", dpi=500, bbox_inches="tight")

# Show the combined plot
# plt.show()
