import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['xtick.labelsize'] = 5
# Assuming your data is in a CSV file named 'exe-time.csv'
# If your data is in a DataFrame, you can skip the read_csv step

w21 = [3.3,6.4,24,97,196, 477]
w22 = [3.4, 6.6,23, 93,194, 518]
df = pd.read_csv('exe-time.csv')

# Convert max_mem to gigabytes
df['max_mem_gb'] = df['max_mem'] / 1024 / 1024

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(20, 5))

# Plot for wall_clock
for panel, group in df.groupby('panel'):
    axes[0].plot(group['tool'], group['wall_clock'], label=f'Panel {panel}')
    
axes[0].set_xticks(df['tool'].unique())
axes[0].set_xlabel('K')
axes[0].set_ylabel('Wall Clock Time (s)')
axes[0].set_title('ksmem computing  Wall Clock Time vs Tool')
axes[0].legend()

# Plot for max_mem
for panel, group in df.groupby('panel'):
    axes[1].plot(group['tool'], group['max_mem_gb'], label=f'Panel {panel}')
    
axes[1].set_xticks(df['tool'].unique())
axes[1].set_xlabel('K')
axes[1].set_ylabel('Max Memory (GB)')
axes[1].set_title('ksmem computing Max Memory')
axes[1].legend()

for panel, group in df.groupby('panel'):
    axes[2].plot(group['tool'], group['w'], label=f'Panel {panel}')

axes[2].set_xticks(df['tool'].unique())  # Set x-ticks to unique tool values
axes[2].set_xlabel('K')
axes[2].set_ylabel('Memory (MB)')
axes[2].set_title('ksmem output size')
axes[2].legend()
#fig.title("ksmem computing")
# Adjust layout for better spacing
plt.tight_layout()

# Save the combined plot to a PDF file
plt.savefig('exe.pdf',dpi=500)

# Show the combined plot
#plt.show()

df = pd.read_csv('make-time.csv')

# Convert max_mem to gigabytes
df['max_mem_gb'] = df['max_mem'] / 1024 / 1024

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(20, 5))

# Plot for wall_clock
for panel, group in df.groupby('panel'):
    axes[0].plot(group['tool'], group['wall_clock'], label=f'Panel {panel}')

axes[0].set_xticks(df['tool'].unique())
axes[0].set_xlabel('K')
axes[0].set_ylabel('Wall Clock Time (s)')
axes[0].set_title('building Wall Clock Time')
axes[0].legend()

# Plot for max_mem
for panel, group in df.groupby('panel'):
    axes[1].plot(group['tool'], group['max_mem_gb'], label=f'Panel {panel}')

axes[1].set_xticks(df['tool'].unique())
axes[1].set_xlabel('K')
axes[1].set_ylabel('Max Memory (GB)')
axes[1].set_title('building Max Memory')
axes[1].legend()

for panel, group in df.groupby('panel'):
    axes[2].plot(group['tool'], group['w'], label=f'Panel {panel}')

axes[2].set_xticks(df['tool'].unique())  # Set x-ticks to unique tool values
axes[2].set_xlabel('K')
axes[2].set_ylabel('Memory (MB)')
axes[2].set_title('index size')
axes[2].legend()
#fig.title("building")
# Adjust layout for better spacing
plt.tight_layout()

# Save the combined plot to a PDF file
plt.savefig('make.pdf', dpi=500)

# Show the combined plot
#plt.show()


