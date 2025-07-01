import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

# Load data
data_4hz = pd.read_csv("LihiriBevan2020-fig1H-4Hz.csv", header=None, names=["x", "y"])
data_20hz = pd.read_csv("LihiriBevan2020-fig1H-20Hz.csv", header=None, names=["x", "y"])

# Invert y (since it came from pixel coordinates)
data_4hz["y"] = data_4hz["y"]
data_20hz["y"] = data_20hz["y"]

# Function to split pre/post
def split_data(df, x_thresh):
    pre = df[df["x"] < x_thresh]["y"].reset_index(drop=True)
    post = df[df["x"] > x_thresh]["y"].reset_index(drop=True)
    return pre, post

# Use appropriate x thresholds
pre_4hz, post_4hz = split_data(data_4hz, 0.5)
pre_20hz, post_20hz = split_data(data_20hz, 4.0)

# Build combined DataFrame
df = pd.DataFrame({
    "Subject": list(range(1, 7)) * 4,
    "Frequency (Hz)": pd.concat([pre_4hz, post_4hz, pre_20hz, post_20hz], ignore_index=True),
    "Time": ['pre']*6 + ['post']*6 + ['pre']*6 + ['post']*6,
    "Stim": ['4 Hz']*12 + ['20 Hz']*12
})

# import pdb
# pdb.set_trace()

# Plotting
sns.set(style='whitegrid')
fig, axes = plt.subplots(1, 2, figsize=(4, 6), sharey=True)
palette = {'pre': 'gray', 'post': 'magenta'}

for i, stim in enumerate(['4 Hz', '20 Hz']):
    ax = axes[i]
    stim_df = df[df['Stim'] == stim]
    
    # Boxplot
    sns.boxplot(data=stim_df, x='Time', y='Frequency (Hz)', palette=palette, ax=ax, width=0.6, showfliers=False)

    # Individual subject lines
    for subj in range(1, 7):
        subj_df = stim_df[stim_df['Subject'] == subj]
        ax.plot(['pre', 'post'], subj_df['Frequency (Hz)'].values, color='gray', marker='o', linewidth=1)

    ax.set_title(stim)
    ax.set_ylim(0, 35)
    ax.set_xlabel('')

    # Add asterisk for significance
    y_max = stim_df['Frequency (Hz)'].max() + 2
    ax.plot([0, 1], [y_max, y_max], color='black', linewidth=1)
    ax.text(0.5, y_max + 1, '*', ha='center', fontsize=16)

plt.tight_layout()
plt.savefig("LahiriBevan2020-figure1H-recreated.png", dpi=300, bbox_inches='tight')
plt.show()
