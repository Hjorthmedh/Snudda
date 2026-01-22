import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.stats import pearsonr

# Read the parameter data
with open('/media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/M1RH-ipsi/M1RH-ipsi_MSN_D1_20Hz_depressing.json-parameters-full.json', 'r') as f:
    data = json.load(f)

# Read the bounds data
with open('/media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/M1RH-ipsi/M1RH-ipsi_MSN_D1_20Hz_depressing.json', 'r') as f:
    bounds_data = json.load(f)

# Parameter names
param_names = ['u_sobol', 'tau_r_sobol', 'tau_f_sobol', 'tau_ratio_sobol', 'cond_sobol']
display_names = ['U', 'τ_R', 'τ_F', 'τ_ratio', 'Conductivity']

# Count parameter sets
n_sets = len([k for k in data.keys() if k != 'meta'])

# Organize parameters and errors
parameters = {name: [] for name in param_names}
errors = []

for i in range(n_sets):
    params = data[str(i)]['parameters']
    for j, name in enumerate(param_names):
        parameters[name].append(params[j])
    errors.append(data[str(i)]['error'])

# Convert to numpy arrays for easier manipulation
errors = np.array(errors)
param_arrays = {name: np.array(parameters[name]) for name in param_names}

# Define threshold for "close to minimum" - let's use bottom 30% of errors
error_threshold_percentile = 30
error_threshold = np.percentile(errors, error_threshold_percentile)
best_indices = errors <= error_threshold

print(f"Analyzing {np.sum(best_indices)} parameter sets with error <= {error_threshold:.4f}")
print(f"(Bottom {error_threshold_percentile}% of all {n_sets} sets)")

# Create scatter plot matrix for best parameters
n_params = len(param_names)
fig, axes = plt.subplots(n_params, n_params, figsize=(16, 16))
fig.suptitle(f'Parameter Correlations (Best {error_threshold_percentile}% by Error)', 
             fontsize=18, fontweight='bold', y=0.995)

for i, (param_i, name_i) in enumerate(zip(param_names, display_names)):
    for j, (param_j, name_j) in enumerate(zip(param_names, display_names)):
        ax = axes[i, j]
        
        if i == j:
            # Diagonal: histogram of parameter values
            ax.hist(param_arrays[param_i][best_indices], bins=15, 
                   color='steelblue', alpha=0.7, edgecolor='black')
            ax.set_ylabel('Count', fontsize=9)
            if i == 0:
                ax.set_title(name_i, fontsize=11, fontweight='bold')
        else:
            # Off-diagonal: scatter plots
            x_data = param_arrays[param_j][best_indices]
            y_data = param_arrays[param_i][best_indices]
            
            # Color by error value
            colors = errors[best_indices]
            scatter = ax.scatter(x_data, y_data, c=colors, 
                               cmap='viridis', alpha=0.6, 
                               s=80, edgecolor='black', linewidth=0.5)
            
            # Calculate correlation
            if len(x_data) > 1:
                corr, p_value = pearsonr(x_data, y_data)
                # Add correlation text
                ax.text(0.05, 0.95, f'r={corr:.3f}', 
                       transform=ax.transAxes, fontsize=9,
                       verticalalignment='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        # Labels
        if j == 0:
            ax.set_ylabel(name_i, fontsize=10, fontweight='bold')
        if i == n_params - 1:
            ax.set_xlabel(name_j, fontsize=10, fontweight='bold')
        
        # Rotate tick labels for readability
        ax.tick_params(axis='both', labelsize=8)
        if i < n_params - 1:
            ax.set_xticklabels([])
        if j > 0:
            ax.set_yticklabels([])
        
        ax.grid(alpha=0.3, linestyle='--')

# Add colorbar
plt.tight_layout(rect=[0, 0.02, 0.98, 0.99])
cbar_ax = fig.add_axes([0.99, 0.15, 0.01, 0.7])
sm = plt.cm.ScalarMappable(cmap='viridis', 
                           norm=plt.Normalize(vmin=errors[best_indices].min(), 
                                            vmax=errors[best_indices].max()))
sm.set_array([])
cbar = fig.colorbar(sm, cax=cbar_ax)
cbar.set_label('Error Value', fontsize=11, fontweight='bold', rotation=270, labelpad=20)

plt.show()

# Print correlation matrix
print("\n" + "=" * 80)
print("Correlation Matrix (Pearson r) for Best Parameter Sets:")
print("=" * 80)
print(f"{'':12s}", end='')
for name in display_names:
    print(f"{name:>12s}", end='')
print()

for i, name_i in enumerate(display_names):
    print(f"{name_i:12s}", end='')
    for j, param_j in enumerate(param_names):
        if i == j:
            print(f"{'1.000':>12s}", end='')
        else:
            param_i = param_names[i]
            x = param_arrays[param_j][best_indices]
            y = param_arrays[param_i][best_indices]
            if len(x) > 1:
                corr, _ = pearsonr(x, y)
                print(f"{corr:>12.3f}", end='')
            else:
                print(f"{'N/A':>12s}", end='')
    print()

# Print summary statistics for best parameters
print("\n" + "=" * 80)
print("Summary Statistics for Best Parameter Sets:")
print("=" * 80)
for param_name, display_name in zip(param_names, display_names):
    values = param_arrays[param_name][best_indices]
    print(f"\n{display_name}:")
    print(f"  Mean:   {np.mean(values):.6e}")
    print(f"  Median: {np.median(values):.6e}")
    print(f"  Std:    {np.std(values):.6e}")
    print(f"  Min:    {np.min(values):.6e}")
    print(f"  Max:    {np.max(values):.6e}")
