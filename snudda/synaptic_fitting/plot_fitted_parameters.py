import matplotlib.pyplot as plt
import numpy as np
import json

# Read the parameter data
with open('/media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/M1RH-ipsi/M1RH-ipsi_MSN_D1_20Hz_depressing.json-parameters-full.json', 'r') as f:
    data = json.load(f)

# Read the bounds data
with open('/media/psf/KTH/2025-11-13-Yvonne-data-teanalysing/Yvonne2019/CategorisedSTP/M1RH-ipsi/M1RH-ipsi_MSN_D1_20Hz_depressing.json', 'r') as f:
    bounds_data = json.load(f)

bounds = bounds_data['model_data']

# Parameter names
param_names = ['u_sobol', 'tau_r_sobol', 'tau_f_sobol', 'tau_ratio_sobol', 'cond_sobol']
display_names = ['U', 'τ_R', 'τ_F', 'τ_ratio', 'Conductivity']
bound_keys = ['U', 'tauR', 'tauF', 'tauRatio', 'cond']

# Count parameter sets (excluding 'meta' key)
n_sets = len([k for k in data.keys() if k != 'meta'])
set_indices = range(n_sets)

# Organize parameters by type
parameters = {name: [] for name in param_names}
errors = []

for i in range(n_sets):
    params = data[str(i)]['parameters']
    for j, name in enumerate(param_names):
        parameters[name].append(params[j])
    errors.append(data[str(i)]['error'])

# Create figure with subplots
fig, axes = plt.subplots(3, 2, figsize=(14, 12))
fig.suptitle('Parameter Set Analysis - M1RH-ipsi MSN D1 20Hz Depressing', fontsize=16, fontweight='bold', y=0.98)

# Flatten axes for easier indexing
axes_flat = axes.flatten()

# Plot each parameter
for idx, (param_name, display_name, bound_key) in enumerate(zip(param_names, display_names, bound_keys)):
    ax = axes_flat[idx]
    
    # Get values and bounds
    values = parameters[param_name]
    bound_min, bound_max = bounds[bound_key]
    
    # Create bar plot with error-based coloring
    colors = plt.cm.viridis(np.array(errors) / max(errors))
    bars = ax.bar(set_indices, values, color=colors, alpha=0.7, edgecolor='black', linewidth=1.5)
    
    # Add horizontal lines for bounds
    ax.axhline(y=bound_min, color='red', linestyle='--', linewidth=2, alpha=0.5, label=f'Min: {bound_min:.2e}')
    ax.axhline(y=bound_max, color='red', linestyle='--', linewidth=2, alpha=0.5, label=f'Max: {bound_max:.2e}')
    
    # Formatting
    ax.set_xlabel('Parameter Set Index', fontsize=11, fontweight='bold')
    ax.set_ylabel('Value', fontsize=11, fontweight='bold')
    ax.set_title(f'{display_name}', fontsize=13, fontweight='bold')
    ax.set_xticks(set_indices)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    ax.legend(fontsize=8, loc='best')
    
    # Use scientific notation for conductivity
    if param_name == 'cond_sobol':
        ax.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))

# Plot error values in the last subplot
ax_error = axes_flat[5]
colors_error = plt.cm.viridis(np.array(errors) / max(errors))
bars = ax_error.bar(set_indices, errors, color=colors_error, alpha=0.7, edgecolor='black', linewidth=1.5)
ax_error.set_xlabel('Parameter Set Index', fontsize=11, fontweight='bold')
ax_error.set_ylabel('Error', fontsize=11, fontweight='bold')
ax_error.set_title('Error Values', fontsize=13, fontweight='bold')
ax_error.set_xticks(set_indices)
ax_error.grid(axis='y', alpha=0.3, linestyle='--')

plt.tight_layout(rect=[0, 0.05, 1, 0.96])
plt.subplots_adjust(hspace=0.35, wspace=0.25)

# Add colorbar for error mapping
sm = plt.cm.ScalarMappable(cmap=plt.cm.viridis, norm=plt.Normalize(vmin=min(errors), vmax=max(errors)))
sm.set_array([])
cbar = fig.colorbar(sm, ax=axes_flat, orientation='horizontal', pad=0.08, fraction=0.03, aspect=40)
cbar.set_label('Error Value', fontsize=11, fontweight='bold')
plt.show()

# Print summary statistics
print("Parameter Set Summary:")
print("=" * 60)
for i in range(n_sets):
    print(f"\nSet {i}: Error = {errors[i]:.6f}")
    for param_name, display_name in zip(param_names, display_names):
        print(f"  {display_name:15s}: {parameters[param_name][i]:.6e}")

print("\n" + "=" * 60)
print(f"Best parameter set (lowest error): Set {errors.index(min(errors))} (error = {min(errors):.6f})")
print(f"Worst parameter set (highest error): Set {errors.index(max(errors))} (error = {max(errors):.6f})")

# Print metadata
print("\n" + "=" * 60)
print("Metadata:")
meta = bounds_data['meta_data']
print(f"  Input region: {meta['input_region']}")
print(f"  Cell type: {meta['cell_type']}")
print(f"  Input frequency: {meta['input_frequency']} Hz")
print(f"  Classification: {meta['classification']}")
print(f"  Number of traces: {len(meta['trace_list'])}")
