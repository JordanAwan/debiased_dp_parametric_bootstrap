import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap, BoundaryNorm

n_samples = [100, 200, 500, 1000, 2000]
epsilons = [0.1, 0.3, 1.0, 3.0, 10.0]

# Coverage data for bootstrap estimator
bootstrap_coverage_data = np.array([
    [0.987, 0.965, 0.977, 0.972, 0.929],
    [0.954, 0.953, 0.939, 0.912, 0.890],
    [0.911, 0.936, 0.892, 0.889, 0.884],
    [0.946, 0.890, 0.895, 0.873, 0.887],
    [0.892, 0.885, 0.875, 0.893, 0.886]
])

# Width data for bootstrap estimator
bootstrap_width_data = np.array([
    [9.506945, 10.142603, 11.0925868, 8.7596316, 3.8108543],
    [6.501721, 6.061232, 3.8878262, 1.5883926, 0.7125495], 
    [3.153237, 2.371824, 1.0570330, 0.6040271, 0.3855435],
    [1.739645, 1.293876, 0.7289715, 0.5067823, 0.3448297],
    [1.591473, 1.145516, 0.6811383, 0.4874412, 0.3386316]
])

# Coverage data for naive estimator
naive_coverage_data = np.array([
    [0.376, 0.301, 0.195, 0.155, 0.130],
    [0.260, 0.209, 0.157, 0.157, 0.178],
    [0.254, 0.285, 0.334, 0.454, 0.541],
    [0.581, 0.622, 0.686, 0.692, 0.693],
    [0.742, 0.731, 0.728, 0.683, 0.710]
])

# Width data for naive estimator
naive_width_data = np.array([
    [3.5096095, 2.8486101, 1.6700156, 0.9719991, 0.5826436],
    [2.3794989, 1.4939593, 0.7714624, 0.4652264, 0.2919145],
    [1.2074834, 0.7899211, 0.4615944, 0.3238381, 0.2276413],
    [0.9532342, 0.6859504, 0.4432774, 0.3183659, 0.2258042],
    [1.1177621, 0.7539813, 0.4631192, 0.3256209, 0.2283101]
])

repro_coverage_data = np.array([
    [0.9400000, 0.9540000, 0.9740000, 0.9790000, 0.9855072],
    [0.9600000, 0.9630000, 0.9739479, 0.9759760, 0.9734423],
    [0.9620000, 0.9720000, 0.9804728, 0.9809810, 0.9817629],
    [0.9700000, 0.9670000, 0.9800000, 0.9859860, 0.9848178],
    [0.9700000, 0.9740000, 0.9832335, 0.9870000, 0.9837563]
])
 

# Transposed data for repro samples width
repro_width_data = np.array([
    [18.585100, 18.750569, 18.053306, 15.480654, 9.747651],
    [18.418816, 15.860604, 10.211779, 3.541586, 1.275519],
    [11.403998, 6.889509, 1.755253, 0.977412, 0.617089],
    [4.451044, 2.005983, 0.964247, 0.764053, 0.539312],
    [2.663767, 1.696158, 1.044317, 0.738545, 0.529185]
])


# Create a figure with 3x2 subplots and two colorbars
fig = plt.figure(figsize=(20, 12))
gs = fig.add_gridspec(2, 4, width_ratios=[1, 1, 1, 0.08])

# Function to create custom colormap for coverage
def create_coverage_cmap():
    # Create a gradient from red to green
    # Using five colors to create a darker green for values between 0.85 and 0.9
    # colors = [(0.8, 0, 0), (0.6, 0.3, 0), (0.2, 0.7, 0), (0, 0.6, 0), (0, 0.9, 0)]  # red -> orange -> yellow-green -> dark green -> bright green
    colors = [
        (0.0, 'red'),       # 0.2
        (0.45, 'red'),      # 0.6 is 45% from 0.2 to 0.99
        (0.55, (0.6, 0.3, 0)),  # 0.7
        (0.70, (0.3, 0.6, 0)),  # 0.8
        (0.85, (0, 0.6, 0)),    # 0.9
        (1.0, (0, 1.0, 0)),     # 0.99
    ]
    # Create a colormap with evenly spaced positions
    cmap = LinearSegmentedColormap.from_list("custom_coverage", colors, N=256)
    cmap.set_over('green')
    return cmap

# Function to create green-only colormap for debiased estimator
def create_debiased_cmap():
    # Create a gradient from dark green to bright green
    colors = [(0, 0.4, 0), (0, 0.7, 0), (0, 0.9, 0)]  # dark green -> medium green -> bright green
    return LinearSegmentedColormap.from_list("debiased_coverage", colors, N=256)

# Function to create custom colormap for colorbar with red below 0.6
def create_colorbar_cmap():
    # Create a colormap that's red below 0.6 and matches the coverage colormap above 0.6
    colors = [(0.8, 0, 0), (0.8, 0, 0), (0.6, 0.3, 0), (0.2, 0.7, 0), (0, 0.6, 0), (0, 0.9, 0)]  # red -> red -> orange -> yellow-green -> dark green -> bright green
    positions = [0, 0.6, 0.6, 0.7, 0.8, 1.0]  # red from 0 to 0.6, then normal gradient
    return LinearSegmentedColormap.from_list("colorbar_coverage", list(zip(positions, colors)), N=256)

# Function to create custom colormap for width plots
def create_width_cmap():
    # Create a gradient from light to dark red
    colors = [(1, 0.8, 0.8), (0.5, 0, 0)]  # light red -> dark red
    return LinearSegmentedColormap.from_list("custom_width", colors, N=256)

# Create a custom normalization that maps values >= 0.9 to 1.0
class CoverageNorm(plt.Normalize):
    def __init__(self, vmin=None, vmax=None, clip=False):
        super().__init__(vmin, vmax, clip)
    
    def __call__(self, value, clip=None):
        # First normalize to 0-1 range
        result = super().__call__(value, clip)
        # Then map values >= 0.9 to 1.0
        if np.isscalar(value):
            if value > 0.893:
                result = 0.9
        else:
            result[value > 0.893] = 0.9
        return result

# Top row - Coverage plots
# Plot 1: Coverage for debiased estimator
ax1 = fig.add_subplot(gs[0, 0])
im1 = sns.heatmap(bootstrap_coverage_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=create_debiased_cmap(),
                  norm=CoverageNorm(vmin=0.1, vmax=1.0),
                  cbar=False)
ax1.set_ylabel('Privacy Parameter (ε)')
ax1.set_title('Coverage of 90% CI\n(Debiased Estimator)')

# Plot 2: Coverage for naive estimator
ax2 = fig.add_subplot(gs[0, 1])
im2 = sns.heatmap(naive_coverage_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=create_coverage_cmap(),
                  norm=CoverageNorm(vmin=0.1, vmax=1.0),
                  cbar=False)
ax2.set_ylabel('Privacy Parameter (ε)')
ax2.set_title('Coverage of 90% CI\n(Naive Estimator)')

# Plot 3: Coverage for repro samples
ax3 = fig.add_subplot(gs[0, 2])
im3 = sns.heatmap(repro_coverage_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=create_coverage_cmap(),
                  norm=CoverageNorm(vmin=0.1, vmax=1.0),
                  cbar=False)
ax3.set_ylabel('Privacy Parameter (ε)')
ax3.set_title('Coverage of 90% CI\n(Repro Samples)')

# Coverage colorbar
cbar_ax1 = fig.add_subplot(gs[0, 3])
# Create a custom colorbar that matches im2 but with red below 0.6
# norm = CoverageNorm(vmin=0.6, vmax=0.99)  # Extend range to show full colorbar
cbar = plt.colorbar(im2.collections[0], cax=cbar_ax1, label='Coverage Rate')
cbar.set_ticks([0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

# Bottom row - Width plots
# Custom colormap for width - dark red for larger values, light red for smaller values
width_cmap = create_width_cmap()  # Using our custom colormap

# Plot 4: Width for debiased estimator
ax4 = fig.add_subplot(gs[1, 0])
im4 = sns.heatmap(bootstrap_width_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=width_cmap,
                  cbar=False)
ax4.set_xlabel('Number of Samples (N)')
ax4.set_ylabel('Privacy Parameter (ε)')
ax4.set_title('Width of 90% CI\n(Debiased Estimator)')

# Plot 5: Width for naive estimator
ax5 = fig.add_subplot(gs[1, 1])
im5 = sns.heatmap(naive_width_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=width_cmap,
                  cbar=False)
ax5.set_xlabel('Number of Samples (N)')
ax5.set_ylabel('Privacy Parameter (ε)')
ax5.set_title('Width of 90% CI\n(Naive Estimator)')

# Plot 6: Width for repro samples
ax6 = fig.add_subplot(gs[1, 2])
im6 = sns.heatmap(repro_width_data,
                  xticklabels=n_samples,
                  yticklabels=epsilons,
                  annot=True,
                  fmt='.2f',
                  cmap=width_cmap,
                  cbar=False)
ax6.set_xlabel('Number of Samples (N)')
ax6.set_ylabel('Privacy Parameter (ε)')
ax6.set_title('Width of 90% CI\n(Repro Samples)')

# Width colorbar
cbar_ax2 = fig.add_subplot(gs[1, 3])
plt.colorbar(im6.collections[0], cax=cbar_ax2, label='CI Width')

plt.tight_layout()
plt.show()




