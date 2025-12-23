# %%
"""Analysis script for sampling kite angle."""
from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import CubicSpline

# Path to GPX used in the original script
GPX_FILE = r"C:\Users\SchwarzN\OneDrive - Université de Fribourg\Private\2026_Greenland\RouteAnalysis\NyquistSnowkite\TrackgesamtRoute.gpx"


def load_process_module() -> object:
    """Dynamically load local process_gpx_profile.py as a module.

    This makes the analysis script robust to being executed from another cwd.
    """
    path_utils = "C:/Users/SchwarzN/OneDrive - Université de Fribourg/Private/2026_Greenland/RouteAnalysis/NyquistSnowkite"
    mod_path = path_utils + "/process_gpx_profile.py"
    if not os.path.exists(mod_path):
        raise FileNotFoundError(f"process_gpx_profile.py not found at {mod_path}")

    spec = importlib.util.spec_from_file_location("process_gpx_profile", str(mod_path))
    module = importlib.util.module_from_spec(spec)
    loader = spec.loader
    assert loader is not None
    loader.exec_module(module)
    return module

# Calculate curvature using the proper formula
# kappa = |y''| / (1 + (y')^2)^(3/2)
def compute_curvature_proper(slope, curvature_raw):
    """Convert raw curvature to proper mathematical curvature."""
    return np.abs(curvature_raw) / (1 + slope**2)**(3/2)

def main(gpx_path: str | None = None, show_plot: bool = True, outdir: str | None = None) -> None:
    mod = load_process_module()
    gpx_path = gpx_path or GPX_FILE
    outdir = outdir or str(Path(gpx_path).resolve().parent / "outputs")

    # Reuse robust extraction and interpolation from process_gpx_profile
    points = mod.extract_points_from_gpx(gpx_path)
    distances, elevations = mod.compute_distances_and_elevations(points)

    delta_elevations = np.diff(elevations)
    avg_delta_elevation = np.mean(np.abs(delta_elevations))
    std_delta_elevation = np.std(np.abs(delta_elevations))
    print(f"Average elevation change between points: {avg_delta_elevation:.2f} m ± {std_delta_elevation:.2f} m")

    x_smooth, y_smooth, slope, curvature, x_uniform, dy_dx_gaussian, d2y_dx2_gaussian, kappa_spline = mod.compute_profile(distances, elevations, n_points=len(points)/1)


    cs = CubicSpline(x_smooth, y_smooth, bc_type="natural")
    x_spline = np.linspace(x_smooth[0], x_smooth[-1], 25000)
    y_spline = cs(x_spline)

    dx_spline = np.diff(x_spline)
    med_dx_spline = np.median(dx_spline)
    mad_dx_spline = np.median(np.abs(dx_spline - med_dx_spline))
    print(f"Median spline sample point distance: {med_dx_spline:.2f} m ± MAD {mad_dx_spline:.2f} m")
    dy_spline = np.diff(y_spline)
    med_dy_spline = np.median(np.abs(dy_spline))
    mad_dy_spline = np.median(np.abs(np.abs(dy_spline) - med_dy_spline))
    print(f"Median spline elevation change between points: {med_dy_spline:.2f} m ± MAD {mad_dy_spline:.2f} m")

    y_spline_smooth = gaussian_filter1d(y_spline, sigma=50)
    delta_y_spline_smooth = np.diff(y_spline_smooth)
    med_delta_y_spline_smooth = np.median(np.abs(delta_y_spline_smooth))
    mad_delta_y_spline_smooth = np.median(np.abs(np.abs(delta_y_spline_smooth) - med_delta_y_spline_smooth))
    print(f"Median smoothed spline elevation change between points: {med_delta_y_spline_smooth:.2f} m ± MAD {mad_delta_y_spline_smooth:.2f} m")

    # make the first and second derivative from the spline (y_spline over x_spline) while takeing a window of 15 points to each side around a point

    delta_y_smooth = np.diff(y_smooth)
    avg_delta_y_smooth = np.mean(np.abs(delta_y_smooth))
    std_delta_y_smooth = np.std(np.abs(delta_y_smooth))
    print(f"Average elevation change between smoothed points: {avg_delta_y_smooth:.2f} m ± {std_delta_y_smooth:.2f} m")

    kappa = compute_curvature_proper(slope, curvature)

    # Save CSV and PNG via the module helpers
    base = Path(gpx_path).stem
    csv_path = os.path.join(outdir, f"{base}_profile.csv")
    png_path = os.path.join(outdir, f"{base}_profile.png")
    mod.save_csv(csv_path, x_smooth, y_smooth, slope, curvature)
    mod.plot_profile(png_path, x_smooth, y_smooth, slope, curvature)

    # Export x_spline and y_spline as a CSV file
    spline_csv_path = os.path.join(outdir, f"{base}_spline.csv")
    np.savetxt(spline_csv_path, np.column_stack((x_spline, y_spline, y_spline_smooth)), delimiter=",", header="x_spline,y_spline,y_spline_smoothed", comments="")

    # Also produce an interactive matplotlib figure similar to original
    fig, axes = plt.subplots(4, 1, figsize=(12, 13))
    axes[0].plot(x_smooth, y_smooth, "b-", label="Height Profile")
    axes[0].scatter(distances, elevations, c="r", s=10, alpha=0.5, label="Original Points")
    axes[0].plot(x_spline, y_spline, 'g--', label='Cubic Spline Fit')
    axes[0].set_ylabel("Elevation (m)")
    axes[0].set_title("Height Profile")
    axes[0].set_xlim(x_smooth[0], x_smooth[-1])
    axes[0].set_ylim(0, 3000)
    axes[0].set_xlim(201000, 202000)
    axes[0].set_ylim(2405, 2420)
    axes[0].legend()
    axes[0].grid()

    axes[1].plot(x_smooth, slope, "g-")
    axes[1].plot(x_uniform, dy_dx_gaussian, "c--", label="Spline Derivative")
    axes[1].set_ylabel("Slope (m/m)")
    axes[1].set_title("First Derivative (Gradient)")
    axes[1].set_xlim(201000, 202000)
    axes[1].set_ylim(-0.05, 0.05)
    axes[1].grid()

    axes[2].plot(x_smooth, curvature, "r-")
    axes[2].plot(x_uniform, d2y_dx2_gaussian, "y--", label="Spline Second Derivative")
    axes[2].set_ylabel("Curvature (1/m)")
    axes[2].set_xlabel("Distance (m)")
    axes[2].set_title("Second Derivative")
    axes[2].grid()

    axes[3].plot(x_smooth, kappa, "m-")
    axes[3].plot(x_uniform, kappa_spline, "k--", label="Spline Curvature")
    axes[3].set_ylabel("Curvature (m⁻¹)")
    axes[3].set_xlabel("Distance (m)")
    axes[3].set_title("Curvature")
    axes[3].grid()

    plt.tight_layout()
    if show_plot:
        plt.show()
    else:
        plt.close(fig)

    total_distance = distances[-1]
    elevation_gain = np.sum(np.diff(elevations)[np.diff(elevations) > 0])
    print(f"Total distance: {total_distance:.2f} m")
    print(f"Total elevation gain: {elevation_gain:.2f} m")

    return x_smooth, slope, curvature, kappa


if __name__ == "__main__":
    x_smooth, slope, curvature, kappa = main()
# %%
# Nyquist analysis
# v/f_s <= R/2 where R = 1/kappa and v/f_s = distance between samples

# mask  kappa to avoid division by zero
kappa_masked = np.where(kappa > 1e-6, kappa, 1e-6)
route_point_distance = 0.5 / kappa_masked  # meters

# Assume a measurement accuracy of 10 meters: > set all distances below that to zero
route_point_distance = np.where(route_point_distance > 10.0, route_point_distance, 0.0)

avg_route_point_distance = np.mean(np.abs(route_point_distance))
std_route_point_distance = np.std(np.abs(route_point_distance))
print(f"Average route point distance between points: {avg_route_point_distance:.2f} m ± {std_route_point_distance:.2f} m")

# interpolate the route point distance with reduced points for better plotting
x_smooth_reduced = np.linspace(x_smooth[0], x_smooth[-1], num=100)
route_point_distance_reduced = np.interp(x_smooth_reduced, x_smooth, route_point_distance)


# # apply a gaussian filter to smooth the route point distance
# route_point_distance_smooth = gaussian_filter1d(route_point_distance, sigma=4)
# route_point_distance = route_point_distance_smooth



# Plot route point distance
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x_smooth_reduced/1e3, route_point_distance_reduced, "b--")
ax.plot(x_smooth/1e3, route_point_distance, "b-", alpha=0.3)
ax.set_ylabel("Sample point distance on route [m]")
ax.set_xlabel("Distance [km]")
ax.set_title("Nyquist Analysis of Route Sampling")
ax.grid()
# ax.set_ylim(0, np.percentile(route_point_distance_reduced, 70))  # limit y-axis for better visibility
# ax.set_ylim(0, 2000)
plt.show()