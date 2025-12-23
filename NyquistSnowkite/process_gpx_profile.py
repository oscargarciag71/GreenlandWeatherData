#!/usr/bin/env python3
"""Process a GPX track into an elevation profile with slope and curvature.

Saves a CSV with distance,elevation,slope,curvature and a PNG with three plots.

Usage example:
  python process_gpx_profile.py \
    -i TrackgesamtRoute.gpx \
    -o outputs \
    --points 2000
"""
from __future__ import annotations

import argparse
import csv
import os
from typing import List, Tuple, Optional

import gpxpy
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter


def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """Return distance in meters between two lat/lon points using haversine."""
    R = 6371000.0  # earth radius in meters
    phi1 = np.radians(lat1)
    phi2 = np.radians(lat2)
    dphi = phi2 - phi1
    dlambda = np.radians(lon2 - lon1)
    a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    return R * c


def extract_points_from_gpx(gpx_path: str) -> List[Tuple[float, float, float]]:
    """Parse GPX and return list of (lat, lon, elevation) for track points or route points.

    Elevation may be None if not present.
    """
    with open(gpx_path, "r", encoding="utf-8") as fh:
        gpx = gpxpy.parse(fh)

    points: List[Tuple[float, float, float]] = []
    # prefer tracks -> segments -> points
    if gpx.tracks:
        for tr in gpx.tracks:
            for seg in tr.segments:
                for p in seg.points:
                    points.append((p.latitude, p.longitude, p.elevation))
    elif gpx.routes:
        for route in gpx.routes:
            for p in route.points:
                points.append((p.latitude, p.longitude, p.elevation))
    else:
        raise ValueError("No tracks or routes found in GPX file")

    if not points:
        raise ValueError("GPX contained no points")

    return points


def compute_distances_and_elevations(points: List[Tuple[float, float, float]]) -> Tuple[np.ndarray, np.ndarray]:
    """Compute cumulative distances (meters) and elevation array (meters).

    - Uses haversine for horizontal distance.
    - Fills missing elevations by 1D interpolation across available points.
    - Removes consecutive duplicate positions (zero delta distance).
    """
    lats = np.array([p[0] for p in points], dtype=float)
    lons = np.array([p[1] for p in points], dtype=float)
    elevs = np.array([np.nan if p[2] is None else float(p[2]) for p in points], dtype=float)

    n = len(lats)
    distances = np.zeros(n, dtype=float)
    for i in range(1, n):
        distances[i] = distances[i - 1] + haversine(lats[i - 1], lons[i - 1], lats[i], lons[i])

    # Remove points where distance didn't increase (duplicates). Keep first occurrence.
    diff = np.diff(distances, prepend=-1.0)
    keep = diff > 0
    # Always keep the first point
    keep[0] = True

    distances = distances[keep]
    elevs = elevs[keep]

    # If all elevations are NaN, we cannot proceed
    if np.isnan(elevs).all():
        raise ValueError("All elevation values are missing in GPX file")

    # Interpolate missing elevation values in index order (safer when distances uneven)
    idx = np.arange(elevs.size)
    good = ~np.isnan(elevs)
    if not np.all(good):
        elevs[np.isnan(elevs)] = np.interp(idx[np.isnan(elevs)], idx[good], elevs[good])

    return distances, elevs


def compute_profile(distances: np.ndarray, elevations: np.ndarray, n_points: int = 1000, outdir: Optional[str] = None, base_name: str = "profile") -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate elevation with cubic spline and compute first/second derivatives.

    Returns x_smooth (m), y_smooth (m), slope (m/m), curvature (m/m^2).
    """
    if distances.size < 2:
        raise ValueError("Need at least two points to interpolate profile")

    
    # distance between sample points — report robust statistics (median + MAD)
    dx = np.diff(distances)
    med_dx = np.median(dx)
    mad_dx = np.median(np.abs(dx - med_dx))
    print(f"Median sample point distance: {med_dx:.2f} m ± MAD {mad_dx:.2f} m")

    # Apply a Hampel (median) filter to remove isolated outliers (GPS spikes)
    def hampel_filter(arr: np.ndarray, window_radius: int = 3, n_sigmas: float = 3.0) -> Tuple[np.ndarray, int]:
        """Replace outliers in `arr` using a Hampel filter.

        - `window_radius` is in samples (radius, not full window length).
        - `n_sigmas` is the number of scaled-MADs to use as threshold.

        Returns filtered array and number of replaced values.
        """
        x = arr.copy()
        n = x.size
        replaced = 0
        for i in range(n):
            left = max(0, i - window_radius)
            right = min(n, i + window_radius + 1)
            window = arr[left:right]
            med = np.median(window)
            mad = np.median(np.abs(window - med))
            if mad == 0:
                threshold = n_sigmas * 1e-9
            else:
                threshold = n_sigmas * 1.4826 * mad
            if abs(arr[i] - med) > threshold:
                x[i] = med
                replaced += 1
        return x, replaced

    # Choose a window radius in meters (for Hampel) converted to samples
    radius_m = 5.0
    window_radius = max(1, int(round(radius_m / 0.2)))
    # Cap radius to something reasonable relative to data length
    window_radius = min(window_radius, max(1, int(distances.size / 4)))
    elevations, n_replaced = hampel_filter(elevations, window_radius=window_radius, n_sigmas=3.0)
    if n_replaced:
        print(f"Hampel filter: replaced {n_replaced} outlier elevation samples (window radius {window_radius} samples)")

    # Resample elevations to an approximately uniform grid for Savitzky-Golay
    # Use approximately the original median spacing so derivatives are in m units
    if distances[-1] <= distances[0]:
        raise ValueError("Invalid distances range for resampling")

    # create uniform grid with step ~= avg_dx (at least 1 meter)
    step = max(1.0, med_dx)
    x_uniform = np.arange(distances[0], distances[-1] + 0.5 * step, step)
    y_uniform = np.interp(x_uniform, distances, elevations)

    # Choose Savitzky-Golay window length in samples (9..15), based on a target smoothing length
    target_window_m = 1500.0
    window_len = int(round(target_window_m / step))
    window_len = min(max(window_len, 9), 15)
    # ensure odd and <= data length
    if window_len % 2 == 0:
        window_len += 1
    if window_len > y_uniform.size:
        window_len = y_uniform.size if (y_uniform.size % 2 == 1) else (y_uniform.size - 1)
        if window_len < 3:
            window_len = 3

    polyorder = 3
    if window_len <= polyorder:
        window_len = polyorder + 2 if (polyorder + 2) % 2 == 1 else polyorder + 3

    dx_uniform = x_uniform[1] - x_uniform[0]
    print(f"Savitzky-Golay: window {window_len} samples (~{window_len*dx_uniform:.1f} m), polyorder {polyorder}, dx {dx_uniform:.2f} m")

    # apply Savitzky-Golay to get smoothed elevation and derivatives
    y_smooth = savgol_filter(y_uniform, window_length=window_len, polyorder=polyorder, mode="interp")
    slope = savgol_filter(y_uniform, window_length=window_len, polyorder=polyorder, deriv=1, delta=dx_uniform, mode="interp")
    curvature = savgol_filter(y_uniform, window_length=window_len, polyorder=polyorder, deriv=2, delta=dx_uniform, mode="interp")

    # Fit a spline through y_smooth over the uniform grid to compute curvature kappa
    cs = CubicSpline(x_uniform, y_smooth, bc_type="natural")
    dy_dx = cs(x_uniform, 1)
    d2y_dx2 = cs(x_uniform, 2)

    # Apply a gaussian filter to derivatives to reduce high-frequency noise
    dy_dx_gaussian = gaussian_filter1d(dy_dx, sigma=5)
    d2y_dx2_gaussian = gaussian_filter1d(d2y_dx2, sigma=5)

    kappa = np.abs(d2y_dx2_gaussian) / (1 + dy_dx_gaussian**2)**(3/2)

    # use the uniform grid as the default output x-axis
    x_smooth = x_uniform

    # --- Compute larger-window median slopes (robust) over 50 m and 100 m ---
    def rolling_median_and_mad(values: np.ndarray, radius_samples: int) -> Tuple[np.ndarray, np.ndarray]:
        n = values.size
        med = np.empty(n, dtype=float)
        mad = np.empty(n, dtype=float)
        for i in range(n):
            left = max(0, i - radius_samples)
            right = min(n, i + radius_samples + 1)
            w = values[left:right]
            m = np.median(w)
            mmad = np.median(np.abs(w - m))
            med[i] = m
            mad[i] = mmad
        return med, mad

    # convert meter windows to sample radii
    radius50 = max(1, int(round(50.0 / dx_uniform)))
    radius100 = max(1, int(round(100.0 / dx_uniform)))
    slope_med50, slope_mad50 = rolling_median_and_mad(slope, radius50)
    slope_med100, slope_mad100 = rolling_median_and_mad(slope, radius100)

    # Print summary robust stats
    print(f"Slope (50 m median) — median: {np.median(slope_med50):.5g}, MAD: {np.median(slope_mad50):.5g}")
    print(f"Slope (100 m median) — median: {np.median(slope_med100):.5g}, MAD: {np.median(slope_mad100):.5g}")


    # Return processed 50 m profile: centers, median elevations, slope, curvature
    return x_smooth, y_smooth, slope, curvature, x_uniform, dy_dx_gaussian, d2y_dx2_gaussian, kappa

    # apply a gaussian filter to smooth dependent on measurement accuracy
    # delta_messurement = 10.0  # meters, assumed GPS accuracy
    # sigma = delta_messurement / (avg_dx*2.355)  # convert FWHM to sigma in samples
    # elevations = gaussian_filter1d(elevations, sigma=sigma)

    # cs = CubicSpline(distances, elevations, bc_type="natural")
    # # Resample to n_points evenly spaced points
    # x_smooth = np.arange(distances[0], distances[-1], 5.0)  # every 5 meters
    # # x_smooth = distances
    # y_smooth = cs(x_smooth)
    # slope = cs(x_smooth, 1)
    # curvature = cs(x_smooth, 2)

    # return x_smooth, y_smooth, slope, curvature


def save_csv(path: str, x: np.ndarray, y: np.ndarray, slope: np.ndarray, curvature: np.ndarray) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["distance_m", "elevation_m", "slope_m_per_m", "curvature_m_per_m2"])
        for a, b, c, d in zip(x, y, slope, curvature):
            writer.writerow([f"{a:.3f}", f"{b:.3f}", f"{c:.6f}", f"{d:.9f}"])


def plot_profile(path_png: str, x: np.ndarray, y: np.ndarray, slope: np.ndarray, curvature: np.ndarray) -> None:
    os.makedirs(os.path.dirname(path_png) or ".", exist_ok=True)
    fig, axs = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    axs[0].plot(x, y, color="tab:blue")
    axs[0].set_ylabel("Elevation (m)")
    axs[0].grid(alpha=0.3)

    axs[1].plot(x, slope * 100.0, color="tab:orange")
    axs[1].set_ylabel("Slope (% = m/ m * 100)")
    axs[1].grid(alpha=0.3)

    axs[2].plot(x, curvature, color="tab:green")
    axs[2].set_ylabel("Curvature (m/m^2)")
    axs[2].set_xlabel("Distance (m)")
    axs[2].grid(alpha=0.3)

    fig.suptitle("Elevation profile and derivatives")
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(path_png, dpi=200)
    plt.close(fig)


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Process GPX to elevation profile, slope, curvature")
    p.add_argument("-i", "--input", required=True, help="Path to GPX file")
    p.add_argument("-o", "--outdir", default="outputs", help="Output directory for CSV and PNG")
    p.add_argument("-n", "--points", type=int, default=1000, help="Number of interpolated points")
    return p


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()

    gpx_path = args.input
    outdir = args.outdir
    n_points = args.points

    base_name = os.path.splitext(os.path.basename(gpx_path))[0]
    csv_path = os.path.join(outdir, f"{base_name}_profile.csv")
    png_path = os.path.join(outdir, f"{base_name}_profile.png")

    print(f"Reading GPX: {gpx_path}")
    points = extract_points_from_gpx(gpx_path)

    distances, elevations = compute_distances_and_elevations(points)
    print(f"Read {distances.size} points, total distance {distances[-1]:.1f} m")

    x, y, slope, curvature = compute_profile(distances, elevations, n_points=n_points, outdir=outdir, base_name=base_name)

    print(f"Saving CSV: {csv_path}")
    save_csv(csv_path, x, y, slope, curvature)

    print(f"Saving PNG: {png_path}")
    plot_profile(png_path, x, y, slope, curvature)

    print("Done.")


if __name__ == "__main__":
    main()
