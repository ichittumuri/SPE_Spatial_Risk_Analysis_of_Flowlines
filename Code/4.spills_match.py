# =============================================================================
# 1. Setup & Imports
# =============================================================================
import os
import warnings
import pandas as pd
import geopandas as gpd
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d
import matplotlib.colors as mcolors
import contextily as ctx
from shapely.geometry import box
from matplotlib.colors import ListedColormap
import matplotlib
from matplotlib.colors import LinearSegmentedColormap


# =============================================================================
# 2. Load Data
# =============================================================================
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
warnings.filterwarnings("ignore", category=RuntimeWarning)

flowlines_gdf = gpd.read_file("interpolated_clean_flowlines.geojson")
spills_gdf = gpd.read_file("clean_spills.geojson")

# =============================================================================
# 3. Match spills to flowlines
# =============================================================================
def match_spills_to_flowlines(spills_gdf, flowlines_gdf, output_file):

    matched_spills      = []
    matched_spill_count = 0
    missing_geom_spill  = 0
    no_op_spill         = 0
    no_match_spill      = 0

    projected_crs = spills_gdf.estimate_utm_crs()

    for idx, spill in spills_gdf.iterrows():
        pt = spill.geometry
        if pt is None or pt.is_empty:
            print(f"[Spill {idx}] Missing geometry – skipping.")
            missing_geom_spill += 1
            continue

        op = spill.get("operator_name", "")
        if pd.isnull(op) or not op.strip():
            print(f"[Spill {idx}] No operator_name – skipping.")
            no_op_spill += 1
            continue

        op_name = op.strip().lower()
        candidates = flowlines_gdf[
            flowlines_gdf["operator_name"]
                .astype(str)
                .str.strip()
                .str.lower()
                .eq(op_name)
        ].copy()
        candidates = candidates[candidates.geometry.notnull()]

        if candidates.empty:
            print(f"[Spill {idx}] No flowlines for operator “{op}”.")
            no_match_spill += 1
            continue

        dists = candidates.geometry.distance(pt).dropna()
        if dists.empty:
            print(f"[Spill {idx}] All distances NaN – skipping.")
            no_match_spill += 1
            continue

        nearest_idx  = dists.idxmin()
        min_dist     = dists.min()
        nearest_line = flowlines_gdf.loc[nearest_idx]

        _, nearest_pt = nearest_points(pt, nearest_line.geometry)

        matched_spill_count += 1
        print(f"[Match {matched_spill_count}] spill {idx} → flowline {nearest_idx} at {min_dist:.2f} m")

        new_spill = spill.copy()
        new_spill["match_point"]          = nearest_pt
        new_spill["match_distance_m"]     = min_dist
        new_spill["matched_flowline_idx"] = nearest_idx
        matched_spills.append(new_spill)

    print("\n=== Spills Matching Complete ===")
    print(f"Total spills:              {len(spills_gdf)}")
    print(f"Matched spills:            {matched_spill_count}")
    print(f"Missing geometry spills:   {missing_geom_spill}")
    print(f"No Operator Name spills:   {no_op_spill}")
    print(f"No line match spills:      {no_match_spill}")

    matched_spills_gdf = gpd.GeoDataFrame(matched_spills)

    matched_spills_gdf = matched_spills_gdf.rename(columns={"geometry": "orig_geometry"})
    matched_spills_gdf["geometry"] = matched_spills_gdf["match_point"]
    matched_spills_gdf = matched_spills_gdf.set_geometry("geometry")

    matched_spills_gdf = matched_spills_gdf.set_crs(projected_crs, allow_override=True)
    matched_spills_gdf = matched_spills_gdf.to_crs(epsg=4326)
    matched_spills_gdf = matched_spills_gdf.drop(columns=["orig_geometry", "match_point"])

    matched_spills_gdf.to_file(output_file, driver="GeoJSON")
    print(f"Wrote {output_file} in EPSG:4326.")

    return matched_spills_gdf


matched_spills_gdf = match_spills_to_flowlines(
    spills_gdf,
    flowlines_gdf,
    "updated_spills.geojson"
)

# =============================================================================
# 4. Plot simple histograms
# =============================================================================
def plot_match_distance_hist(matched_spills_gdf, bins=30, zoom_max=None, figsize=(10, 6), title=None): 
    distances = matched_spills_gdf["match_distance_m"].dropna()

    plt.figure(figsize=figsize)
    plt.hist(distances, bins=bins, edgecolor="black")
    if title:
        plt.title(title)
    plt.xlabel("Distance (meters)")
    plt.ylabel("Frequency")
    if zoom_max is not None:
        plt.xlim(0, zoom_max)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

plot_match_distance_hist(matched_spills_gdf, bins=40) #title="Distance Between Original and Matched Spill Points"
plot_match_distance_hist(matched_spills_gdf, bins=6000, zoom_max=1000) # title="Distance Between Original and Matched Spill Points (Zoomed In)"

# =============================================================================
# 5. Find threshold distance between modes
# =============================================================================
def find_threshold_distance(
    gdf,
    column="match_distance_m",
    bins=50,
    sigma=1.5,
    prominence=5,
    count_cutoff=30
):
    distances = gdf[column].dropna()
    distances = distances[distances > 0]
    if distances.empty:
        return None, {"reason": "no_positive_distances"}

    mean_dist = float(distances.mean())
    median_dist = float(distances.median())

    min_val = float(distances.min())
    max_val = float(distances.max())
    if min_val <= 0 or max_val <= 0 or min_val == max_val:
        bins_array = bins
    else:
        bins_array = np.logspace(np.log10(min_val), np.log10(max_val), bins)

    counts, bin_edges = np.histogram(distances, bins=bins_array)
    counts_smooth = gaussian_filter1d(counts.astype(float), sigma=sigma)

    peak_idx, _ = find_peaks(counts_smooth, prominence=prominence)
    peak_idx = np.sort(peak_idx)

    split_index = None
    last_idx = None
    threshold_value = None

    if len(peak_idx) >= 2:
        first_mode, second_mode = int(peak_idx[0]), int(peak_idx[1])
        valley_slice = slice(first_mode, second_mode)
        valley_rel = int(np.argmin(counts_smooth[valley_slice]))
        split_index = first_mode + valley_rel

        for i in range(split_index - 1, -1, -1):
            if counts[i] > count_cutoff:
                last_idx = i
                threshold_value = float(bin_edges[last_idx + 1])
                break

    details = {
        "counts": counts,
        "bin_edges": bin_edges,
        "counts_smooth": counts_smooth,
        "peak_idx": peak_idx,
        "mean": mean_dist,
        "median": median_dist,
        "min_val": min_val,
        "max_val": max_val,
        "split_index": split_index,
        "last_idx": last_idx
    }
    return threshold_value, details

def plot_distance_histogram(
    gdf,
    details,
    threshold_value=None,
    column="match_distance_m",
    figsize=(10, 6),
    color="skyblue",
    count_cutoff=None
):
    distances = gdf[column].dropna()
    distances = distances[distances > 0]

    plt.figure(figsize=figsize)
    plt.hist(distances, bins=details["bin_edges"], edgecolor="black", color=color)

    plt.axvline(details["mean"], color='red', linestyle='--', label=f'Mean: {details["mean"]:.1f} m')
    plt.axvline(details["median"], color='green', linestyle='--', label=f'Median: {details["median"]:.1f} m')

    if threshold_value is not None:
        plt.axvline(threshold_value, color='blue', linestyle='--', linewidth=2,
                    label=f'Threshold: {threshold_value:.2f} m')

    if details["min_val"] > 0 and details["min_val"] != details["max_val"]:
        plt.xscale("log")

    plt.xlabel("Distance (meters, log scale)")
    plt.ylabel("Frequency")
    plt.grid(axis='y')
    plt.legend()

    title_parts = ["Log Distance Between Original and Matched Spill Points"]
    plt.title(" ".join(title_parts))
    plt.tight_layout()
    plt.show()

# Threshold
threshold, details = find_threshold_distance(
    matched_spills_gdf,
    bins=50,
    sigma=1.5,
    prominence=5,
    count_cutoff=30
)

if threshold is not None:
    print(f"Threshold distance: {threshold:.2f} m")
else:
    print("No threshold found before the second mode.")

# Plot log histogram
plot_distance_histogram(
    matched_spills_gdf,
    details,
    threshold_value=threshold,
    count_cutoff=30
)

# =============================================================================
# 6. Save spills under threshold
# =============================================================================
def save_spills_under_threshold(gdf, threshold, output_file):
    filtered_gdf = gdf[gdf["match_distance_m"] < threshold].copy()
    filtered_gdf = filtered_gdf.to_crs(epsg=4326)
    filtered_gdf.to_file(output_file, driver="GeoJSON")
    print(f"Wrote {output_file} in EPSG:4326.")
    return filtered_gdf

spills_under_85m = save_spills_under_threshold(
    matched_spills_gdf,
    threshold=85.07,
    output_file="updated_spills_under_85m.geojson"
)

# =============================================================================
# 7. De-duplicate spills_under_85m by coordinate (keep first)
# =============================================================================
spills_under_85m = gpd.read_file("updated_spills_under_85m.geojson")

# keep only valid Point geometries
spills_under_85m = spills_under_85m[
    spills_under_85m.geometry.notna()
    & (~spills_under_85m.geometry.is_empty)
    & (spills_under_85m.geometry.type == "Point")
].copy()

# coord key with rounding to tame float noise (~1e-6 degrees)
spills_under_85m["coord_key"] = spills_under_85m.geometry.apply(
    lambda p: (round(p.x, 6), round(p.y, 6))
)

before = len(spills_under_85m)
spills_under_85m = spills_under_85m.drop_duplicates(subset=["coord_key"], keep="first").drop(columns=["coord_key"])
after = len(spills_under_85m)

print(f"Spills ≤85.07 m before dedup: {before}")
print(f"Spills ≤85.07 m after  dedup: {after}")

spills_under_85m.to_file("updated_spills_under_85m.geojson", driver="GeoJSON")

# =============================================================================
# 8. Dataset summary
# =============================================================================
def summarize_spill_datasets(spills_file, matched_file, matched_under_file):
    spills_raw = gpd.read_file(spills_file)
    matched_spills = gpd.read_file(matched_file)
    matched_under = gpd.read_file(matched_under_file)

    n_raw = len(spills_raw)
    n_matched = len(matched_spills)
    n_under = len(matched_under)

    print(f"Original spill points ({spills_file}):               {n_raw}")
    print(f"Matched spills ({matched_file}):                     {n_matched}")
    print(f"Matched spills under threshold ({matched_under_file}): {n_under}")
    print(f"Retention rate under threshold:                      {n_under / n_raw:.2%}")

    return {
        "n_raw": n_raw,
        "n_matched": n_matched,
        "n_under": n_under,
        "retention_rate": n_under / n_raw
    }

summary = summarize_spill_datasets(
    "spills.geojson",
    "updated_spills.geojson",
    "updated_spills_under_85m.geojson"
)

# =============================================================================
# 9. Join spills with flowline attributes (prefer spill columns on conflicts)
# =============================================================================
joined_rows = []

match_idx_col = "matched_flowline_idx"
if match_idx_col not in spills_under_85m.columns:
    match_idx_col = "matched_crudeoil_idx"

for i, spill in spills_under_85m.iterrows():
    line_idx = spill.get(match_idx_col)
    if pd.isna(line_idx) or line_idx not in flowlines_gdf.index:
        print(f"[Skip spill {i}] No valid match index: {line_idx}")
        continue

    # Flowline attributes (no geometry)
    flowline_row = flowlines_gdf.loc[line_idx]
    flow_dict = flowline_row.drop(labels=["geometry"]).to_dict()

    # Spill attributes to keep (no prefix), keep its own match_distance_m
    spill_keep_cols = [c for c in spill.index if c not in {"geometry", match_idx_col}]
    spill_dict = spill[spill_keep_cols].to_dict()

    # Merge: flowline first, then spill overwrites any name collisions
    merged = {**flow_dict, **spill_dict}
    merged["geometry"] = spill.geometry  # always use spill point geometry

    joined_rows.append(merged)

joined_gdf = gpd.GeoDataFrame(joined_rows, geometry="geometry", crs=spills_under_85m.crs)
if joined_gdf.crs is None or joined_gdf.crs.to_epsg() != 4326:
    joined_gdf = joined_gdf.to_crs(epsg=4326)

joined_gdf["X"] = joined_gdf.geometry.x
joined_gdf["Y"] = joined_gdf.geometry.y
joined_gdf["coords"] = list(zip(joined_gdf["X"], joined_gdf["Y"]))

cols_to_drop = ["lon", "lat", "X", "Y", "coords"]
joined_gdf = joined_gdf.drop(columns=[c for c in cols_to_drop if c in joined_gdf.columns])

desired_order = [
    "unique_id", "operator_name", "operator_number",
    "flowline_id", "location_id", "status",
    "flowline_action", "location_type", "fluid",
    "material", "diameter_in", "length_ft",
    "max_operating_pressure", "line_age_yr", "construct_date",
    "match_distance_m", "incident_date", "root_cause", "risk",
    "geometry"
]
current_cols = joined_gdf.columns.tolist()
ordered_cols = [col for col in desired_order if col in current_cols] + \
               [col for col in current_cols if col not in desired_order]

joined_gdf = joined_gdf[ordered_cols]

joined_gdf.to_file("spills_w_flowline_attributes.geojson", driver="GeoJSON")
print(f"Saved {len(joined_gdf)} spill-points with flowline attributes (spill columns take precedence).")