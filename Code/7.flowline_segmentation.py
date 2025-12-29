# =============================================================================
# 1. Setup & Imports
# =============================================================================
import os
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import MultiLineString, LineString, Point
from shapely.ops import nearest_points, transform
import matplotlib.pyplot as plt
import seaborn as sns
import pyproj

# =============================================================================
# 2. Load & Prepare Data
# =============================================================================
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

matched_flowlines_gdf = gpd.read_file('combined_flowlines_final.geojson')
matched_flowlines_gdf = matched_flowlines_gdf.to_crs(epsg=4326)

# =============================================================================
# 3. Compute Global Min/Max BEFORE Filtering
# =============================================================================
geod = pyproj.Geod(ellps="WGS84")

def compute_geod_length(geom):
    if isinstance(geom, MultiLineString):
        return sum(geod.geometry_length(line) for line in geom.geoms)
    elif isinstance(geom, LineString):
        return geod.geometry_length(geom)
    else:
        return 0.0

matched_flowlines_gdf['geod_length_m'] = matched_flowlines_gdf['geometry'].apply(compute_geod_length)

original_min = matched_flowlines_gdf['geod_length_m'].min()
original_max = matched_flowlines_gdf['geod_length_m'].max()
print(f"Original global min length (m): {original_min:.4f}")
print(f"Original global max length (m): {original_max:.2f}")

# =============================================================================
# 4. Filter for Flowlines ≥ 1 Meter
# =============================================================================
min_valid_length = 1.0
filtered_flowlines_gdf = matched_flowlines_gdf[
    matched_flowlines_gdf['geod_length_m'] >= min_valid_length
].copy()

n_original = len(matched_flowlines_gdf)
n_filtered = len(filtered_flowlines_gdf)
print(f"Remaining flowlines ≥ {min_valid_length} m: {n_filtered}")
print(f"Filtered out: {n_original - n_filtered}")

filtered_min = filtered_flowlines_gdf['geod_length_m'].min()
filtered_max = filtered_flowlines_gdf['geod_length_m'].max()
print(f"Filtered global min length (m): {filtered_min:.4f}")
print(f"Filtered global max length (m): {filtered_max:.2f}")

# =============================================================================
# 5. Build Length Array for Plotting
# =============================================================================
line_lengths = []
for geom in filtered_flowlines_gdf['geometry']:
    if isinstance(geom, MultiLineString):
        for line in geom.geoms:
            line_lengths.append(geod.geometry_length(line))
    elif isinstance(geom, LineString):
        line_lengths.append(geod.geometry_length(geom))

line_lengths = np.array(line_lengths)
line_lengths = line_lengths[line_lengths >= 1.0]

# =============================================================================
# 6. Plot Histograms
# =============================================================================

# Histogram 1: Full distribution
plt.figure(figsize=(10, 6))
sns.histplot(line_lengths, bins=50)
plt.xlabel('Length (meters)')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# Histogram 2: Truncated at 500 m
max_length_to_plot = 500
truncated_lengths = line_lengths[line_lengths <= max_length_to_plot]

plt.figure(figsize=(10, 6))
sns.histplot(truncated_lengths, bins=50)
plt.xlabel('Length (meters)')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# Histogram 3: Log-scaled X-axis
mean_len = np.mean(line_lengths)
median_len = np.median(line_lengths)

plt.figure(figsize=(10, 6))
sns.histplot(line_lengths, bins=50, kde=False, log_scale=(True, False), color='skyblue')

plt.axvline(mean_len, color='red', linestyle='--', label=f'Mean: {mean_len:.1f} m')
plt.axvline(median_len, color='green', linestyle='--', label=f'Median: {median_len:.1f} m')

plt.xlabel('Length (meters, log scale)')
plt.ylabel('Frequency')
plt.grid(True, which='major', axis='x')
plt.legend()
plt.show()

# =============================================================================
# 7. Generate Points Every 50 Meters Along Flowlines
# =============================================================================
og_length = len(matched_flowlines_gdf)
print(f"Original number of flowlines: {og_length}")

# Project to meters for spacing calculations
gdf_m = matched_flowlines_gdf.to_crs(matched_flowlines_gdf.estimate_utm_crs())

points_list = []
for idx, row in gdf_m.iterrows():
    geom = row.geometry
    if geom.is_empty:
        continue

    # Merge MultiLineString parts into one continuous line if needed
    if geom.geom_type == "MultiLineString":
        from shapely.ops import linemerge
        geom = linemerge(geom)

    length = geom.length
    distances = np.arange(0, length + 1e-9, 50)  # every 50 m

    for d in distances:
        pt = geom.interpolate(d)
        new_row = row.drop(labels="geometry")
        new_row["geometry"] = pt
        points_list.append(new_row)

# Create points GeoDataFrame and reproject to original CRS
points_gdf = gpd.GeoDataFrame(points_list, geometry="geometry", crs=gdf_m.crs)
points_gdf = points_gdf.to_crs(matched_flowlines_gdf.crs)

# Save intermediate point file
points_gdf.to_file("flowline_points_50m.geojson", driver="GeoJSON")
print(f"Segmented into {len(points_gdf)} point features.")

# =============================================================================
# 8. Validate & Prepare Points
# =============================================================================
points_gdf = gpd.read_file("flowline_points_50m.geojson")

# Keep only valid Point geometries
valid = (
    points_gdf.geometry.notna()
    & (~points_gdf.geometry.is_empty)
    & (points_gdf.geometry.type == "Point")
)
points_gdf = points_gdf.loc[valid].copy()

# Coordinate key (rounded to avoid float precision noise)
points_gdf["coord_key"] = points_gdf.geometry.apply(lambda p: (round(p.x, 6), round(p.y, 6)))

print(f"Total valid points: {len(points_gdf)}")
print(f"Unique coordinate locations: {points_gdf['coord_key'].nunique()}")

# =============================================================================
# 9. Identify Unique & Duplicate Points
# =============================================================================
vc = points_gdf["coord_key"].value_counts()

# Unique: occurs once
unique_mask = points_gdf["coord_key"].map(vc).eq(1)
unique_50_points = points_gdf.loc[unique_mask].copy()

# Duplicate: occurs ≥ 2 times
dup_mask_all = points_gdf["coord_key"].map(vc).ge(2)
dup_50_points = points_gdf.loc[dup_mask_all].copy()

print(f"Unique (non-dup) rows: {len(unique_50_points)}")
print(f"Duplicate rows (incl. first): {len(dup_50_points)}")

# =============================================================================
# 10. Keep Duplicates with Spill Matches
# =============================================================================
spills_gdf = gpd.read_file("spills_w_flowline_attributes.geojson")  # contains 'unique_id'
spills_ids = set(spills_gdf["unique_id"].dropna())

# Duplicates that match a spill's unique_id
dup_in_spills = dup_50_points[dup_50_points["unique_id"].isin(spills_ids)].copy()

# Keep only first occurrence per coord_key
dup_in_spills_kept = dup_in_spills.sort_index().drop_duplicates(subset=["coord_key"], keep="first")

print(f"Duplicate rows with unique_id in spills: {len(dup_in_spills)}")
print(f"Kept from that subset: {len(dup_in_spills_kept)}")

# =============================================================================
# 11. Remove Spill-Matched Duplicates from Remaining
# =============================================================================
kept_keys = set(dup_in_spills_kept["coord_key"])
dup_remaining = dup_50_points[~dup_50_points["coord_key"].isin(kept_keys)].copy()

print(f"Remaining duplicate rows after spills-keep: {len(dup_remaining)}")

# =============================================================================
# 12. Randomly Keep One Row per Remaining Duplicate
# =============================================================================
dup_remaining_keep = (
    dup_remaining.groupby("coord_key", group_keys=False)
                 .apply(lambda df: df.sample(n=1, random_state=42))
)

print(f"Random-kept from remaining duplicate groups: {len(dup_remaining_keep)}")

# =============================================================================
# 13. Combine All Kept Points
# =============================================================================
final_points = gpd.GeoDataFrame(
    pd.concat([unique_50_points, dup_in_spills_kept, dup_remaining_keep], ignore_index=True),
    geometry="geometry",
    crs=points_gdf.crs
)

# Drop helper column
final_points = final_points.drop(columns=["coord_key"])

print(f"Final points after de-duplication: {len(final_points)}")

# =============================================================================
# 14. Validation
# =============================================================================
expected_unique_coords = points_gdf["coord_key"].nunique()
got = len(final_points)

print(f"Expected unique coords: {expected_unique_coords}")
print(f"Got final rows:         {got}")
print("PASS ✅" if got == expected_unique_coords else "MISMATCH ❌")

# =============================================================================
# 15. Check for NAs and Zeros
# =============================================================================
def check_na_zero(df):
    results = []
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            na_count = df[col].isna().sum()
            zero_count = (df[col] == 0).sum()
            results.append({
                'column': col,
                'na_count': na_count,
                'zero_count': zero_count
            })
    return pd.DataFrame(results)

print("\n--- NA and Zero Counts: final_points ---")
print(check_na_zero(final_points))

cols_to_check = ['diameter_in', 'max_operating_pressure', 'length_ft']
for col in cols_to_check:
    if col in final_points.columns:
        final_points = final_points[final_points[col] != 0]

final_points = final_points.reset_index(drop=True)
print(f"Dataset size after removing rows with zeros in {cols_to_check}: {len(final_points)}")
print(check_na_zero(final_points))

# =============================================================================
# 16. Summary Statistics
# =============================================================================
print("\n===== Pipeline Summary =====")
print(f"Original MultiLineString flowlines: {len(matched_flowlines_gdf)}")
print(f"Final de-duplicated points:         {len(final_points)}")
print("================================\n")

final_points.to_file("flowline_points_50m_dedup.geojson", driver="GeoJSON")
