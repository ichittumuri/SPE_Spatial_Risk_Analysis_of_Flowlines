# =============================================================================
# 1. Setup & Imports
# =============================================================================
import os
import pandas as pd
import geopandas as gpd
import rasterio
from pyproj import CRS

# =============================================================================
# 2. Load Datasets & Ensure Same CRS
# =============================================================================
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')

combined = gpd.read_file("combined_flowlines_spills.geojson")

# =============================================================================
# 5. Sample Elevation from DEM
# =============================================================================
dem = rasterio.open('output_USGS30m.tif')
dem_crs = CRS(dem.crs)

# Reproject combined to DEM CRS if needed
if combined.crs != dem_crs:
    combined = combined.to_crs(dem_crs)

def get_elevation(point, dem):
    try:
        val = list(dem.sample([(point.x, point.y)]))[0][0]
        if dem.nodata is not None and val == dem.nodata:
            return None
        return float(val)
    except Exception:
        return None

combined['elevation'] = combined.geometry.apply(lambda pt: get_elevation(pt, dem))
combined = combined.drop(columns=['index_right'], errors='ignore')

# =============================================================================
# 6. Join Census Tract Population Density (spatial join)
# =============================================================================
pop_density = gpd.read_file('Population_Density_(Census_Tracts)').to_crs(combined.crs)

print('Summary of Census Tract Data:')
print(pop_density.info())
print('\nFirst few rows of the data:')
print(pop_density.head())

joined = gpd.sjoin(combined, pop_density, how='left', predicate='within')
combined['avg_population'] = joined['Populati_1']

# =============================================================================
# 7. Check for NAs and Zeros
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

print("\n--- NA and Zero Counts: combined ---")
print(check_na_zero(combined))

# =============================================================================
# 8 Fill missing values via spatial IDW from nearby points
# =============================================================================
import numpy as np
from sklearn.neighbors import BallTree

def idw_fill(gdf, col, k=8, power=2.0, max_dist=None):
    # coords (n, 2)
    xs = gdf.geometry.x.to_numpy()
    ys = gdf.geometry.y.to_numpy()
    coords = np.column_stack([xs, ys])

    known = gdf[col].notna().to_numpy()
    missing = ~known

    # nothing to do
    if known.sum() == 0 or missing.sum() == 0:
        return gdf

    # fit on knowns
    tree = BallTree(coords[known], metric='euclidean')
    k_eff = min(k, known.sum())  # can't query more neighbors than we have

    # query neighbors for missing rows
    dists, nbr_idx = tree.query(coords[missing], k=k_eff)  # (m, k)
    vals = gdf.loc[known, col].to_numpy()[nbr_idx]         # (m, k)

    # optionally mask out far neighbors
    if max_dist is not None:
        far = dists > max_dist
        vals = np.where(far, np.nan, vals)
        dists = np.where(far, np.nan, dists)

    # weights: 1 / d^power  (avoid div-by-zero; if d==0, weight is huge)
    dists = np.where(dists <= 0, 1e-9, dists)
    w = 1.0 / (dists ** power)

    # handle all-NaN neighbor sets (can happen with max_dist)
    w = np.where(np.isnan(vals), 0.0, w)
    vals = np.where(np.isnan(vals), 0.0, vals)

    num = (w * vals).sum(axis=1)
    den = w.sum(axis=1)

    est = np.full(missing.sum(), np.nan, dtype=float)
    good = den > 0
    est[good] = num[good] / den[good]

    # fallback: if no valid neighbors (den==0), use nearest neighbor (k=1) regardless of distance
    if (~good).any():
        d1, i1 = tree.query(coords[missing][~good], k=1)
        v1 = gdf.loc[known, col].to_numpy()[i1[:, 0]]
        est[~good] = v1

    # write back
    gdf.loc[missing, col] = est
    return gdf

# Run IDW fill for the two columns with NaNs.
# (Your GeoDataFrame is already in the DEM CRS here, which is perfect.)
combined = idw_fill(combined, 'elevation', k=8, power=2.0, max_dist=1500)   # ~1.5 km cap (tweak as you like)
combined = idw_fill(combined, 'avg_population', k=8, power=2.0, max_dist=5000) # tracts can be larger; 5 km is reasonable

# Quick check
print("\n--- After IDW fill ---")
print(combined[['elevation','avg_population']].isna().sum())

# =============================================================================
# 9. Column Ordering, Finalize, and Save (Full Dataset)
# =============================================================================
desired = [
    'unique_id', 'operator_name', 'operator_number', 'flowline_id', 'location_id',
    'status', 'flowline_action', 'location_type', 'fluid', 'material',
    'diameter_in', 'length_ft', 'max_operating_pressure', 'line_age_yr', 'elevation', 'avg_population',
    'construct_date', 'incident_date', 'geod_length_m', 'match_distance_m',
    'root_cause', 'risk', 'coord_key', 'geometry'
]

ordered_cols = [c for c in desired if c in combined.columns] + \
               [c for c in combined.columns if c not in desired]

combined = combined[ordered_cols]
combined = combined.set_geometry('geometry')  # re-affirm geometry column
combined.to_file("final_dataset_full.geojson", driver="GeoJSON")

# =============================================================================
# 10. Export Selected Subset
# =============================================================================
desired_subset = [
    'unique_id', 'status', 'flowline_action', 'location_type', 'fluid', 'material',
    'diameter_in', 'length_ft', 'max_operating_pressure', 'line_age_yr', 'elevation',
    'match_distance_m', 'risk', 'geometry' 
]

final_subset = combined[[c for c in desired_subset if c in combined.columns]]
final_subset.to_file("final_dataset_subset.geojson", driver="GeoJSON")

# =============================================================================
# 11. Check for NAs and Zeros
# =============================================================================

print("\n--- NA and Zero Counts: combined ---")
print(check_na_zero(combined))

print("\n--- NA and Zero Counts: final_subset ---")
print(check_na_zero(final_subset))

# =============================================================================
# 12. Dataset Summary
# =============================================================================
print(f"Number of rows: {len(final_subset)}")
print(f"Number of columns: {final_subset.shape[1]}")

if 'risk' in final_subset.columns:
    print(f"Rows with risk = 1: {(final_subset['risk'] == 1).sum()}")
    print(f"Rows with risk = 0: {(final_subset['risk'] == 0).sum()}")





