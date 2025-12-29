# --- Start of Script 1: 1.excel_to_gbd.ipynb ---
# This code transforms the original data formats into GeoJSON files, removing unnecessary attributes along the way.
# The gbd files are in CRS : EPSG-26913

import os
import pandas as pd
import openpyxl
import geopandas as gpd
from shapely.geometry import LineString

os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
pd.options.display.max_columns = None

def drop_rows_with_word(df, word="gathering"):
    df = df.copy()
    geom_col = df.geometry.name if hasattr(df, "geometry") else None
    mask = pd.Series(False, index=df.index)
    for c in df.columns:
        if c == geom_col:
            continue
        if pd.api.types.is_object_dtype(df[c]):
            mask |= df[c].astype(str).str.contains(word, case=False, na=False)
    return df.loc[~mask].copy()

# 1) load datasets
gdf0 = gpd.read_file("ECMC_Flowline_Data_Access/COGCC_Form44_Off_Location_Flowlines_Approved_CONFIDENTIAL.gdb")
gdf1 = gpd.read_file("ECMC_Flowline_Data_Access/COGCC_Form44_Crude_Oil_Produced_Water_Transfer_Flowlines_Approved_CONFIDENTIAL.gdb")
flowlines = pd.read_excel('FlowlineSpreadsheet_Mines.xlsx')
spills = pd.read_excel('Flowline-Related Spills (through 2024).xlsx')

# 2) drop 'gathering' everywhere
gdf0 = drop_rows_with_word(gdf0, "gathering")
gdf1 = drop_rows_with_word(gdf1, "gathering")
flowlines = drop_rows_with_word(flowlines, "gathering")
spills = drop_rows_with_word(spills, "gathering")

# 3) process spills
columns_list = spills.columns.tolist()
print("Columns in spills DataFrame:", columns_list)

selected_columns = ['trkg_num', 'Operator Name', 'facility_type', 'Spill_Desc', 'Spill Type', 'Root Cause', 'Preventative Measure',
                    'Root Cause Type', 'Detailed Root Cause Type', 'Long', 'Lat','facility_status', 'Gathering?', 'Metallic?', 'incident_date']
spills_selected = spills[selected_columns]

num_rows = len(spills_selected)
print("Number of rows:", num_rows)

spills_cleaned = spills_selected
only_nan_columns = spills_cleaned.isna().all()
columns_with_all_nan = only_nan_columns[only_nan_columns].index.tolist()
print("Columns with all NaN values:", columns_with_all_nan)

spills_filtered = spills_cleaned

num_rows = len(spills_filtered)
print("New number of rows:", num_rows)

spl_gdf = gpd.GeoDataFrame(spills_filtered, geometry=gpd.points_from_xy(spills_filtered.Long, spills_filtered.Lat), crs='EPSG:4326')

# 4) combine gdbs; build flowlines GeoDataFrame
if gdf0.crs != gdf1.crs:
    gdf1 = gdf1.to_crs(gdf0.crs)

gdf = pd.concat([gdf0, gdf1], ignore_index=True)
if 'Doc_Num' in gdf.columns:
    gdf.drop(columns=['Doc_Num'], inplace=True)

columns_list = flowlines.columns.tolist()
print("Columns in spills DataFrame:", columns_list)

selected_columns = ['LOCATION_ID', 'FLOWLINEID', 'STARTLOCATIONID', 'FLOWLINEACTION', 'ENTIRELINEREMOVED', 'ACTIONDESCRIPTION', 'RECEIVE_DATE', 'OPERATOR_NUM', 'COMPANY_NAME', 'LOCATIONTYPE', 'ENDLAT', 'ENDLONG', 'STARTLAT', 'STARTLONG',
                    'PIPEMATERIAL', 'BEDDINGMATERIAL', 'TYPEOFFLUIDTRANS', 'MAXOPPRESSURE', 'CONSTRUCTDATE']
flowlines_selected = flowlines[selected_columns]

if 'geometry' not in flowlines_selected.columns:
    flowlines_selected['geometry'] = ''

for index, row in flowlines_selected.iterrows():
    geom = LineString([(row['STARTLONG'],row['STARTLAT']),(row['ENDLONG'],row['ENDLAT'])])
    flowlines_selected.at[index,'geometry'] = geom

fl_gdf = gpd.GeoDataFrame(flowlines_selected, geometry='geometry', crs='EPSG:4326')

if fl_gdf.crs != gdf.crs:
    fl_gdf.to_crs(gdf.crs, inplace=True)
    print('Change fl crs to gdf crs')

# 5) align spills CRS and save bases
if spl_gdf.crs != gdf.crs:
    spl_gdf.to_crs(gdf.crs, inplace=True)
    print('Change spl crs to gdf crs')

gdf.to_file('crudeoil_offlocation.geojson', driver='GeoJSON')
fl_gdf.to_file('flowlines.geojson', driver='GeoJSON')
spl_gdf.to_file('spills.geojson', driver='GeoJSON')

# 6) match crude→flowlines
flowlines_gdf = gpd.read_file('flowlines.geojson')
crudeoil_gdf = gpd.read_file('crudeoil_offlocation.geojson')

if flowlines_gdf.crs != crudeoil_gdf.crs:
    flowlines_gdf = flowlines_gdf.to_crs(crudeoil_gdf.crs)

matched_rows, matched_count, skipped_count = [], 0, 0

for idx, crude in crudeoil_gdf.iterrows():
    geom = crude.geometry
    if geom is None:
        print(f"[{idx}] Missing geometry – skipping.")
        skipped_count += 1
        continue

    op = crude.get("Operator", "")
    if pd.isnull(op) or not op.strip():
        print(f"[{idx}] No Operator – skipping.")
        skipped_count += 1
        continue

    op_name = op.strip().lower()
    candidates = flowlines_gdf[
        flowlines_gdf["COMPANY_NAME"].str.strip().str.lower().eq(op_name)
    ].copy()
    candidates = candidates[candidates.geometry.notnull()]

    if candidates.empty:
        print(f"[{idx}] No flowlines for operator “{op}”.")
        skipped_count += 1
        continue

    dists = candidates.geometry.distance(geom).dropna()
    if dists.empty:
        print(f"[{idx}] All distances NaN – skipping.")
        skipped_count += 1
        continue

    nearest_idx = dists.idxmin()
    min_dist    = dists.min()

    matched_count += 1
    print(f"[Match {matched_count}] crudeoil index {idx} → flowline index {nearest_idx} at {min_dist:.2f} m")

    new_row = crude.copy()
    for col in flowlines_gdf.columns:
        if col == flowlines_gdf.geometry.name:
            continue
        new_row[col] = flowlines_gdf.loc[nearest_idx, col]
    new_row["flowline_match_distance_m"] = min_dist
    matched_rows.append(new_row)

total = len(crudeoil_gdf)
print("\n=== Complete ===")
print(f" Total crude-oil features:         {total}")
print(f" Successfully matched:             {matched_count}")
print(f" Skipped (no match or missing data): {skipped_count}")

matched_crudeoil_gdf = gpd.GeoDataFrame(matched_rows, crs=crudeoil_gdf.crs)
matched_crudeoil_gdf["unique_id"] = range(1, len(matched_crudeoil_gdf) + 1)
matched_crudeoil_gdf = drop_rows_with_word(matched_crudeoil_gdf, "gathering")
matched_crudeoil_gdf.to_file("combined_flowlines.geojson", driver="GeoJSON")

# 7) sizes
print("\n=== Saved dataset sizes ===")
print("crudeoil_offlocation.geojson rows:", len(gdf))
print("flowlines.geojson rows:", len(fl_gdf))
print("spills.geojson rows:", len(spl_gdf))
print("combined_flowlines.geojson rows:", len(matched_crudeoil_gdf))

# === Saved dataset sizes ===
""" 
crudeoil_offlocation.geojson rows: 259977
flowlines.geojson rows: 21211
spills.geojson rows: 1115
full_length_flowlines.geojson rows: 159048 
"""