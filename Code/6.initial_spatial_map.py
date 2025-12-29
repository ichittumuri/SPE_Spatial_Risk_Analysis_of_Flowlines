# plot_spills_flowlines.py
import os
import warnings
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import contextily as ctx
from shapely.geometry import box

# -----------------------------
# 1) Paths + basic setup
# -----------------------------
os.chdir('/Users/ichittumuri/Desktop/MINES/COGCC-Risk-Analysis/Data')
warnings.filterwarnings("ignore", category=UserWarning)

FLOWLINES_FP = "combined_flowlines_final.geojson"
SPILLS_FP    = "updated_spills_under_85m.geojson" 
# "spills_w_flowline_attributes.geojson"
DIST_COL     = "match_distance_m"

# -----------------------------
# 2) Load
# -----------------------------
flowlines = gpd.read_file(FLOWLINES_FP)
spills = gpd.read_file(SPILLS_FP)
spills = spills[spills.geometry.type.isin(["Point", "MultiPoint"])].copy()

# -----------------------------
# 3) Full map
# -----------------------------
f_3857 = flowlines.to_crs(3857)
s_3857 = spills.to_crs(3857)

fig, ax = plt.subplots(figsize=(10, 10))
f_3857.plot(ax=ax, color="#2C7BB6", linewidth=1.3, alpha=0.8, zorder=2) # 0B50B8

d = s_3857[DIST_COL]
norm = mcolors.Normalize(vmin=float(d.min()), vmax=float(d.max()))
s_3857.plot(
    ax=ax,
    column=DIST_COL,
    cmap="YlOrRd",
    norm=norm,
    markersize=12,
    edgecolor="black",
    linewidth=0.2,
    legend=True,
    zorder=3,
    legend_kwds={"label": "Distance to Flowline (m)", "shrink": 0.6},
)

ctx.add_basemap(ax, source=ctx.providers.CartoDB.Positron, zorder=1)
# ax.set_title("Spills ≤ 85.07 m — Full Extent", fontsize=12)
ax.set_axis_off()
plt.tight_layout()
plt.savefig("spills_flowlines_full.png", dpi=200)
plt.show()

# -----------------------------
# 4) Cropped map (Boulder–Greeley) — match R bounds
# -----------------------------
# ------------------------------------------------------------
 
# ------------------------------------------------------------
from pyproj import Transformer
import numpy as np

# -----------------------------
# 1) R bounds (EXACT)
# -----------------------------
# --- auto-trim controls ---
BOTTOM_DROP_PCT   = 0.01   # drop bottom 1% of flowlines by centroid Y
BOTTOM_EXTRA_PAD_M = 80    # extra meters above that cutoff (avoid edge clutter)

# padding around final flowline box (after trimming)
PX, PY = 120, 120          # set to 0 for razor-tight

# ============================================================
# 1) Exact R bbox (lon/lat) -> Web Mercator
# ============================================================
min_lon, min_lat = -105.85, 39.56
max_lon, max_lat = -104.59, 40.68

tf = Transformer.from_crs(4326, 3857, always_xy=True)
xmin_b, ymin_b = tf.transform(min_lon, min_lat)
xmax_b, ymax_b = tf.transform(max_lon, max_lat)

# ============================================================
# 2) Clip vectors to R bbox
# ============================================================
clip_geom = gpd.GeoDataFrame(geometry=[box(min_lon, min_lat, max_lon, max_lat)],
                             crs=4326).to_crs(3857).geometry.iloc[0]
f_clip = f_3857.clip(clip_geom)
s_clip = s_3857.clip(clip_geom)

if f_clip.empty:
    raise ValueError("No flowlines inside the bbox; cannot compute tight extent.")

# ============================================================
# 3) AUTO-trim very bottom flowlines, then recompute tight box
# ============================================================
# centroid Y for each flowline geometry (already in meters)
flow_cy = f_clip.geometry.centroid.y
y_cut = flow_cy.quantile(BOTTOM_DROP_PCT) + BOTTOM_EXTRA_PAD_M

# keep only features above the cutoff
f_clip = f_clip[flow_cy > y_cut]

# (optional) keep spills above same cutoff so none show below the map frame
if not s_clip.empty:
    s_clip = s_clip[s_clip.geometry.centroid.y > y_cut]

# tight bounds from remaining flowlines
fminx, fminy, fmaxx, fmaxy = f_clip.total_bounds

# clamp to R bbox and add gentle padding
xmin = max(fminx - PX, xmin_b)
xmax = min(fmaxx + PX, xmax_b)
ymin = max(fminy - PY, ymin_b)
ymax = min(fmaxy + PY, ymax_b)

# tiny epsilon to hide tile seams
eps = 1.0
xmin += eps; ymin += eps; xmax -= eps; ymax -= eps

# ============================================================
# 4) Fetch basemap exactly for this extent
# ============================================================
img, ext = ctx.bounds2img(xmin, ymin, xmax, ymax, source=ctx.providers.CartoDB.Positron)

# ============================================================
# 5) Plot
# ============================================================
fig, ax = plt.subplots(figsize=(10, 12))

# Basemap locked to exact extent
ax.imshow(img, extent=ext, origin="upper", zorder=1)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_aspect("equal")
ax.autoscale(False)

# Flowlines
f_clip.plot(ax=ax, color="#2C7BB6", linewidth=1, alpha=0.95, zorder=2)

# Spills colored by distance
if not s_clip.empty:
    d2 = s_clip[DIST_COL]
    norm2 = mcolors.Normalize(vmin=float(d2.min()), vmax=float(d2.max()))
    s_clip.plot(
        ax=ax, column=DIST_COL, cmap="YlOrRd", norm=norm2,
        markersize=12, edgecolor="black", linewidth=0.2,
        legend=True, zorder=3,
        legend_kwds={"label": "Distance to Flowline (m)", "shrink": 0.6},
    )

# Finish
ax.margins(0)
ax.set_axis_off()
plt.tight_layout()
plt.savefig("spills_flowlines_TIGHT_FLOWLINES_TRIMMED.png", dpi=220,
            bbox_inches="tight", pad_inches=0)
plt.show()