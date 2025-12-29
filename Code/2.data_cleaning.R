library(sf)           
library(lubridate)
library(dplyr)
library(stringr)
library(lubridate)
library(VIM)
library(zoo)

setwd("~/Desktop/MINES/COGCC-Risk-Analysis/Data")

flowlines <- st_read("combined_flowlines.geojson")
colnames(flowlines)

# Line Age
flowlines$CONSTRUCTDATE <- ymd(flowlines$CONSTRUCTDATE)
today <- Sys.Date()
flowlines$line_age_yr <- as.numeric(difftime(today, flowlines$CONSTRUCTDATE, units = "days")) / 365.25

# Operator name cleanup
op_name_fix <- c(
  "KINDER MORGAN CO2 CO LP" = "KINDER MORGAN CO2 CO LLC",
  "BEEMAN OIL & GAS INC"    = "BEEMAN OIL & GAS LLC"
)
if ("Operator" %in% colnames(flowlines)) {
  flowlines$Operator <- dplyr::recode(flowlines$Operator, !!!op_name_fix)
}

# drop unwanted columns
columns_to_remove <- c(
  "SHAPE_Length", "STARTLOCATIONID", "ENTIRELINEREMOVED", "ACTIONDESCRIPTION", "RECEIVE_DATE",
  "COMPANY_NAME", "ENDLAT", "ENDLONG", "STARTLAT", "STARTLONG", "PIPEMATERIAL",
  "BEDDINGMATERIAL", "TYPEOFFLUIDTRANS", "flowline_match_distance_m")
cols_to_drop <- intersect(columns_to_remove, colnames(flowlines))
cat("Dropping these existing columns:", paste(cols_to_drop, collapse = ", "), "\n")
flowlines <- flowlines %>% select(-all_of(cols_to_drop))

# add risk column
flowlines$risk <- 0

# new order
new_order <- c(
  "unique_id", "Operator", "OPERATOR_NUM", "FLOWLINEID", "LOCATION_ID", "Status", "FLOWLINEACTION",
  "LOCATIONTYPE", "Fluid", "Material", "Diam_in", "Length_ft", "MAXOPPRESSURE",
  "line_age_yr", "CONSTRUCTDATE", "risk", "geometry")
flowlines <- flowlines %>% select(all_of(new_order))

# check
columns_to_check <- c("Status", "FLOWLINEACTION", "LOCATIONTYPE", "Fluid", "Material")

for (column in columns_to_check) {
  if (column %in% colnames(flowlines)) {
    vals <- unique(flowlines[[column]])
    cat("Unique values in", column, ":", paste(vals, collapse = ", "), "\n")
  } else {
    cat(column, ": Column not found in DataFrame.\n")
  }
}

# Mapping
status_mapping <- c(
  "Active" = "Active", "ACTIVE" = "Active", "Actove" = "Active", "Avtive" = "Active",
  "Actve" = "Active", "active" = "Active", "Actuve" = "Active", "Flowline" = "Active",
  "Out of Service" = "Out of Service", "OOS" = "Out of Service", "OutofService" = "Out of Service",
  "Out-of-Service" = "Out of Service", "Out Of Service" = "Out of Service",
  "Out of service" = "Out of Service", "OUT OF SERVICE" = "Out of Service",
  "OOSLAT" = "Out of Service",
  "Abandoned" = "Abandoned", "abandoned" = "Abandoned", "Abandoned in Place" = "Abandoned",
  "ABANDONED" = "Abandoned", "Abandon" = "Abandoned", "Abadnon" = "Abandoned",
  "TA" = "Abandoned", "Abandonded" = "Abandoned", "Plugged" = "Abandoned", "AB" = "Abandoned",
  "Inactive" = "Inactive", "InActive" = "Inactive", "INACTIVE" = "Inactive", "IDLE" = "Inactive",
  "Pending Analysis" = "Pending Analysis", "PA" = "Pending Analysis",
  "ABiP" = "Pending Analysis", "ABIP" = "Pending Analysis", "AIP" = "Pending Analysis",
  "Shut In" = "Shut In", "shut in" = "Shut In", "Shut in" = "Shut In", "SI" = "Shut In",
  "Pre-Abandonment" = "Pre-Abandonment", "Pre Abandonment" = "Pre-Abandonment",
  "PreAbandonment" = "Pre-Abandonment", "Pre-Abandoned" = "Pre-Abandonment",
  "Removed" = "Removed", "REMOVED" = "Removed", "REMOVED DUPLICATE" = "Removed",
  "Partial Removed see comment" = "Removed",
  "New Construction" = "New Construction",
  "PreCommission" = "Pre-Commissioned", "Pre-commissioned" = "Pre-Commissioned",
  "Pre-Commissioned" = "Pre-Commissioned",
  "Flushed and capped" = "Decommissioned",
  "Flow filled and capped" = "Decommissioned",
  "Discontinued" = "Decommissioned",
  "Injection" = "Injection", "Injecting" = "Injection",
  "Future" = "Planned", "Planned" = "Planned",
  "Nitrogen" = "DROP_ROW",
  "Air Line" = "DROP_ROW",
  "3rd Party Line" = "DROP_ROW",
  "Turned into Air Line" = "DROP_ROW",
  "Unknown" = NA, "Status" = NA, "unk" = NA, "0" = NA
)
flowlines <- flowlines %>%
  mutate(Status = recode(Status, !!!status_mapping)) %>%
  filter(Status != "DROP_ROW")
print(unique(flowlines$Status))

flowlineaction_mapping <- c(
  "Out of Service" = "Out of Service",
  "Removed From Service" = "Out of Service",
  "Pre-Abandonment Notice" = "Pre-Abandonment Notice",
  "Abandonment Verification" = "Abandonment",
  "Realignment" = "Realignment",
  "Registration" = "Registration",
  "Abandonment" = "Abandonment"
)
flowlines <- flowlines %>%
  mutate(FLOWLINEACTION = recode(FLOWLINEACTION, !!!flowlineaction_mapping))
print(unique(flowlines$FLOWLINEACTION))

locationtype_mapping <- c(
  "Production Facilities" = "Production Facilities",
  "Well Site" = "Well Site",
  "Manifold" = "Manifold",
  "Compressor Station" = "Compressor Station",
  "Gathering Line" = "Gathering Line",
  "Crude Oil Transfer Line" = "Crude Oil Transfer Line",
  "Produced Water Transfer System" = "Produced Water Transfer System"
)
flowlines <- flowlines %>%
  mutate(LOCATIONTYPE = recode(LOCATIONTYPE, !!!locationtype_mapping))
print(unique(flowlines$LOCATIONTYPE))

fluid_mapping <- c(
  "Natual Gas" = "Natural Gas", "Natural Gas Production" = "Natural Gas",
  "Natural Gas Lift" = "Natural Gas", "Natuarl Gas" = "Natural Gas",
  "Natural Gas High Pressure" = "Natural Gas", "Natural Gas Supply" = "Natural Gas",
  "Natural Gas Sales" = "Natural Gas", "Supply Gas" = "Natural Gas",
  "Fuel Gas" = "Natural Gas",
  "Oil" = "Crude Oil", "Crude Oil" = "Crude Oil",
  "Crude Oil Unprocessed" = "Crude Oil",
  "Crude Oil Emulsion" = "Crude Oil Emulsion",
  "Crude Oil Emmulsion, Water And Oil" = "Crude Oil Emulsion",
  "Crude Oill Emulsion" = "Crude Oil Emulsion",
  "Crude Oil And Water Emulsion" = "Crude Oil Emulsion",
  "Oil Water Emulsion" = "Crude Oil Emulsion",
  "Oil/Water" = "Crude Oil Emulsion", "Oil Water" = "Crude Oil Emulsion",
  "Oil And Water" = "Crude Oil Emulsion", "Emulsion" = "Crude Oil Emulsion",
  "Crude Oil / Produced Water" = "Crude Oil Emulsion",
  "Produced Water" = "Produced Water", "Water" = "Produced Water",
  "Saltwater" = "Produced Water", "Injection Produced Water" = "Produced Water",
  "Produced  Water" = "Produced Water", "Produced/Waste Water" = "Produced Water",
  "Waste Water/Produced Water/Formation Water" = "Produced Water",
  "Fresh Water" = "Produced Water", "(Other) Treated Produced Water" = "Produced Water",
  "Co2" = "Co2/Produced Water", "C02/Prod Water" = "Co2/Produced Water",
  "Co2/Prod Water" = "Co2/Produced Water", "Co2Produced Water" = "Co2/Produced Water",
  "Co2/Produced Wtaer" = "Co2/Produced Water", "Co2 Produced Water" = "Co2/Produced Water",
  "Co2/Produced" = "Co2/Produced Water",
  "3 Phase" = "Multiphase", "Multiphase" = "Multiphase", "Multi-Phase" = "Multiphase",
  "Mulitphase" = "Multiphase", "Multi Phase" = "Multiphase", "Mulit Phase" = "Multiphase",
  "Mutliphase" = "Multiphase", "Oil, Water, Gas" = "Multiphase", "Oil-Gas-Water" = "Multiphase",
  "Gas, Oil And Water" = "Full Well Stream", "Oil /Water/Gas" = "Full Well Stream",
  "Oil/Gas/Water" = "Full Well Stream", "Oil, Gas, Water" = "Full Well Stream",
  "Gas,  Oil And Water" = "Full Well Stream", "Natural Gas/Condensate/Produced Water" = "Full Well Stream",
  "Condensate" = "Condensate",
  "Poly" = "Polymer Fluids", "Polymer fluids" = "Polymer Fluids",
  "Liquid" = "Other", "Liquids (Wtr/Cond)" = "Other",
  "Unprocessed Production Fluids" = "Other", "Production Fluids" = "Other",
  "Produced Fluids" = "Other", "Other" = "Other",
  "Other (Natural Gas Liquid)" = "Other", "Other- Nitrogen/ Natural Gas" = "Other",
  "Service Line" = "Other", "Air" = "Other", "Cox V2" = "Other",
  "Steel" = "Other", "Fluid" = "Other", "Vent" = "Other",
  "Unk" = NA, "Unknown" = NA
)
flowlines <- flowlines %>%
  mutate(
    Fluid = str_to_title(str_trim(Fluid)),
    Fluid = recode(Fluid, !!!fluid_mapping)
  )
print(sort(unique(na.omit(flowlines$Fluid))))

material_mapping <- c(
  "Carbon Steel" = "Carbon Steel", "Carbonsteel" = "Carbon Steel", "Carbon  Steel" = "Carbon Steel",
  "Carbon Steel Sch 80" = "Carbon Steel", "Carbon Steel - Hdpe" = "Carbon Steel/HDPE",
  "Carbon Steel-Hdpe" = "Carbon Steel/HDPE", "Carbon Steel And Hdpe" = "Carbon Steel/HDPE",
  "Carbon Steel/Hdpe" = "Carbon Steel/HDPE", "Carbon Steel/HDPE" = "Carbon Steel/HDPE",
  "Hdpe/Steel" = "Carbon Steel/HDPE", "HDPE/Steel" = "Carbon Steel/HDPE",
  "Hdpe Lined Steel" = "Carbon Steel/HDPE", "Hdpe/Steel, Flexsteel" = "Carbon Steel/HDPE/Flexsteel",
  "Carbon Steel, Hdpe,Stainless Steel" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel, Hdpe, Stainless Steel" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel/Stainless Steel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel/Hdpe/Stainless" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel/Hdpe/Stainless Steel" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel/Stainless/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless/Carbon Steel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless/ Carbon Steel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless Steel/Carbon Steel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless/Carbonsteel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Satinless/Carbon Steel/Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless/Carbon Steel/Hspe" = "Carbon Steel/HDPE/Stainless Steel",
  "Stainless/Carbon Steel/ Hdpe" = "Carbon Steel/HDPE/Stainless Steel",
  "Carbon Steel/Fiberglass" = "Fiberglass/Carbon Steel",
  "Carbon Steel Mixed With Fiberglass" = "Fiberglass/Carbon Steel",
  "Fiberglass" = "Fiberglass", "Fibergalss" = "Fiberglass", "Fiberspar" = "Fiberglass",
  "Fiber Glass" = "Fiberglass", "Flexspar" = "Fiberglass", "Fiberoptic" = "Fiberglass",
  "Fiberglass & Fiberspar" = "Fiberglass", "Fiberglass And Carbon Steel" = "Fiberglass/Carbon Steel",
  "Fiberglass Sleaved W/ Hdpe" = "Fiberglass/HDPE", "Fiberglass And Hdpe" = "Fiberglass/HDPE",
  "Hdpe/Fiberglass" = "Fiberglass/HDPE",
  "Steel" = "Steel", "Stainless" = "Steel", "Stainless Steel" = "Steel",
  "Coated Steel" = "Steel", "Lined Steel" = "Steel", "Flexpipe" = "Steel",
  "Flex Steel" = "Steel", "Flexsteel" = "Steel",
  "Other (Flexsteel)" = "Steel", "Other (Flex Steel)" = "Steel", "Other (Flex Pipe)" = "Steel",
  "Hdpe" = "HDPE", "Hdpe Poly" = "HDPE", "Hdpi Poly" = "HDPE",
  "Hdpe Line Sdr 7" = "HDPE", "Hdpe Sdr7" = "HDPE", "Hdpe Sdr11" = "HDPE",
  "Composite Hdpe" = "HDPE", "High-Density Polyethylene (Hdpe)" = "HDPE",
  "Poly" = "Polycarbonate", "Polyline" = "Polycarbonate", "Polycarbonate" = "Polycarbonate",
  "Other (Poly)" = "Polycarbonate", "Sdr 7 Poly" = "Polycarbonate",
  "Poly & Steel" = "Polycarbonate/Steel", "Steel/Poly" = "Polycarbonate/Steel",
  "Poly/Steel" = "Polycarbonate/Steel", "Polycarbonate/Steel" = "Polycarbonate/Steel",
  "Poly Sleeved Steel" = "Polycarbonate/Steel", "Hdpe Poly Sdr 11" = "Polycarbonate/HDPE",
  "Sdr7 Polyethelyne" = "Polyethylene", "Sdr 11 Poly Pipe" = "Polyethylene",
  "Sdr 11 Poly" = "Polyethylene", "Poly Pipe" = "Polyethylene",
  "Sdr_Poly" = "Polyethylene", "Sdr-11" = "Polyethylene", "Sdr-11 Poly" = "Polyethylene",
  "Poly Sdr 7" = "Polypropylene", "Poly Sdr-7" = "Polypropylene",
  "Pvc" = "PVC",
  "Duplex" = "Duplex",
  "Zaplock" = "Other",
  "Plastic" = "Other", "Polypipe" = "Other", "Core Linepipe" = "Other", "Flex Pipe" = "Other",
  "Shawcor Fp150" = "Other", "Zapock" = "Other", "Other (Hdpe And Tubing)" = "Other",
  "Other (Please Specify)" = "Other", "Other (Stainless Steel)" = "Other",
  "Fplp" = "Other", "Oil" = "Other", "Gas" = "Other", "Co2/Produced Water" = "Other",
  "Other" = "Other",
  "Unknown" = NA, "Other (Unknown)" = NA, "Unk" = NA, "0" = NA,
  "Material" = NA, "Sdr" = NA, "Flowline" = NA
)
flowlines <- flowlines %>%
  mutate(
    Material = str_to_title(str_trim(Material)),
    Material = recode(Material, !!!material_mapping)
  )
print(unique(flowlines$Material))

# rename columns
flowlines <- flowlines %>%
  rename(
    unique_id = unique_id,
    operator_name = Operator,
    operator_number = OPERATOR_NUM,
    flowline_id = FLOWLINEID,
    location_id = LOCATION_ID,
    status = Status,
    flowline_action = FLOWLINEACTION,
    location_type = LOCATIONTYPE,
    fluid = Fluid,
    material = Material,
    diameter_in = Diam_in,
    length_ft = Length_ft,
    max_operating_pressure = MAXOPPRESSURE,
    line_age_yr = line_age_yr,
    construct_date = CONSTRUCTDATE,
    geometry = geometry
  )

#------SPILLS-------
spills <- st_read("spills.geojson")
colnames(spills)

# Operator name cleanup
if ("Operator" %in% colnames(spills)) {
  spills$Operator <- dplyr::recode(spills$Operator, !!!op_name_fix)
}

# drop columns
columns_to_remove <- c(
  "trkg_num", "facility_type", "Spill_Desc", "Spill.Type", "Root.Cause",
  "Preventative.Measure", "Detailed.Root.Cause.Type", "facility_status", "Gathering.", "Metallic.")
cols_to_drop <- intersect(columns_to_remove, colnames(spills))
cat("Dropping these existing columns:", paste(cols_to_drop, collapse = ", "), "\n")
spills <- spills %>% select(-all_of(cols_to_drop))

# add risk column
spills$risk <- 1

# new order
new_order <- c(
  "Operator.Name", "Root.Cause.Type", "incident_date", "risk", "Long", "Lat", "geometry")
spills <- spills %>% select(all_of(new_order))

# rename columns
spills <- spills %>%
  rename(
    operator_name = Operator.Name,
    root_cause = Root.Cause.Type,
    lon = Long,
    lat = Lat
  )

# check:
root_cause_mapping <- c(
  "Equipment failure" = "Equipment Failure"
)
spills <- spills %>%
  mutate(root_cause = recode(root_cause, !!!root_cause_mapping))
print(unique(spills$root_cause))

# missing geometry 
empties <- st_is_empty(spills)
if (!any(empties)) stop("No empty geometries to fill.")

df <- st_drop_geometry(spills[empties, ])
df$lon <- suppressWarnings(as.numeric(df$lon))
df$lat <- suppressWarnings(as.numeric(df$lat))

swap <- which((abs(df$lat) > 90 & abs(df$lon) <= 180) |
                (abs(df$lon) > 180 & abs(df$lat) <= 90))
if (length(swap)) {
  tmp <- df$lon[swap]; df$lon[swap] <- df$lat[swap]; df$lat[swap] <- tmp
}

ok <- is.finite(df$lon) & is.finite(df$lat) &
  df$lon >= -180 & df$lon <= 180 &
  df$lat >=  -90 & df$lat <=  90
df <- df[ok, , drop = FALSE]
fill_idx <- which(empties)[ok]

pts <- st_as_sf(df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
pts <- st_transform(pts, st_crs(spills))
st_geometry(spills)[fill_idx] <- st_geometry(pts)

sum(st_is_empty(spills))

# drop NA's None's 0's
spills <- na.omit(spills)

# Save as GeoJSON
st_write(spills, "clean_spills.geojson", driver = "GeoJSON", delete_dsn = TRUE)
st_write(flowlines, "clean_combined_flowlines.geojson", driver = "GeoJSON", delete_dsn = TRUE)
cat("# flowlines rows:", nrow(flowlines), "\n")
cat("# spills rows:", nrow(spills), "\n")
