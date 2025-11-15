library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

PATH = "C:/Users/utente/Nextcloud/01_work/01_research/01_projects/02_SPATIAL_COX/package/"


# 1) Read your London boroughs shapefile
london <- st_read(file.path(PATH, "London_Borough_Excluding_MHW.shp"), quiet = TRUE) |>
  st_make_valid()
# 2) Get the UK outline as sf
uk <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf")

# 3) Match CRSs
london <- st_transform(london, st_crs(uk))

# 4) Dissolve boroughs to a single London polygon (optional but nicer)
london_u <- st_union(st_geometry(london))

# 5) Plot: UK base in light gray, London filled/accented
p <- ggplot() +
  geom_sf(data = uk, fill = "grey95", color = "grey70", linewidth = 0.3) +
  geom_sf(data = london_u, fill = "#E41A1C", color = "#E41A1C", linewidth = 0.2) +
  coord_sf(xlim = st_bbox(uk)[c("xmin", "xmax")],
           ylim = st_bbox(uk)[c("ymin", "ymax")],
           expand = FALSE) +
  theme_void()

p

# 6) Save to an image if you want
ggsave("uk_london_highlight.png", p, width = 7, height = 8, dpi = 300)
