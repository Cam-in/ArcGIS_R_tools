####
## Original script provided by Pablo Guerrero (BASE, UDEC, IEB)
###

# Cargar paquetes
library(terra)
library(ggplot2)

# Rutas a capas ambientales
path_2000 <- "/Users/pabloguerrero/Library/CloudStorage/Dropbox/0 paper Chionis/0 revision 1/layers correlation/2000"
path_2010 <- "/Users/pabloguerrero/Library/CloudStorage/Dropbox/0 paper Chionis/0 revision 1/layers correlation/2010"

# Rutas a archivos de ocurrencia
occ_minor_path <- "/Users/pabloguerrero/Library/CloudStorage/Dropbox/0 paper Chionis/0 revision 1/layers correlation/OCC_C_minor_FilteredwYear.csv"
occ_albus_path <- "/Users/pabloguerrero/Library/CloudStorage/Dropbox/0 paper Chionis/0 revision 1/layers correlation/OCC_C_albus_FilteredwYear.csv"

# Cargar coordenadas
occ_minor <- read.csv(occ_minor_path)
occ_albus <- read.csv(occ_albus_path)

# Confirmar nombres correctos de columnas
coord_names <- c("Longitude", "Latitude")

# Función para cargar la primera capa de cada archivo .nc
load_nc_layer <- function(path, file) {
  rast(file.path(path, file))[[1]]
}

# Variables
vars <- c("phyc", "so", "sws", "thetao")

# Obtener archivos .nc
files_2000 <- list.files(path_2000, pattern = "\\.nc$", full.names = FALSE)
files_2010 <- list.files(path_2010, pattern = "\\.nc$", full.names = FALSE)
names(files_2000) <- sapply(strsplit(files_2000, "_"), `[`, 1)
names(files_2010) <- sapply(strsplit(files_2010, "_"), `[`, 1)

# Correlaciones y ANOVA
cor_results <- list()
anova_results <- list()
all_values <- data.frame()  # para almacenar valores extraídos

# Función para extraer valores por especie
extract_species <- function(occ_data, species_name, r2000, r2010, varname) {
  if (!all(coord_names %in% names(occ_data))) {
    stop(paste("ERROR: Las columnas 'Longitude' y 'Latitude' no están presentes en", species_name))
  }
  pts <- vect(cbind(occ_data$Longitude, occ_data$Latitude), type = "points", crs = crs(r2000))
  val_2000 <- extract(r2000, pts)[,2]
  val_2010 <- extract(r2010, pts)[,2]
  data.frame(
    species = species_name,
    variable = varname,
    decade = rep(c("2000", "2010"), each = length(val_2000)),
    value = c(val_2000, val_2010)
  )
}

# Loop por variable ambiental
for (v in vars) {
  r2000 <- load_nc_layer(path_2000, files_2000[v])
  r2010 <- load_nc_layer(path_2010, files_2010[v])
  
  # Alinear espacialmente
  common_ext <- intersect(ext(r2000), ext(r2010))
  r2000_crop <- crop(r2000, common_ext)
  r2010_crop <- crop(r2010, common_ext)
  if (!compareGeom(r2000_crop, r2010_crop, stopOnError = FALSE)) {
    r2010_crop <- resample(r2010_crop, r2000_crop)
  }
  
  # Correlación global
  vals_2000 <- values(r2000_crop)
  vals_2010 <- values(r2010_crop)
  idx <- complete.cases(vals_2000, vals_2010)
  cor_val <- cor(vals_2000[idx], vals_2010[idx], method = "pearson")
  cor_results[[v]] <- cor_val
  
  # Extraer valores y hacer ANOVA por especie
  df_minor <- extract_species(occ_minor, "C_minor", r2000_crop, r2010_crop, v)
  df_albus <- extract_species(occ_albus, "C_albus", r2000_crop, r2010_crop, v)
  df_combined <- rbind(df_minor, df_albus)
  all_values <- rbind(all_values, df_combined)  # guardar para graficar
  
  for (sp in unique(df_combined$species)) {
    df_sp <- subset(df_combined, species == sp)
    res <- aov(value ~ decade, data = df_sp)
    anova_results[[paste(v, sp, sep = "_")]] <- summary(res)
  }
}

# Mostrar correlaciones globales
print(cor_results)

# Mostrar resultados ANOVA
for (k in names(anova_results)) {
  cat("\n---", k, "---\n")
  print(anova_results[[k]])
}

# Generar gráfico comparativo
ggplot(all_values, aes(x = decade, y = value, fill = decade)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_grid(variable ~ species, scales = "free_y") +
  theme_minimal(base_size = 12) +
  labs(
    title = "Environmental values at species occurrences by decade",
    x = "Decade",
    y = "Environmental value"
  ) +
  theme(legend.position = "none")
