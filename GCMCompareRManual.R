library(sf)
library(raster)
library(tidyverse)
library(sp)
library(rgdal)
library(fasterize)
library(maps)
library(maptools)
library(ggrepel)

#Carpeta de capas

Table_Layers <- read_csv("envidatS3paths.txt", col_names = FALSE) %>% 
  magrittr::set_names("Link") %>% 
  mutate(model = str_remove_all(Link, "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/cmip5/2061-2080/bio/CHELSA_bio_mon_")) %>% 
  tidyr::separate(model, sep = "_", into = c("Model", "RCP", "Coso", "Coso2", "Bio", "Year", "Version")) %>% 
  dplyr::select(Link, Model, RCP, Year, Bio) %>% 
  dplyr::filter(RCP == "rcp85")




Layers <- Table_Layers %>% 
  mutate(Bio = as.numeric(Bio)) %>% 
  dplyr::group_split(Model) %>% 
  purrr::map(~arrange(.x, Bio)) %>% 
  purrr::map(~mutate(.x, Bio = formatC(Bio, digits = 1, flag = "0", mode = "integer")))

Poliogono <- getData(name = "GADM", country = "GRL", level =0)

STACK <- list()
for(j in 1:nrow(Layers[[i]])){
  download.file(Layers[[i]]$Link[j], destfile = paste0("Temp", Layers[[i]]$Bio[j], ".tif"))
  STACK[[j]] <- raster(paste0("Temp", Layers[[i]]$Bio[j], ".tif")) %>%
  crop(Poliogono) %>% magrittr::set_names(paste0("bio", Layers[[i]]$Bio[j]))
  message(paste(j, "of", 19))
}

Mask <- fasterize::fasterize(st_as_sf(Poliogono), STACK[[1]])

STACK <- STACK %>% reduce(stack)

STACK <- STACK*Mask

writeRaster(STACK, paste(unique(Layers[[i]]$Model), unique(Layers[[i]]$RCP),".tif", sep = "_"))

To_erase <- list.files(pattern = "Temp", full.names = T)

file.remove(To_erase)
#Datos de GCMs

table_GCMs <- read_csv("/home/derek/Documents/GCMCompareR2/data/GCMs_details.csv") %>% 
  dplyr::select(-"Actual name")

# Poligono para cortar




## Seleccionar años, uno de "2021.2040","2041.2060","2061.2080", "2081.2100"

year_type <- "2081.2100"

# Seleccionar escenario  "ssp126", "ssp245", "ssp370" "ssp585"

scenario_type <-  "ssp585"

# Resolución 

res_sel <-"10m"

## Tipo de comparación entre "bio_bio" y "bio_several"

compare_type <- "bio_bio"

## Elegir las bios bio1

bio_to_plot <- c("bio1", "bio12")

scenarios <- list.files(paste0(Layers,"clim_data/worldclimCMIP6"), pattern = res_sel) %>%
  dplyr::as_tibble() %>%
  magrittr::set_names("GCM") %>%
  tidyr::separate(col = GCM, into = c("wc_", "resolution_", "bio_", "gcm_", "scenario_", "year_"), sep = "_") %>% dplyr::select(-wc_, -bio_) %>%
  slice(-1)

available_gcms <- scenarios %>%
  filter(resolution_ == res_sel) %>%
  filter(year_ == year_type) %>%
  # filter(scenario_ == input$scenario_type %>% str_replace("rcp", "")) %>%
  filter(scenario_ == scenario_type) %>% #str_replace("rcp", "")) %>%
  pull(gcm_)


clim_vars_files <- paste0("wc2.1_",
                          res_sel,
                          "_bioc_",
                          available_gcms %>% gsub("-", "\\.", .), 
                          "_", scenario_type, "_",
                          year_type %>% gsub("-", "\\.", .)) %>%
  purrr::map(~list.files(paste0(Layers,"clim_data/worldclimCMIP6/", .x),
                         # pattern = .x,
                         full.names = T)) %>%
  purrr::map(~ raster::stack(.x)) %>%
  purrr::map(~setNames(.x, names(.x) %>% gsub(".*\\.bio","",.) %>% paste0("bio",.)))

rclim_vars_complete <- clim_vars_files %>% 
  purrr::map(~raster::subset(.x, bio_to_plot)) %>% 
  magrittr::set_names(available_gcms) %>% 
 purrr::map(~crop(.x, Poliogono)) %>% 
  purrr::map(~mask(.x, Poliogono)) 

clim_baseline_complete <- paste0(Layers,"clim_data/worldclimCMIP6/wc2.1_", res_sel) %>% 
  list.files(full.names = T) %>% stack %>% 
  setNames(names(.) %>% gsub(".*_bio","bio",.)) %>% 
  subset(bio_to_plot) %>% 
  crop(Poliogono) %>% 
  mask(Poliogono)

clim_delta <- rclim_vars_complete %>%
  purrr::map(~ .x - clim_baseline_complete)


clim_deltaperc <- rclim_vars_complete %>%
  purrr::map(~ 100 * (.x - clim_baseline_complete) / clim_baseline_complete)

clim_ens <- rclim_vars_complete %>%
  purrr::reduce(`+`) %>%                                  # Sum all layers
  raster::calc(fun = function(x){x / length(rclim_vars_complete)})  # Divide by the number of layers

clim_delta_ensemble <- clim_ens - clim_baseline_complete
  
clim_diff <- rclim_vars_complete %>%
    purrr::map(~.x - clim_ens)



table_diff_scaled <- list()
for(i in names(rclim_vars_complete[[1]])){
  temp_table <- rclim_vars_complete %>% purrr::map(~ .x[[i]]) %>%
    purrr::map_dfc(~ raster::values(.x)) %>% 
    filter_all(all_vars(!is.na(.))) %>% t() %>% 
    scale() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("GCM") %>%
    tibble::as_tibble()
  table_diff_scaled[[length(table_diff_scaled) + 1]] <- temp_table
  names(table_diff_scaled)[length(table_diff_scaled)] <- i
}

table_diff_scaled <- table_diff_scaled %>%
  purrr::map_dfc(~ rowMeans(.x[,2:ncol(.x)], na.rm = T)) %>%
  dplyr::bind_cols(table_diff_scaled[[1]][,1])
# Combine temperature and precipitation variables separatedly
table_scaled <- tibble::tibble(GCM = table_diff_scaled %>%
                                 dplyr::select(GCM) %>% dplyr::pull(GCM),
                               x_axis = table_diff_scaled %>%
                                 dplyr::select(bio_to_plot[[1]]) %>%
                                 rowMeans(na.rm = T),
                               y_axis = table_diff_scaled %>%
                                 dplyr::select(bio_to_plot[[2]]) %>%
                                 rowMeans(na.rm = T)) %>%
  dplyr::mutate(Distance = raster::pointDistance(c(0, 0),
                                                 .[,2:3],
                                                 lonlat = FALSE)) %>%
  dplyr::arrange(Distance)

std_error <- mean(table_scaled$Distance) + 2*sd(table_scaled$Distance)
table_scaled <- table_scaled %>%
  dplyr::mutate(Within_circle = ifelse(Distance <= std_error, TRUE, FALSE))


#### unscaled table

rclim_vars_complete[[length(rclim_vars_complete) + 1]] <- clim_baseline_complete
rclim_vars_complete[[length(rclim_vars_complete) + 1]] <- clim_ens
names(rclim_vars_complete)[9:10] <- c("Presente","Ensamble")

table_diff_unscaled <- list()
for(i in names(rclim_vars_complete[[1]])){
  temp_table <- rclim_vars_complete %>% purrr::map(~ .x[[i]]) %>%
    purrr::map_dfc(~ raster::values(.x)) %>% 
    filter_all(all_vars(!is.na(.))) %>% t() %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("GCM") %>%
    tibble::as_tibble()
  table_diff_unscaled[[length(table_diff_unscaled) + 1]] <- temp_table
  names(table_diff_unscaled)[length(table_diff_unscaled)] <- i
}

table_diff_unscaled <- table_diff_unscaled %>%
  purrr::map_dfc(~ rowMeans(.x[,2:ncol(.x)], na.rm = T)) %>%
  dplyr::bind_cols(table_diff_unscaled[[1]][,1])
# Combine temperature and precipitation variables separatedly
table_unscaled <- tibble::tibble(GCM = table_diff_unscaled %>%
                                 dplyr::select(GCM) %>% dplyr::pull(GCM),
                               x_axis = table_diff_unscaled %>%
                                 dplyr::select(bio_to_plot[[1]]) %>%
                                 rowMeans(na.rm = T),
                               y_axis = table_diff_unscaled %>%
                                 dplyr::select(bio_to_plot[[2]]) %>%
                                 rowMeans(na.rm = T)) %>% mutate(Tipo = ifelse(GCM == "Presente", "Presente", ifelse(GCM == "Ensamble", "Ensamble", "GCM")))

#### Otra cosa

circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){r = diameter / 2
tt <- seq(0,2*pi,length.out = npoints)
xx <- center[1] + r * cos(tt)
yy <- center[2] + r * sin(tt)
return(data.frame(x = xx, y = yy))}

circle <- circleFun(center = c(0,0), diameter = std_error * 2, npoints = 300)



Carpeta_results <- "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/"

saveRDS(rclim_vars_complete, paste0(Carpeta_results,"GCMs.rds"))

saveRDS(clim_baseline_complete, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Baseline.rds")

saveRDS(clim_delta, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Delta.rds")

saveRDS(clim_deltaperc, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Delta_porc.rds")

saveRDS(clim_ens, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Ensamble.rds")

saveRDS(clim_delta_ensemble, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Delta_Ensamble.rds")

saveRDS(clim_diff, "/home/derek/Documents/Pelambres_Cambio_climatico/Archivos_Pelambres/Delta_diff_ensemble.rds")

saveRDS(circle, paste0(Carpeta_results,"Circle.rds"))

saveRDS(table_scaled, paste0(Carpeta_results,"table_scaled.rds"))

saveRDS(table_unscaled, paste0(Carpeta_results,"table_unscaled.rds"))

