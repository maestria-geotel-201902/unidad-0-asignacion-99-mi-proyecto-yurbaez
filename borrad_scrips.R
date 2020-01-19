library(sf)
library(tidyverse)
library(gstat)
library(stars)
library(tmap)
library(ez)
library(RColorBrewer)
library (sp)
library(spdep) 
library(lmtest)
library(spData)
source('lisaclusters.R')

## cargar los datos desde el gpkg ....
(datos <- st_read('paramsoutlet_orden1.gpkg', crs = 32619))
(datos <- datos %>% st_difference())
(pol1 <- st_read(dsn = 'r_stream_basins_1.geojson', crs = 32619))
pol2 <- st_read(dsn = 'r_stream_basins_2.geojson', crs = 32619)
pol3 <- st_read(dsn = 'r_stream_basins_3.geojson', crs = 32619)
pol4 <- st_read(dsn = 'r_stream_basins_4.geojson', crs = 32619)


datos %>% dplyr::filter(Max_order_Strahler==1)

## union espacial

VARSEL <- datos %>% 
  dplyr::select(
    DD = Drainage_Density_km_over_km2,
    SF = Shape_Factor,
    ER = Elongation_Ratio,
    TS = Total_Stream_Length_km,
    MS = Mean_Slope
  )

Varselpol1 <- pol1 %>% st_join(left = F, VARSEL)
Varselpol2 <- Varselpol1 %>% 
  mutate(logDD = log(DD))

xy <- Varselpol2 %>%
  st_centroid() %>% 
  mutate(x=unlist(map(geometry,1)),
         y=unlist(map(geometry,2))) %>% 
  st_drop_geometry() %>% 
  select(fid, x, y)
Varselpol3 <- Varselpol2 %>% 
  inner_join(xy) %>% kable
Varselpol3


datosnum <- datos %>%
  st_drop_geometry() %>% 
  select_if(is.numeric) %>% 
  select_if(~ sum(!is.na(.))>0) %>% 
  select_if(function(x) var(x, na.rm=T)!=0)

datosnum %>% ezCor(r_size_lims = 2:3, label_size = 2)
datosnum %>% cor