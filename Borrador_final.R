library(sf)
library(tidyverse)
library(tmap)
library(ez)
library(RColorBrewer)
library (sp)
library(spdep) #Para crear vecinos y pesos
source('lisaclusters.R')
library(gstat)

(datos <- st_read('paramsoutlet_orden1.gpkg', crs = 32619))
(datos <- datos %>% st_difference())
(pol1 <- st_read(dsn = 'r_stream_basins_1.geojson', crs = 32619))
pol2 <- st_read(dsn = 'r_stream_basins_2.geojson', crs = 32619)
pol3 <- st_read(dsn = 'r_stream_basins_3.geojson', crs = 32619)
pol4 <- st_read(dsn = 'r_stream_basins_4.geojson', crs = 32619)

datos %>% dplyr::filter(Max_order_Strahler==1)

datos %>%
  select_if(is.numeric) %>%
  gather(variable, valor, -geom) %>%
  st_drop_geometry() %>% 
  group_by(variable) %>% 
  summarise(m=mean(valor, na.rm=T))

datos %>%
  select_if(is.numeric) %>%
  gather(variable, valor, -geom) %>% 
  tm_shape() + tm_dots(col = 'valor') + tm_facets(by='variable', free.coords = F, free.scales = T)

#Tabla cols numericas, con varianza y sin NA
datosnum <- datos %>%
  st_drop_geometry() %>% 
  select_if(is.numeric) %>% 
  select_if(~ sum(!is.na(.))>0) %>% 
  select_if(function(x) var(x, na.rm=T)!=0)

datosnum %>% ezCor(r_size_lims = 2:3, label_size = 2)

datosnum %>% cor

# Union espacial
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
  inner_join(xy)
Varselpol3

# VECINDAD



Vecxcont <- poly2nb(Varselpol3)
summary(Vecxcont)

Varselpol3.sp <- as_Spatial(Varselpol3)
coords <- coordinates(Varselpol3.sp)
VecxK <- knn2nb(knearneigh(coords, k=5))

PesoW <- nb2listw(VecxK)

datos <- datos %>% st_difference()
coords <- coordinates(as_Spatial(datos))
nb <- knn2nb(knearneigh(coords, k = 5))
summary(nb)
wW <- nb2listw(nb)

#ESDA

p1 <- tm_shape(Varselpol3) +
  tm_fill(col = "DD", style = 'jenks', palette = brewer.pal(9, name = 'Reds')) +
  tm_borders(lwd = 0.5)
p1

p2 <- tm_shape(Varselpol3) +
  tm_fill(col = "logDD", style = 'jenks',
          palette = brewer.pal(9, name = 'Reds'), midpoint = NA) +
  tm_borders(lwd = 0.5)
tmap_arrange(p1, p2)

Varselpol3 %>% st_drop_geometry() %>%
  gather(variable, valor, -(fid:label)) %>%
  ggplot() + aes(sample=valor) +
  stat_qq() + stat_qq_line() + theme_bw() +
  theme(text = element_text(size = 14)) +
  facet_wrap(~variable, scales = 'free')

Varselpol3 %>% st_drop_geometry() %>%
  gather(variable, valor, -(fid:label)) %>% group_by(variable) %>%
  summarise(prueba_normalidad=shapiro.test(valor)$p.value)


gmoranw <- moran.test(x = log1p(datos$Drainage_Density_km_over_km2), listw = wW)
gmoranw


lisamap(objesp = Varselpol3,
        var ='logDD',
        pesos = PesoW,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 14,
        fuentedatos = 'ENHOGAR 2017',
        titulomapa = paste0('Clusters LISA de DD'))



