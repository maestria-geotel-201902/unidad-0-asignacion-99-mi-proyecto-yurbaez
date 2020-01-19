## Universidad Autonoma de Santo Domingo
## Maestría Teledeteccion y Ciencias de la Informacion Geografica
## Materia: Analisis Espacial
## Profesor: José Ramón Martínez
## USO DE VARIABLES MORFOMETRICAS EN EL ANALISIS DE LA DENSIDAD DE DRENAJE
## DE LAS MICROCUENCAS DE ORDEN 1 (STRAHLER) DEL RIO OCOA 
## EN LA REPUBLICA DOMINICANA
## Maestrantes: Alba Cadete, Mirel Volcán, Yoenny Urbáez

## LIBRERIAS A UTILIZAR
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

## IMPORTACION, ORGANIZACION DE DATOS E INTEROPERABILIDAD

# Cargar datos de variables 
(datos <- st_read('paramsoutlet_orden1.gpkg', crs = 32619))
(datos <- datos %>% st_difference())
(pol1 <- st_read(dsn = 'r_stream_basins_1.geojson', crs = 32619))
pol2 <- st_read(dsn = 'r_stream_basins_2.geojson', crs = 32619)
pol3 <- st_read(dsn = 'r_stream_basins_3.geojson', crs = 32619)
pol4 <- st_read(dsn = 'r_stream_basins_4.geojson', crs = 32619)

# Orden de Red Cuencas 1, clasificacion de Strahler
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

# Tabla cols numericas, con varianza 
datosnum <- datos %>%
  st_drop_geometry() %>% 
  select_if(is.numeric) %>% 
  select_if(~ sum(!is.na(.))>0) %>% 
  select_if(function(x) var(x, na.rm=T)!=0)

# Evaluacion de correlacion entre las variables como criterio de seleccion
datosnum %>% ezCor(r_size_lims = 2:3, label_size = 2)
datosnum %>% cor

# Union espacial de las variables seleccionadas
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

## Creación de objeto XY con atributos del objeto Varselpol2 
## mediante el centroide de los polígonos
xy <- Varselpol2 %>%
  st_centroid() %>% 
  mutate(x=unlist(map(geometry,1)),
         y=unlist(map(geometry,2))) %>% 
  st_drop_geometry() %>% 
  select(fid, x, y)

## Creación del objeto Varselpol3 mediante unión de XY y Varselpol2
Varselpol3 <- Varselpol2 %>% 
  inner_join(xy)
Varselpol3

# VECINDAD

##  Analisis de vecindad por contiguidad
Vecxcont <- poly2nb(Varselpol3)
summary(Vecxcont)

## Análisis de Vecindad por cantidad de los 5 vecinos más cercanos
Varselpol3.sp <- as_Spatial(Varselpol3)
coords <- coordinates(Varselpol3.sp)
VecxK <- knn2nb(knearneigh(coords, k=5))

## Análisis de Vecindad por peso de observaciones vecinas en Varselpol3
PesoW <- nb2listw(VecxK)
PesoW

PesowB <- nb2listw(VecxK, style = 'B')
PesowB

## Análisis de Vecindad por peso de observaciones vecinas en la data completa
datos <- datos %>% st_difference()
coords <- coordinates(as_Spatial(datos))
nb <- knn2nb(knearneigh(coords, k = 5))
summary(nb)


# ESDA

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

lisamap(objesp = Varselpol3,
        var ='logDD',
        pesos = PesoW,
        tituloleyenda = 'Significancia\n("x-y", léase\ncomo "x"\nrodeado de "y"',
        leyenda = T,
        anchuratitulo = 1000,
        tamanotitulo = 16,
        fuentedatos = 'SRTM',
        titulomapa = paste0('Clusters LISA de logDD'))

## MODELIZACION

#Para la modelización se requiere cargar la siguiente librería de tidyverse, sf, spdep, lmtest, 
## además de cargar los objetos generados en los análisis de Vencindad y autocorrelación.

#Varselpol3, sf con los datos fuente.
#PesoW, pesos estandarizados por filas (estilo “W”).
#PesoB, pesos binarios (estilo “B”).

#Variable seleccionada Densidad de Drenaje 

Varselpctlog <- Varselpol3 %>% mutate_each(
  funs(PCT=round(./DD,4)*100,
       PCTLOG=log(round(./DD,4)*100)),
  -1, -2, -geometry, -label)

Varselpctlog 

# Comprobando autocorrelación mediante la prueba moran global

moran.plot(Varselpol3$logDD, PesoW)

(gmoranw <- moran.test(na.action = na.exclude, zero.policy = T, x = log1p(datos$Drainage_Density_km_over_km2), listw = PesoW))

(gmoranb <- moran.test(na.action = na.exclude, zero.policy = T, x = log1p(datos$Drainage_Density_km_over_km2), listw = PesowB))

gmoranb <- moran.test(na.action = na.exclude, zero.policy = T, x = Varselpctlog$logDD_PCT, listw = PesowB)
gmoranb

gmoranwl <- moran.test(na.action = na.exclude, zero.policy = T, x = Varselpctlog$logDD_PCTLOG, listw = PesoW)
gmoranwl

gmoranbl <- moran.test( na.action = na.exclude, zero.policy = Tx = Varselpctlog$logDD_PCTLOG, listw = PesowB)
gmoranbl

(gmoranwale<-moran.test(na.action = na.exclude, zero.policy = T, x=rnorm(3027),listw = PesoW))

#Si el valor de p es inferior al nivel de significancia (comúnmente fijado en 0.05 o 0.01), se rechaza la hipótesis nula “No hay autocorrelación espacial”. Por lo tanto, concluimos que hay, a priori, autocorrelación espacial, tanto para la variable original (sufijo _PCT) como la versión transformada (sufijo _PCTLOG).
#Evaluemos si el supuesto de normalidad se cumple:

shapiro.test(Varselpctlog$logDD_PCT)
shapiro.test(Varselpctlog$logDD_PCTLOG)


# Modelo lineal común, utilizando las versiones transformadas de las variables, evaluamos homocedasticidad

# Modelo lineal, utilizando todas las variables

modlin <- Varselpol3 %>%
  select(logDD, TS, MS, SF, ER) %>%
  st_drop_geometry() %>%
  lm(logDD ~ ., .)
modlin %>% summary 

modlinc <- Varselpctlog %>%
  select(contains('_PCTLOG')) %>%
  st_drop_geometry() %>%
  lm(logDD_PCTLOG ~ ., .)
modlinc %>% summary

modlinc %>% bptest

sar <- Varselpctlog %>% select(contains('_PCTLOG')) %>%
  st_drop_geometry() %>%
  spautolm(
    formula = logDD_PCTLOG ~ .,
    data = .,
    listw =PesoW)
summary(sar)

sar2 <- Varselpctlog %>% select(contains('_PCTLOG')) %>%
  st_drop_geometry() %>%
  spautolm(
    formula = logDD_PCTLOG ~ TS_PCTLOG + MS_PCTLOG + SF_PCTLOG + ER_PCTLOG,
    data = .,
    listw = PesoW)
summary(sar2)

Sar3 <- Varselpctlog %>% select(contains('_PCTLOG')) %>%
  st_drop_geometry() %>%
  spautolm(
    formula = logDD_PCTLOG ~ TS_PCTLOG + MS_PCTLOG + SF_PCTLOG,
    data = .,
    listw = PesoW)
summary(Sar3)

## AUTOCORRELACION ESPACIAL LOCAL

Varselpol_lomo <- localmoran(Varselpctlog$'logDD', listw = PesoW)
summary(Varselpol_lomo)

Varselpctlog$sVarselpctlogDD <- scale (Varselpctlog$'logDD') %>% as.vector()

Varselpctlog$laglogDD <-lag.listw(PesoW, Varselpctlog$'logDD')

summary(Varselpctlog$sVarselpctlogDD)

summary(Varselpctlog$laglogDD)

puntz <- Varselpctlog$sVarselpctlogDD
rezag <- Varselpctlog$laglogDD
df <- data.frame(puntz, rezag)

moran.plot(puntz, PesoW)

# Diagrama de dispersión de Moran en ggplot

ggplot(df, aes(puntz, rezag)) +
  geom_point() + geom_smooth(method = 'lm', se = F) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')


# Variable nueva sobre significancia de la correlación local, rellena con NAs
Varselpctlog$quad_sig <- NA


# Cuadrante high-high quadrant
Varselpctlog[(Varselpctlog$sVarselpctlog >= 0 &
                Varselpctlog$laglogDD >= 0) &
               (Varselpol_lomo[, 4] <= 0.05), "quad_sig"] <- "high-high"


# Cuadrante low-low
Varselpctlog[(Varselpctlog$sVarselpctlog <= 0 &
                Varselpctlog$laglogDD >= 0) & 
               (Varselpol_lomo[, 4] <= 0.05), "quad_sig"] <- "low-low"

# Cuadrante high-low
Varselpctlog[(Varselpctlog$sVarselpctlog >= 0 & 
                Varselpctlog$laglogDD <= 0) & 
               (Varselpol_lomo[, 4] <= 0.05), "quad_sig"] <- "high-low"

# Cuadrante low-high
Varselpctlog[(Varselpctlog$sVarselpctlog <= 0 & 
                Varselpctlog$laglogDD >= 0) & 
               (Varselpol_lomo[, 4] <= 0.05), "quad_sig"] <- "low-high"

# No significativas
Varselpctlog[(Varselpol_lomo[, 4] > 0.05), "quad_sig"] <- "not signif."

#Convertir a factorial
Varselpctlog$quad_sig <- as.factor(Varselpctlog$quad_sig)

# Mapa Significancia
Varselpctlog %>% 
  ggplot() +
  aes(fill = quad_sig) + 
  geom_sf(color = "white", size = .05)  +
  theme_void() + scale_fill_brewer(palette = "Set1")




## KRIGING ORDINARIO

## Sistema de Coordenadas WGS84 UTM Zona 19
## EPSG: 32619

crsdestino <- 32619

## Los variogramas muestrales 
## A partir del variograma muestral, generamos un variograma modelo que será el que utlizará
## la función krige para realizar la interpolación

## Creamos el objeto orden1logdd que es el resultado de la variable transformada logDD 

orden1logdd <- datos %>% mutate(logDD=log(Drainage_Density_km_over_km2)) %>% select(logDD)

## Creamos el objeto v para representar el variograma v
v <- variogram(logDD~1, orden1logdd)
plot(v)

## Variograma Modelo se utlizará la función krige para realizar la interpolación
## Se necesita un variograma que crezca de inmediato, por esto vamos a utilizar el variograma modelo 

## Variograma modelo Esferico con rango de 1000 metros

v_m <- fit.variogram(v, vgm(model = "Sph", range = 1000))
v_m
plot(v, v_m)

## model     psill     range
## Sph    0.1658994   988.8336


## Variograma modelo Exponencial con rango de 1000 metros

v_m2 <- fit.variogram(v, vgm(model = "Exp", range = 1000))
v_m2

##   model     psill    range
##   Exp    0.1658997  219.8573

plot(v, v_m2)


## Variograma modelo Gauseano

v_m3 <- fit.variogram(v, vgm(model = "Gau", range = 1000))
v_m3

##   model     psill    range
##    Gau    0.165893 315.4503

plot(v, v_m3, plot.numbers = T)
plot(v, v_m3)

attr(v_m, 'SSErr')
## [1] 2.496814e-07

attr(v_m2, 'SSErr') 
## [1] 2.505845e-07

attr(v_m3, 'SSErr')
## [1] 6.398799e-07

grd <- st_bbox(orden1logdd) %>%
  st_as_stars(dx = 100) %>% #1000 metros=1km de resolución espacial
  st_set_crs(crsdestino)
grd

plot(grd)

## 
k <- krige(formula = logDD~1, locations = orden1logdd, newdata = grd, model = v_m)

plot(k)


## Utilicemos ggplot para representar el objeto stars.

ggplot() +
  geom_stars(data = k, aes(fill = var1.pred, x = x, y = y)) + 
  scale_fill_gradient(low="#deebf7", high="#3182bd") +
  geom_sf(data = st_cast(Varselpol3, "MULTILINESTRING")) +
  geom_sf(data = orden1logdd) +
  geom_sf_text(data = Varselpol3, aes(label=''), check_overlap = T, size = 1) +
  theme_bw()

ggplot() +
  geom_stars(data = exp(k), aes(fill = var1.pred, x = x, y = y)) + 
  scale_fill_gradient(low="#deebf7", high="#3182bd", trans = 'log10') +
  geom_sf(data = st_cast(Varselpol3, "MULTILINESTRING")) +
  geom_sf(data = orden1logdd) +
  geom_sf_text(data = Varselpol3, aes(label=''), check_overlap = T, size = 1) +
  theme_bw()
