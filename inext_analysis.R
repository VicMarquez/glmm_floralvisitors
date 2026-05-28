## Análisis de Diversidad con iNext
#Biblio de referencia: https://cran.r-project.org/web/packages/iNEXT/vignettes/Introduction.html
## install iNEXT package from CRAN
install.packages("iNEXT")
## import packages
library(iNEXT)
library(ggplot2)
library(tidyverse)
library(readODS)
library(dplyr)

# Cargar datos
#datos <- read_ods("nuevos_analisis_polinizadores.ods")
datos <- read_ods(file.choose())

# Ver nombres problemáticos
datos$`Nombre Cientifico`[duplicated(datos$`Nombre Cientifico`)]
# Ver filas con nombres problemáticos
datos %>% 
  filter(`Nombre Cientifico` %in% c("genero1", "genero2", "genero3", "-")) %>%
  select(Orden, Familia, Genero, `Nombre Cientifico`)
#uno familia y genero para evitar los duplicados
datos <- datos %>%
  mutate(`Nombre Cientifico` = paste(Familia, `Nombre Cientifico`, sep = "_"))
# Ver la fila problemática
datos %>% 
  filter(`Nombre Cientifico` == "-_-")
# Renombrar manualmente
datos$`Nombre Cientifico`[datos$Orden == "Diptera" & datos$`Nombre Cientifico` == "-_-"] <- "Diptera_sp"
datos$`Nombre Cientifico`[datos$Orden == "Coleoptera" & datos$`Nombre Cientifico` == "-_-"] <- "Coleoptera_sp"
# Verificar
sum(duplicated(datos$`Nombre Cientifico`))

# Crear matriz
pol <- datos %>%
  select(`Nombre Cientifico`, `ABRA 1`:`CAN 12`) %>%
  column_to_rownames("Nombre Cientifico")

colnames(pol) <- c("ABRA1","ABRA2","ABRA3",
                   "GOY4","GOY5","GOY6",
                   "EMB7","EMB8","EMB9",
                   "CAN10","CAN11","CAN12")

write.table(pol, "Polinizadores_Diversidad.txt",
            sep = "\t", quote = FALSE, col.names = NA)

head(pol)

# Verificar
head(pol)
dim(pol)

#Transponer para que las filas sean sitios y columnas especies
pol_t<- as.data.frame(t(pol))

# Asignar nombres correctos a columnas
colnames(pol_t) <- paste0("sp", seq_len(ncol(pol_t)))  # sp1, sp2, ..., spN
pol_t$site <- rownames(pol_t)
pol_t

#Asignar grupo según el nombre del sitio
pol_t$group <- case_when(
  grepl("^ABRA|^EMB", pol_t$site) ~ "P",
  grepl("^GOY|^CAN", pol_t$site) ~ "S",
  TRUE ~ NA_character_
)

# Separar por grupo
grupo_P <- pol_t %>% filter(group == "P") %>% select(starts_with("sp"))
grupo_S <- pol_t %>% filter(group == "S") %>% select(starts_with("sp"))

# Convertir a presencia-ausencia (0/1)
grupo_P_pa <- (grupo_P > 0) * 1
grupo_S_pa <- (grupo_S > 0) * 1

# Calcular frecuencia de ocurrencia por especie (número de sitios donde ocurre)
freq_P <- colSums(grupo_P_pa)
freq_S <- colSums(grupo_S_pa)

# Número de sitios
n_sites_P <- nrow(grupo_P_pa)
n_sites_S <- nrow(grupo_S_pa)

# Armar vectores como espera iNEXT
iNEXT_P <- c(n_sites_P, freq_P[freq_P > 0])
iNEXT_S <- c(n_sites_S, freq_S[freq_S > 0])

# Crear lista
data_iNEXT <- list(P = iNEXT_P, S = iNEXT_S)

# Correr iNEXT
out <- iNEXT(
  data_iNEXT,
  q = 0,
  datatype = "incidence_freq",
  endpoint = max(n_sites_P, n_sites_S) * 2
)

# Graficar
ggiNEXT(out, type = 1) +
  labs(title = "Curvas de rarefacción y extrapolación (riqueza)",
       x = "Número de unidades de muestreo (sitios)",
       y = "Riqueza de especies (q = 0)") +
  theme_minimal()


#############
# Separar por grupo
grupo_P <- pol_t %>% filter(group == "P") %>% select(starts_with("sp"))
grupo_S <- pol_t %>% filter(group == "S") %>% select(starts_with("sp"))

# Sumar abundancias por especie en cada grupo
abund_P <- colSums(grupo_P)
abund_S <- colSums(grupo_S)

# Crear lista
data_iNEXT_abund <- list(P = abund_P[abund_P > 0], S = abund_S[abund_S > 0])

# Correr iNEXT para abundancias
out_abund <- iNEXT(
  data_iNEXT_abund,
  q = 0,  # q = 0 es riqueza de especies
  datatype = "abundance",
  endpoint = 2 * max(sapply(data_iNEXT_abund, sum))  # extrapola al doble del total de individuos
)

# Graficar
ggiNEXT(out_abund, type = 1) +
  labs(title = "Curvas de rarefacción y extrapolación (abundancia)",
       x = "Número de individuos",
       y = "Riqueza de especies (q = 0)") +
  theme_minimal()
##########################
# completitud de muestreo ###
# Para abundancia
DataInfo(data_iNEXT_abund, datatype = "abundance")

# Para incidencia
DataInfo(data_iNEXT, datatype = "incidence_freq")

out_abund$iNextEst$coverage_based

out$iNextEst$coverage_based 
library(viridisLite)
library(viridis)
ggiNEXT(out_abund, type = 3) +  # type = 3 es cobertura
  labs(title = "Sample coverage vs. species richness",
       x = "Sample coverage",
       y = "Species richness") +
  scale_color_viridis_d(name = "Grupo", option = "D") +  # opción D es la clásica
  scale_fill_viridis_d(name = "Grupo", option = "D") +
  theme_minimal()


ggiNEXT(out, type = 3) +  # type = 3 es cobertura
  labs(title = "Sample coverage vs. species richness",
       x = "Sample coverage",
       y = "Species richness") +
  scale_color_viridis_d(name = "Grupo", option = "D") +  # opción D es la clásica
  scale_fill_viridis_d(name = "Grupo", option = "D") +
  theme_minimal()


