####################################################################
# ANALISIS DE COMUNIDADES DE VISITANTES FLORALES
# Fecha: 2026
# Descripcion: Analisis de comunidades de polinizadores en sistemas
# ganaderos usando GLMMs con distribucion binomial negativa
####################################################################
# PAQUETES 
library(readODS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(glmmTMB)
library(DHARMa)
library(vegan)

# BASE DE DATOS TRATAMIENTO
# IMPORTANTE: el orden de los sitios debe coincidir con las columnas
# de las planillas ODS (ABRA 1-3, GOY 4-6, EMB 7-9, CAN 10-12)

tratamiento <- data.frame(
  plots    = c("ABRA 1","ABRA 2","ABRA 3","GOY 4","GOY 5","GOY 6",
                "EMB 7","EMB 8","EMB 9","CAN 10","CAN 11","CAN 12"),
  treatment = factor(c("P","P","P","S","S","S","P","P","P","S","S","S")),
  sitio    = factor(c("ABRA","ABRA","ABRA","GOY","GOY","GOY",
                       "EMB","EMB","EMB","CAN","CAN","CAN")))
####################################################################
# 1. CARGAR Y PREPARAR DATOS
####################################################################

# CARGAR PLANILLAS 
# IMPORTANTE: colocar los archivos ODS en el directorio de trabajo
datos <- read_ods("nuevos_analisis_polinizadores.ods")
head(datos)
ls() # objetos cargados 
rasgos <- read_ods("rasgosfuncionales_apiformes.ods")

# MATRIZ DE ORDENES 
abund_orden <- datos %>%
group_by(Orden) %>%
summarise(across(`ABRA 1`:`CAN 12`, sum, .names = "{.col}"))

nombres_ordenes <- abund_orden$Orden
abund_matrix <- as.matrix(abund_orden[, -1])
rownames(abund_matrix) <- nombres_ordenes

# Transponer los sitios en filas, ordenes en columnas
abund_t <- t(abund_matrix)
head(abund_t)

#abundancia total de cada orden.  
abundancia_total_orden <- datos %>%
  group_by(Orden) %>%
  summarise(Total = sum(across(`ABRA 1`:`CAN 12`), na.rm = TRUE)) %>%
  arrange(desc(Total))

# Ver el resultado
print(abundancia_total_orden)

# MATRIZ DE GRUPOS FUNCIONALES ABEJAS 
hym <- datos %>%
  filter(Orden %in% c("Hymenoptera", "Otros Hymenoptera")) %>%
  filter(Familia != "Vespidae") %>%
  mutate(grupo = case_when(
    `Nombre Cientifico` == "Apis mellifera" ~ "Honeybee",
    `Nombre Cientifico` %in% c("Plebeia catamarcensis", "Plebeia molesta") ~ "Stingless bees",
    Familia == "Halictidae" ~ "Halictidae",
    Familia %in% c("Apidae","Andrenidae","Colletidae","Megachilidae") ~ "Other bees"
  ))

hym_grupos <- hym %>%
  group_by(grupo) %>%
  summarise(across(`ABRA 1`:`CAN 12`, sum)) %>%
  as.data.frame()

nombres_grupos <- hym_grupos$grupo
hym_grupos$grupo <- NULL
hym_matrix <- t(as.matrix(hym_grupos))
colnames(hym_matrix) <- nombres_grupos

# DIVERSIDAD Y ABUNDANCIA TOTAL 
tratamiento$shannon         <- diversity(abund_t, index = "shannon")
tratamiento$abundancia_total <- rowSums(abund_t)
tratamiento$riqueza_sp      <- specnumber(abund_t)

####################################################################
# 2. MODELOS (Shannon, Abundancia, Riqueza)
####################################################################
m_shannon <- glmmTMB(shannon ~ treatment + (1|sitio), data = tratamiento, family = gaussian)
summary(m_shannon)
residuos_shannon <- simulateResiduals(fittedModel = m_shannon)
plot(residuos_shannon)

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  0.76819    0.04902   15.67   <2e-16 ***
#  treatmentS   0.15667    0.06933    2.26   0.0238 * 
#Conway-Maxwell-Poisson
m_riqueza <- glmmTMB(riqueza_sp ~ treatment + (1|sitio), data = tratamiento, family = compois)
summary(m_riqueza) # este modelo tiene problemas de convergencia
#el intercepto tiene NaN en el error estándar. Esto ocurre porque el efecto aleatorio de bloque está absorbiendo toda la varianza del intercepto.
plot(simulateResiduals(m_riqueza))

m_riqueza_simple <- glmmTMB(riqueza_sp ~ treatment, data = tratamiento, family = compois)
summary(m_riqueza_simple)
plot(simulateResiduals(m_riqueza_simple))

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  1.57554    0.03147   50.07   <2e-16 ***
#  treatmentS  -0.07146    0.05520   -1.29    0.196  

m_abundancia <- glmmTMB(abundancia_total ~ treatment + (1|sitio), data = tratamiento, family =  nbinom2 )
summary(m_abundancia)
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   6.2288     0.1611   38.65  < 2e-16 ***
#  treatmentS   -1.3575     0.2300   -5.90 3.57e-09 ***

plot(simulateResiduals(m_abundancia))

####################################################################
# 4. MODELOS POR ORDEN (glmmTMB)
####################################################################
ord_df <- as.data.frame(abund_t)
ord_df$treatment <- tratamiento$treatment
ord_df$sitio    <- tratamiento$sitio

# Coleoptera
m_col      <- glmmTMB(Coleoptera ~ treatment + (1|sitio), data = ord_df, family = nbinom2)
summary(m_col)

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   5.7505     0.1864  30.856  < 2e-16 ***
#treatmentS   -1.6180     0.2676  -6.046 1.48e-09 ***

plot(simulateResiduals(m_col)) # los supuestos esta bien 
# El tratamiento es sig y la abundancia de coleoptera disminuye significativamente en silvopasturas
#Z=-6.04

# Hymenoptera
m_hym      <- glmmTMB(Hymenoptera ~ treatment + (1|sitio), data = ord_df, family = nbinom2)
summary(m_hym)

#AIC       BIC    logLik -2*log(L)  df.resid 
#133.8     135.7     -62.9     125.8         8 
#Random effects:
#Conditional model:
#Groups Name        Variance  Std.Dev. 
#bloque (Intercept) 6.154e-10 2.481e-05
#Number of obs: 12, groups:  bloque, 4

#Dispersion parameter for nbinom2 family (): 4.39
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   5.1591     0.1973  26.154  < 2e-16 ***
#  treatmentS   -1.0759     0.2823  -3.812 0.000138 ***

plot(simulateResiduals(m_hym)) # 

# La misma tendencia, los hymenopteros disminuyen sig en las silvopasturas Z=-3.81

# Diptera
m_dip <- glmmTMB(Diptera ~ treatment + (1|sitio), data = ord_df, family = nbinom2)
summary(m_dip)
#Conditional model:
#Estimate Std. Error z value Pr(>|z|)    3
#(Intercept)   1.5404     0.2076   7.418 1.18e-13 ***
#treatmentS   -0.3878     0.3212  -1.207    0.227  

plot(simulateResiduals(m_dip)) # estan ok
# el tipo de manejo no es significativo 

# Lepidoptera
m_lep      <- glmmTMB(Lepidoptera ~ treatment + (1|sitio), data = ord_df, family = nbinom2)
summary(m_lep)
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   1.8971     0.3218   5.896 3.73e-09 ***
#  treatmentS   -0.5534     0.4749  -1.165    0.244  
plot(simulateResiduals(m_lep))
# el tipo de manejo no es significativo

####################################################################
# 5. MODELOS GRUPOS FUNCIONALES ABEJAS (glmmTMB)
####################################################################
hym_df<- as.data.frame(hym_matrix)
hym_df<- as.data.frame(hym_matrix)
hym_df$treatment <- tratamiento$treatment
hym_df$sitio <- tratamiento$sitio
# Verificar
names(hym_df)
# Total de individuos por grupo funcional
hym %>%
  group_by(grupo) %>%
  summarise(total = sum(across(`ABRA 1`:`CAN 12`))) %>%
  arrange(desc(total))

#STINGLESS BEES
m_stingless <- glmmTMB(`Stingless bees` ~ treatment + (1|sitio), data = hym_df, family = nbinom2)
summary(m_stingless)

#AIC       BIC    logLik -2*log(L)  df.resid 
#117.9     119.8     -54.9     109.9         8 
#Random effects:
#  Conditional model:
#  Groups Name        Variance Std.Dev.
#sitio (Intercept) 0.3354   0.5791  
#Number of obs: 12, groups:  bloque, 4
#Dispersion parameter for nbinom2 family (): 0.965 
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   4.8686     0.5960   8.169  3.1e-16 ***
#  treatmentS   -2.8553     0.9903  -2.883  0.00394 ** 
plot(simulateResiduals(m_stingless)) #ok
# el tipo de manejo afecta sig a la abundancia de stingless bees

# OTHER BEES
m_other <- glmmTMB(`Other bees` ~ treatment + (1|sitio), data = hym_df, family = nbinom2)
summary(m_other)
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   2.8409     0.4100   6.930 4.22e-12 ***
#  treatmentS   -0.1828     0.5848  -0.313    0.755 

# el tipo de manejo no afecta sig a la abundancia de other bees
plot(simulateResiduals(m_other))

# HALICTIDAE 
m_halictidae <- glmmTMB(Halictidae ~ treatment + (1|sitio), data = hym_df, family = nbinom2)
summary(m_halictidae)
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  3.44843    0.27358  12.605   <2e-16 ***
#  treatmentS  -0.09114    0.38765  -0.235    0.814  
plot(simulateResiduals(m_halictidae)) #ok
# el tipo de uso del suelo no afecta sig a halictidae

# HONEYBEE 
m_honeybee <- glmmTMB(Honeybee ~ treatment + (1|sitio), data = hym_df, family = nbinom2)
summary(m_honeybee) 
#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)
#(Intercept)  -0.4055     0.5592  -0.725    0.468
#treatmentS   -0.6931     0.9356  -0.741    0.459
plot(simulateResiduals(m_honeybee))
# apis no se ve sig afectada por el tipo de manejo
testZeroInflation(simulateResiduals(m_honeybee)) # no esta inflado en 0
####################################################################
# 6. MODELOS RASGOS FUNCIONALES (glmmTMB).
####################################################################
# OBJETO hym_long_traits
hym_long_traits <- hym %>%
  pivot_longer(cols = `ABRA 1`:`CAN 12`,
               names_to = "sitio",
               values_to = "abundancia") %>%
  mutate(
    treatment = case_when(
      sitio %in% c("ABRA 1","ABRA 2","ABRA 3") ~ "P",
      sitio %in% c("GOY 4","GOY 5","GOY 6")    ~ "S",
      sitio %in% c("EMB 7","EMB 8","EMB 9")    ~ "P",
      sitio %in% c("CAN 10","CAN 11","CAN 12") ~ "S"
    ),
   sitio = case_when(
      sitio %in% c("ABRA 1","ABRA 2","ABRA 3") ~ "ABRA",
      sitio %in% c("GOY 4","GOY 5","GOY 6")    ~ "GOY",
      sitio %in% c("EMB 7","EMB 8","EMB 9")    ~ "EMB",
      sitio %in% c("CAN 10","CAN 11","CAN 12") ~ "CAN"
    )
  ) %>%
  left_join(rasgos %>% select(Specie, Behavior, `Nest location`, Diet),
            by = c("Nombre Cientifico" = "Specie"))

# Verificar
head(hym_long_traits)

# ?????? BEHAVIOR ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
beh_data <- hym_long_traits %>%
  filter(!is.na(Behavior), Behavior != "cleptoparasitic") %>%
  group_by(sitio, treatment, sitio, Behavior) %>%
  summarise(abun = sum(abundancia), .groups = "drop") %>%
  group_by(sitio) %>%
  mutate(total = sum(abun)) %>%
  ungroup()

# ?????? NEST LOCATION ???????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
nest_data <- hym_long_traits %>%
  filter(!is.na(`Nest location`)) %>%
  group_by(sitio, treatment, sitio, `Nest location`) %>%
  summarise(abun = sum(abundancia), .groups = "drop") %>%
  group_by(sitio) %>%
  mutate(total = sum(abun)) %>%
  ungroup()

# ?????? DIET ??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????
diet_data <- hym_long_traits %>%
  filter(!is.na(Diet), Diet != "NA") %>%
  group_by(sitio, treatment, sitio, Diet) %>%
  summarise(abun = sum(abundancia), .groups = "drop") %>%
  group_by(sitio) %>%
  mutate(total = sum(abun)) %>%
  ungroup()

# Especies e individuos por rasgo de sociabilidad
hym_long_traits %>%
  filter(!is.na(Behavior), Behavior != "cleptoparasitic") %>%
  group_by(Behavior) %>%
  summarise(
    n_especies   = n_distinct(`Nombre Cientifico`),
    n_individuos = sum(abundancia)
  )

#Behavior n_especies n_individuos
#<chr>         <int>        <dbl>
#  1 social            6          824
#2 solitary         12          133

# Especies e individuos por nest location
hym_long_traits %>%
  filter(!is.na(`Nest location`)) %>%
  group_by(`Nest location`) %>%
  summarise(
    n_especies   = n_distinct(`Nombre Cientifico`),
    n_individuos = sum(abundancia)
  )
#Nest location` n_especies n_individuos
#<chr>                <int>        <dbl>
#  1 above ground             6          822
#2 below ground            13          137
# Especies e individuos por diet
hym_long_traits %>%
  filter(!is.na(Diet), Diet != "NA") %>%
  group_by(Diet) %>%
  summarise(
    n_especies   = n_distinct(`Nombre Cientifico`),
    n_individuos = sum(abundancia)
  )
#Diet        n_especies n_individuos
#<chr>            <int>        <dbl>
#  1 oligolectic          3           28
#2 polylectic          10          917
##############

# BEHAVIOR/SOCIALITY 
m_beh <- glmmTMB(abun ~ treatment * Behavior + (1|sitio), data = beh_data,family = nbinom2)
summary(m_beh)
plot(simulateResiduals(m_beh))

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                   4.8416     0.3290  14.716  < 2e-16 ***
#  treatmentS                   -2.4744     0.4804  -5.151 2.60e-07 ***
#  Behaviorsolitary             -2.5224     0.4812  -5.242 1.59e-07 ***
#  treatmentS:Behaviorsolitary   2.6402     0.6892   3.831 0.000128 ***


# NEST LOCATION
m_nest <- glmmTMB(abun ~ treatment * `Nest location` + (1|stio),
                  data = nest_data,
                  family = nbinom2)
summary(m_nest)
plot(simulateResiduals(m_nest))

#AIC       BIC    logLik -2*log(L)  df.resid 
#205.0     212.1     -96.5     193.0        18 

#Random effects:
  
#  Conditional model:
#  Groups Name        Variance  Std.Dev. 
#sitio (Intercept) 6.023e-10 2.454e-05
#Number of obs: 24, groups:  bloque, 4

#Dispersion parameter for nbinom2 family (): 1.54 

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                              4.8363     0.3315  14.589  < 2e-16 ***
#  treatmentS                              -2.4384     0.4833  -5.045 4.54e-07 ***
#  `Nest location`below ground             -2.4537     0.4836  -5.074 3.90e-07 ***
#  treatmentS:`Nest location`below ground   2.5407     0.6928   3.667 0.000245 ***

# DIET 
m_diet <- glmmTMB(abun ~ treatment * Diet + (1|sitio),
                  data = diet_data,
                  family = nbinom2)
summary(m_diet)
plot(simulateResiduals(m_diet))

#AIC       BIC    logLik -2*log(L)  df.resid 
#178       185       -83       166        18 

#Random effects:
#  
#  Conditional model:
#  Groups Name        Variance  Std.Dev. 
#sitio (Intercept) 8.434e-10 2.904e-05
#Number of obs: 24, groups:  bloque, 4

#Dispersion parameter for nbinom2 family (): 1.55 

#Conditional model:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)                 1.0415     0.4081   2.552   0.0107 *  
#  treatmentS                 -0.4353     0.6042  -0.720   0.4713    
#Dietpolylectic              3.8501     0.5248   7.336  2.2e-13 ***
#  treatmentS:Dietpolylectic  -1.4774     0.7683  -1.923   0.0545 .

####################################################################
# DF RESIDUAL DE CADA MODELO
####################################################################
df.residual(m_abundancia)
df.residual(m_shannon)
df.residual(m_riqueza)
df.residual(m_col)
df.residual(m_hym)
df.residual(m_dip)
df.residual(m_lep)
df.residual(m_other)
df.residual(m_nest)
df.residual(m_beh)
df.residual(m_diet)

####################################################################
# GRAFICO PARA RELACIONES ENTRE TRAITS
####################################################################
library(ggalluvial)

# Preparar datos
traits_data <- rasgos %>%
  filter(Behavior != "cleptoparasitic",
         !is.na(`Nest location`),
         !is.na(Diet),
         Diet != "NA") %>%
  group_by(Behavior, `Nest location`, Diet) %>%
  summarise(n = n(), .groups = "drop")

# Gráfico alluvial
ggplot(traits_data,
       aes(axis1 = Behavior,
           axis2 = `Nest location`,
           axis3 = Diet,
           y = n)) +
  geom_alluvium(aes(fill = Behavior), alpha = 0.7) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("social" = "#6A0DAD", "solitary" = "#FFD700")) +
  scale_x_discrete(limits = c("Sociality", "Nest Location", "Diet"),
                   expand = c(0.1, 0.1)) +
  labs(title = "Trait Associations in the Floral Visitor Community",
       y = "Number of Species",
       fill = "Sociality") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom")
