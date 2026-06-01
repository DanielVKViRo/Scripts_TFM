library(MASS)
library(deSolve)
library(lpSolve)
library(tidyverse)
library(writexl)
library(readxl)

########################## Ponemos el path donde están los documentos  
base_dir <- "C:/Users/danie/Documents/Master_Bioinformatica_y_Biologia_Computacional/TFM/scripts de Abel para Daniel/scripts de Abel para Daniel/Comunidades_20esp_por_com"
setwd(base_dir)
##########################

########################## Tipos de combinaciones
all_groups   <- c("SD","SC","GD","GC")
get_combinations <- function(vec) {
  unlist(lapply(seq_along(vec),
                function(k) combn(vec, k, simplify = FALSE)),
         recursive = FALSE)
}
all_combs  <- get_combinations(all_groups)
comb_names <- sapply(all_combs, paste, collapse = "+")
##########################

###### Inicializamos donde guardar los resultados de las invasiones para la primera simulacion
results_invasions <- data.frame(
  type_community = character(),
  initial_richness = integer(),
  l = numeric(),
  kappa = integer(),
  case_letter = character(),
  resistance = numeric(),
  resistance_SD = integer(),
  resistance_GD = integer(),
  resistance_GC = integer(),
  resistance_SC = integer(),
  displacement = numeric(),
  displacement_SD = integer(),
  displacement_GD = integer(),
  displacement_GC = integer(),
  displacement_SC = integer(),
  disruption = numeric(),
  disruption_SD = integer(),
  disruption_GD = integer(),
  disruption_GC = integer(),
  disruption_SC = integer(),
  augmentation = numeric(),
  augmentation_SD = integer(),
  augmentation_GD = integer(),
  augmentation_GC = integer(),
  augmentation_SC = integer()
)

# Definimos los distintos tamaños de recursos
sizes_res <- c("S", "M", "L")

# Cargamos las especies
Pool <- readRDS("Pool_200.rds")

## Iteramos por cada tipo de comunidad  
for (current_comb in comb_names) {
  
  ## Iteramos por cada tamaño de comunidad
  for (size_r in sizes_res) {
    
    ## Cargamos los archivos
    file_results <- paste0("Results_matrices_", current_comb, "_20rep_", size_r, "_invasions.xlsx")
    assembly <- read_excel(file_results)
    
    ## Iteramos por cada comunidad 
    for (com_num in seq(1, 151200, 840)) {
      current_community <- assembly[com_num:(com_num + 839), ]
      
      ## Inicializamos los valores de los tipos de invasión
      augmentation <- 0
      augmentation_SD <- 0
      augmentation_GD <- 0
      augmentation_GC <- 0
      augmentation_SC <- 0
      
      resistance <- 0
      resistance_SD <- 0
      resistance_GD <- 0
      resistance_GC <- 0
      resistance_SC <- 0
      
      displacement <- 0
      displacement_SD <- 0
      displacement_GD <- 0
      displacement_GC <- 0
      displacement_SC <- 0
      
      disruption <- 0
      disruption_SD <- 0
      disruption_GD <- 0
      disruption_GC <- 0
      disruption_SC <- 0
      
      ## Analizamos las 40 invasiones a la comunidad
      
      for (inv_num in seq(1, 840, 21)) {
        
        # Tomamos la invasión inicial 
        current_invasion <- current_community[inv_num:(inv_num + 20), ]
        
        # Tomamos el tipo del invasor
        invasor_type <- substr(current_invasion[21,]$species_id, 1, 2)
        
        # Calculamos la riqueza antes y después de la invasión
        abundances_init <- as.numeric(current_invasion$abundance)
        abundances <- current_invasion$abundance_invasion
        
        species_init <- sum(abundances_init[1:20] > 0)
        species_invas <- sum(abundances[1:20] > 0)
        
        # Vemos cuál es el tipo de invasión
        if (abundances[21] > 0) {
          if (abs(species_init - species_invas)>0){
            displacement <- displacement + 1
            if (invasor_type == "SD") {displacement_SD <- displacement_SD + 1}
            else if (invasor_type == "GD") {displacement_GD <- displacement_GD + 1}
            else if (invasor_type == "GC") {displacement_GC <- displacement_GC + 1}
            else {displacement_SC <- displacement_SC + 1}
          } else {
            augmentation <- augmentation + 1
            if (invasor_type == "SD") {augmentation_SD <- augmentation_SD + 1}
            else if (invasor_type == "GD") {augmentation_GD <- augmentation_GD + 1}
            else if (invasor_type == "GC") {augmentation_GC <- augmentation_GC + 1}
            else {augmentation_SC <- augmentation_SC + 1}
          }
        } else {
          if (abs(species_init - species_invas)>0){
            disruption <- disruption+ 1
            if (invasor_type == "SD") {disruption_SD <- disruption_SD + 1}
            else if (invasor_type == "GD") {disruption_GD <- disruption_GD + 1}
            else if (invasor_type == "GC") {disruption_GC <- disruption_GC + 1}
            else {disruption_SC <- disruption_SC + 1}
          } else {
            resistance <- resistance + 1
            if (invasor_type == "SD") {resistance_SD <- resistance_SD + 1}
            else if (invasor_type == "GD") {resistance_GD <- resistance_GD + 1}
            else if (invasor_type == "GC") {resistance_GC <- resistance_GC + 1}
            else {resistance_SC <- resistance_SC + 1}
          }
        }
        
      }
      
      # Porcentaje sobre total de invasiones
      resistance <- resistance / 40 * 100
      resistance_SD <- resistance_SD / 40 * 100
      resistance_GD <- resistance_GD / 40 * 100
      resistance_GC <- resistance_GC / 40 * 100
      resistance_SC <- resistance_SC / 40 * 100
      
      disruption <- disruption / 40 * 100
      disruption_SD <- disruption_SD / 40 * 100
      disruption_GD <- disruption_GD / 40 * 100
      disruption_GC <- disruption_GC / 40 * 100
      disruption_SC <- disruption_SC / 40 * 100
      
      augmentation <- augmentation / 40 * 100
      augmentation_SD <- augmentation_SD / 40 * 100
      augmentation_GD <- augmentation_GD / 40 * 100
      augmentation_GC <- augmentation_GC / 40 * 100
      augmentation_SC <- augmentation_SC / 40 * 100
      
      displacement <- displacement / 40 * 100
      displacement_SD <- displacement_SD / 40 * 100
      displacement_GD <- displacement_GD / 40 * 100
      displacement_GC <- displacement_GC / 40 * 100
      displacement_SC <- displacement_SC / 40 * 100
      
      # Resto de parámetros
      l <- current_invasion$l[1]
      kappa <- current_invasion$kappa[1]
      case_letter <- current_invasion$case_letter[1]
      
      results_analysis <- data.frame(
        type_community = current_comb,
        initial_richness = species_init,
        l = l,
        kappa = kappa,
        case_letter = case_letter,
        resistance = resistance,
        resistance_SD = resistance_SD,
        resistance_GD = resistance_GD,
        resistance_GC = resistance_GC,
        resistance_SC = resistance_SC,
        displacement = displacement,
        displacement_SD = displacement_SD,
        displacement_GD = displacement_GD,
        displacement_GC = displacement_GC,
        displacement_SC = displacement_SC,
        disruption = disruption,
        disruption_SD = disruption_SD,
        disruption_GD = disruption_GD,
        disruption_GC = disruption_GC,
        disruption_SC = disruption_SC,
        augmentation = augmentation,
        augmentation_SD = augmentation_SD,
        augmentation_GD = augmentation_GD,
        augmentation_GC = augmentation_GC,
        augmentation_SC = augmentation_SC
      )
      
      results_invasions <- rbind(
        results_invasions,
        results_analysis
      )
      
      
    }
  }
}

# Quitamos las comunidades con riqueza inicial nula
results_invasions <- results_invasions[results_invasions$initial_richness > 0,]
# Convertir kappa a valor numerico
results_invasions$kappa <- as.numeric(results_invasions$kappa)

####################################################################################################
# Graficar los resultados
####################################################################################################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(rstatix)


#------------------------------------------------------------------------------


# =========================
# Datos
# =========================

df1 <- results_invasions[,c("case_letter", "kappa", "augmentation")]
df2 <- results_invasions[,c("case_letter", "kappa", "displacement")]
df3 <- results_invasions[,c("case_letter", "kappa", "disruption")]
df4 <- results_invasions[,c("case_letter", "kappa", "resistance")]

df1 <- df1 %>% rename(invasion = augmentation)
df2 <- df2 %>% rename(invasion = displacement)
df3 <- df3 %>% rename(invasion = disruption)
df4 <- df4 %>% rename(invasion = resistance)

df1$grupo <- "Aumento"
df2$grupo <- "Desplazamiento"
df3$grupo <- "Disrupción"
df4$grupo <- "Resistencia"

df_all <- bind_rows(df1, df2, df3, df4)

df_all$kappa <- recode(df_all$kappa,
                       "10" = "κ bajo",
                       "60" = "κ alto",
                       "100" = "κ alto")

df_all$grupo <- factor(df_all$grupo,
                       levels = c("Resistencia",
                                  "Desplazamiento",
                                  "Aumento",
                                  "Disrupción"))

df_all$kappa <- factor(df_all$kappa,
                       levels = c("κ bajo", "κ alto"))

# =========================
# Resumen
# =========================

df_sum <- df_all %>%
  group_by(case_letter, grupo, kappa) %>%
  summarise(
    media = mean(invasion, na.rm = TRUE),
    sd = sd(invasion, na.rm = TRUE),
    .groups = "drop"
  )

# =========================
# Tests (κ bajo vs κ alto)
# =========================

comparaciones <- df_all %>%
  group_by(case_letter, grupo) %>%
  summarise(
    p = t.test(invasion ~ kappa)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

y_pos <- df_sum %>%
  group_by(case_letter, grupo) %>%
  summarise(
    y = max(media + sd) + 5,
    .groups = "drop"
  )

comparaciones <- left_join(comparaciones, y_pos,
                           by = c("case_letter", "grupo"))

comparaciones$label <- case_when(
  comparaciones$p_adj < 0.0001 ~ "****",
  comparaciones$p_adj < 0.001 ~ "***",
  comparaciones$p_adj < 0.01  ~ "**",
  comparaciones$p_adj < 0.05  ~ "*"
)

comparaciones$group1 <- "κ bajo"
comparaciones$group2 <- "κ alto"

# =========================
# Plot
# =========================

pd <- position_dodge(width = 0.8)

ggplot(df_sum, aes(x = grupo, y = media, fill = kappa)) +
  
  geom_col(position = pd, width = 0.7) +
  
  geom_errorbar(
    aes(ymin = media - sd, ymax = media + sd),
    position = pd,
    width = 0.2,
    linewidth = 0.8
  ) +
  
  facet_wrap(~ case_letter, nrow = 1) +
  
  stat_pvalue_manual(
    comparaciones,
    label = "label",
    xmin = "grupo",
    xmax = "grupo",
    y.position = "y",
    tip.length = 0.01,
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "Tipo de invasión",
    y = "Porcentaje de resultado de invasión",
    fill = expression(kappa)
  ) +
  
  theme_minimal(base_size = 16) +
  
  theme(
    
    axis.title = element_text(size = 16),
    
    axis.text.x = element_text(
      angle = 40,
      hjust = 1,
      vjust = 1,
      size = 13
    ),
    
    axis.text.y = element_text(size = 13),
    
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )



# Gráfica de resistencia según la especie invasora

# ----------------------------
# Datos en formato largo
# ----------------------------
df <- results_invasions[, c("kappa", "case_letter",
                            "resistance_GD", "resistance_GC",
                            "resistance_SD", "resistance_SC")]

df_long <- df %>%
  pivot_longer(cols = c(resistance_GD, resistance_GC, resistance_SD, resistance_SC),
               names_to = "tipo_invasion",
               values_to = "porcentaje") %>%
  mutate(
    tipo_invasion = recode(tipo_invasion,
                           "resistance_GD" = "GD",
                           "resistance_GC" = "GC",
                           "resistance_SD" = "SD",
                           "resistance_SC" = "SC"),
    kappa = recode(kappa,
                   "10" = "κ bajo",
                   "60" = "κ alto",
                   "100" = "κ alto")
  )

# ----------------------------
# Estadística (comparaciones por pares)
# ----------------------------
stat.test <- df_long %>%
  group_by(kappa, case_letter) %>%
  pairwise_t_test(
    porcentaje ~ tipo_invasion,
    p.adjust.method = "BH"
  ) %>%
  arrange(p.adj) %>%
  add_xy_position(x = "tipo_invasion") %>%
  group_by(kappa, case_letter) %>%
  mutate(
    y.position = max(df_long$porcentaje, na.rm = TRUE) +
      row_number() * 5
  ) %>%
  ungroup() %>%
  filter(p.adj.signif != "ns")

# ----------------------------
# Resumen para barras
# ----------------------------
df_sum <- df_long %>%
  group_by(kappa, case_letter, tipo_invasion) %>%
  summarise(
    media = mean(porcentaje, na.rm = TRUE),
    sd = sd(porcentaje, na.rm = TRUE),
    .groups = "drop"
  )

pd <- position_dodge(0.7)

# ----------------------------
# Gráfico final
# ----------------------------
p <- ggplot(df_sum, aes(x = tipo_invasion, y = media, fill = tipo_invasion)) +
  
  geom_bar(stat = "identity", position = pd, width = 0.65, color = "black") +
  
  geom_errorbar(aes(ymin = media - sd, ymax = media + sd),
                width = 0.2, position = pd) +
  
  facet_grid(kappa ~ case_letter) +
  
  labs(
    x = "Tipo de especie",
    y = "Porcentaje de resultado de invasión"
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

# ----------------------------
# Significancia estadística
# ----------------------------
p + stat_pvalue_manual(
  stat.test,
  label = "p.adj.signif",
  tip.length = 0.01
)




# Gráfica según especie invasora y tipo de comunidad

# -----------------------------
# 1. Selección de datos
# -----------------------------
df <- results_invasions[, c(
  "type_community",
  "kappa",
  "case_letter",
  "resistance_SD",
  "resistance_SC",
  "resistance_GD",
  "resistance_GC"
)]

# -----------------------------
# 2. Agrupar comunidades
# -----------------------------
df$type_community <- ifelse(
  df$type_community %in% c("SC", "SD", "SD+SC"),
  "solo especialistas",
  "con generalistas"
)

# -----------------------------
# 3. Formato largo
# -----------------------------
df_long <- df %>%
  pivot_longer(
    cols = c(resistance_GD, resistance_GC, resistance_SD, resistance_SC),
    names_to = "tipo_invasion",
    values_to = "porcentaje"
  )

# -----------------------------
# 4. Recodificar invasores
# -----------------------------
df_long$tipo_invasion <- recode(
  df_long$tipo_invasion,
  "resistance_GD" = "generalista",
  "resistance_GC" = "generalista",
  "resistance_SD" = "especialista",
  "resistance_SC" = "especialista"
)

df_long$tipo_invasion <- factor(df_long$tipo_invasion)

# -----------------------------
# 5. Recodificar kappa
# -----------------------------
df_long$kappa <- recode(
  as.character(df_long$kappa),
  "10" = "κ bajo",
  "60" = "κ alto",
  "100" = "κ alto"
)

# -----------------------------
# 6. Resumen para barras
# -----------------------------
df_sum <- df_long %>%
  group_by(type_community, kappa, case_letter, tipo_invasion) %>%
  summarise(
    media_resistencia = mean(porcentaje, na.rm = TRUE),
    sd_resistencia = sd(porcentaje, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# 7. Test estadístico
# -----------------------------
stat_test <- df_long %>%
  group_by(type_community, kappa, case_letter) %>%
  t_test(porcentaje ~ tipo_invasion) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  add_xy_position(x = "type_community", dodge = 0.7)

# -----------------------------
# 8. Gráfico con barras
# -----------------------------
p <- ggplot(
  df_sum,
  aes(
    x = type_community,
    y = media_resistencia,
    fill = tipo_invasion
  )
) +
  
  geom_col(
    position = position_dodge(0.7),
    width = 0.6
  ) +
  
  geom_errorbar(
    aes(
      ymin = media_resistencia - sd_resistencia,
      ymax = media_resistencia + sd_resistencia
    ),
    position = position_dodge(0.7),
    width = 0.2
  ) +
  
  facet_grid(kappa ~ case_letter) +
  
  theme_minimal() +
  
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  
  labs(
    x = "Tipo de comunidad",
    y = "Porcentaje de resistencia",
    fill = "Tipo de invasor"
  )

# -----------------------------
# 9. Añadir significancia
# -----------------------------
p +
  stat_pvalue_manual(
    stat_test,
    label = "p.adj.signif",
    hide.ns = TRUE
  )

# Gráfica de tipo de invasión según riqueza inicial
df1 <- results_invasions[,c("initial_richness", "augmentation")]
df2 <- results_invasions[,c("initial_richness", "displacement")]
df3 <- results_invasions[,c("initial_richness", "disruption")]
df4 <- results_invasions[,c("initial_richness", "resistance")]

df1 <- df1 %>% rename(invasion = augmentation)
df2 <- df2 %>% rename(invasion = displacement)
df3 <- df3 %>% rename(invasion = disruption)
df4 <- df4 %>% rename(invasion = resistance)

df1$group <- "Aumento"
df2$group <- "Desplazamiento"
df3$group <- "Disrupción"
df4$group <- "Resistencia"

df_all <- bind_rows(df1, df2, df3, df4)

df_sum <- df_all %>%
  group_by(initial_richness, group) %>%
  summarise(
    mean_inv = mean(invasion, na.rm = TRUE),
    sd_inv = sd(invasion, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_sum, aes(x = initial_richness, y = mean_inv,
                   color = group, fill = group, group = group)) +
  geom_ribbon(aes(ymin = mean_inv - sd_inv,
                  ymax = mean_inv + sd_inv),
              alpha = 0.25, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  labs(
    x = "Riqueza inicial de la comunidad",
    y = "Porcentaje de resultados de invasión",
    color = "Tipo de invasión",
    fill = "Tipo de invasión"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title = element_text(face = "bold")
  )