##############################################
# 1) GENERACIÓN DEL POOL DE ESPECIES (200)
##############################################

# Limpieza y directorio de trabajo
rm(list = ls())
setwd("C:/Users/danie/Documents/Scripts_TFM")

# Parámetros
n_sp         <- 200
n_res        <- 20    # número de recursos
alpha_high   <- 20    # concentración alta Dirichlet
alpha_low    <- 0.1   # concentración baja Dirichlet
n_per_group  <- 50    # especies por cada grupo ED, EC, GD, GC
thresh       <- 0.01  # umbral mínimo para cortar pesos

reps_per_comb <- 20

# Librerías
library(MCMCpack)
library(tidyverse)

# Función para umbral y renormalización
threshold_vector <- function(v, thresh = 0.01) {
  v[v < thresh] <- 0
  if (sum(v) > 0) v <- v / sum(v)
  v
}

# Creadores controlados de cada tipo
create_SD_controlled <- function(id, thresh = 0.01) {
  cons_res <- sample(1:n_res, 1)
  sec_res  <- sample(setdiff(1:n_res, cons_res), 1)
  a1 <- rep(alpha_low, n_res); a1[cons_res] <- alpha_high
  a2 <- rep(alpha_low, n_res); a2[sec_res]  <- alpha_high
  list(
    id          = id,
    group       = "SD",
    label       = paste0("S", cons_res, "_D", sec_res),
    consumption = threshold_vector(rdirichlet(1, a1), thresh),
    secretion   = threshold_vector(rdirichlet(1, a2), thresh)
  )
}

create_GD_controlled <- function(id, thresh = 0.01) {
  cons_idx <- sort(sample(1:n_res, 10))
  sec_res  <- sample(setdiff(1:n_res, cons_idx), 1)
  a1 <- rep(alpha_low, n_res); a1[cons_idx] <- alpha_high
  a2 <- rep(alpha_low, n_res); a2[sec_res]  <- alpha_high
  list(
    id          = id,
    group       = "GD",
    label       = paste0("G", paste(cons_idx, collapse = "."), "_D", sec_res),
    consumption = threshold_vector(rdirichlet(1, a1), thresh),
    secretion   = threshold_vector(rdirichlet(1, a2), thresh)
  )
}

create_SC_controlled <- function(id, thresh = 0.01) {
  cons_res <- sample(1:n_res, 1)
  sec_idx  <- sort(sample(setdiff(1:n_res, cons_res), 10))
  a1 <- rep(alpha_low, n_res); a1[cons_res] <- alpha_high
  a2 <- rep(alpha_low, n_res); a2[sec_idx]  <- alpha_high
  list(
    id          = id,
    group       = "SC",
    label       = paste0("S", cons_res, "_D", paste(sec_idx, collapse = ".")),
    consumption = threshold_vector(rdirichlet(1, a1), thresh),
    secretion   = threshold_vector(rdirichlet(1, a2), thresh)
  )
}

create_GC_controlled <- function(id, thresh = 0.01) {
  cons_idx <- sort(sample(1:n_res, 10))
  sec_idx  <- setdiff(1:n_res, cons_idx)
  a1 <- rep(alpha_low, n_res); a1[cons_idx] <- alpha_high
  sec_partial <- rdirichlet(1, rep(alpha_high, length(sec_idx)))
  sec_vec <- numeric(n_res); sec_vec[sec_idx] <- sec_partial
  list(
    id          = id,
    group       = "GC",
    label       = paste0("G", paste(cons_idx, collapse = "."), "_C"),
    consumption = threshold_vector(rdirichlet(1, a1), thresh),
    secretion   = threshold_vector(sec_vec, thresh)
  )
}

# Generar listas de especies
ED_list <- map(seq_len(n_per_group), ~ create_SD_controlled(paste0("SD_", .x), thresh))
GD_list <- map(seq_len(n_per_group), ~ create_GD_controlled(paste0("GD_", .x), thresh))
EC_list <- map(seq_len(n_per_group), ~ create_SC_controlled(paste0("SC_", .x), thresh))
GC_list <- map(seq_len(n_per_group), ~ create_GC_controlled(paste0("GC_", .x), thresh))
pool_list <- c(ED_list, GD_list, EC_list, GC_list)

# Construir data.frame base
pool_df <- tibble(
  id    = map_chr(pool_list, "id"),
  group = map_chr(pool_list, "group"),
  label = map_chr(pool_list, "label")
)

# Construir matrices de consumo y secreción sin advertencia
consumption_matrix <- do.call(rbind, lapply(pool_list, `[[`, "consumption"))
colnames(consumption_matrix) <- paste0("C", 1:n_res)
cons_mat <- as_tibble(consumption_matrix)

secretion_matrix <- do.call(rbind, lapply(pool_list, `[[`, "secretion"))
colnames(secretion_matrix) <- paste0("S", 1:n_res)
sec_mat <- as_tibble(secretion_matrix)

# Unir todo
pool_df_final <- bind_cols(pool_df, cons_mat, sec_mat)

# Nombre dinámico para el pool
pool_filename <- paste0("Pool_", n_per_group * 4, ".rds")

# Guardar pool
saveRDS(pool_df_final, pool_filename)
cat("Pool generado con", n_sp, "especies y guardado en", pool_filename, "\n\n")



##############################################
# 2) SIMULACIÓN DE COMUNIDADES
##############################################

pool_filename <- "Pool_200.rds"

# Tomamos 20 especies por comunidad
reps_per_comb <- 20

# Parámetros de simulación
require_full_consumption <- TRUE
max_attempts             <- 1000

# Cargar pool
pool_df_final <- readRDS(pool_filename)

# Función para todas las combinaciones no vacías
all_groups   <- c("SD","SC","GD","GC")
get_combinations <- function(vec) {
  unlist(lapply(seq_along(vec),
                function(k) combn(vec, k, simplify = FALSE)),
         recursive = FALSE)
}
all_combs  <- get_combinations(all_groups)
comb_names <- sapply(all_combs, paste, collapse = "+")

# Patrones rotativos para 3 grupos
species_patterns <- list(
  c(7,7,6), c(6,7,7), c(7,6,7)
)

# Bucle de simulación
for(i_comb in seq_along(all_combs)) {
  current_comb <- all_combs[[i_comb]]
  comb_name    <- comb_names[i_comb]
  cat("Procesando combinación:", comb_name, "\n")
  
  num_groups  <- length(current_comb)
  subset_pool <- filter(pool_df_final, group %in% current_comb)
  
  matrices_C <- vector("list", reps_per_comb)
  matrices_S <- vector("list", reps_per_comb)
  
  for(rep_i in seq_len(reps_per_comb)) {
    # Determinar número de especies por grupo
    species_per_group <- switch(
      as.character(num_groups),
      "1" = rep(20,1),
      "2" = rep(10,2),
      "3" = species_patterns[[ (rep_i - 1) %% 3 + 1 ]],
      "4" = rep(5,4)
    )
    
    attempt <- 1
    repeat {
      # Muestreo sin reemplazo y única ocurrencia de cada ID
      idxs <- unlist(mapply(function(grp, n_sp) {
        sample(filter(subset_pool, group == grp)$id, n_sp, replace = FALSE)
      }, current_comb, species_per_group, SIMPLIFY = FALSE))
      if (length(unique(idxs)) != length(idxs)) {
        stop("IDs repetidos en réplica ", rep_i, " de ", comb_name)
      }
      
      sel   <- filter(subset_pool, id %in% idxs)
      mat_C <- as.matrix(sel[, grep("^C\\d+", colnames(sel))])
      mat_S <- as.matrix(sel[, grep("^S\\d+", colnames(sel))])
      rownames(mat_C) <- rownames(mat_S) <- sel$id
      
      ok_cons <- !require_full_consumption || all(colSums(mat_C) > 0)
      if (ok_cons) break
      
      attempt <- attempt + 1
      if (attempt > max_attempts) {
        stop("Cobertura no lograda tras ", max_attempts,
             " intentos en réplica ", rep_i, " de ", comb_name)
      }
    }
    
    matrices_C[[rep_i]] <- mat_C
    matrices_S[[rep_i]] <- mat_S
  }
  
  # Guardar resultados
  saveRDS(
    list(matrices_C = matrices_C, matrices_S = matrices_S),
    paste0("matrices_", comb_name, "_", reps_per_comb, "rep.rds")
  )
  rm(matrices_C, matrices_S)
  gc()
}

cat("Simulación completada.\n")

