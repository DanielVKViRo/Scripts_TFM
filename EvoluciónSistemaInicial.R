# Limpieza total y carga de funciones
rm(list = ls())

library(MASS)
library(deSolve)
library(lpSolve)
library(tidyverse)
library(writexl)

base_dir <- "C:/Users/danie/Documents/Master_Bioinformatica_y_Biologia_Computacional/TFM/scripts de Abel para Daniel/scripts de Abel para Daniel"

source(file.path(base_dir, "convergence_param_func2.R"))
source(file.path(base_dir, "dynamic_integrator.R"))
source(file.path(base_dir, "event_func2.R"))

n_spec     <- 20
n_res      <- 20
tau        <- 0.1
chi        <- rep(1, n_res)
abundancia <- rep(5, n_spec)

# valores que uso de kappa
kappa_list <- c(10, 60, 100)

# l en el rango que uso
l_list     <- c(0.1, 0.5, 0.9)

# tiempos de simulación
times <- seq(0, 10000, by = 2)

# pool con g y mu
g_mu_pool <- readRDS(file.path(base_dir, "Pool_200.rds"))

# archivos de matrices
mat_files <- list.files(base_dir, pattern = "^matrices.*\\.rds$", full.names = TRUE)
if(length(mat_files)==0) stop("No hay archivos matrices .rds")

for(mat_path in mat_files){
  base_name <- tools::file_path_sans_ext(basename(mat_path))
  mats      <- readRDS(mat_path)
  S_list    <- mats$matrices_S
  C_list    <- mats$matrices_C
  n_sets    <- length(S_list)
  
  # resultados por caso (ahora L, M, S (5 recursos) y XS (2 recursos))
  results_cases <- list(L = tibble(), M = tibble(), S = tibble())
  results_cases_res <- list(L = tibble(), M = tibble(), S = tibble())
  
  for(set_id in seq_len(n_sets)){
    Smat <- S_list[[set_id]]
    Cmat <- C_list[[set_id]]
    species_ids <- rownames(Smat)
    if(is.null(species_ids)) stop("Sin nombres en filas de Smat")
    
    # Selecciones de recursos (una por set)
    sel_half <- sample.int(n_res, n_res %/% 2)  # mitad
    sel_five <- sample.int(n_res, 5)            # S: 5 recursos
    
    # plantillas (0/1) para cada caso: L, M, S(=5), XS(=2)
    kappa_cases <- list(
      L  = rep(1, n_res),   # todos
      M  = numeric(n_res),  # mitad
      S  = numeric(n_res)  # cinco recursos
    )
    kappa_cases$M[sel_half] <- 1
    kappa_cases$S[sel_five] <- 1
    
    for(lB in l_list){
      for(kB in kappa_list){
        # condición inicial común (nivel basal en recursos = kB * tau)
        init_common <- c(abundancia, rep(kB * tau, n_res))
        names(init_common) <- c(species_ids, paste0("R", 1:n_res))
        
        for(case_letter in names(kappa_cases)){
          # construir vector kappa para la dinámica (vector length n_res)
          kappa_vec <- kappa_cases[[case_letter]] * kB
          
          # REINICIALIZAR la comprobación de convergencia para ESTA integración
          convergence_param_func2(init_common, times)
          
          # parámetros para ode (kappa como vector)
          parms <- list(
            Cmat      = Cmat,
            Smat      = Smat,
            l         = rep(lB, n_spec),
            g         = 100,
            w         = rep(1, n_res),
            chi       = chi,
            mu        = 0.5,
            tau       = tau,
            kappa     = kappa_vec,
            n_spec    = n_spec,
            n_res     = n_res,
            threshold = 1e6
          )
          
          # por seguridad (si tu event_func2 lo usa)
          parms$t_inside <- times
          
          # ejecutar integración con manejo de errores (no detener todo el loop)
          out <- tryCatch(
            {
              ode(y = init_common,
                  times = times,
                  func = CRm_Modsecret,
                  parms = parms,
                  rootfunc = event_func2)
            },
            error = function(e){
              warning(sprintf("ode fallo en set %d l=%g k=%g case=%s : %s",
                              set_id, lB, kB, case_letter, e$message))
              return(NULL)
            }
          )
          
          # extraer resultado final (o NA si falló)
          if(is.null(out)){
            abunds <- rep(NA_real_, length(species_ids))
          } else {
            final <- tail(as.data.frame(out), 1)
            final[] <- lapply(final, function(col) if(is.numeric(col)) replace(col, col < 1e-3, 0) else col)
            abunds <- as.numeric(final[1, species_ids])
          }
          
          res <- as.numeric(final[1, paste0("R", seq(1, 20, 1))])
          
          # guardar fila por especie
          results_cases[[case_letter]] <- bind_rows(
            results_cases[[case_letter]],
            tibble(
              pair_id    = set_id,
              l          = lB,
              kappa      = kB,
              species_id = species_ids,
              abundance  = abunds,
              res_abundance = res,
              case_letter = case_letter,
              supplied_idx = paste0(which(kappa_vec > 0), collapse = ",")
            )
          )
          
        } # end for case_letter
      } # end for kB
    } # end for lB
  } # end for set_id
  
  # guardar un archivo por caso (sufijos _L _M _S)
  suffix_map <- c(L = "_L", M = "_M", S = "_S")
  for(case_letter in names(results_cases)){
    out_name <- file.path(base_dir, paste0("Results_", base_name, suffix_map[case_letter], ".xlsx"))
    write_xlsx(results_cases[[case_letter]], out_name)
    cat("Guardado:", basename(out_name), "->", nrow(results_cases[[case_letter]]), "filas\n")
  }
  
}

