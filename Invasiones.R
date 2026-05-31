
library(MASS)
library(deSolve)
library(lpSolve)
library(tidyverse)
library(writexl)
library(readxl)

base_dir <- "C:/Users/danie/Documents/Master_Bioinformatica_y_Biologia_Computacional/TFM/scripts de Abel para Daniel/scripts de Abel para Daniel"
setwd(base_dir)

source(file.path(base_dir, "convergence_param_func2.R"))
source(file.path(base_dir, "dynamic_integrator.R"))
source(file.path(base_dir, "event_func2.R"))

n_spec <- 20
n_res <- 20
tau <- 0.1

times <- seq(0, 10000, by = 2)

# Cargamos las especies
Pool <- readRDS("Pool_200.rds")

# Función para todas las combinaciones no vacías
all_groups   <- c("SD","SC","GD","GC")
get_combinations <- function(vec) {
  unlist(lapply(seq_along(vec),
                function(k) combn(vec, k, simplify = FALSE)),
         recursive = FALSE)
}
all_combs  <- get_combinations(all_groups)
comb_names <- sapply(all_combs, paste, collapse = "+")

# Definimos los distintos tamaños de recursos
sizes_res <- c("S", "M", "L")

# Recorremos cada comunidad
for (current_comb in comb_names) {
  file_community <- paste0("matrices_", current_comb, "_20rep.rds")
  community <- readRDS(file_community)
  
  # Recorremos cada tamaño
  for (size_r in sizes_res) {
    file_results <- paste0("Results_matrices_", current_comb, "_20rep_", size_r, ".xlsx")
    
    assembly <- read_excel(file_results)
    
    results_cases <- data.frame(
      pair_id = integer(),
      l = numeric(),
      kappa = integer(),
      species_id = character(),
      abundance = numeric(),
      case_letter = character(),
      supplied_idx = character(),
      abundance_invasion = numeric()
    )
    
    for (i in seq(1, nrow(assembly_res), by=20)) {
      current_assembly <- assembly[i:(i+19),]
      
      # Este comando sirve para obtener un vector de los recursos necesarios
      ids_res <- as.numeric(str_split(current_assembly$supplied_idx[1], ",")[[1]])
      resources <- rep(0, n_res)
      resources[ids_res] <- 1
      
      res_abundance <- current_assembly$res_abundance
      
      # Especies de la comunidad y 
      species_ids <- current_assembly$species_id
      
      # Id de la comunidad
      community_id <- current_assembly$pair_id[1]
      
      # Parámetros del modelo
      kappa <- current_assembly$kappa[1]
      l <- current_assembly$l
      
      for (spec_type in all_groups) {
        # Tomamos 10 especies del tipo spec_type que no estén en la comunidad
        invasors <- sample(setdiff(Pool[Pool$group==spec_type,]$id, species_ids), 10)
        
        for (invasor in invasors) {
          # Creamos el data frame donde guardaremos la información
          results_invasion <- rbind(current_assembly, c(community_id, l[1], kappa, invasor, 0.01, size_r, current_assembly$supplied_idx[1]))
          
          #Tomamos el vector de consumo y secreción del invasor
          invasor_Cvec <- as.matrix(Pool[Pool$id==invasor,][4:23])
          invasor_Svec <- as.matrix(Pool[Pool$id==invasor,][24:43])
          
          # Creamos la matriz resultante y reinicializamos las anteriores
          Cmat <- community$matrices_C[[community_id]]
          Smat <- community$matrices_S[[community_id]]
          
          Cmat <- rbind(Cmat, invasor_Cvec)
          Smat <- rbind(Smat, invasor_Svec)
          
          init_ab <- c(current_assembly$abundance, 0.01, res_abundance)
          names(init_ab) <- c(species_ids, invasor, paste0("R",1:n_res))
          
          parms <- list(
            Cmat      = Cmat,
            Smat      = Smat,
            l         = rep(l[1], n_spec + 1),
            g         = 100,
            w         = rep(1, n_res),
            chi       = rep(1, n_res),
            mu        = rep(0.5, n_spec + 1),
            tau       = 0.1,
            kappa     = kappa*resources,
            n_spec    = n_spec + 1,     # número de especies
            n_res     = n_res,      # número de recursos
            threshold = 1e6   # umbral de abundancia
          )
          
          convergence_param_func2(init_ab, times)
          
          out <- ode(
            y        = init_ab,
            times    = times,
            func     = CRm_Modsecret,
            parms    = parms,
            rootfunc = event_func2
          )
          
          # extraer ultima
          final <- tail(as.data.frame(out),1)
          # sin filtro de umbral
          # extraer abundancias de especies
          final[] <- lapply(final, function(col) {
            if (is.numeric(col)) replace(col, col < 1e-3, 0) else col
          })
          
          # Tomamos las abundancias después de la invasión
          abunds_invasion <- as.numeric(final[1, c(species_ids, invasor)])
          
          # Guardamos los datos en el data frame
          results_invasion$abundance_invasion <- abunds_invasion
          
          # Añadimos los resultados a results_cases
          results_cases <- rbind(results_cases, results_invasion)
        }
      }
    }
    
    out_name <- paste0("Results_matrices_", current_comb, "_20rep_", size_r, "_invasions.xlsx")
    write_xlsx(results_cases, out_name)
    cat("Guardado:", out_name, "->", nrow(results_cases), "filas\n")
    
  }
}