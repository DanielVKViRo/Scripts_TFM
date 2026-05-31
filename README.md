# Código del Trabajo Fin de Máster

Este repositorio contiene el código utilizado en el Trabajo Fin de Máster titulado:

**Efectos de las interacciones metabólicas en la diversidad y estabilidad de las comunidades microbianas**

## Contenido

El repositorio incluye los scripts empleados para:

* Generar el pool inicial de especies y construir comunidades a partir de él.
* Obtener el estado de equilibrio de las comunidades iniciales.
* Realizar los experimentos de invasión.
* Analizar los resultados obtenidos y generar las figuras incluidas en la memoria.

Las figuras finales utilizadas en el trabajo se encuentran incluidas en el repositorio. Los resultados completos de las simulaciones no se han incluido debido a su volumen, pero pueden reproducirse mediante los scripts proporcionados.

## Descripción de los scripts

* `GeneracionEspecies.R`: genera el pool inicial de especies y construye las comunidades.
* `EvolucionSistemaInicial.R`: evoluciona las comunidades iniciales hasta alcanzar el equilibrio. Para ello utiliza los scripts auxiliares `dynamic_integrator.R`, `event_func2.R` y `convergence_param_func2.R`.
* `Invasiones.R`: realiza los experimentos de invasión sobre las comunidades en equilibrio.
* `analisis_invasion.R`: analiza los resultados de las invasiones y genera las figuras correspondientes.

## Autor

Daniel Vidales Rosa

Máster en Bioinformática y Biología Computacional

Universidad Autónoma de Madrid
