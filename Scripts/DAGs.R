## ---------------------------
##
## Script name: 
##
## Author: Dr. Joan Dudney
##
## Date Created: 2023-05-30
##
## Copyright (c) Joan Dudney, 2023
## Email: dudney@ucsb.edu
##
## ---------------------------
##
## Notes: creating a dag following
##   https://cran.r-project.org/web/packages/ggdag/vignettes/intro-to-dags.html
##
## ---------------------------

librarian::shelf(ggdag, ggplot2)
theme_set(theme_dag())

## simple
dagify(y ~ x) %>%
  ggdag()

## nonlinear
dagify(y ~ ~x) %>%
  ggdag()

## bidirectional
dagify(y ~ ~x) %>%
  ggdag_canonical()

## not a dag becuase it's cyclic
dagify(
  y ~ x,
  x ~ a
) %>%
  ggdag()


whitebark <- dagify(growth ~ precip + temp + soils + density,
                         density ~ precip + temp + soils,
                         soils ~ precip + temp,
                         labels = c(
                           "growth" = "Tree\n Growth",
                           "precip" = "Precipitation",
                           "temp" = "Temperature",
                           "soils" = "Soils",
                           "density" = "Stand\n Density"
                         ),
                         latent = "soils",
                         exposure = "precip",
                         outcome = "growth"
)

whitebark <- dagify(growth ~ precip + temp + soils + density,
                    density ~ precip + temp + soils,
                    soils ~ precip + temp,
                    latent = "soils",
                    exposure = "precip",
                    outcome = "growth"
)

ggdag(whitebark, text = T)


concept <- dagify(A ~ B + C,
                    B ~ C,
                    outcome = "A"
)

ggdag(concept, text = T)

