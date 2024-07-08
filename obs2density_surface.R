# Density surface modelling
# based on: https://distancesampling.org/R/vignettes/mexico-analysis.html
library(dsm)
library(ggplot2)

# plotting options
gg.opts <- theme(panel.grid.major=element_blank(),
                panel.grid.minor=element_blank(),
                panel.background=element_blank())

# make the results reproducible
set.seed(11123)

# load data 
data <- sf::st_read('D:\\WMU\\survey_data\\501_moose_locations.shp')
print(data)
plot(data)
