library(RLdbRDA)
library(vegan)

devtools::load_all()

data(varespec)
data(varechem)

out <- rldbrda(varespec, varechem)
out

plot_data <- prepare_plot_data(out)
plot_data

g <- plot_dbrda(plot_data)
g

