rm(list = ls())   # Remove everything from the Global environment

library(tidyverse)
library(rmutil)   # For generating from pareto distribution
source("distributionsL.R")

## Plot the true valuation distributions.
x.vec <- seq(0, 20, by = 0.01)

F.unif <- punif(x.vec, min = 1, max = 20)
F.piecewise.unif <- sapply(x.vec, function(x){piecewise.unif.cdf(x)} )
F.pareto <- ppareto(x.vec, m = 3, s = 100)
F.gamma <- pgamma(x.vec, shape = 10, rate = 2)
F.beta <- pbeta(x.vec, shape1 = 2, shape2 = 2)


F.tibble <- tibble(F.unif, F.piecewise.unif, F.pareto, F.gamma, F.beta)
colnames(F.tibble) <- c("uniform", "piecewise_uniform", "pareto", "gamma", "beta")

# ## color scheme
# legend.vec <- c("Uniform" = "darkred", 
#                 "Piecewise Uniform" = "darkorange", 
#                 "Pareto" = "darkviolet", 
#                 "Gamma" = "green", 
#                 "Beta" = "steelblue")
# 
# true.F.plots <- ggplot(data = F.tibble, aes(x = x.vec)) +
#   geom_line(mapping = aes(y = uniform, color = "Uniform")) + 
#   geom_line(mapping = aes(y = piecewise_uniform, color = "Piecewise Uniform")) +
#   geom_line(mapping = aes(y = pareto, color = "Pareto")) +
#   geom_line(mapping = aes(y = gamma, color = "Gamma")) +
#   geom_line(mapping = aes(y = beta, color = "Beta")) +
#   labs(x = "x", y = "F(x)", color = "True F: ") + 
#   scale_color_manual(values = legend.vec) +
#   theme(legend.position = "top")

## Line type scheme
lt.vec    <- c("Uniform" = "longdash", "Piecewise Uniform" = "dashed",
               "Pareto" = "dotted", "Gamma" = "dotdash", "Beta" = "solid")
shape.vec <- c("Uniform" = 2, "Piecewise Uniform" = 1,
               "Pareto" = 4, "Gamma" = 0, "Beta" = 8)

F.long <- F.tibble %>%
  mutate(x = x.vec) %>%
  pivot_longer(cols = c(uniform, piecewise_uniform, pareto, gamma, beta),
               names_to = "dist", values_to = "F") %>%
  mutate(dist = recode(dist,
                       uniform = "Uniform",
                       piecewise_uniform = "Piecewise Uniform",
                       pareto = "Pareto",
                       gamma = "Gamma",
                       beta = "Beta"),
         dist = factor(dist, levels = names(lt.vec)))

# One marker per 1 unit of x (100 grid steps), staggered by 0.2 per curve
F.points <- F.long %>%
  group_by(dist) %>%
  filter((row_number() - 1 - 20 * (as.integer(dist) - 1)) %% 50 == 0) %>%
  ungroup()

true.F.plots <- ggplot(F.long, aes(x = x, y = F,
                                   linetype = dist, shape = dist)) +
  geom_line() +
  geom_point(data = F.points, size = 2) +
  labs(x = "x", y = "F(x)", linetype = "True F: ", shape = "True F: ") +
  scale_linetype_manual(values = lt.vec) +
  scale_shape_manual(values = shape.vec) +
  theme(legend.position = "inside",
        legend.position.inside = c(1, 0),
        legend.justification = c(1, 0),
        legend.background = element_rect(color = "grey50"),
        legend.key.width = unit(1.5, "cm"))
true.F.plots

ggsave(filename = "Plots/TrueF_plots_BW.png", plot = true.F.plots, height = 6, width = 8)