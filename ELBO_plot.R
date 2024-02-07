
# generate the ELBO plot
ELBO_df <- readRDS(file = "./ELBO_df.rds")

ggplot(data = ELBO_df, 
       aes(x = iter, y = ELBO, group = rep)) + 
  geom_line(alpha = 0.1) + 
  geom_point(alpha = 0.1, size=0.2) +
  facet_wrap(~ scene, scale = "free") + 
  labs(x = "Iteration", y = "ELBO") + 
  scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())