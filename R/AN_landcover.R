
lakes_lc[, c(-1, -2)] <- (lakes_lc[, c(-1, -2)] / rowSums(lakes_lc[, c(-1, -2)]))

lakes_lc <- lakes_lc %>% 
  mutate(Artificial_surface =
           `Tissu urbain discontinu   (ha)` +
           `Tissu urbain continu` +
           `Equipements sportifs et de loisirs  (ha)` +
           `Zones industrielles ou commerciales et installations publiques` +
           `Réseaux routier et ferroviaire et espaces associés` +
           Aéroports +
           `Extraction de matériaux` +
           Décharges +
           Chantiers +
           `Espaces verts urbains`,
         Agricultural_surface =
           Vignobles +
           `Vergers et petits fruits` +
           `Prairie et autres surfaces toujours en herbe à usage agricole  (ha)` +
           `Systèmes culturaux complexes   (ha)` +
           `Terres arables hors périmètres d'irrigation   (ha)` +
           `Espaces agricoles interrompus par espace naturels  (ha)`,
         Natural_surface =
           `Glaciers et neiges éternelles  (ha)` +
           `Marais intérieurs (ha)` +
           `Tourbières  (ha)` +
           `Plans d'eau` +
           `Pelouses et pâturages naturels  (ha)` +
           `Forêt de conifères (ha)` +
           `Forêt de Feuillus  (ha)` +
           `Forêt mélangée  (ha)` +
           `Forêt et végétation arbustive en mutation   (ha)` +
           `Landes et broussailles  (ha)` +
           `Végétation clairsemée  (ha)` +
           `Roches Nues  (ha)` +
           `Sable, plage et dune   (ha)` +
           `Végétation sclérophylle` +
           `Marais maritimes`)

lakes_lc <- lakes_lc %>%
  mutate(Anthropised_surface = Artificial_surface + Agricultural_surface,
         Natural_Artificial_ratio = (Natural_surface + 1) / (Artificial_surface + 1),
         Natural_Anthropised_ratio = (Natural_surface + 1) / (Anthropised_surface + 1),
         Altitude_category = factor(ifelse(`Altitude (m)` > 1400, "High", "Low")))


lc_plot1 <- ggplot(lakes_lc) +
  geom_boxplot(aes(fct_shift(Altitude_category), Natural_surface)) +
  xlab("Altitude") + ylab("Natural surface") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic()

lc_plot2 <- ggplot(lakes_lc) +
  geom_boxplot(aes(fct_shift(Altitude_category), Agricultural_surface)) +
  xlab("Altitude") + ylab("Agricultural surface") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic()

lc_plot3 <- ggplot(lakes_lc) +
  geom_boxplot(aes(fct_shift(Altitude_category), Artificial_surface)) +
  xlab("Altitude") + ylab("Artificial surface") +
  scale_y_continuous(labels=scales::percent) +
  theme_classic()


lakes_pop <- lakes_pop %>% 
  left_join(lakes_bv) %>% 
  mutate(Altitude_category = factor(ifelse(ALTITUDE > 1400, "High", "Low")))

pop_plot <- ggplot(lakes_pop) +
  geom_boxplot(aes(fct_shift(Altitude_category), POP20km)) +
  scale_y_continuous(trans = scales::log10_trans(),
                     breaks = c(1000, 10000, 100000),
                     labels = scales::math_format(expr = 10^.x, format = force)(3:5)) +
  xlab("Altitude") + ylab("Population") +
  theme_classic()

lc_plots <- cowplot::plot_grid(lc_plot1, lc_plot2, lc_plot3, pop_plot,
                   labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2)

rm(lc_plot1, lc_plot2, lc_plot3, pop_plot)