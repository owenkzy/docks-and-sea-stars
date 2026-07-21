#useful R code

#reading in csv----

straight1 <- read_csv("Data/ch_a_c1_master.csv")

#making factors and reorganizing levels----

%>%
  mutate(type = factor(type, levels = c("dock", "inter-dock", "no-dock")))
#or 
%>% 
  mutate(type = as.factor(type))

#DHARMa

simulateResiduals(model, plot = TRUE)

#dot whisker plot examples below----

ggplot(data = data_long,
       aes(x = movie_type, y = score,
           colour = movie_type)) +
  
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.3),
    alpha = 0.4,
    size = 2) +
  
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  labs(x = "Type of movie",
       y = "score") +
  pub_theme() +
  theme(legend.position = "none")
ggsave("ssuno.png",  width = 6, height = 4, units = "in", dpi = 300)

#another dot whisker example

ggplot(data = data_long,
       aes(x = scoring_type, y = score,
           colour = movie_type, group = movie_type)) +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1,
    dodge.width = 0.3),
    alpha = 0.4,
    size = 2) +
  stat_summary(fun = mean,
               geom = "point",
               position = position_dodge(width = 0.3),
               size = 3) +
  stat_summary(fun.data = mean_cl_normal,
               geom = "errorbar",
               position = position_dodge(width = 0.3),
               width = 0.2) +
  facet_wrap_paginate(~movie) +
  pub_theme()


#pub theme below----

pub_theme1 <- function() {
  theme_classic() +
    theme(
      # axis
      axis.text = element_text(size = 11, colour = "black"),
      axis.title = element_text(size = 11, colour = "black"),
      axis.line = element_line(colour = "black", linewidth = 0.1),
      axis.ticks = element_line(color = "black", linewidth = 0.2),
      axis.ticks.length = unit(0.1, "lines"),
      # legend
      legend.title = element_text(size = 11, colour = "black"),
      legend.text = element_text(size = 11, colour = "black"),
      # panel border
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.5),
      panel.spacing = unit(0.1, "lines")
    ) 
}

