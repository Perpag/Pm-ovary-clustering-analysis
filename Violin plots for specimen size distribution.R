library(ggplot2)

#Fig.S2a

ggplot(Supplementary_Data_2,
       aes(x = "oogonia in germ cell nests",
           y = `measured diameter (μm)`)) +
  
  # violin shape
  geom_violin(fill = "grey90", color = "black", adjust = 0.5) +
  
  # individual points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  # mean point
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               size = 0.8) +
  
  ylab("diameter (µm)") +
  xlab("") +
  theme_classic(base_size = 28) 

ggsave("name.tiff", units="in", width=8, height=9, dpi=300, compression = 'lzw')

#Fig.S2b
ggplot(Supplementary_Data_2,
       aes(x = "oogonia in germ cell nests",
           y = `estimated volume [μm^3]`)) +
  
  # violin shape
  geom_violin(fill = "grey90", color = "black", adjust = 0.5) +
  
  # individual points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  # mean point
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               size = 0.8) +
  
  ylab("Volume [μm^3]") +
  xlab("") +
  theme_classic(base_size = 28)

ggsave("name.tiff", units="in", width=8, height=9, dpi=300, compression = 'lzw')


#Fig.S2c
ggplot(Supplementary_Data_2,
       aes(x = "smallerst nanos 1+ cell adjacent to germ cell nests",
           y = `measured diameter (μm)`)) +
  
  # violin shape
  geom_violin(fill = "grey90", color = "black", adjust = 0.5) +
  
  # individual points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  # mean point
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               size = 0.8) +
  
  ylab("diameter (µm)") +
  xlab("") +
  theme_classic(base_size = 28) 

ggsave("name.tiff", units="in", width=8, height=9, dpi=300, compression = 'lzw')

#Fig.S2d
ggplot(Supplementary_Data_2,
       aes(x = "smallerst nanos 1+ cell adjacent to germ cell nests",
           y = `estimated volume [μm^3]`)) +
  
  # violin shape
  geom_violin(fill = "grey90", color = "black", adjust = 0.5) +
  
  # individual points
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  
  # mean point
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.1,
               size = 0.8) +
  
  
  ylab("Volume [μm^3]") +
  xlab("") +
  theme_classic(base_size = 28)

ggsave("name.tiff", units="in", width=8, height=9, dpi=300, compression = 'lzw')

#Fig.S2e

ggplot(Supplementary_Data_2,
       aes(x = Material,  
           y = `Volume [μm^3]`)) +
  
  # violin
  geom_violin(fill = "grey90",
              color = "black",
              adjust = 0.5,
              width = 0.8,
              trim = FALSE) +
  
  # jittered points
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.8) +
  
  # mean point (diamond)
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.15,
               size = 0.8) +
  
  ylab("Volume [μm^3]") +
  xlab("") +
  theme_classic(base_size = 28)

ggsave("name.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')

#Fig.S2f

ggplot(Supplementary_Data_2,
       aes(x = Material,  
           y = `Diameter [μm] (estimated)`)) +
  
  # violin
  geom_violin(fill = "grey90",
              color = "black",
              adjust = 0.5,
              width = 0.8,
              trim = FALSE) +
  
  # jittered points
  geom_jitter(width = 0.15,
              size = 2,
              alpha = 0.8) +
  
  # mean point (diamond)
  stat_summary(fun = mean,
               geom = "point",
               shape = 23,
               size = 4,
               fill = "black") +
  
  # mean ± SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.15,
               size = 0.8) +
  
  ylab("diameter (µm)") +
  xlab("") +
  theme_classic(base_size = 28) 

ggsave("name.tiff", units="in", width=10, height=9, dpi=300, compression = 'lzw')
