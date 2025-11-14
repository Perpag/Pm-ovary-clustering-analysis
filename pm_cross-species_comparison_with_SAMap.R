library(ggalluvial)
library(networkD3)
library(dplyr)
library(ggplot2)
library(reshape2)
library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(stringr)
library(ggnewscale)


#In R Studio

pm<-pm_ovary_integrated_umap
sp<-sp_ovary
spl<-sp_3dpf_integrated_umap
dm<-dm_ovary
dr<-dr_ovary
hs<-hs_ovary
mm<-mm_ovary
mm_pg<-pg_integrated

pm[["RNA3"]] <- as(object = pm[["RNA"]], Class = "Assay")
DefaultAssay(pm) <- "RNA3"
pm[["RNA"]] <- NULL
pm <- RenameAssays(object = pm, RNA3 = 'RNA')
SaveH5Seurat(pm, filename = "pm_ovary.h5Seurat")
Convert("pm_ovary.h5Seurat", dest = "h5ad")

sp[["RNA3"]] <- as(object = sp[["RNA"]], Class = "Assay")
DefaultAssay(sp) <- "RNA3"
sp[["RNA"]] <- NULL
sp <- RenameAssays(object = sp, RNA3 = 'RNA')
SaveH5Seurat(sp, filename = "sp_adult_ovary.h5Seurat")
Convert("sp_adult_ovary.h5Seurat", dest = "h5ad")

spl[["RNA3"]] <- as(object = spl[["RNA"]], Class = "Assay")
DefaultAssay(spl) <- "RNA3"
spl[["RNA"]] <- NULL
spl <- RenameAssays(object = spl, RNA3 = 'RNA')
SaveH5Seurat(spl, filename = "sp_3df_larva.h5Seurat")
Convert("sp_3df_larva.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
mm <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_adult_ovary.h5Seurat")
Convert("mm_adult_ovary.h5Seurat", dest = "h5ad")

dr[["RNA3"]] <- as(object = dr[["RNA"]], Class = "Assay")
DefaultAssay(dr) <- "RNA3"
dr[["RNA"]] <- NULL
dr <- RenameAssays(object = dr, RNA3 = 'RNA')
SaveH5Seurat(dr, filename = "dr_adult_ovary.h5Seurat")
Convert("dr_adult_ovary.h5Seurat", dest = "h5ad")

dm[["RNA3"]] <- as(object = dm[["RNA"]], Class = "Assay")
DefaultAssay(dm) <- "RNA3"
dm[["RNA"]] <- NULL
dm <- RenameAssays(object = dm, RNA3 = 'RNA')
SaveH5Seurat(dm, filename = "dm_adult_ovary.h5Seurat")
Convert("dm_adult_ovary.h5Seurat", dest = "h5ad")

hs[["RNA3"]] <- as(object = hs[["RNA"]], Class = "Assay")
DefaultAssay(hs) <- "RNA3"
hs[["RNA"]] <- NULL
hs <- RenameAssays(object = hs, RNA3 = 'RNA')
SaveH5Seurat(hs, filename = "hs_fetal_ovary.h5Seurat")
Convert("hs_fetal_ovary.h5Seurat", dest = "h5ad")

mm[["RNA3"]] <- as(object = mm[["RNA"]], Class = "Assay")
DefaultAssay(mm) <- "RNA3"
mm[["RNA"]] <- NULL
hs <- RenameAssays(object = mm, RNA3 = 'RNA')
SaveH5Seurat(mm, filename = "mm_pituitary_gland.h5Seurat")
Convert("mm_pituitary_gland.h5Seurat", dest = "h5ad")

#In python

cd samap_directory

python3

from samap.mapping import SAMAP
from samap.analysis import (get_mapping_scores, GenePairFinder,
                            sankey_plot, chord_plot, CellTypeTriangles, 
                            ParalogSubstitutions, FunctionalEnrichment,
                            convert_eggnog_to_homologs, GeneTriangles)
from samalg import SAM
import pandas as pd

% bash map_genes.sh --tr1 /directory of P.miniata proteome fasta file --t1 prot --n1 pm --tr2 //directory of species tested proteome fasta file --t2 prot --n2 xx #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

fn1 = "directory of pm_ovary.h5ad"
fn2 = "directory of species compared to"

filenames = {'pm':fn1,'xx':fn2} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison

sm = SAMAP(
  filenames,
  f_maps = 'data/maps/',
  save_processed=True #if False, do not save the processed results to `*_pr.h5ad`
)

sm.run(pairwise=True)
samap = sm.samap # SAM object with 2 species stitched together
keys = {'pm':'seurat_clusters','xx':'seurat_clusters'} #xx= spl or dm or dr or hs or mm or sp or mm_pg depending on the comparison
D,MappingTable = get_mapping_scores(sm,keys,n_top = 0)
D.head()
df = MappingTable

print(df)

df.to_csv('mapping_table_sea_star_vs_x_species.csv')

#In R Studio

#Sankey plots for all apart from for human fetal ovary

df <- load mapping scores
data <- melt(df)
data$X <- as.factor(data$X)
dat <- data %>%
  transmute(
    species_left  = as.character(X),
    species_right = as.character(variable),
    value         = suppressWarnings(as.numeric(value))
  ) 


pm_left  <- mean(grepl("_pm", dat$species_left),  na.rm = TRUE)
pm_right <- mean(grepl("_pm", dat$species_right), na.rm = TRUE)

if (is.finite(pm_left) && is.finite(pm_right) && pm_left < pm_right) {
  tmp <- dat$species_left
  dat$species_left  <- dat$species_right
  dat$species_right <- tmp
}


left_num <- as.integer(stringr::str_extract(dat$species_left, "\\d+$"))

dat <- dat %>%
  mutate(
    cluster_left = case_when(
      left_num %in% 1:4   ~ "epithelium",
      left_num %in% 5:8   ~ "follicle I",
      left_num %in% 9:12  ~ "follicle II",
      left_num %in% 13:14 ~ "germ cell nest",
      left_num %in% 15:19 ~ "oocytes",
      left_num == 20      ~ "neurons",
      left_num %in% 21:22 ~ "muscles",
      left_num %in% 23:25 ~ "follicle III",
      left_num %in% 26:27 ~ "immune",
      TRUE ~ "Unknown"
    )
  )


cluster_cols <- c(
  "epithelium"     = "#63A8FF",
  "follicle I"     = "#F27676",
  "follicle II"    = "#F4A300",  
  "germ cell cyst" = "#21FF00",  
  "oocytes"        = "#009F72",
  "neurons"        = "#FFEF00",
  "muscles"        = "#9B6BFF",
  "follicle III"   = "#FFB6C1",
  "immune"         = "#ADD8E6",
  "unknown"        = "grey70"
)

# Left boxes use a per-stratum palette derived from cluster_cols
stratum_colors <- dat %>%
  distinct(species_left, cluster_left) %>%
  mutate(color = cluster_cols[as.character(cluster_left)]) %>%
  { setNames(.$color, .$species_left) }


dat$flow_id <- seq_len(nrow(dat))

# Ordering vectors
left_levels  <- unique(dat$species_left[order(left_num)]) 
right_levels <- unique(dat$species_right)

# Common factor levels so layers combine cleanly
all_levels <- c(left_levels, setdiff(right_levels, left_levels))

lodes <- bind_rows(
  transmute(dat, side = "Left",  stratum = species_left,  flow_id, value, cluster_left),
  transmute(dat, side = "Right", stratum = species_right, flow_id, value, cluster_left)
) %>%
  mutate(
    x_num   = ifelse(side == "Left", 1L, 2L),                         # 1 = Left, 2 = Right
    stratum = factor(stratum, levels = all_levels, ordered = TRUE)
  )


ggplot(lodes, aes(x = x_num, stratum = stratum, alluvium = flow_id, y = value)) +
  # Colored connections (by LEFT cluster)
  geom_alluvium(aes(fill = cluster_left), width = 0.25, alpha = 0.85) +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  
  ggnewscale::new_scale("fill") +
  
  # RIGHT boxes (white with black borders)
  geom_stratum(
    data = dplyr::filter(lodes, x_num == 2L),
    fill = "white", color = "black", width = 0.25
  ) +
  # LEFT boxes (colored by their own label using the derived palette)
  geom_stratum(
    data = dplyr::filter(lodes, x_num == 1L),
    aes(fill = after_stat(stratum)), color = "black", width = 0.25
  ) +
  geom_text(
    data = dplyr::filter(lodes, x_num %in% c(1L, 2L)),
    stat = "stratum", aes(label = after_stat(stratum)), size = 3
  ) +
  scale_fill_manual(values = stratum_colors, guide = "none") +
  
  # Lock axis positions and labels
  scale_x_continuous(breaks = c(1, 2), labels = c("Left", "Right"), expand = c(0.15, 0.05)) +
  labs(y = "value", x = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )
ggsave("Sankey_plot_pm_ovary_vs_x species.tiff", units="in", width=12, height=13
       , dpi=300, compression = 'lzw')


#Sankey plot for human fetal ovary


data <- melt(df)
data$X <- as.factor(data$X)
dat <- data %>%
  transmute(
    species_left  = as.character(X),
    species_right = as.character(variable),
    value         = suppressWarnings(as.numeric(value))
  ) 

pm_left  <- mean(grepl("_pm", dat$species_left),  na.rm = TRUE)
pm_right <- mean(grepl("_pm", dat$species_right), na.rm = TRUE)

if (is.finite(pm_left) && is.finite(pm_right) && pm_left < pm_right) {
  tmp <- dat$species_left
  dat$species_left  <- dat$species_right
  dat$species_right <- tmp
}

left_num <- as.integer(stringr::str_extract(dat$species_left, "\\d+$"))

dat <- dat %>%
  mutate(
    cluster_left = case_when(
      left_num %in% 1:4   ~ "epithelium",
      left_num %in% 5:8   ~ "follicle I",
      left_num %in% 9:12  ~ "follicle II",
      left_num %in% 13:14 ~ "germ cell cyst",
      left_num %in% 15:19 ~ "oocytes",
      left_num == 20      ~ "neurons",
      left_num %in% 21:22 ~ "muscles",
      left_num %in% 23:25 ~ "follicle III",
      left_num %in% 26:27 ~ "immune",
      TRUE ~ "Unknown"
    )
  )

cluster_cols <- c(
  "epithelium"     = "#63A8FF",
  "follicle I"     = "#F27676",
  "follicle II"    = "#F4A300",  
  "germ cell cyst" = "#21FF00",  
  "oocytes"        = "#009F72",
  "neurons"        = "#FFEF00",
  "muscles"        = "#9B6BFF",
  "follicle III"   = "#FFB6C1",
  "immune"         = "#ADD8E6",
  "unknown"        = "grey70"
)

# Left boxes use a per-stratum palette derived from cluster_cols
stratum_colors <- dat %>%
  distinct(species_left, cluster_left) %>%
  mutate(color = cluster_cols[as.character(cluster_left)]) %>%
  { setNames(.$color, .$species_left) }

dat$flow_id <- seq_len(nrow(dat))

# Ordering vectors
left_levels  <- unique(dat$species_left[order(left_num)])  # pm1..pm27 order by suffix
right_levels <- unique(dat$species_right)

# Common factor levels so layers combine cleanly
all_levels <- c(left_levels, setdiff(right_levels, left_levels))

lodes <- bind_rows(
  transmute(dat, side = "Left",  stratum = species_left,  flow_id, value, cluster_left),
  transmute(dat, side = "Right", stratum = species_right, flow_id, value, cluster_left)
) %>%
  mutate(
    x_num   = ifelse(side == "Left", 1L, 2L),                         # 1 = Left, 2 = Right
    stratum = factor(stratum, levels = all_levels, ordered = TRUE)
  )

# ---- Sankey plot ----
ggplot(lodes, aes(x = x_num, stratum = stratum, alluvium = flow_id, y = value)) +
  # Colored connections (by LEFT cluster)
  geom_alluvium(aes(fill = cluster_left), width = 0.25, alpha = 0.85) +
  scale_fill_manual(values = cluster_cols, name = "Cluster") +
  
  ggnewscale::new_scale("fill") +
  
  # RIGHT boxes (white with black borders)
  geom_stratum(
    data = dplyr::filter(lodes, x_num == 2L),
    fill = "white", color = "black", width = 0.25
  ) +
  # LEFT boxes (colored by their own label using the derived palette)
  geom_stratum(
    data = dplyr::filter(lodes, x_num == 1L),
    aes(fill = after_stat(stratum)), color = "black", width = 0.25
  ) +
  
  # ---- LABELS ----
# Left labels (inside)
geom_text(
  data = dplyr::filter(lodes, x_num == 1L),
  stat = "stratum", aes(label = after_stat(stratum)),
  size = 3
) +
  # Right labels (outside with repel)
  geom_text_repel(
    data = dplyr::filter(lodes, x_num == 2L),
    stat = "stratum", aes(label = after_stat(stratum)),
    size = 3, nudge_x = 0.3, hjust = 0, direction = "y",
    segment.size = 0.2
  ) +
  
  scale_fill_manual(values = stratum_colors, guide = "none") +
  
  # Lock axis positions and labels
  scale_x_continuous(breaks = c(1, 2),
                     labels = c("Left", "Right"),
                     expand = c(0.15, 0.4)) +
  labs(y = "value", x = NULL) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )
