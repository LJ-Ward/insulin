############################################################
# Box + Jitter
############################################################

library(tidyverse)
library(ggpubr)

############################
# Load data
############################
df <- read.csv(
  "df1.csv",
  sep=",", dec=".",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

df <- df[!(df$Primary_ID %in% c("mzmed","rtmed")), ]
df <- df[df$ClassID %in% c("insulin","control"), ]

############################
# Features + names
############################
features <- c("m1876","m2659","m1530","m1927","m2073")

names_lookup <- c(
  m1876="DiAcSpm",
  m1530="C6:0",
  m2659="SAH",
  m1927="Inosine",
  m2073="EPA"
)

############################
# Long format + log2
############################
plot_df <- df %>%
  select(ClassID, all_of(features)) %>%
  pivot_longer(-ClassID,
               names_to="feature",
               values_to="intensity") %>%
  mutate(
    intensity_log2 = log2(as.numeric(intensity) + 1),
    metabolite = names_lookup[feature],
    ClassID = factor(ClassID, levels=c("control","insulin"))
  ) %>%
  filter(is.finite(intensity_log2))  

############################
# Stats
############################
stats_df <- plot_df %>%
  group_by(metabolite) %>%
  summarise(
    log2FC = mean(intensity_log2[ClassID=="insulin"], na.rm=TRUE) -
      mean(intensity_log2[ClassID=="control"], na.rm=TRUE),
    p = wilcox.test(intensity_log2 ~ ClassID)$p.value,
    .groups="drop"
  ) %>%
  mutate(
    label = paste0(
      metabolite,
      "\nlog2FC = ", sprintf("%.2f", log2FC),
      "\np = ", format(p, digits=2, scientific=TRUE)
    )
  )

lab_lookup <- setNames(stats_df$label, stats_df$metabolite)

############################
# Manual order
############################
plot_df$metabolite <- factor(
  plot_df$metabolite,
  levels=c("DiAcSpm","C6:0","SAH","Inosine","EPA")
)

############################
# Plot
############################
ggplot(plot_df, aes(ClassID, intensity_log2, colour=ClassID)) +
  
  geom_boxplot(width=0.7, outlier.shape=NA, fill="white") +
  geom_jitter(width=0.12, alpha=0.25, size=1.8) +
  
  facet_wrap(
    ~ metabolite,
    nrow=1,
    labeller = labeller(metabolite = lab_lookup)
  ) +
  
  scale_colour_manual(values=c(
    control="#1f78b4",
    insulin="#ff7f00"
  )) +
  
  labs(x=NULL, y="Log2 intensity") +
  
  theme_classic(base_size=14) +
  theme(
    legend.position="none",
    strip.text = element_text(
      face="bold",
      size=13,
      lineheight=1.1
    )
  )
