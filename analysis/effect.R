############################################################
# FOREST PLOT — Cliff's delta with CI
############################################################
library(tidyverse)
library(effsize)

############################################################
# LOAD DATA 
############################################################
df <- read.csv(
  "df1.csv",
  sep=",", dec=".",
  stringsAsFactors=FALSE,
  check.names=FALSE
)

df <- df[df$ClassID %in% c("insulin","control"), ]

############################################################
# FEATURE → NAME + CLASS LOOKUP
############################################################
lookup <- tibble(
  feature=c(
    "m0052","m0067","m0209","m0222","m0290","m0320","m0462","m0680",
    "m1330","m1382","m1530","m1705","m1795","m1859","m1876","m1927",
    "m2044","m2073","m2081","m2106","m2432","m2438","m2522","m2554",
    "m2613","m2655","m2659","m2684","m2690","m2777","m3117","m4180"
  ),
  name=c(
    "1-Pyr2C","2-Hexenal","1-Methylnicotinamide","L-Valine",
    "MetO","Methylhypoxanthine","7-Methylguanine[M+H]","7-Methylguanine[M+Na]",
    "Lysylproline","C4:0-OH[M+H]","C6:0","C7:0","Ribothymidine","C4:0-OH[M+K]",
    "DiAcSpm","Inosine","ATRA","EPA[M+H]","EPA[M+H](13C1)","C7:0-OH",
    "C12:0[M+H]","C12:0[M+H](13C1)","C12:0-OH","C14:2-OH","C12:0-DC",
    "C15:0-OH","SAH","C14:0-OH[M+H]","C14:0-OH[M+H](13C1)","S1P","C16:1",
    "SM(34:2)"
  ),
  class=c(
    "Other","Fatty acids / lipids","Other","Amino acid / peptide",
    "Other","Purine / nucleoside","Purine / nucleoside","Purine / nucleoside",
    "Amino acid / peptide","Acylcarnitines","Acylcarnitines","Acylcarnitines",
    "Purine / nucleoside","Acylcarnitines","Other","Purine / nucleoside",
    "Fatty acids / lipids","Fatty acids / lipids","Fatty acids / lipids",
    "Acylcarnitines","Acylcarnitines","Acylcarnitines","Acylcarnitines",
    "Acylcarnitines","Acylcarnitines","Acylcarnitines","Other",
    "Acylcarnitines","Acylcarnitines","Sphingolipids","Acylcarnitines","Sphingolipids"
  )
)

############################################################
# CLASS COLOURS (background shading)
############################################################
class_cols <- c(
  "Acylcarnitines"        = "#6A5ACD",
  "Amino acid / peptide"  = "#E64B35",
  "Fatty acids / lipids"  = "#1B9E77",
  "Purine / nucleoside"   = "#4DBBD5",
  "Sphingolipids"         = "#CC79A7",
  "Other"                 = "grey40"
)

############################################################
# LONG FORMAT + log2 transform
############################################################
plot_df <- df %>%
  select(ClassID, all_of(lookup$feature)) %>%
  pivot_longer(-ClassID, names_to="feature", values_to="intensity") %>%
  mutate(intensity_log2 = log2(as.numeric(intensity) + 1)) %>%
  left_join(lookup, by="feature")

############################################################
# CLIFF'S DELTA + CI
############################################################
eff_df <- plot_df %>%
  group_by(name, class) %>%
  group_modify(~{
    cd <- cliff.delta(.x$intensity_log2 ~ .x$ClassID)
    tibble(
      delta   = as.numeric(cd$estimate),
      ci_low  = cd$conf.int[1],
      ci_high = cd$conf.int[2]
    )
  }) %>%
  ungroup()

############################################################
# ORDERING
############################################################
class_levels <- unique(lookup$class)

eff_df <- eff_df %>%
  mutate(
    class = factor(class, levels=class_levels),
    up = delta > 0
  ) %>%
  arrange(class, desc(up), desc(abs(delta))) %>%
  mutate(name = factor(name, levels=name))

############################################################
# BACKGROUND SHADING BLOCKS
############################################################
blocks <- eff_df %>%
  mutate(y = as.numeric(name)) %>%
  group_by(class) %>%
  summarise(
    ymin=min(y)-0.5,
    ymax=max(y)+0.5
  )

############################################################
# PLOT
############################################################
ggplot(eff_df, aes(x=delta, y=name)) +
  
  # class shading
  geom_rect(data=blocks,
            aes(xmin=-Inf, xmax=Inf, ymin=ymin, ymax=ymax, fill=class),
            inherit.aes=FALSE,
            alpha=0.25) +
  
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  
  # CI bars (colour-mapped)
  geom_errorbarh(aes(xmin=ci_low, xmax=ci_high, colour=up),
                 height=0.15,
                 linewidth=1.3) +
  
  # points
  geom_point(aes(colour=up), size=3) +
  
  scale_fill_manual(values=class_cols, guide="none") +
  
  scale_colour_manual(values=c(
    "FALSE"="#FF7F0E",   # inuslin higher
    "TRUE"="#1f78b4"   # control higher
  ), guide="none") +
  
  scale_x_continuous(limits=c(-0.6,0.6)) +
  
  labs(
    x="Effect size (Cliff’s delta)",
    y=NULL
  ) +
  
  theme_classic(base_size=14)
