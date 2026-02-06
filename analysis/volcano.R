############################
# Volcano Plot
# Colored by metabolite class
# Labels = FDR < 0.05 only
############################
library(tidyverse)
library(ggrepel)

# Load data
df <- read.csv("df1.csv",
               sep=",", dec=".",
               stringsAsFactors=FALSE,
               check.names=FALSE)

df <- df[!(df$Primary_ID %in% c("mzmed","rtmed")), ]
df <- df[df$ClassID %in% c("insulin","control"), ]

# Matrix
mat <- df[, !(names(df) %in% c("Primary_ID","ClassID"))]
mat <- as.data.frame(lapply(mat, as.numeric))
mat <- as.matrix(mat)

classes <- df$ClassID
mat_log <- log2(mat + 1)

# Stats
results <- data.frame(feature = colnames(mat_log))

# Calculate log2 fold change
results$log2FC <- apply(mat_log, 2, function(x)
  mean(x[classes=="insulin"], na.rm=TRUE) -
    mean(x[classes=="control"], na.rm=TRUE)
)

# Compute p-values using t-test
results$pvalue <- apply(mat_log, 2, function(x){
  a <- x[classes=="insulin"]; b <- x[classes=="control"]
  if(length(a)<2 || length(b)<2) return(1)
  t.test(a,b)$p.value
})

# Adjust for FDR
results$FDR <- p.adjust(results$pvalue, "BH")
results$neglog10FDR <- -log10(results$FDR)

# Remove any rows with NA values in log2FC or neglog10FDR
results <- results[complete.cases(results$log2FC, results$neglog10FDR), ]

# Annotation lookup (manually provided metabolite details)
lookup <- data.frame(
  feature=c(
    "m0052","m0067","m0209","m0222","m0290","m0320","m0462","m0680",
    "m1330","m1382","m1530","m1705","m1795","m1859","m1876","m1927",
    "m2044","m2073","m2081","m2106","m2432","m2438","m2522","m2554",
    "m2613","m2655","m2659","m2684","m2690","m2777","m3117","m4180"),
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
    "Acylcarnitines","Acylcarnitines","Sphingolipids","Acylcarnitines","Sphingolipids"),
  stringsAsFactors=FALSE
)

results <- left_join(results, lookup, by="feature")

# Label ALL significant hits (FDR < 0.05)
labels_df <- results %>%
  filter(!is.na(name)) %>%
  arrange(FDR) %>%
  slice(1:10)

# Colours (updated for better readability)
class_cols <- c(
  "Acylcarnitines"        = "#6A5ACD",   # purple-blue
  "Amino acid / peptide"  = "#E64B35",   # red-orange
  "Fatty acids / lipids"  = "#1B9E77",   # dark green
  "Purine / nucleoside"   = "#4DBBD5",   # cyan
  "Sphingolipids"         = "#CC79A7",   # magenta
  "Other"                 = "grey40"
)

# Plot (with larger circles and adjusted text labels)
ggplot(results, aes(log2FC, neglog10FDR, colour=class)) +
  geom_point(size=3.5, alpha=0.9) +
  geom_text_repel(
    data=labels_df,
    aes(label=name),
    size=5,
    max.overlaps=Inf,
    show.legend = FALSE
  ) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  geom_vline(xintercept=c(-1,1), linetype="dotted") +
  scale_colour_manual(
    values = class_cols,
    na.value = "grey80",
    breaks = names(class_cols)
  ) +
  labs(
    x="log2 fold change (insulin / control)",
    y="-log10(FDR)",
    colour="Metabolite class"
  ) +
  theme_classic(base_size = 14) +
  theme(panel.background = element_blank(),   # Remove the background
        panel.grid.major = element_blank(),    # Remove major grid lines
        panel.grid.minor = element_blank())    # Remove minor grid lines
