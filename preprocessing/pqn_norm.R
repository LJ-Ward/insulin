###############################
##  PQN normalisation script ##
###############################

## --------- SETTINGS (EDIT THESE) ---------

# Path to your xcms output file (xlsx or csv)
input_file  <- "df1.xlsx"

# Sheet name or index (only used for Excel files)
sheet       <- 1

# Output file for PQN-normalised data
output_file <- "df1_norm.csv"

# Number of sample-level metadata rows at the TOP
# (for your example: "group" and "set" -> 2 rows)
n_sample_meta_rows <- 2

# Name of the FIRST sample-intensity column
# (for your example: "QC01_RPpos.mzML")
first_sample_col_name <- "QC01.mzML"

# Row label (first column) used to define the reference group
# and the value in that row to use as reference
# -> use all samples where set == "training"
reference_meta_row_label <- "set"
reference_meta_value     <- "training"

# Small value added to intensities before PQN (normally 0)
epsilon <- 0

## --------- END OF SETTINGS ---------


#############################################
##           Load / Install Packages       ##
#############################################

if (grepl("\\.xlsx?$", input_file, ignore.case = TRUE)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    install.packages("readxl")
  }
  library(readxl)
}


#############################################
##                 Read Data               ##
#############################################

if (grepl("\\.xlsx?$", input_file, ignore.case = TRUE)) {
  raw_data <- readxl::read_excel(input_file, sheet = sheet)
} else {
  raw_data <- utils::read.csv(input_file, check.names = FALSE)
}

raw_data <- as.data.frame(raw_data, check.names = FALSE)

# Locate first sample column
first_sample_col <- match(first_sample_col_name, colnames(raw_data))
if (is.na(first_sample_col)) {
  stop("Could not find column named '", first_sample_col_name,
       "' in the input header.")
}

# Columns before intensities = feature metadata
feature_meta_cols <- seq_len(first_sample_col - 1)

# Split metadata rows (top block) and feature rows
if (n_sample_meta_rows > 0) {
  sample_meta  <- raw_data[seq_len(n_sample_meta_rows), , drop = FALSE]
  feature_data <- raw_data[-seq_len(n_sample_meta_rows), , drop = FALSE]
} else {
  sample_meta  <- NULL
  feature_data <- raw_data
}

feature_meta <- feature_data[, feature_meta_cols, drop = FALSE]
intensities  <- feature_data[, first_sample_col:ncol(feature_data), drop = FALSE]


#############################################
##        Convert Intensities to Numeric   ##
#############################################

intensities_num <- as.data.frame(
  lapply(intensities, function(x) as.numeric(as.character(x))),
  check.names = FALSE
)

if (!is.null(epsilon) && epsilon != 0) {
  intensities_num <- intensities_num + epsilon
}


#############################################
##    Build Reference Matrix (set=training)##
#############################################

if (is.null(sample_meta)) {
  stop("No sample_meta; cannot select reference by metadata.")
}

first_meta_col <- feature_meta_cols[1]

meta_row_idx <- which(sample_meta[, first_meta_col] == reference_meta_row_label)

if (length(meta_row_idx) == 0) {
  stop("Could not find metadata row labeled '",
       reference_meta_row_label, "'.")
}
if (length(meta_row_idx) > 1) {
  stop("More than one metadata row labeled '",
       reference_meta_row_label, "'. Please check file.")
}

meta_vec <- as.character(
  sample_meta[meta_row_idx, first_sample_col:ncol(sample_meta)]
)

keep <- meta_vec == reference_meta_value

if (!any(keep)) {
  stop("No samples with ", reference_meta_row_label, " = '",
       reference_meta_value, "' were found.")
}

message("Using ", sum(keep),
        " samples with ", reference_meta_row_label, " = '",
        reference_meta_value, "' as PQN reference.")

ref_mat <- intensities_num[, keep, drop = FALSE]

# Reference spectrum = median intensity across reference samples
reference_spectrum <- apply(ref_mat, 1, median, na.rm = TRUE)
reference_spectrum[reference_spectrum == 0] <- NA


#############################################
##      Compute PQN Dilution Factors       ##
#############################################

quotients <- sweep(intensities_num, 1, reference_spectrum, "/")
dilution_factors <- apply(quotients, 2, median, na.rm = TRUE)


#############################################
##        Apply PQN Normalisation          ##
#############################################

pqn_intensities <- sweep(intensities_num, 2, dilution_factors, "/")


#############################################
##         Reassemble & Write Output       ##
#############################################

pqn_intensities_out <- as.data.frame(
  lapply(pqn_intensities, function(x) signif(x, 6)),
  check.names = FALSE
)

# Add PQN_factor row for reference
pqn_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(raw_data)))
colnames(pqn_row) <- colnames(raw_data)
pqn_row[1, feature_meta_cols[1]] <- "PQN_factor"
pqn_row[1, first_sample_col:ncol(raw_data)] <-
  dilution_factors[colnames(intensities_num)]

final_output <- rbind(
  sample_meta,
  pqn_row,
  cbind(feature_meta, pqn_intensities_out)
)

utils::write.csv(final_output, output_file, row.names = FALSE)
message("Done! PQN-normalised data written to: ", output_file)
