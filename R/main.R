
#' Encode DNA Sequences
#'
#' @description Converts DNA sequences into categorical
#'   nucleotide features for machine learning prediction.
#' @param dna_strings A character vector of equal-length
#'   DNA sequences containing A, T, C and G.
#' @return A data frame with factor columns representing
#'   nucleotide positions.
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "GAACA"))

dna_encoding <- function(dna_strings){
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}


#' Predict m6A Sites for Multiple Samples
#'
#' @description Uses a trained random forest model to
#'   predict m6A modification probabilities and statuses
#'   for multiple samples.
#' @param ml_fit A trained random forest model.
#' @param feature_df A data frame containing gc_content,
#'   RNA_type, RNA_region, exon_length,
#'   distance_to_junction, evolutionary_conservation
#'   and DNA_5mer.
#' @param positive_threshold Numeric probability
#'   threshold for classifying positive predictions.
#'   The default is 0.5.
#' @return A data frame containing the original features
#'   and two additional columns: predicted_m6A_prob
#'   and predicted_m6A_status. Probabilities strictly
#'   greater than the threshold are classified Positive.
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction"))
#' prediction_multiple(ml_fit, head(example_df, 3))
#' @import randomForest
#' @importFrom stats predict
#' @export

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
                  "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) # Check errors if incorrect column names of input data.frame
  feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  infered_prob <- predict(ml_fit, newdata = cbind(feature_df, dna_encoding(feature_df$DNA_5mer)), type = "prob")
  feature_df$predicted_m6A_prob <- infered_prob[,"Positive"]
  feature_df$predicted_m6A_status <- "Negative"
  feature_df$predicted_m6A_status[feature_df$predicted_m6A_prob > positive_threshold] <- "Positive"
  return(feature_df) # Return a data.frame with predicted m6A prob and predicted m6A status
}


#' Predict m6A Modification for a Single Sample
#'
#' @description Predicts the m6A modification probability
#'   and status of one sample using a trained random
#'   forest model.
#' @param ml_fit A trained random forest model.
#' @param gc_content Numeric GC content.
#' @param RNA_type Character value describing RNA type.
#' @param RNA_region Character value describing RNA region.
#' @param exon_length Numeric log2-transformed exon length.
#' @param distance_to_junction Numeric log2-transformed
#'   distance to the nearest splice junction.
#' @param evolutionary_conservation Numeric conservation
#'   score.
#' @param DNA_5mer A five-nucleotide DNA sequence.
#' @param positive_threshold Numeric probability threshold.
#'   The default is 0.5.
#' @return A named character vector containing
#'   predicted_m6A_prob and predicted_m6A_status.
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds", package = "m6APrediction"))
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction"))
#' prediction_single(
#'   ml_fit = ml_fit,
#'   gc_content = example_df$gc_content[1],
#'   RNA_type = example_df$RNA_type[1],
#'   RNA_region = example_df$RNA_region[1],
#'   exon_length = example_df$exon_length[1],
#'   distance_to_junction =
#'     example_df$distance_to_junction[1],
#'   evolutionary_conservation =
#'     example_df$evolutionary_conservation[1],
#'   DNA_5mer = example_df$DNA_5mer[1]
#' )
#' @export

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  pred_outcome <- prediction_multiple(ml_fit = ml_fit,
                                      data.frame(gc_content = gc_content,
                                                 RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
                                                 RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
                                                 exon_length = exon_length,
                                                 distance_to_junction = distance_to_junction,
                                                 evolutionary_conservation = evolutionary_conservation,
                                                 DNA_5mer = DNA_5mer),
                                      positive_threshold = positive_threshold)
  returned_vector <- c(pred_outcome$predicted_m6A_prob, pred_outcome$predicted_m6A_status)
  names(returned_vector) <- c("predicted_m6A_prob", "predicted_m6A_status")
  return(returned_vector) # Return a named vector with predicted m6A prob and predicted m6A status
}


