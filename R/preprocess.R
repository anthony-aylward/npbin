#!/usr/bin/Rscript
#===============================================================================
# preprocess_counts.R
#===============================================================================

# Preprocess allele count data




# Function definitions ---------------------------------------------------------

#' Get the indices of required fields in a header.
#'
#' Returns a list providing the index of each required field.
#'
#' \code{get_header_indices} accepts a character vector, which could be the
#' entire header in one element separated by whitespace or a vector with
#' several elements, one for each field. It uses \code{which} to determine
#' the inidces of the required fields.
#'
#' @param header A character vector
#' @return A list mapping required fields to their respective indices
get_header_indices <- function(header) {
  if (length(header) == 1) {
    split_header <- strsplit(header, "\\s+")[[1]]
  } else {
    split_header <- header
  }
  list(
    chr_index = which(split_header == "chr"),
    location_index = which(split_header == "location"),
    total_index = which(split_header == "m"),
    ref_count_index = which(split_header == "xm")
  )
}

#' Convert a logical value to a character.
#'
#' Converts TRUE to "P", FALSE to "M".
#'
#' This is a utility for construction of the "winning.chip" field required by
#' NPBin functions.
#'
#' @param logical_val A logical value
#' @return "P" if input was TRUE, "M" if it was FALSE.
logical_to_winning_chip <- function(logical_val) {
  if (logical_val) {'P'} else {'M'}
}

#' Preprocess allele count data
#'
#' Perform preprocessing steps to prepare data for NPBin analysis
#'
#' This function extracts the required fields from the input data and returns
#' a data frame formatted for NPBin analysis.
#'
#' @param data_frame A data frame containing at least the four required fields
#' @param chr_index Index of the chromosome field
#' @param chr_index Index of the location field
#' @param chr_index Index of the total read count field
#' @param chr_index Index of the reference allele count field
#' @return A data frame 
#' @export
preprocess_counts <- function(
  data_frame,
  chr_index = 1,
  location_index = 2,
  total_index = 3,
  ref_count_index = 4
) {
  m <- data_frame[[total_index]]
  xm <- data_frame[[ref_count_index]]
  winning_chip_logical <- xm < m / 2
  data.frame(
    chr = data_frame[[chr_index]],
    location = data_frame[[location_index]],
    m = m,
    xm = xm,
    winning.chip = sapply(winning_chip_logical, logical_to_winning_chip)
  )
}

#' Preprocess allele count data, handling a header if necessary.
#'
#' This function wraps \code{preprocess_counts} with header handling.
#'
#' Checks if the input data includes a header by checking column types of
#' the input data frame. If a header is present, the indices of the required
#' fields are noted and used with \code{preprocess_counts}
#'
#' @param data_frame A data frame containing required fields for NPBin analysis
#' @return A data frame formatted for NPBin analysis
#' @export
preprocess_file <- function(file_path) {
  header_indices <- list(
    chr_index = 1,
    location_index = 2,
    total_index = 3,
    ref_count_index = 4
  )
  data_frame <- read.table(file_path, stringsAsFactors=FALSE)
  
  if (!(TRUE %in% lapply(data_frame, is.integer))) {
    data_frame <- read.table(args[[1]], header=TRUE)
    header_indices <- get_header_indices(names(data_frame))
  }
  
  preprocess_counts(
    data_frame,
    chr_index = header_indices[["chr_index"]],
    location_index = header_indices[["location_index"]],
    total_index = header_indices[["total_index"]],
    ref_count_index = header_indices[["ref_count_index"]]
  )
}