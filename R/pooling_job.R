cacheEnv <- new.env()

cacheEnv$output_semi_std_full     <- NULL
cacheEnv$output_semi_overlap_full <- NULL
cacheEnv$output_non_std_full      <- NULL
cacheEnv$output_non_overlap_full  <- NULL


cacheEnv$groupings_std     <- c("n", "p", "d", "lambda", "sign")
cacheEnv$groupings_overlap <- c("n", "d", "s", "c")


#' Full output generation
#'
#' @param param_mode "semi" or "non"
#' @param group_mode "std" or "overlap"
#' @return Full output
#' @export
get_output_full <- function(param_mode, group_mode){

  path_output <- sprintf("output_%s_%s", param_mode, group_mode)
  output_list <- list.files(path = path_output, pattern = "output", full.names = TRUE,
                            recursive = TRUE)

  # Rare event that there are empty files
  info <- file.info(output_list)
  output_list <- rownames(info[info$size != 0, ])


  # Granular output
  output_full <- sapply(output_list, outputParser, group_mode, simplify = FALSE)
  output_full <- dplyr::bind_rows(output_full)

  indx.processed <- as.integer(stringi::stri_extract_first_regex(basename(output_list), "[0-9]+"))
  output_full    <- cbind(graphNum = indx.processed, output_full)

  col_types  <- sapply(output_full, typeof)
  col_double <- colnames(output_full)[col_types == "double"]

  output_full <- dplyr::mutate_at(output_full, col_double, dplyr::funs(round(.,2)))


  # Remove duplicates
  output_full <- unique(output_full)

  # Update global variables
  assign(sprintf("output_%s_%s_full", param_mode, group_mode), output_full, envir = cacheEnv)

  return(output_full)

}




#' All observations with meeting columns
#'
#' @param param_mode "semi" or "non"
#' @param group_mode "std" or "overlap"
#' @return Non-summarized output (with meeting columns)
#' @export
#' @importFrom magrittr %>%
get_output_meet_full <- function(param_mode, group_mode){

  output_full <- eval(parse(text = sprintf("output_%s_%s_full", param_mode, group_mode)),
                      envir = cacheEnv)

  groupings <- eval(parse(text = sprintf("groupings_%s", group_mode)),
                      envir = cacheEnv)

  if (group_mode == "std") {
    stats      <- c("overlap", "post_var", "MSE")
    stats_meet <- c("overlap", "var", "MSE")
  } else {
    stats      <- c("overlap", "MSE")
    stats_meet <- stats
  }

  # Select variables
  vars_meet <- c(groupings, stats)
  output_full_meet <- dplyr::select(output_full, vars_meet)

  # Aesthetic column name changes
  colnames(output_full_meet) <- replace(colnames(output_full_meet),
                                        which(colnames(output_full_meet) %in% stats),
                                        stats_meet)

  # Order
  output_full_meet <- output_full_meet[ , c(groupings, stats_meet)]

  # Write to disk
  utils::write.table(output_full_meet, file = sprintf("mcmc_peixoto_%s_%s.txt", param_mode, group_mode),
                     sep = "\t", row.names = FALSE)

  return(output_full_meet)

}



#' Summary of indices processed
#'
#' @param param_mode "semi" or "non"
#' @param group_mode "std" or "overlap"
#' @return Tally of indices procced
#' @export
#' @importFrom magrittr %>%
get_samps_processed <- function(param_mode, group_mode) {

  output_full <- eval(parse(text = sprintf("output_%s_%s_full", param_mode, group_mode)),
                      envir = cacheEnv)

  if (group_mode == "std"){

    samps_processed <- output_full %>%
      group_by(n, d, p, lambda, sign) %>%
      mutate(samps = n(),
             indx_beg = min(graphNum)) %>%
      select(n, d, p, lambda, sign, samps, indx_beg) %>%
      distinct()

    parameters <- dplyr::select(parameters, -c(a, b, c))
    samps_processed <- dplyr::right_join(samps_processed, parameters,
                                         by = c("n", "d", "p", "lambda", "sign"))


  } else {

    samps_processed <- output_full %>%
      group_by(n, d, s, c) %>%
      mutate(samps = n(),
             indx_beg = min(graphNum)) %>%
      select(n, d, s, c, samps, indx_beg) %>%
      distinct()

    parameters_overlap_temp <- parameters_overlap %>%
      dplyr::mutate_at("c", dplyr::funs(round(.,2)))

    samps_processed <- dplyr::right_join(samps_processed, parameters_overlap_temp,
                                         by = c("n", "d", "s", "c"))

  }

  return(samps_processed)

}



#' Summary of full output
#'
#' @param param_mode "semi" or "non"
#' @param group_mode "std" or "overlap"
#' @return Full summarized output
#' @export
#' @importFrom magrittr %>%
get_output_sum <- function(param_mode, group_mode){


  output_full <- eval(parse(text = sprintf("output_%s_%s_full", param_mode, group_mode)),
                      envir = cacheEnv)

  col_types <- sapply(output_full, typeof)
  col_num   <- colnames(output_full)[col_types != "character"]

  if (group_mode == "std") {

    output_full_sum <- output_full %>%
      dplyr::select(c(col_num, "sign")) %>%
      dplyr::group_by(n, p, d, lambda, sign) %>%
      dplyr::summarize_all(mean, na.rm = TRUE)

  } else {

    output_full_sum <- output_full %>%
      dplyr::select(col_num) %>%
      dplyr::group_by(n, d, s, c) %>%
      dplyr::summarize_all(mean, na.rm = TRUE)

  }

  return(output_sum_full)

}


#' Summary of full output for meeting
#'
#' @param output_full_meet full results (from cluster)
#' @param param_mode "semi" or "non"
#' @param group_mode "std" or "overlap"
#' @return Full summarized output
#' @export
#' @importFrom magrittr %>%
get_output_sum_meet <- function(output_full_meet, param_mode, group_mode){

  col_types <- sapply(output_full_meet, typeof)
  col_num   <- colnames(output_full_meet)[col_types != "character"]

  if (group_mode == "std") {

    output_sum_meet <- output_full_meet %>%
      dplyr::select(c(col_num, "sign")) %>%
      dplyr::group_by(n, p, d, lambda, sign) %>%
      dplyr::summarize_all(mean, na.rm = TRUE)

  } else {

    output_sum_meet <- output_full_meet %>%
      dplyr::select(col_num) %>%
      dplyr::group_by(n, d, s, c) %>%
      dplyr::summarize_all(mean, na.rm = TRUE)

  }

  # Write to disk
  utils::write.table(output_sum_meet, file = sprintf("mcmc_peixoto_%s_%s_sum.txt", param_mode, group_mode),
                     sep = "\t", row.names = FALSE)

  return(output_sum_meet)

}
