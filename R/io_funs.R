
#' Parse single-row peixoto output
#'
#' @param filePath path to output file
#' @return A one-row data frame
#' @export
outputParser <- function(filePath, group_mode) {

  output <- utils::read.csv(filePath, header = TRUE, stringsAsFactors = FALSE)

  return(output)

}


#' Produce (single-row) algorithm output as txt file
#'
#' @param file_in file name of graph
#' @param pyJobOut raw output from python
#' @export
outputGen <- function(file_in, pyJobOut) {

  # Build path and file name
  params <- unlist(regmatches(file_in, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",file_in)))
  sign <- regmatches(file_in, regexpr("(pos|neg)", file_in))

  file_out <- gsub("graph", "output", file_in)
  file_out <- gsub("csv", "txt", file_out)
  dir_out  <- paste("output", params["n"], params["d"], params["p"], params["lambda"], sign, sep = "/")
  path_out <- paste(dir_out, file_out, sep = "/")


  # Parse python output
  stats     <- strsplit(pyJobOut, "\a")[[1]][2]
  stats     <- unlist(regmatches(stats, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",stats)))

  results <- data.frame(n = params["n"], p = params["p"],  d = params["d"],
                        lambda = params["lambda"], overlap = stats[1], time_tot = stats[2],
                        time_minSBM = stats[3], time_eq = stats[4], time_harvest = stats[5],
                        countZero = stats[6], entropy = stats[7])

  utils::write.csv(results, file = path_out, row.names = FALSE)

}
