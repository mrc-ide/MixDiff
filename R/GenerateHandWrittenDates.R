#######################################
### functions to generate a plot of a handwritten date ###
#######################################

#' Plots a handwritten date
#'
#' @param date A date
#' @param d an object of class mnist containing a database of handwritten digits (see \code{rmnist::load_mnist}).
#' @return NULL
#' @import rmnist
#' @export
#' @examples
#' plot_handwritten_date(Sys.Date()) # print today's date
plot_handwritten_date <- function(date = as.Date("01/01/2017", format = "%d/%m/%Y"), d = load_mnist(download_if_missing = TRUE)) # date has to be a date
{
  # convert date to string with "/" separator
  date <- as.character(date, format = "%d/%m/%Y")

  # create a "dot" image to separate day, month and year
  width <- 3
  dot_mat <- matrix(0, 28, 28)
  for (i in 20 + (1:width))
  {
    dot_mat[i, (28 - width) / 2 + (1:width)] <- 1
  }

  separator <- t(dot_mat)
  class(separator) <- c("mnist_digit", "matrix")
  attr(separator, "label") <- 0 ### this S3 class forces to have an integer label...
  attr(separator, "label") <- "mnist_digit"
  attr(separator, "data") <- "matrix"

  # create a panel of 10 images for dd/mm/yyyy
  par(mfrow = c(1, 10), mar = c(0, 0, 0, 0))

  # generate a random sample of handwritten images for each character in the date string, forcing repeated numbers to have the same handwriting every time they appear
  date_idx <- strsplit(date, NULL)[[1]]
  date_idx_unique <- unique(date_idx)
  chosen_image_unique <- lapply(date_idx_unique, function(k) {
    if (k == "/") chosen_image <- separator else chosen_image <- d[[sample(which(d$label %in% as.numeric(k)), 1)]]
    return(chosen_image)
  })

  # do the plotting
  for (idx in date_idx)
  {
    plot(chosen_image_unique[[match(idx, date_idx_unique)]], box = FALSE)
  }
}
