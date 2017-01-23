# Load the MNIST digit recognition dataset into R
# http://yann.lecun.com/exdb/mnist/
# assume you have all 4 files and gunzip'd them
# creates train$n, train$x, train$y  and test$n, test$x, test$y
# e.g. train$x is a 60000 x 784 matrix, each row is one digit (28x28)
# call:  show_digit(train$x[5,])   to see a digit.
# brendan o'connor - gist.github.com/39760 - anyall.org

load_mnist <- function() {
  load_image_file <- function(filename) {
    ret = list()
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    ret$n = readBin(f,'integer',n=1,size=4,endian='big')
    nrow = readBin(f,'integer',n=1,size=4,endian='big')
    ncol = readBin(f,'integer',n=1,size=4,endian='big')
    x = readBin(f,'integer',n=ret$n*nrow*ncol,size=1,signed=F)
    ret$x = matrix(x, ncol=nrow*ncol, byrow=T)
    close(f)
    ret
  }
  load_label_file <- function(filename) {
    f = file(filename,'rb')
    readBin(f,'integer',n=1,size=4,endian='big')
    n = readBin(f,'integer',n=1,size=4,endian='big')
    y = readBin(f,'integer',n=n,size=1,signed=F)
    close(f)
    y
  }
  train <<- load_image_file('mnist/train-images-idx3-ubyte')
  test <<- load_image_file('mnist/t10k-images-idx3-ubyte')
  
  train$y <<- load_label_file('mnist/train-labels-idx1-ubyte')
  test$y <<- load_label_file('mnist/t10k-labels-idx1-ubyte')  
}

show_digit <- function(arr784, col=gray(12:1/12), ...) {
  image(matrix(arr784, nrow=28)[,28:1], col=col, ...)
}

print_handwritten_date <- function(date=as.Date("01/01/2017", format="%d/%m/%Y")) # date has to be a date
{
  # convert date to string with "/" separator
  date <- as.character(date, format="%d/%m/%Y")
  
  # create a "dot" image to separate day, month and year
  width <- 3
  dot_vect <- rep(0, 28^2)
  for(i in 20+(1:width))
  {
    dot_vect[28*i+(28-width)/2+(1:width)] <- 1
  }
  separator <- dot_vect
  
  # create a panel of 10 images for dd/mm/yyyy
  par(mfrow=c(1, 10), mar=c(0, 0, 0, 0))
  
  # generate a random sample of handwritten images for each character in the date string, forcing repeated numbers to have the same handwriting every time they appear 
  date_idx <- strsplit(date, NULL)[[1]]
  date_idx_unique <- unique(date_idx)
  chosen_image_unique <- lapply(date_idx_unique, function(k) {
    if (k == "/") chosen_image <- separator else chosen_image <- train$x[sample(which(train$y %in% k), 1), ]; return(chosen_image)
  })
  
  # do the plotting
  for(idx in date_idx)
  {
    show_digit(chosen_image_unique[[match(idx, date_idx_unique)]], axes=FALSE)
  }
  
}

### To print a handwritten version of a date, use the following
# load_mnist()
# date <- Sys.Date() # this has to be in Date format
# print_handwritten_date(date)