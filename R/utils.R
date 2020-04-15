# internal function checking if it is possible to use fread.
# fread is the fastest way to load data but if dataset are too big
# fread can generate error (because of lack of memory)
# this function compare the size of a file to the max allowed size
# in Mb (default value is 800Mb) allowed to use the fread function
can_use_fread <- function (filePath, maxSize = 800) {
  fileSize <- file.info(filePath)$size/1024^2
  if (fileSize <= maxSize) {
    return(TRUE)
  }
  return(FALSE)
}

