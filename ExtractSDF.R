# Extractng from SDF

sdfname <- "UniqueBasetest1.sdf"
outsdfname <- "uniqueBasetest6.sdf"
numfname <- "computed.txt"

nummatr <- read.csv(numfname, header=FALSE)
numlist <- list(nummatr[,1])[[1]]

f <- file(outsdfname, "w")

lines <- readLines(sdfname)
nlines <- length(lines)
imol <- 0
istart <- 1
for (i in 1:nlines) {
  if (any(grep("\\$\\$\\$\\$", lines[i]))) {
    imol <- imol + 1
    if (imol %in% numlist) {
      for (j in istart:i)
        cat(lines[j], file=f, sep="\n")
    }
    istart <- i + 1
  }
}

close(f)
