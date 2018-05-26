# Adding of properties to SDF file

# Input sdf file without properties
sdfname <- "uniquebasetest6.sdf"

# Input text file with properties
propsfname <- "uniquebasetest5.txt"

# Output sdf file with inserted properties
outsdfname <- "uniquebasetest7.sdf"

# Read in file with properties
props <- read.table(propsfname, header=TRUE)

# Read sdf file without properties
lines <- readLines(sdfname)
nlines <- length(lines)

# Write new sdf file with inserted ptoperties
f <- file(outsdfname, "w")
imol <- 0
for (i in 1:nlines) {
  cat(lines[i], file = f, sep = "\n")
  if (any(grep("\\$\\$\\$\\$", lines[i+1]))) {
    imol <- imol + 1
    for (prop in colnames(props)) {
      cat(sprintf(">  <%s>\n", prop), file=f)
      cat((props[imol, prop]), file=f, sep="\n")
      cat("\n", file=f)
    }
  }
}
close(f)
