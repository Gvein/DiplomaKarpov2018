source("cinf-sdf.R")
source("cinf-mol.R")
source("cinf-isomorph.R")

fnameimp <- "str92dextracted.sdf"
SymmList <- list()
NoSymmList <- list()
SymmNumberList <- list()
UnSymmNumberList <- list()
nsymm <- 0
nnosymm <- 0
mdb <- read_sdf(fnameimp)
nmol <- length(mdb)
for (imol in 1:nmol) {
  mol <- mdb[[imol]]
  natoms <- length(mol$atoms)
  atomlist <- list()
  for (atom in mol$atoms) {
    atomlist <- c(atomlist, atom$el)
  }
  ct <- mol_get_ct(mol, bond_orders = 1)
  isomorph <- find_substr_isomorph(atomlist, ct, atomlist, ct)
  nperm <- length(isomorph)
  atomn1 <- 0
  atomn2 <- 0
  for (iatom in 1:length(mol$atoms)) {
    atom <- mol$atoms[[iatom]]
    if (atom$el == "N" && !atomn1) {
      atomn1 <- iatom
      next
    }
    if (atom$el == "N" && !atomn2) {
      atomn2 <- iatom
      break
    }
  }
  cat("Molecule=", imol, " atoms=", natoms, " Permutations=", nperm)
  cat(" N1=", atomn1, " N2=", atomn2)
  symm <- FALSE
  for (perm in isomorph) {
    if (perm[atomn1] == atomn2) 
    {
      cat(" Symmerty")
      symm <- TRUE
      break
    }
  }
  if (symm) {
    nsymm <- nsymm + 1
    SymmNumberList[[nsymm]] <- imol
    SymmList[[nsymm]] <- mol
  } else {
    nnosymm <- nnosymm + 1
    UnSymmNumberList[[nnosymm]] <- imol
    
    NoSymmList[[nnosymm]] <- mol
  }
  cat("\n")
}
#cat(SymmList)
write_sdf(SymmList, "NoMetalSymmBase.sdf")
write_sdf(NoSymmList, "NoMetalUnSymmBase.sdf")
#write_sdf(SymmNumberList, "SymmNumberBase.sdf")
#write_sdf(UnSymmNumberList, "UnSymmNumberBase.sdf")
#file.create("SymmNumberBase.txt")
#file.create("UnSymmNumberBase.txt")
#for (i in 1:length(SymmNumberList)){
#  write.table(c(SymmNumberList[[i]], ""),
#              "SymmNumberBase.txt", row.names = F, col.names = F, quote = F, append = T)
#}
#for (i in 1:length(UnSymmNumberList)){
#  write.table(c(UnSymmNumberList[[i]], ""),
#              "UnSymmNumberBase.txt", row.names = F, col.names = F, quote = F, append = T)
#}
f <- file("SymmNumberBase92dex.txt", "w")
for (i in SymmNumberList) cat(i, file=f, sep="\n")
close(f)
f <- file("UnSymmNumberBase92dex.txt", "w")
for (i in UnSymmNumberList) cat(i, file=f, sep="\n")
close(f)
