# Detection of duplicate structures

source("cinf-sdf.R")
source("cinf-mol.R")
source("cinf-ptable.R")
source("cinf-isomorph.R")

input_file <- "UniqueBaseTest.sdf"
Unique_list_file <- "UniqueBaseTest.txt"

ct_list <- list()
atomlist_list <- list()
natoms_list <- list()
inv_list <- list()
unique_list <- list()
nunique <- 0
inv_list <- list()

mdb <- read_sdf(input_file)
nmol <- length(mdb)
for (imol in 1:nmol) {
  cat("Molecule=", imol)
  mol <- mdb[[imol]]
  natoms <- length(mol$atoms)
  atomlist <- list()
  for (atom in mol$atoms) {
    atomlist <- c(atomlist, atom$el)
  }
  ct <- mol_get_ct(mol, bond_orders = 1)
  
  # form matrix gmatr for computing invariant
  gmatr <- ct
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
    numel <- PT$NumEl[[atom$el]]
    gmatr[iatom,iatom] <- numel + 0.1 * atom$nh
  }
  for (iatom1 in 1:(natoms-1))
    for (iatom2 in (iatom1+1):natoms) {
      gmatr[iatom1,iatom2] <- -gmatr[iatom1,iatom2]
      gmatr[iatom2,iatom1] <- gmatr[iatom2,iatom2]
    }
  
  # compute invariant
  inv <- log(det(gmatr))
  cat(" invariant=", inv)
  
  # search in the list of previously found invariants
  nfound <- 0
  if (nunique > 0)
    for (i in 1:nunique)
      if (abs(inv - inv_list[[i]]) < 0.00001) {
        if (natoms != natoms_list[[i]]) next
        isomorph <- find_substr_isomorph(atomlist, ct, atomlist_list[[i]], ct_list[[i]])
        perm <- length(isomorph)
        if (perm > 0) {
          nfound <- unique_list[[i]]
          break
        }
      }
  
  # if not found than expand the list of unique structures
  if (nfound == 0) {
    nunique <- nunique + 1
    unique_list[[nunique]] <- imol
    inv_list[[nunique]] <- inv
    ct_list[[nunique]] <- ct
    atomlist_list[[nunique]] <- atomlist
    natoms_list[[nunique]] <- natoms
  } else {
    cat(" identical to ", nfound)
  }
  
  # search in the list of previously found invariants
  cat("\n")
}
cat ("Number of unique structures =", nunique)

# Write numbers to file
f <- file(Unique_list_file, "w")
for (num in unique_list) cat(num, file=f, sep="\n")
close(f)

