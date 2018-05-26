# Reads and writes files in Sybyl mol2-file

# Read molecular database from file filename
read_mol2 <- function(filename) {

  # read all lines into memory
  lines <- readLines(filename)
  nlines <- length(lines)

  moldbase <- list()

  # scan the MOLECULE blocks
  mol_block_start <- integer()
  imol <- 0
  for (i in 1:nlines) {
    if (any(grep("@<TRIPOS>MOLECULE", lines[i]))) {
	imol <- imol + 1
      mol_block_start[imol] <- i
    }
  }
  nmols <- imol

  # scan the ATOM blocks
  atom_block_start <- integer()
  imol <- 0
  for (i in 1:nlines) {
    if (any(grep("@<TRIPOS>ATOM", lines[i]))) {
	imol <- imol + 1
      atom_block_start[imol] <- i
    }
  }

  # scan the BOND blocks
  bond_block_start <- integer()
  imol <- 0
  for (i in 1:nlines) {
    if (any(grep("@<TRIPOS>BOND", lines[i]))) {
	imol <- imol + 1
      bond_block_start[imol] <- i
    }
  }

  # loop over molecules
  for (imol in 1:nmols) {
    hline <- lines[mol_block_start[imol]+2]
    res <- gregexpr("[0-9]+", hline)
    starts <- res[[1]]
    natoms <- as.integer(substr(hline, starts[1], starts[2] - 1))
    nbonds <- as.integer(substr(hline, starts[2], starts[3] - 1))
    molecule <- list()

    # loop over atoms
    atoms <- list()
    for (iatom in 1:natoms) {
      aline <- lines[atom_block_start[imol] + iatom]
      res <- gregexpr("[^ \t]+", aline)
      starts <- res[[1]]
      x_str <- substr(aline, starts[3], starts[4]-1)
      x <- as.double(x_str)
      y_str <- substr(aline, starts[4], starts[5]-1)
      y <- as.double(y_str)
      z_str <- substr(aline, starts[5], starts[6]-1)
      z <- as.double(z_str)
      syb <- substr(aline, starts[6], starts[7]-1)
      syb <- sub('[[:space:]]+$', '', syb)
      el <- sub('[.].*', '', syb)
      pch_str <- substring(aline, starts[9])
      pch <- as.double(pch_str)
      if (el == "CL") el <- "Cl"
      if (el == "BR") el <- "Br"
      atom <- list()
      atom[["el"]] <- el
      atom[["syb"]] <- syb
      if (syb == "C.cat" || syb == "N.pl3")
	    atom[["ch"]] <- 1
	  else
	    atom[["ch"]] <- 0
      atom[["nh"]] <- 0
      atom[["pch"]] <- pch
      atom[["x"]] <- x
      atom[["y"]] <- y
      atom[["z"]] <- z
      atoms[[iatom]] <- atom
    }
    molecule[["atoms"]] <- atoms

    # loop over bonds
    bonds <- list()
    if (nbonds > 0) {
      for (ibond in 1:nbonds) {
        bline <- lines[bond_block_start[imol] + ibond]
        res <- gregexpr("[^ \t]+", bline)
        starts <- res[[1]]
        at1_str <- substr(bline, starts[2], starts[3]-1)
        at1 <- as.integer(at1_str)
        at2_str <- substr(bline, starts[3], starts[4]-1)
        at2 <- as.integer(at2_str)
		if (is.na(starts[5])) {
          bo_str <- substring(bline, starts[4])
		} else {
          bo_str <- substring(bline, starts[4], starts[5]-1)	
		}
        if (any(grep("a", bo_str))) {
          if (any(grep("ar", bo_str))) {
            bo <- 4
          } else {
            bo <- 1
          }
		} else if (any(grep("un", bo_str))) {
		  bo <- 8
        } else {
          bo <- as.integer(bo_str)
        }
        bond <- list()
        bond[["at1"]] <- at1
        bond[["at2"]] <- at2
        bond[["bo"]] <- bo
		bond[["syb"]] <- bo_str
        bonds[[ibond]] <- bond
      }
    }
    molecule[["bonds"]] <- bonds
    moldbase[[imol]] <- molecule
  }
  moldbase
}

# Writes molecular database mdb to file fname in Sybyl mol2 format
write_mol2 <- function(mdb, fname) {
  of <- file(fname, "w")
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    mol <- mdb[[imol]]
	cat("@<TRIPOS>MOLECULE\n", file=of)
	cat("*****\n", file=of)
	natoms <- length(mol$atoms)
	nbonds <- length(mol$bonds)
	cat(sprintf("%5d%5d%5d%5d%5d\n", natoms, nbonds, 0, 0, 0), file=of)
	cat("SMALL\n", file=of)
	cat("GASTEIGER\n", file=of)
	cat("\n", file=of)
	cat("@<TRIPOS>ATOM\n", file=of)
	for (iatom in 1:natoms) {
	  atom <- mol$atoms[[iatom]]
	  aname <- sprintf("%s%d", atom$el, iatom)
	  cat(sprintf("%7d %-10s %10.6f %10.6f %10.6f %-10s 1 LIG     %10.6f\n", 
	    iatom, aname, atom$x, atom$y, atom$z, atom$syb, atom$pch), file=of)
	}
	cat("@<TRIPOS>BOND\n", file=of)
	for (ibond in 1:nbonds) {
	  bond <- mol$bonds[[ibond]]
	  if (bond$bo == 4) {bname <- "ar"}
	  cat(sprintf("  %5d %5d %5d %-s\n", ibond, bond$at1, bond$at2, bond$syb), file=of)
	}
  }
  close(of)
}

