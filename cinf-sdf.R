# Reads sdf-file

source("cinf-hydrogens.R")
source("cinf-molattribs.R")

read_sdf <- function (filename, coord=TRUE, prop=TRUE, 
                      to_numeric=FALSE, delete_expl_hydr=FALSE) {

  trim <- function(s) {
    if (substr(s, 2, 2) == " ") {
	  substr(s, 1, 1)
	} else {
	  substr(s, 2, 2) <- tolower(substr(s, 2, 2))
	  s
	}
  }
  
  lines <- readLines(filename)
  nlines <- length(lines)
  moldbase <- list()
  imol <- 0
  last_line <- -1
  for (i in 1:nlines) {
    if (any(grep("\\$\\$\\$\\$", lines[i]))) {
	  imol <- imol + 1
        first_line <- last_line + 2
	  last_line <- i - 1
	  expl_hydr <- FALSE
	  
      # Read line with the number of atoms and bonds
	  natoms <- as.integer(substr(lines[first_line+3], 1, 3))
	  nbonds <- as.integer(substr(lines[first_line+3], 4, 6))
	  molecule <- list()
	  
	  # Read the name of the molecule
	  name <- lines[first_line]
	  name <- gsub(" $", "", name)
	  if (name != "") molecule$name <- name

	  # Read atoms
	  first_atom_block <- first_line + 4
	  last_atom_block <- first_atom_block + natoms - 1
	  atoms <- list()
	  for (j in first_atom_block:last_atom_block) {
		atom <- list()
		atom[["el"]] <- trim(substr(lines[j], 32, 33))
		if (atom[["el"]] == "H") expl_hydr <- TRUE
		atom[["nh"]] <- as.integer(substr(lines[j], 40, 42))
		ch <- as.integer(substr(lines[j], 37, 39))
		if (ch > 0) ch <- 4 - ch
		atom[["ch"]] <- ch
		if (coord) {
		  atom[["x"]] <- as.numeric(substr(lines[j], 1, 10))
		  atom[["y"]] <- as.numeric(substr(lines[j], 11, 20))
		  atom[["z"]] <- as.numeric(substr(lines[j], 21, 30))
		}
		atoms[[j - first_atom_block + 1]] <- atom
	  }
	  molecule[["atoms"]] <- atoms
	  
	  # Read bonds
	  if (nbonds > 0) {
	    first_bond_block <- last_atom_block + 1
	    last_bond_block <- first_bond_block + nbonds - 1
	    bonds <- list()
	    for (j in first_bond_block:last_bond_block) {
	      bond <- list()
		  bond[["at1"]] <- as.integer(substr(lines[j], 1, 3))
		  bond[["at2"]] <- as.integer(substr(lines[j], 4, 6))
		  bond[["bo"]]  <- as.integer(substr(lines[j], 7, 9))
		  bonds[[j - first_bond_block + 1]] <- bond
	    }
	    molecule[["bonds"]] <- bonds
	  } else {
	    last_bond_block <- last_atom_block
	  }

        # Read additional information
        j = last_bond_block + 1
        while (!any(grep("M  END", lines[j]))) {
          if (any(grep("^V  ", lines[j]))) {
            iat <- as.integer(substr(lines[j], 4, 6))
            lab <- substr(lines[j], 8, 8)
            atom <- atoms[[iat]]
            atom[["lab"]] <- lab
            atoms[[iat]] <- atom
          }
          j <- j + 1
        }
	  molecule[["atoms"]] <- atoms
	  
	  # Read properties
	  if (prop) {
		first_prop_block <- j + 1
		props <- list()
		if (first_prop_block < last_line) {
		  for (j in seq(first_prop_block, last_line - 2, 3)) {
		    r <- regexpr("<.+>", lines[j])
			first_propname <- r + 1
			last_propname <- r + attr(r, "match.length") - 2
			propname <- substr(lines[j], first_propname, last_propname)
			propval <- lines[j+1]
			if (to_numeric) propval <- as.numeric(propval)
			if (!is.na(propval)) props[[propname]] <- propval
		  }
		}
		molecule[["props"]] <- props
	  }
	  
	  if (delete_expl_hydr) molecule <- del_expl_hydr(molecule)
      molecule <- add_impl_hydr(molecule)	  
	  
	  moldbase[[imol]] <- molecule
	}
  }
  add_mol_attribs(moldbase)
}

# Writes sdf-file
write_sdf <- function(mdb, filename) {
  f <- file(filename, "w")
  nmols <- length(mdb)
  for (imol in 1:nmols) {
    # Write header of molecule
	if (!is.null(mdb[[imol]]$name)) {
	  name <- mdb[[imol]]$name
	} else {
	  name <- ""
	}
	cat(name, file=f)
	cat("\n", file=f)
    cat("   - CINF - \n", "\n", file=f)
    mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
    nbonds <- length(mol$bonds)
    nprops <- length(mol$props)
    # Count labels
    nlab <- 0
    atoms <- mol[["atoms"]]
    for (iatom in 1:natoms) {
      if (!is.null(atoms[[iatom]]$lab)) nlab <- nlab + 1
    }
    cat(sprintf("%3d%3d  0  0  0  0  0  0  0  0%3d V2000", natoms, nbonds, nlab+1), file=f, sep="\n")

    # Write atoms
    for (iatom in 1:natoms) {
      atom <- mol$atoms[[iatom]]
      x <- atom$x
      y <- atom$y
      z <- atom$z
      el <- atom$el     
      ch <- atom$ch
      if (ch!=0) ch <- 4 - atom$ch
      nh <- atom$nh
      cat(sprintf("%10.4f%10.4f%10.4f %-2s  0%3d%3d  0  0  0", x, y, z, el, ch, nh), file=f, sep="\n")
    }

    # Write bonds
    if (nbonds > 0) {
      for (ibond in 1:nbonds) {
        bond <- mol$bonds[[ibond]]
        at1 <- bond$at1
        at2 <- bond$at2
        bo <- bond$bo
        cat(sprintf("%3d%3d%3d  0  0  0", at1, at2, bo), file=f, sep="\n")
      }
    }

    # Write labels
    if (nlab > 0) {
      for (iatom in 1:natoms) {
        if (!is.null(atoms[[iatom]]$lab)) {
          cat(sprintf("V  %3d %s", iatom, atoms[[iatom]]$lab), file=f, sep="\n")
        }
      }
    }

    cat("M  END", file=f, sep="\n")

    # Write properties
    if (nprops > 0) {
      prop_names <- names(mol$props)
      for (iprop in 1:nprops) {
        prop = mol$props[[iprop]]
        cat(sprintf(">  <%s> (%d)", prop_names[iprop], imol), file=f, sep="\n")
        if (is.numeric(prop)) {
          cat(sprintf("%g", prop), "", file=f, sep="\n")
        } else {
          cat(sprintf("%s", prop), "", file=f, sep="\n")
        }
      }
    }
    cat("$$$$", file=f, sep="\n")
  }
  close(f)
}

