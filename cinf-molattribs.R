# Computes additional molecular attributes for molecules

add_mol_attribs <- function(moldbase, base=TRUE) {

  # Adds to atoms:
  #   vd_ - vertex degree
  #   va_ - valence
  #   pi_ - the number of pi-electrons
  #   ar_ - aromaticity
  #   ne_ - vector of neighbours
  #   bo_ - vector of bond orders with neighbours  
  add_base_mol_attribs <- function(mdb) {
    for (imol in 1:nmol) {
	  mol <- mdb[[imol]]
	  natoms <- length(mol$atoms)
      nbonds <- length(mol$bonds)
      if (natoms > 0) {
	    for (iatom in 1:natoms) {
	      atom <- mol$atoms[[iatom]]
	      atom$vd_ <- 0
	      atom$va_ <- atom$nh + abs(atom$ch)
		  atom$pi_ <- 0
		  atom$ar_ <- FALSE
		  atom$ne_ <- integer()
		  atom$bo_ <- integer()
		  mdb[[imol]]$atoms[[iatom]] <- atom
	    }
	  }
	  if (nbonds > 0) {
	    for (ibond in 1:nbonds) {
	      bond <- mdb[[imol]]$bonds[[ibond]]
		  atom1 <- mdb[[imol]]$atoms[[bond$at1]]
		  atom2 <- mdb[[imol]]$atoms[[bond$at2]]
		  atom1$vd_ <- atom1$vd_ + 1
		  atom2$vd_ <- atom2$vd_ + 1
		  if (bond$bo < 4) {
		    atom1$va_ <- atom1$va_ + bond$bo
		    atom2$va_ <- atom2$va_ + bond$bo
		    atom1$pi_ <- atom1$pi_ + bond$bo - 1
		    atom2$pi_ <- atom2$pi_ + bond$bo - 1
		  } else if (bond$bo == 4) {
		    atom1$va_ <- atom1$va_ + 1.5
		    atom2$va_ <- atom2$va_ + 1.5
		    atom1$pi_ <- 1
		    atom2$pi_ <- 1
			atom1$ar_ <- TRUE
			atom2$ar_ <- TRUE
		  }
		  atom1$ne_[atom1$vd_] <- bond$at2
		  atom2$ne_[atom2$vd_] <- bond$at1
		  atom1$bo_[atom1$vd_] <- bond$bo
		  atom2$bo_[atom2$vd_] <- bond$bo
		  mdb[[imol]]$atoms[[bond$at1]] <- atom1
		  mdb[[imol]]$atoms[[bond$at2]] <- atom2
	    }
		for (iatom in 1:natoms) {
          atom <<- mol$atoms[[iatom]]
		  if (atom$va_ == 4.5) atom$va_ <- 4
		}
	  }
	}
	mdb
  }

  nmol <- length(moldbase)
  moldbase1 <- moldbase
  if (base) {
    moldbase1 <- add_base_mol_attribs(moldbase1)
  }
  moldbase1
}
