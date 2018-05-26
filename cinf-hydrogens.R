# Handling hydrogens in molecules

source("cinf-ptable.R")

# Deletes explicit hydrogens
# if add_impl==TRUE, than the appropriate implicit hydrogens are added upon each deletion
del_expl_hydr <- function(mol, add_impl=TRUE, recalc_attribs=FALSE) {
  to_delete <- integer()
  nh <- 0
  natoms <- length(mol$atoms)
  nbonds <- length(mol$bonds)
  for (iatom in 1:natoms) {
    if (mol$atoms[[iatom]]$el == "H") {
	  nh <- nh + 1
	  to_delete[nh] <- iatom
	}
  }
  if (nh > 0) {
    ind <- 1:natoms
	natoms1 <- natoms
	for (i in nh:1) {
	  hatom <- to_delete[i]
	  if (hatom == 1) {
	    ind <- ind[2:natoms1]
	  } else if (hatom == natoms1) {
	    ind <- ind[1:(natoms1-1)]
	  } else {
	    ind <- ind[c(1:(hatom-1), (hatom+1):natoms1)]
	  }
	  natoms1 <- natoms1 - 1
	}
	mol$atoms <- mol$atoms[ind]
	new_ia <- integer()
	for (i in 1:natoms) new_ia[i] <- 0
	for (i in 1:natoms1) new_ia[ind[i]] <- i
	natoms <- natoms - nh
	if (recalc_attribs) for (iatom in 1:natoms) {
	  atom <- mol$atoms[[iatom]]
	  if (atom$vd_ > 0) {
	    ind <- integer()
		lind <- 0
	    for (i in 1:atom$vd_) {
		  atom$ne_[i] <- new_ia[atom$ne_[i]]
		  if (atom$ne_[i] > 0) {
		    lind <- lind + 1
			ind[lind] = i
		  }
		}
		if (lind != atom$vd_) {
		  atom$ne_ <- atom$ne_[ind]
		  atom$vd_ <- lind
        }		
	  }
	  mol$atoms[[iatom]] <- atom
	}
	
	if (nbonds > 0 && natoms > 0) {
	  ind <- logical()
	  for (ibond in 1:nbonds) {
	    bond <- mol$bonds[[ibond]]
	    ind[ibond] <- TRUE
	    bond$at1 <- new_ia[bond$at1]
	    bond$at2 <- new_ia[bond$at2]
	    if (bond$at1 == 0) {
		  ind[ibond] <- FALSE
		  if (add_impl) mol$atoms[[bond$at2]]$nh <- mol$atoms[[bond$at2]]$nh + 1 
	    }
	    if (bond$at2 == 0) {
	      ind[ibond] <- FALSE
		  if (add_impl) mol$atoms[[bond$at1]]$nh <- mol$atoms[[bond$at1]]$nh + 1 
	    }
	    mol$bonds[[ibond]] <- bond
	  }
	  mol$bonds <- mol$bonds[ind]
	}
  }
  mol
}

# Adds implicit hydrogens
add_impl_hydr <- function(mol) {
  natoms <- length(mol$atoms)
  nbonds <- length(mol$bonds)
  val <- integer()
  for (iatom in 1:natoms) {
    val[iatom] <- abs(mol$atoms[[iatom]]$ch)
  }
  if (nbonds > 0) for (ibond in 1:nbonds) {
    bond <- mol$bonds[[ibond]]
	at1 <- bond$at1
	at2 <- bond$at2
	bo <- bond$bo
	val[at1] <- val[at1] + bo
	val[at2] <- val[at2] + bo
  }  
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
	nh <- 0
	if (val[iatom] <= PT$ComVal[[atom$el]]) {
	  if (atom$ch <= 0) {
	    nh <- PT$ComVal[[atom$el]] - val[iatom]
	  } else {
	    nh <- PT$MaxVal[[atom$el]] - val[iatom]
	  }
	} else {
	  if (val[iatom] <= PT$MaxVal[[atom$el]]) {
	    nh <- (val[iatom] - PT$ComVal[[atom$el]]) %/% 2
	  } else {
	    nh <- 0
	  }
	}
	atom$nh <- nh
    mol$atoms[[iatom]] <- atom	
  }
  mol
}

