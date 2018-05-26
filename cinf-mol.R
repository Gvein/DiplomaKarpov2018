# Different procedures concerning molecules

# Returns array of chemical element labels
mol_get_chelabs <- function(mol) {
  atoms <- mol[["atoms"]]
  natoms <- length(atoms)
  chelabs <- array(dim=natoms)
  for (iatom in 1:natoms) {
    chelabs[iatom] <- atoms[[iatom]]$el
  }
  chelabs
}

# Returns connection table for a molecule
mol_get_ct <- function(mol, bond_orders=0) {
  atoms <- mol[["atoms"]]
  natoms <- length(atoms)
  bonds <- mol[["bonds"]]
  nbonds <- length(bonds)
  ct <- matrix(0, nrow=natoms, ncol=natoms)
  for (ibond in 1:nbonds) {
    bond <- bonds[[ibond]]
	if (bond_orders) {
      ct[bond$at1, bond$at2] <- bond$bo
      ct[bond$at2, bond$at1] <- bond$bo
	} else {
      ct[bond$at1, bond$at2] <- 1
      ct[bond$at2, bond$at1] <- 1
	}
  }
  ct
}

# Extracts substructure from a molecule
substruct <- function(mol, oldnum) {
  na_small <- length(oldnum)
  na_big <- length(mol$atoms)
  nb_big <- length(mol$bonds)
  newnum <- integer(na_big)
  for (i in 1:na_small) newnum[oldnum[i]] <- i
  newmol <- list()
  
  # Copy atoms
  atoms <- list()
  for (ia in 1:na_small) atoms[[ia]] <- mol$atoms[[oldnum[ia]]]
  
  # Copy bonds
  bonds <- list()
  ib1 <- 0
  for (ib in 1:nb_big) {
    bond <- mol$bonds[[ib]]
	if ((newnum[bond$at1] > 0) && (newnum[bond$at2] > 0)) {
	  bond$at1 <- newnum[bond$at1]
	  bond$at2 <- newnum[bond$at2]
	  ib1 <- ib1 + 1
	  bonds[[ib1]] <- bond
	}
  }
  
  # Assemble substructure
  newmol$atoms <- atoms
  newmol$bonds <- bonds
  newmol  
}

# Extracts substructure from a molecule using mask
substr_mask <- function(mol, mask) {
  mask_size <- length(mask)
  new_size <- 0
  for (i in 1:mask_size) {
    if (mask[i]) {
	  new_size <- new_size + 1
	}
  }
  oldnum <- integer(new_size)
  j <- 0
  for (i in 1:mask_size) {
    if (mask[i]) {
	  j <- j + 1
	  oldnum[j] <- i
	}
  }
  substruct(mol, oldnum)
}

# Computes distance matrix
calc_distance_matrix <- function(connTable) {
  n <- dim(connTable)[1]
  distMatrix <- matrix(0, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
	  if (connTable[i,j]) {
	    distMatrix[i,j] <- 1
	  } else {
	    distMatrix[i,j] <- 0
	  }
	}
  }
  repeat {
    nc <- 0
	for (i in 1:(n-1)) {
	  for (j in (i+1):n) {
	    if (!distMatrix[i,j]) {
		  md <- 10000
		  for (k in 1:n) {
		    if (distMatrix[i,k] * distMatrix[j,k]) {
			  s <- distMatrix[i,k] + distMatrix[j,k]
			  if (s < md) {md <- s}
			}
		  }
		  if (md < 10000) {
		    distMatrix[i,j] <- md
			distMatrix[j,i] <- md
			nc <- nc + 1
		  }
		}
	  }
	}
    if (!nc) break
  }
  distMatrix
}

