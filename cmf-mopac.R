# Interface to the MOPAC program

source("cinf-mol2.R")

# Returns charge of molecule
get_mol_charge <- function(mol) {
  natoms <- length(mol$atoms)
  charge <- 0
  for (i in 1:natoms) charge <- charge + mol$atoms[[i]]$pch
  charge
}

# Writes MOPAC inpout file with Cartesian coordinates
write_mopac_input_file <- function(mol, fname) {
  charge <- get_mol_charge(mol)
  natoms <- length(mol$atoms)
  of <- file(fname, "w")
  cat(sprintf("AUX LARGE CHARGE=%g SPARKLE SUPER PRECISE PM7\n", round(charge)), file=of)
  cat("\n", file=of)
  cat("\n", file=of)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
    cat(sprintf("%-2s%9.5f 1%9.5f 1%9.5f 1\n", atom$el, atom$x, atom$y, atom$z), file=of)
  }
  close(of)
}

run_mopac <- function(fname)
{
  shell(paste("mopac2016.exe", "aaa.mop"))
}

read_mopac_aux_file <- function(fname) {

  # proceed to next line
  next_line <- function() {
    iline <<- iline + 1
	aline <<- aux_lines[iline]
  }
  
  # read ncount string tokens
  read_strings <- function(count) {
	res <- character(count)
    icount <- 0
	while (icount < count) {
	  next_line()
	  tokens <- strsplit(aline, " +", perl=TRUE)[[1]]
	  ntokens <- length(tokens)
	  for (itoken in 2:ntokens) {
	    token <- tokens[itoken]
		res[icount+itoken-1] <- token
	  }
	  icount <- icount + ntokens - 1
	}
	res
  }

  # read ncount integers
  read_integers <- function(count) {
    string_tokens <- read_strings(count)
	as.integer(string_tokens)
  }
  
  # read ncount real numbers
  read_reals <- function(count) {
    string_tokens <- read_strings(count)
	as.numeric(string_tokens)
  }
  
  # read coordinates of atoms in molecule
  read_coordinates <- function(natoms) {
    coord <- matrix(nrow=natoms, ncol=3)
	for (iatom in 1:natoms) {
	  next_line()
	  coord[iatom,1] <- as.numeric(substr(aline, 1, 10)) 
	  coord[iatom,2] <- as.numeric(substr(aline, 11, 20))
	  coord[iatom,3] <- as.numeric(substr(aline, 21, 30))
	}
	coord
  }
  
  # read lower half triangle
  read_lower_half_triangle <- function(size) {
    count <- size * (size + 1) / 2
    matr <- matrix(nrow=size, ncol=size)
	next_line()
	values <- read_reals(count)
	# copy array values to lower triangle
	ival <- 0
	for (i in 1:size)
	  for (j in 1:i) {
	    ival <- ival + 1
		matr[i,j] <- values[ival]
		matr[j,i] <- matr[i,j]
	  }
	matr
  }
  
  # read rectangular matrix
  read_matrix <- function(nrows, ncols) {
    count <- nrows * ncols
	matr <- matrix(nrow=nrows, ncol=ncols)
	values <- read_reals(count)
	ival <- 0
	for (i in 1:nrows)
	  for (j in 1:ncols) {
	    ival <- ival + 1
		matr[i,j] <- values[ival]
	  }
	matr
  }
  
  aux_lines <- readLines(fname)
  nlines <- length(aux_lines)
  iline <- 0
  while (iline < nlines) {
    next_line()
    
    if (substr(aline,2,6)=="ERROR") {
      return(list())
    }
	
	# Read scalar values
	r <-  regexpr("([A-Z_]+):[A-Z/]+=([0-9+-.D]+)$", aline, perl=TRUE)
	if (r > 0) {
	  start1 <- attr(r, "capture.start")[1]
	  stop1 <- start1 + attr(r, "capture.length")[1] - 1
	  name <- substr(aline, start1, stop1)
	  start2 <- attr(r, "capture.start")[2]
	  stop2 <- start2 + attr(r, "capture.length")[2] - 1
      valstr1 <- substr(aline, start2, stop2)
	  valstr2 <- sub("D", "E", valstr1)
	  value <- as.numeric(valstr2)
	  if (name == "HEAT_OF_FORMATION") {
	    heat_of_formation <- value
		next
	  }
	  if (name == "TOTAL_ENERGY") {
	    total_energy <- value
		next
	  }
	  if (name == "DIPOLE") {
	    dipole <- value
		next
	  }
	}
	
	# Read vectors and matrices
	r <-  regexpr("([A-Z_:]+).([0-9]+).= *$", aline, perl=TRUE)
	if (r > 0) {
	  start1 <- attr(r, "capture.start")[1]
	  stop1 <- start1 + attr(r, "capture.length")[1] - 1
	  name <- substr(aline, start1, stop1)
	  start2 <- attr(r, "capture.start")[2]
	  stop2 <- start2 + attr(r, "capture.length")[2] - 1
	  count <- as.integer(substr(aline, start2, stop2))
	  if (name == "ATOM_EL") {
	    natoms <- count
		atom_el <- read_strings(count)
		next
	  }
	  if (name == "ATOM_CORE") {
	    atom_core <- read_integers(count)
		next
	  }
	  if (name == "ATOM_X:ANGSTROMS") {
	    atom_x <- read_coordinates(natoms)
		next
	  }
	  if (name == "AO_ATOMINDEX") {
	    naorbitals <- count
	    ao_atomindex <- read_integers(count)
		next
	  }
	  if (name == "ATOM_SYMTYPE") {
	    atom_symtype <- read_strings(count)
		next
	  }
	  if (name == "AO_ZETA") {
	    ao_zeta <- read_reals(count)
		next
	  }
	  if (name == "ATOM_PQN") {
	    atom_pqn <- read_integers(count)
		next
	  }
	  if (name == "ATOM_X_OPT:ANGSTROMS") {
	    atom_x_opt <- read_coordinates(natoms)
		next
	  }
	  if (name == "ATOM_CHARGES") {
	    atom_charges <- read_reals(count)
		next
	  }
	  if (name == "OVERLAP_MATRIX") {
	    overlap_matrix <- read_lower_half_triangle(naorbitals)
		next
	  }
	  if (name == "EIGENVECTORS") {
	    nmorbitals <- count / naorbitals
        eigenvectors <- read_matrix(nmorbitals, naorbitals)
        next 		
	  }
	  if (name == "TOTAL_DENSITY_MATRIX") {
	    total_density_matrix <- read_lower_half_triangle(naorbitals)
		next
	  }
	  if (name == "EIGENVALUES") {
	    eigenvalues <- read_reals(count)
		next
	  }
	  if (name == "MOLECULAR_ORBITAL_OCCUPANCIES") {
	    molecular_orbital_occupancies <- read_reals(count)
		next
	  }
	}
  }
  list(
    natoms = natoms,
	naorbitals = naorbitals,
	nmorbitals = nmorbitals,
    atom_el = atom_el, 
	atom_core = atom_core, 
	atom_x = atom_x, 
	ao_atomindex = ao_atomindex,
	atom_symtype = atom_symtype,
	ao_zeta = ao_zeta,
	atom_pqn = atom_pqn,
	atom_x_opt = atom_x_opt,
	atom_charges = atom_charges,
	overlap_matrix = overlap_matrix,
	eigenvectors = eigenvectors,
	total_density_matrix = total_density_matrix,
	eigenvalues = eigenvalues,
	molecular_orbital_occupancies = molecular_orbital_occupancies,
	heat_of_formation = heat_of_formation,
	total_energy = total_energy,
	dipole = dipole
  )
}

read_mopac_out_file <- function(fname, mol) {
 
  # proceed to next line
  next_line <- function() {
    iline <<- iline + 1
	aline <<- arc_lines[iline]
  }
  
  get_num_heavy_atoms <- function(mol) {
	iheavy <- 0	
	for (atom in mol$atoms) {
	  if (atom$el != "H") iheavy <- iheavy + 1
	}
	return(iheavy)
  }
  
  arc_lines <- readLines(fname)
  nlines <- length(arc_lines)
  iline <- 0
  
  # Skip to superdelocalizabilities
  while (iline < nlines) {
    next_line()
    r <- regexpr("           SUPERDELOCALIZABILITIES", aline)
	if (r > 0) break
  }
  
  # Read several scalar values
  next_line()
  next_line()
  mulliken_electronegativity <- as.numeric(substr(aline, 34, 44))
  next_line()
  parr_pople_absolute_hardness <- as.numeric(substr(aline, 34, 44))
  next_line()
  schuurmann_mo_shift_alpha <- as.numeric(substr(aline, 34, 44))
  next_line()
  next_line()
  ehomo <- as.numeric(substr(aline, 34, 44))
  next_line()
  elumo <- as.numeric(substr(aline, 34, 44))
  for (i in 1:4) next_line()
  
  # Read arrays with superdelocalizabilities
  num_heavy_atoms <- get_num_heavy_atoms(mol)
  heavy_atom_index <- integer(num_heavy_atoms)
  Dn <- double(num_heavy_atoms)
  De <- double(num_heavy_atoms)
  qZ <- double(num_heavy_atoms)
  piS <- double(num_heavy_atoms)
  chomo <- double(num_heavy_atoms)
  clumo <- double(num_heavy_atoms)
  # Read Dn(r), De(r) and q(r)-Z(r)
  for (i in 1:num_heavy_atoms) {
    next_line()
    heavy_atom_index[i] <- as.integer(substr(aline, 5, 7))
	Dn[i] <- as.numeric(substr(aline, 11, 20))
	De[i] <- as.numeric(substr(aline, 24, 33))
	qZ[i] <- as.numeric(substr(aline, 37, 46))
  }
  for (i in 1:5) next_line()
  # Read piS
  for (i in 1:num_heavy_atoms) {
    next_line()
	piS[i] <- as.numeric(substr(aline, 11, 20))
  }
  for (i in 1:5) next_line()
  # Read homo and lumo
  for (i in 1:num_heavy_atoms) {
    next_line()
	chomo[i] <- as.numeric(substr(aline, 24, 33))
	clumo[i] <- as.numeric(substr(aline, 37, 46))
  }
  
  list(
    num_heavy_atoms = num_heavy_atoms,
	heavy_atom_index = heavy_atom_index,
	Dn = Dn,
	De = De,
	qZ = qZ,
	piS = piS,
	chomo = chomo,
	clumo = clumo,
	mulliken_electronegativity = mulliken_electronegativity,
	parr_pople_absolute_hardness = parr_pople_absolute_hardness,
	schuurmann_mo_shift_alpha = schuurmann_mo_shift_alpha,
	ehomo = ehomo,
	elumo = elumo
  )
}

calc_fukui <- function(mol, mopac_aux) {

  get_num_heavy_atoms <- function(mol) {
	iheavy <- 0	
	for (atom in mol$atoms) {
	  if (atom$el != "H") iheavy <- iheavy + 1
	}
	return(iheavy)
  }

  # Initialization
  num_heavy_atoms <- get_num_heavy_atoms(mol)
  Dn <- double(num_heavy_atoms)
  De <- double(num_heavy_atoms)
  piS <- double(num_heavy_atoms)
  chomo <- double(num_heavy_atoms)
  clumo <- double(num_heavy_atoms)
  
  # Find heavy atoms
  natoms <- length(mol$atoms)
  heavy_atom_index <- integer(num_heavy_atoms)
  num_heavy_atom <- integer(natoms)
  iheavy <- 0
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
	if (atom$el != "H") {
	  iheavy <- iheavy + 1
	  heavy_atom_index[iheavy] <- iatom
	  num_heavy_atom[iatom] <- iheavy
	}
  }
  
  # Compute global values
  nelectrons <- sum(mopac_aux$molecular_orbital_occupancies)
  ihomo <- nelectrons / 2
  ilumo <- ihomo + 1
  ehomo <- mopac_aux$eigenvalues[ihomo]
  elumo <- mopac_aux$eigenvalues[ilumo]
  alpha <- (ehomo + elumo) / 2
  mulliken_electronegativity <- -alpha
  parr_pople_absolute_hardness <- (ehomo - elumo) / 2
  schuurmann_mo_shift_alpha <- alpha
  
  # Compute local values
  # Find indexes for atomic orbitals
  first_ao <- integer(natoms)
  last_ao <- integer(natoms)
  for (iaorb in 1:mopac_aux$naorbitals) {
    iatom <- mopac_aux$ao_atomindex[iaorb]
    if (first_ao[iatom] == 0) first_ao[iatom] <- iaorb
    last_ao[iatom] <- iaorb
  }
  # Find indexes for atomic orbitals of heavy atoms
  first_ao_heavy <- integer(num_heavy_atoms)
  last_ao_heavy <- integer(num_heavy_atoms)
  for (iaorb in 1:mopac_aux$naorbitals) {
    iatom <- mopac_aux$ao_atomindex[iaorb]
	ihatom <- num_heavy_atom[iatom]
	if (ihatom > 0) {
	  if (first_ao_heavy[ihatom] == 0) first_ao_heavy[ihatom] <- iaorb
	  last_ao_heavy[ihatom] <- iaorb
	}
  }
  # Compute Fukui indexes
  for (ihatom in 1:num_heavy_atoms) {
	  if (first_ao_heavy[ihatom]) {
    for (iaorb in first_ao_heavy[ihatom]:last_ao_heavy[ihatom]) {
      # Compute chomo and clumo
	  chomo[ihatom] <- chomo[ihatom] + 2 * mopac_aux$eigenvectors[ihomo, iaorb]^2
	  clumo[ihatom] <- clumo[ihatom] + 2 * mopac_aux$eigenvectors[ilumo, iaorb]^2
	  # Compute Dn - nucleophilic delocalizabilities
	  for (ivac in ilumo:mopac_aux$nmorbitals) {
	    Dn[ihatom] <- Dn[ihatom] + 
		2 * mopac_aux$eigenvectors[ivac, iaorb]^2 / (alpha - mopac_aux$eigenvalues[ivac])
	  }
	  # Compute De - electrophilic delocalizabilities
	  for (iocc in 1:ihomo) {
	    De[ihatom] <- De[ihatom] + 
		2 * mopac_aux$eigenvectors[iocc, iaorb]^2 / (mopac_aux$eigenvalues[iocc] - alpha)
	  }
	  # Compute piS - self-polarizabilities
	  for (iocc in 1:ihomo) {
	    for (ivac in ilumo:mopac_aux$nmorbitals) {
		  piS[ihatom] <- piS[ihatom] +
		  4 * mopac_aux$eigenvectors[ivac, iaorb]^2 * mopac_aux$eigenvectors[iocc, iaorb]^2 /  
		  (mopac_aux$eigenvalues[ivac] - mopac_aux$eigenvalues[iocc])
		}
	  }
    }
  }
  }

  list(
    num_heavy_atoms = num_heavy_atoms,
	num_heavy_atom = num_heavy_atom,
	heavy_atom_index = heavy_atom_index,
	first_ao = first_ao,
	last_ao = last_ao,
	first_ao_heavy = first_ao_heavy,
	last_ao_heavy = last_ao_heavy,
	ihomo = ihomo,
	Dn = Dn,
	De = De,
	piS = piS,
	chomo = chomo,
	clumo = clumo,
	mulliken_electronegativity = mulliken_electronegativity,
	parr_pople_absolute_hardness = parr_pople_absolute_hardness,
	schuurmann_mo_shift_alpha = schuurmann_mo_shift_alpha,
	ehomo = ehomo,
	elumo = elumo
  )
}

calc_mol_mopac <- function(mol) {
  write_mopac_input_file(mol, "aaa.mop")
  run_mopac("aaa.mop")
  mopac_aux <- read_mopac_aux_file("aaa.aux")
  if (length(mopac_aux) > 0) {
# mopac_out <- read_mopac_out_file("aaa.out", mol)
    mopac_fukui <- calc_fukui(mol, mopac_aux)
  } else {
    mopac_fukui <- list()
  }
  c(mopac_aux, mopac_fukui)
}

cmf_calc_mdb_mopac_mem <- function(mdb, print_mopac=TRUE, ...) {
  mopac_mdb_res <- list()
  ncomp <- length(mdb)
  nerrors <- 0
  for (imol in 1:ncomp) {
##  for (imol in 213:217) {
    mol <- mdb[[imol]]
	  mopac_mol_res <- calc_mol_mopac(mol)
	  mopac_mdb_res[[imol]] <- mopac_mol_res
	  if (length(mopac_mol_res) > 0) {
	    if (print_mopac) {
	      cat(paste("Structure", imol, "of", ncomp,
        " natoms=",  mopac_mol_res$natoms,
        " heat=", mopac_mol_res$heat_of_formation,
	      " E(HOMO)=", mopac_mol_res$ehomo, 
		    " E(LUMO)=", mopac_mol_res$elumo, 
		    "\n"))
	    }
	  } else {
	    cat(paste("Structure", imol, "not computed", "\n"))
	    nerrors <- nerrors + 1
	  }
	  flush.console()
  }
  if (nerrors > 0) {
    cat(paste("Computations failed on", nerrors, "structures\n"))
    flush.console()
  }
  mopac_mdb_res
}

cmf_calc_mdb_mopac <- function
(
  mdb_fname = "ligands-train.mol2",
  mopac_res_fname = "ligands-mopac-res-train.RData",
  ...
)
{
   mdb <- read_mol2(mdb_fname)
   mopac_mdb_res <- cmf_calc_mdb_mopac_mem(mdb, ...)
   save(mopac_mdb_res, file=mopac_res_fname)
}


#mdb <- read_mol2("EuSymm.mol2")
#mol <- mdb[[1]]
#res <- calc_mol_mopac(mol)
