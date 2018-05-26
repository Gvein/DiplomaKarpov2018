r_data_name <- "uniquebasetest5.RData"
metal <- "Eu"
ncoord <- 4


load(r_data_name)
nmol <- length(mopac_mdb_res)
dist_mat <- matrix(0,nmol,ncoord)
iat_mat <- matrix(0,nmol,ncoord)
elem_mat <- matrix("",nmol,ncoord)
x_mat <- matrix(0,nmol,ncoord)
y_mat <- matrix(0,nmol,ncoord)
z_mat <- matrix(0,nmol,ncoord)
dist_mat_sort <- matrix(0,nmol,ncoord)
iat_mat_sort <- matrix(0,nmol,ncoord)
elem_mat_sort <- matrix("",nmol,ncoord)
x_mat_sort <- matrix(0,nmol,ncoord)
y_mat_sort <- matrix(0,nmol,ncoord)
z_mat_sort <- matrix(0,nmol,ncoord)
for(imol in 1:nmol) {
  mopac_mol_res <- mopac_mdb_res[[imol]]
  natoms <- mopac_mol_res$natoms
  nat_metal <- 0
  for (iat in 1:natoms) {
    if (mopac_mol_res$atom_el[iat] == metal) {
      nat_metal <- iat
      break
    }
  }
  
  dist_array <- numeric(natoms)
  met_x <- mopac_mol_res$atom_x_opt[nat_metal,1]
  met_y <- mopac_mol_res$atom_x_opt[nat_metal,2]
  met_z <- mopac_mol_res$atom_x_opt[nat_metal,3]
  for (iat in 1:natoms) {
    at_x <- mopac_mol_res$atom_x_opt[iat,1]
    at_y <- mopac_mol_res$atom_x_opt[iat,2]
    at_z <- mopac_mol_res$atom_x_opt[iat,3]
    dist <- sqrt((at_x-met_x)^2 + (at_y-met_y)^2 + (at_z-met_z)^2)
    dist_array[iat] <- dist
  }
  dist_index <- sort(dist_array, index.return=TRUE)
  
  for (j in 1:ncoord) {
    dist_mat[imol,j] <- dist_index$x[j+1]
    iat_mat[imol,j] <- dist_index$ix[j+1]
    elem_mat[imol,j] <- mopac_mol_res$atom_el[iat_mat[imol,j]]
    x_mat[imol,j] <- mopac_mol_res$atom_x_opt[iat,1]
    y_mat[imol,j] <- mopac_mol_res$atom_x_opt[iat,2]
    z_mat[imol,j] <- mopac_mol_res$atom_x_opt[iat,3]
  }
  
  index_sort <- sort(iat_mat[imol,], index.return=TRUE)$ix
  dist_mat_sort[imol,] <- dist_mat[imol,index_sort]
  iat_mat_sort[imol,] <- iat_mat[imol,index_sort]
  elem_mat_sort[imol,] <- elem_mat[imol,index_sort]
  x_mat_sort[imol,] <- x_mat[imol,index_sort]
  y_mat_sort[imol,] <- y_mat[imol,index_sort]
  z_mat_sort[imol,] <- z_mat[imol,index_sort]
}


colnames_list <- list()
for (icol in 1:ncoord) {
  colnames_list[[icol]] <- sprintf("%s-%s%d", metal, elem_mat_sort[1,icol], icol)
}
df <- data.frame(dist_mat_sort)
colnames(df) <- colnames_list
write.table(df, file="SymmStr9Eu3dexDist.txt")



