r_data_name <- "uniquebasetest5.RData"
load(r_data_name)
nmol <- length(mopac_mdb_res)
heat_of_formation <- numeric(nmol)
total_energy <- numeric(nmol)
for (imol in 1:nmol) {
  mopac_mol_res <- mopac_mdb_res[[imol]]
  heat_of_formation[imol] <- mopac_mol_res$heat_of_formation
  total_energy[imol] <- mopac_mol_res$total_energy
}
df <- data.frame(hform=heat_of_formation, tenergy=total_energy)
write.table(df, file="uniquebasetest5.txt")