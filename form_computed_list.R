# Form list of structures successfully processed by MOPAC
# and forms new rdata file with these structures

mopac_res_fname <- "uniquebasetest3.RData"
computed_list_fname <- "computed.txt"
new_mopac_res_fname <- "uniquebasetest5.RData"

computed_list <- list()
new_mopac_mdb_res <- list()
load(mopac_res_fname)
nmol <- length(mopac_mdb_res)
icomp <- 0
for (imol in 1:nmol) {
  if (length(mopac_mdb_res[[imol]]) > 0) {
    icomp <- icomp + 1
    computed_list[[icomp]] <- imol
    new_mopac_mdb_res[[icomp]] <- mopac_mdb_res[[imol]]
  }
}
mopac_mdb_res <- new_mopac_mdb_res
f <- file(computed_list_fname, "w")
for (num in computed_list) cat(num, file=f, sep="\n")
close(f)
save(mopac_mdb_res, file=new_mopac_res_fname)
