# Performing quantum-chemical calculations with MOPAC2012
# for the combines (training+test) set

source("cmf-mopac.R")

# File name of molecular database
mdb_fname <- "uniquebasetest1.mol2"

# File name for the results of quantum-chemical calculations
mopac_res_fname <- "uniquebasetest3.RData"

cmf_calc_mdb_mopac(
  mdb_fname = mdb_fname,
  mopac_res_fname = mopac_res_fname
)
