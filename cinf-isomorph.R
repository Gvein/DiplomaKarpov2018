# Searching isomorphisms

# Finds all substructure isomorphisms
find_substr_isomorph <- function (substr_lab, substr_ct, str_lab, str_ct) {
  
  isomorph <- function(stp) {
    for (i in 1:str_size)
      if (!used[i] && substr_lab[[stp]]==str_lab[[i]]) {
        cc[stp] <<- i
        if (stp > 1)
          for (j in 1:(stp-1))
            if ((substr_ct[stp,j]!=0) && (substr_ct[stp,j]!=str_ct[i,cc[j]])) { 
              to_exit <- TRUE
              next
            }
        if (to_exit) {
          to_exit <- FALSE
          next
        }
        if (stp == substr_size) {
          num_matches <<- num_matches + 1
          isom_list[[num_matches]] <<- cc
        }
        else {
          used[i] <<- TRUE
          isomorph(stp+1)
          used[i] <<- FALSE
        }
      }
  }

  substr_size <- length(substr_lab)
  str_size <- length(str_lab)
  used <- logical(str_size)
  cc <- integer(substr_size)
  to_exit <- FALSE
  num_matches <- 0
  isom_list <- list()
  isomorph(1)
  isom_list
} 

# Finds mol2 in mol1
test_isomorph <- function() {
  # mol1
  mol1_lab <- c("C", "C", "C", "C", "C", "N", "O")
  mol1_ct <- matrix(0, nrow=7, ncol=7)
  mol1_ct[1,2] <- mol1_ct[2,1] <- 1
  mol1_ct[2,3] <- mol1_ct[3,2] <- 1
  mol1_ct[3,4] <- mol1_ct[4,3] <- 1
  mol1_ct[4,5] <- mol1_ct[5,4] <- 1
  mol1_ct[5,6] <- mol1_ct[6,5] <- 2
  mol1_ct[1,6] <- mol1_ct[6,1] <- 1
  mol1_ct[4,7] <- mol1_ct[7,4] <- 2

  # mol2
  mol2_lab <- c("O", "C", "C", "N")
  mol2_ct <- matrix(0, nrow=4, ncol=4)
  mol2_ct[1,2] <- mol2_ct[2,1] <- 2
  mol2_ct[2,3] <- mol2_ct[3,2] <- 1
  mol2_ct[3,4] <- mol2_ct[4,3] <- 2

  isomorph <- find_substr_isomorph(mol2_lab, mol2_ct, mol1_lab, mol1_ct)
}
