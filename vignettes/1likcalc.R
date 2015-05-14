###################################################
## Calculate the profile likelihoods over a grid ##
###################################################

library(spatialtalk)
data(anom1962)

#############################################################
##
## Great Circle Distance Matrix
##
#############################################################

## Calculate the distance matrix and save it for later use

# fields implementation
d1 = rdist.earth(loc)

# New implementation
d = rdist_earth1(loc)

## If you don't have enough memory to calculate the distance matrix
## all at once, try this:

##d = matrix(NA, n, n)
##index = round(seq(1, n, length = 5))
##for(i in 1:4){
##  d[index[i]:index[i+1],] <- rdist.earth(loc[index[i]:index[i+1],], loc)
##}

# Add d to the data file
save(z, loc, n, d, file = "anom1962.RData")



#############################################################
##
## nll over grid
##
#############################################################


## Calculate the profile nll's over a grid

# Time required for a single iteration of nll (no tapering)

#- Single runs

# Old function
org.process = system.time({
    org.out = nll(20)
})

# New function
nll.process = system.time({
      n.out = nll_arma(20, d, z, n)
})

#- Multiple runs

rho.seq = seq(20, 80, length = 30)

# Calculate using sapply
org.time = system.time({
  nll.seq.org = sapply(rho.seq, nll)
})

## Parallelized version

# Figure out the total number of cores and use 1 less. (better stability)
n.cores = parallel::detectCores() - 1

# Parallelized with parallel package
library(parallel)

# Create snow Cluster
parallel.time.snow = rep(NA,n.cores)

for(core in 1:n.cores){
  #Create cluster
  cl = makeCluster(core)
  
  # Parallelize
  parallel.time.snow[i] = system.time({
    # Run parallelization
    v = parSapply(cl = cl, X = rho.seq, FUN = nll2, d=d, z=z, n=n)
  })[3]
  
  # Kill the cluster
  stopCluster(cl)
}


# Parallelized with OpenMP in C++
parallel.time.cpp = rep(NA,n.cores)

for(core in 1:n.cores){
 parallel.time.cpp[i] = system.time({
   # Call the C++ parallelized version of nll
   nll.seq = nll_parallel(rho.seq, d, z, n, core)
  })[3]
}

#############################################################
##
## Tapering Setup Functions
##
#############################################################

## You can change the taper range (gamma in the paper) here

#- Original function
system.time({ 
  setup = make.tapersetup(d, wendland2.1, taprange = 50)
})

#- New function
system.time({ 
  setup.cpp = make_tapersetup_R(d, taprange = 50)
})

#- Modified function for Eigen setup
system.time({ 
  setup.eigen = make_tapersetup_eigen(d,taprange = 50)
})

# Free up some memory
rm(d); gc() 



#############################################################
##
## Tapering Method I and II functions
##
#############################################################

rho.seq = seq(20, 80, length = 30)

# Old tapering functions:
nll.1taper.seq = sapply(rho.seq, nll.1taper, setup = setup)

nll.2taper.seq = sapply(rho.seq, nll.2taper, setup = setup)

# Old tapering functions parallelized:

# Create storage for 1taper timings
parallel.time.snow.1taper = rep(NA,n.cores)

for(core in 1:n.cores){
  #Create cluster
  cl = makeCluster(core)
  
  # Parallelize
  parallel.time.snow.1taper[i] = system.time({
    # Run parallelization
    v = parSapply(cl = cl, X = rho.seq, FUN = nll.1taper, setup = setup)
  })[3]
  
  # Kill the cluster
  stopCluster(cl)
}

# Create storage for 2taper timings
parallel.time.snow.2taper = rep(NA,n.cores)

for(core in 1:n.cores){
  #Create cluster
  cl = makeCluster(core)
  
  # Parallelize
  parallel.time.snow.2taper[i] = system.time({
    # Run parallelization
    v = parSapply(cl = cl, X = rho.seq, FUN = nll.2taper, setup = setup)
  })[3]
  
  # Kill the cluster
  stopCluster(cl)
}


# Create storage for 2taper timings
parallel.time.nll1 = rep(NA,n.cores)

# Create storage for 2taper timings
parallel.time.nll2 = rep(NA,n.cores)

# New tapering functions:
for(core in 1:n.cores){
  parallel.time.nll1[core] = system.time({
  nll1.val = nll_1taper_parallel(rho.seq,
                                 setup.eigen$n,
                                 setup.eigen$good.dists,
                                 setup.eigen$taps,
                                 setup.eigen$ia,
                                 setup.eigen$ja,
                                 z,
                                 setup.eigen$rescol,
                                 core)
  })[3]

}

for(core in 1:n.cores){
  parallel.time.nll2[core] = system.time({
    nll2.val = nll_2taper_parallel(
                          rho.seq,
                          setup.eigen$n,
                          setup.eigen$good.dists,
                          setup.eigen$taps,
                          setup.eigen$ia,
                          setup.eigen$ja,
                          z,
                          setup.eigen$rescol, core)
  })[3]
}

# Save results to file to avoid recomputing
save(rho.seq, nll.seq, nll.1taper.seq, nll.2taper.seq,
     file = "likcalc.RData")
