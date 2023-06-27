# Load the Methods.R script
## Required Packages ##
# install.packages("dplyr")
# install.packages("pbapply")
# install.packages("lme4")
# install.packages("robustlmm")
# install.packages("MASS")
# install.packages("reshape2")
# install.packages("rsvddpd")
source("Methods.R")

# Example of simulating one dataset (seed 1, 2% contamination for all types of contamination) with and without contamination

shift_values = c(4, 7, 10) # Shift values k = 4, 7 or 10
varinflated_values = 0.2 # Variance inflated values c = .2
pointmass_values = 10 # Pointmass values c = 10
percent_values = 0.02 # Percentage equals 2%
#percent_values = c(0.02, 0.05, 0.1)
seed.values = 1
#seeds.values = 1:100

cl = NULL # If you want, you cant add a cluster with the help of the package 'parallel' to speedup simulations
# require(parallel)
# cl = makeCluster(getOption("cl.cores", 2))
# clusterEvalQ(cl, library(robustlmm))
# clusterExport(cl, varlist = ls())

list.tmp = list()

for(i in 1:length(seeds.values)){
  set.seed(seeds.values[i])
  print(paste0("Seed: ", seeds.values[i]))
  data.tmp = sim.amb(seed = seeds.values[i],
                      Ngen = 100, Nenv = 8, Ncomp = 2,
                      effectGlobal = c(mean = 15, sd = sqrt(3)),
                      effectGen = c(mean = 5, sd = 1),
                      effectEnv = c(mean = 8, sd = sqrt(2)),
                      k = c(28, 15))
  data.rep.tmp = Create.Replications(data.tmp, reps = 2, sig = 1)

  for(j in 1:length(percent_values)){
    print(paste0("Percentage: ", percent_values[j]))
    for(k in 1:length(shift_values)){
      print(paste0("Shift: ", shift_values[k]))
      data.cont.tmp = Data_Contamination(data = data.rep.tmp,
                                          percentage = percent_values[j],
                                          seed = seeds.values[i],
                                          type = "shift",
                                          k = shift_values[k])
      errvar.tmp = suppressMessages(Error_Var(data.cont.tmp$data.cont,
                                              cluster = cl))
      w.tmp = R.Weights(data.cont.tmp$data.cont, errvar.tmp)
      
      list.tmp[[paste0("shift_perc", percent_values[j]*100, "_k",shift_values[k])]] = list(data = data.cont.tmp,
                       errorvariances = errvar.tmp,
                       weights = w.tmp)
      #list.tmp = list(data = data.cont.tmp,
      #                errorvariances = errvar.tmp,
      #                weights = w.tmp)
      #name.tmp = paste0("Seed", seeds.values[i],
      #                   "_perc", percent_values[j] * 100,
      #                   "_shift", shift_values[k],
      #                   ".RData")
      #saveRDS(list.tmp, file = name.tmp)
    }
    for(k in 1:length(varinflated_values)){
      print(paste0("HighVariance: ", varinflated_values[k]))
      data.cont.tmp = Data_Contamination(data = data.rep.tmp,
                                          percentage = percent_values[j],
                                          seed = seeds.values[i],
                                          type = "varinflated",
                                          c = varinflated_values[k])
      errvar.tmp = suppressMessages(Error_Var(data.cont.tmp$data.cont,
                                              cluster = cl))
      w.tmp = R.Weights(data.cont.tmp$data.cont, errvar.tmp)
      
      list.tmp[[paste0("varinflated_perc", percent_values[j]*100)]] = list(data =data.cont.tmp,
                       errorvariances = errvar.tmp,
                       weights = w.tmp)
      
     # list.tmp = list(data =data.cont.tmp,
     #                 errorvariances = errvar.tmp,
     #                 weights = w.tmp) 
     # name.tmp = paste0("Seed", seeds.values[i],
     #                    "_perc", percent_values[j] * 100,
     #                    "_hvar", k,
     #                    ".RData")
     # saveRDS(list.tmp, file = name.tmp)
    }
    for(k in 1:length(pointmass_values)){
      print(paste0("PointMass: ", pointmass_values[k]))
      data.cont.tmp = Data_Contamination(data = data.rep.tmp,
                                          percentage = percent_values[j],
                                          seed = seeds.values[i],
                                          type = "pointmass",
                                          c = pointmass_values[k])
      errvar.tmp = suppressMessages(Error_Var(data.cont.tmp$data.cont,
                                              cluster = cl))
      w.tmp = R.Weights(data.cont.tmp$data.cont, errvar.tmp)
      
      list.tmp[[paste0("pointmass_perc", percent_values[j]*100)]] = list(data =data.cont.tmp,
                       errorvariances = errvar.tmp,
                       weights = w.tmp)
      
      #list.tmp = list(data = data.cont.tmp,
      #                errorvariances = errvar.tmp,
      #                weights = w.tmp)
      #name.tmp = paste0("Seed", seeds.values[i],
      #                   "_perc", percent_values[j] * 100,
      #                   "_pmass", pointmass_values[k],
      #                   ".RData")
      #saveRDS(list.tmp, file = name.tmp)
    }
  }
  print(paste0("No contamination"))
  errvar.tmp = suppressMessages(Error_Var(data.rep.tmp,
                                          cluster = cl))
  w.tmp = R.Weights(data.rep.tmp, errvar.tmp)
  
  list.tmp[["nocontamination"]] = list(data = data.rep.tmp,
                   errorvariances = errvar.tmp,
                   weights = w.tmp)
  
  #name.tmp = paste0("Seed", seeds.values[i], ".RData")
  #saveRDS(list.tmp, file = name.tmp)
  
}

list.tmp$shift_perc2_k4 # data with shift k = 4
list.tmp$shift_perc2_k7 # data with shift k = 7
list.tmp$shift_perc2_k10 # data with shift k = 10
list.tmp$varinflated_perc2 # data with variance inflated
list.tmp$pointmass_perc2 # data with pointmass
list.tmp$nocontamination # data without contamination
