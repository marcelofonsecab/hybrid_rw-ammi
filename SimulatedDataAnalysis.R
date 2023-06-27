# Read the data from the specified RData file
# In this example its a data with 2 percent of contamination and Shift with k = 4
alldata = readRDS("Simulated Data\\SimulatedData_percent2_shift4.RData")

# All possible datas are:
# "Simulated Data\\SimulatedData_notcontaminated.RData" (not contaminated data)
# "Simulated Data\\SimulatedData_percent2_shift4.RData" (2% and shift k = 4)
# "Simulated Data\\SimulatedData_percent5_shift4.RData" (5% and shift k = 4)
# "Simulated Data\\SimulatedData_percent10_shift4.RData" (10% and shift k = 4)
# "Simulated Data\\SimulatedData_percent2_shift7.RData" (2% and shift k = 7)
# "Simulated Data\\SimulatedData_percent5_shift7.RData" (5% and shift k = 7)
# "Simulated Data\\SimulatedData_percent10_shift7.RData" (10% and shift k = 7)
# "Simulated Data\\SimulatedData_percent2_shift10.RData" (2% and shift k = 10)
# "Simulated Data\\SimulatedData_percent5_shift10.RData" (5% and shift k = 10)
# "Simulated Data\\SimulatedData_percent10_shift10.RData" (10% and shift k = 10)
# "Simulated Data\\SimulatedData_percent2_pointmass.RData" (2% and shift c = 10)
# "Simulated Data\\SimulatedData_percent5_pointmass.RData" (5% and shift c = 10)
# "Simulated Data\\SimulatedData_percent10_pointmass.RData" (10% and shift c = 10)
# "Simulated Data\\SimulatedData_percent2_varinflated.RData" (2% and shift c = .2)
# "Simulated Data\\SimulatedData_percent5_varinflated.RData" (5% and shift c = .2)
# "Simulated Data\\SimulatedData_percent10_varinflated.RData" (10% and shift c = .2)


# Load the Methods.R script
source("Methods.R")

# Attach the first dataset in alldata
attach(alldata[[1]])

# Perform all SVDs from all models on the data
allsvds = All_SVDS(data$data.not_contaminated,
                   data$data.contaminated,
                   weights = weights)

# Calculate GEI metrics for the obtained SVD results
.allMetrics = GEI.Metrics(allsvds)

# Iterate over the remaining datasets in alldata
for (i in 2:5) {
  attach(alldata[[i]])  # Attach the current dataset
  
  # Perform all SVDs from all models on the data of the current dataset
  allsvds = All_SVDS(data$data.not_contaminated,
                     data$data.contaminated,
                     weights = weights)
  
  # Calculate GEI metrics for the current SVD results
  .tmpMetrics = GEI.Metrics(allsvds)
  
  # Append the metrics of the current dataset to the overall metrics
  .allMetrics = appendList(.allMetrics, .tmpMetrics)
}

# The names of the models are defined as a character vector representing the names of each model.
names = c("AMMI",
          "W-AMMI{1}", "W-AMMI{2}", "W-AMMI{3}",
          "W-AMMI{4}", "W-AMMI{5}"," W-AMMI{6}",
          "W-AMMI{7}", "W-AMMI{8}", "W-AMMI{9}",
          "R-AMMI{1}", "R-AMMI{2}", "R-AMMI{3}",
          "R-AMMI{4}", "R-AMMI{5}", "R-AMMI{6}",
          "R-AMMI{7}", "R-AMMI{8}", "R-AMMI{9}",
          "RW-AMMI{1}", "RW-AMMI{2}", "RW-AMMI{3}",
          "RW-AMMI{4}", "RW-AMMI{5}", "RW-AMMI{6}",
          "RW-AMMI{7}", "RW-AMMI{8}", "RW-AMMI{9}")

# Metricsapply function is used to calculate the "MPEV" metric for all models in .allMetrics. 
# hypothesis.testing is set to TRUE, indicating that hypothesis testing is performed.
# h0 is set to 100, which represents the null hypothesis value.
# side is set to "greater", indicating that the alternative hypothesis is greater than h0.
Values = Metricsapply(.allMetrics, "MPEV",
                      hypothesis.testing = T,
                      h0 = 100, side = "greater")

# The names of the models are assigned to the corresponding values in Values$Metric and Values$pvalues.
names(Values$Metric) = names(Values$pvalues) = names

# It needs to do it for each metric, for every contamination, including the non-contaminated scenario
# And bind all values into one matrix to make the Tables which are showed in the paper

#### Biplots
# For the biplots we got a random seed using the sample() function, to select one simulated data to analysis the biplot
set.seed(768622) # the number 768622 was generated from a uniform distribution from -999999 to 999999
dataseed = sample(1:100, 1)

alldata = readRDS("Simulated Data\\SimulatedData_percent5_shift7.RData") # for contamination with 5% and shift k = 7
attach(alldata[[dataseed]])

# For biplots without contamination we used the following
allsvds.notcontaminated = All_SVDS(data$data.not_contaminated,
                                   data$data.not_contaminated,
                                   weights = weights)
BiplotCreation(allsvds.notcontaminated$WAMMI$WAMMI_lmm.Weights.Env$SVD)


# Perform all SVDs from all models on the data from the selected seed
allsvds = All_SVDS(data$data.not_contaminated,
                   data$data.contaminated,
                   weights = weights)
