## code to prepare `SimulatedPointCounts` dataset goes here
# A few different protocols for distance binning
distance_protocols <- list(p1 = c(0.5,1,Inf))

# A few different protocols for time binning
time_protocols <- list(p1 = seq(1,10))

survey_data_frame <-
  gen_survey_data_frame(sim_rep = 5,
                        tau = 0.2,#seq(0.2,2,0.4),
                        phi = 0.1,#c(0.1,0.5,2.5),
                        Density = 0.1,#c(0.1,0.5,2.5),
                        distance_protocols = distance_protocols,
                        time_protocols = time_protocols) |> 
  dplyr::mutate(seed = 1:1 ) 

SimulatedPointCounts <- simulate_point_counts(survey_data_frame)

usethis::use_data(SimulatedPointCounts, overwrite = TRUE)
