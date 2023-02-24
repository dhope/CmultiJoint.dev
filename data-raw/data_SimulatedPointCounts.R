## code to prepare `SimulatedPointCounts` dataset goes here
# A few different protocols for distance binning
distance_protocols <- list(p1 = c(0.5,1,Inf),
                           p2 = c(0.5,1,2),
                           p3 = c(0.5,1),
                           p4 = c(0.5,Inf),
                           p5 = Inf,
                           p6 = 1)

# A few different protocols for time binning
time_protocols <- list(p1 = c(3,5,10),
                       p2 = 3,
                       p3 = 5,
                       p4 = 10,
                       p5 = seq(1,10,1))

survey_data_frame <-
  gen_survey_data_frame(sim_rep = 5,
                        tau = 0.2,#seq(0.2,2,0.4),
                        phi = 0.1,#c(0.1,0.5,2.5),
                        Density = 0.1,#c(0.1,0.5,2.5),
                        distance_protocols = distance_protocols,
                        time_protocols = time_protocols) 
  survey_data_frame$seed = 1:nrow(survey_data_frame) 

SimulatedPointCounts <- simulate_point_counts(survey_data_frame,
                                              ncores = floor(parallel::detectCores()/2) )

usethis::use_data(SimulatedPointCounts, overwrite = TRUE)
