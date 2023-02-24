simulate_single_point_count <- function(tau, phi,
                                        Density,
                                        distance_protocols,
                                        time_protocols,
                                        seed,
                                        dim = 10,# landscape x and y dimensions (100 metre increments)
                                        ...){
  if(!requireNamespace("dplyr", quietly = TRUE)){
    stop("This function requires dplyr. Please install it first.")
  }
  # Parameters for this survey
  tau_true <- tau
  phi_true <- phi
  Density_true <- Density*1000
  
  # ------------------------------------
  # Place birds on landscape around observer (centred on landscape)
  # ------------------------------------
  N_landscape <- round(Density_true*dim^2) # Number of birds to place on landscape
  withr::with_seed(seed,{
  birds <- tibble::tibble(x = runif(N_landscape,-dim/2,dim/2),
                      y = runif(N_landscape,-dim/2,dim/2)) |> 
    dplyr::mutate(dist=sqrt(x^2 + y^2)) |>   # Distances to observer
    # Remove birds outside maximum distance
    dplyr::filter( dist <= (dim/2)) |> 
    dplyr::arrange(dist) |> 
    dplyr::mutate(bird_id = dplyr::row_number())
  })
  
  N <-  nrow(birds)
  
  
  # ------------------------------------
  # Simulate bird cues, based on phi_true
  # ------------------------------------
  

  withr::with_seed(seed, {
  cues <- tidyr::expand_grid(bird_id = 1:N, cue_number = 1:50) |> 
    dplyr::arrange(bird_id,cue_number) |> 
    dplyr::group_by(bird_id) |> 
    dplyr::mutate(time = cumsum(rexp(50, phi_true))) |> 
    dplyr::ungroup() |> 
    dplyr::left_join(birds, by = "bird_id") |> 
    # ------------------------------------
  # Determine which cues are detected, based on tau_true
  # ---------
    dplyr::mutate(
      p = exp(-(dist/tau_true)^2),  # Probability each cue is detected
      detected = rbinom(N*50, 1, p)) # binary variable: was cue actually detected?
    
  })
   
      
  # ------------------------------------
  # Isolate first detected cue for each bird
  # ------------------------------------
  
  dat <- dplyr::filter(cues,detected == 1) |> 
    dplyr::group_by(bird_id) |> 
    dplyr::slice_min(order_by = 
                       cue_number,
                     n=1, with_ties = F) |> 
    dplyr::ungroup()
  
  # ------------------------------------
  # Transcription: distance and time bins
  # ------------------------------------
  
  rint <- distance_protocols
  tint <- time_protocols
  # nrint <- length(rint)
  # ntint <- length(tint)
  # browser()
  # Separate into distance and time bins
  dat$rint <- cut(dat$dist,c(0,rint))
  dat$tint <- cut(dat$time,c(0,tint))
  dat <- na.omit(dat)
  
  Y <- table(dat[,c("rint","tint")]) |> as.matrix()
   # Data to analyze
  
  
  
  return( list(Y = Y,rint = rint,tint = tint))
  
  
}


gen_survey_data_frame <- function(
  sim_rep,
  tau,
  phi,
  Density,
  distance_protocols,
  time_protocols
   ){
  tidyr::expand_grid(
    tau,
    phi,
    Density,
    distance_protocols,
    time_protocols,
    nsurveys =1:sim_rep
  ) 
}


simulate_point_counts <- function(survey_data_frame){
  if(!all(c("tau", "phi", "Density","distance_protocols", 
            "time_protocols", 
            "seed") %in% names(survey_data_frame))) 
            stop("simualate_point_counts requires the following columns:
                                          tau, phi, Density,distance_protocols, time_protocols, seed")
  
  
  # # Maximum number of distance bins
  # mdbin <- sapply(distance_protocols,function(x) length(x)) %>% max()
  # 
  # # Maximum number of time bins
  # mtbin <- sapply(time_protocols,function(x) length(x)) %>% max()
  # 
  # # -------------------------------------------------
  # # Arrays to store point count data
  # # -------------------------------------------------
  # 
  # Yarray <- array(NA,dim=c(nrow(survey_data_frame),mdbin,mtbin))
  # rarray <- array(NA,dim=c(nrow(survey_data_frame),mdbin))
  # tarray <- array(NA,dim=c(nrow(survey_data_frame),mtbin))
  
  
  surveys <- 
    purrr::pmap(survey_data_frame,
                simulate_single_point_count ) |> 
    purrr::transpose()
  browser()
  # -------------------------------------------------
  # Arrays to store point count data
  # -------------------------------------------------
  Yarray <- abind::abind(surveys$Y, along = 3) |> 
    aperm(c(3,1,2))
  rarray <-  surveys$rint |> 
    abind::abind(along = 2) |> 
    aperm(c(2,1))
  tarray <- surveys$tint |> 
    abind::abind(along = 2) |> 
    aperm(c(2,1))
  
  return(list(Yarray=Yarray, rarray=rarray, tarray=tarray))
  
}