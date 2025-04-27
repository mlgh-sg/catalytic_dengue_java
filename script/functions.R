#### function to simulate dengue cases
#### input:
#### nT, nA, lamH, lam, rho, gamma, pop, amin, amax
#### nT: number of years of cases simulations (based on the last nT of lam)
#### nA: number of age groups
#### lamH: a vector of historical lambda - 1/4 total FoI
#### lam: lambda over time during observation period - a vector - 1/4 total FoI
#### rho: reporting rate of dengue incidence
#### gamma: relative reporting of primary infections compared to secondary
#### pop: population by age group
#### amin: lower bound of age groups (minimum of 0)
#### amax: upper bound of age groups (maximum of 99)
#### this approach edit slightly the way the algorithm works in o'driscoll et al
#### 1 more flexible historical lambda
#### 2 the way progression of immunity worked, 1 year means already exposed for a year previously

lam <- rep(0.03,100)

simcases <- function(nT, nA, lamH, lam, rho, gamma, pop, amin, amax){
  
  # with slight tweaks from o'driscoll et al
  # allowing a more dynamic historical lamH
  
  # immune profiles at beginning of time period
  age <- seq(0,99)
  age_idx <- 1:length(age)
  age_rev_idx <- length(age):1
  
  susc <- mono <- matrix(NA, nrow=nT+1, ncol=100)
  inc1 <- inc2 <- matrix(NA, nrow=nT, ncol=100)
  
  # initial susceptibility & monoinfection proportions
  for (i in 1:100){
    
    if (i > length(lamH)) {
      sum_lamH <- sum(lamH[1:length(lamH)])
    } else {
      sum_lamH <- sum(lamH[1:i])
    }
    
    susc[1,i] <- exp(-4*sum_lamH)
    mono[1,i] <- 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH))
      
  }
  
  # loop the catalytic model through time of observations
  for (t in 2:(nT+1)){
    
    susc[t,1] <- exp(-4*lam[t-1])
    mono[t,1] <- 4*exp(-3*lam[t-1])*(1-exp(-lam[t-1]))
    susc[t,2:100] <- susc[t-1,1:99] - 4*lam[t-1]*susc[t-1,1:99]
    mono[t,2:100] <- mono[t-1,1:99] + 4*lam[t-1]*susc[t-1,1:99] - 3*lam[t-1]*mono[t-1,1:99]
    inc1[t-1,1] <- 4*lam[t-1]*1
    inc2[t-1,1] <- 3*lam[t-1]*0
    inc1[t-1,2:100] <- 4*lam[t-1]*susc[t-1,1:99]
    inc2[t-1,2:100] <- 3*lam[t-1]*mono[t-1,1:99]
    
  }
  
  # reported cases
  cases_S <- matrix(NA, nrow=nT, ncol=nA)
  cases_PS <- matrix(NA, nrow=nT, ncol=nA)
  prop_secondary <- matrix(NA, nrow=nT, ncol=nA)
  
  for(t in 1:nT) for(a in 1:nA){
    age_index_min <- which(age==amin[a])
    age_index_max <- which(age==amax[a])
    
    cases_S[t,a] = round(rho*(mean(inc2[t,age_index_min:age_index_max]))*pop[t,a],0)
    cases_PS[t,a] = round(rho*(mean(inc2[t,age_index_min:age_index_max]) + gamma*mean(inc1[t,age_index_min:age_index_max]))*pop[t,a],0)
    prop_secondary[t,a] =  round(rho*mean(inc2[t,age_index_min:age_index_max])*pop[t,a],0)/cases_PS[t,a]
  }
  
  totalcases_S <- rowSums(cases_S)
  totalcases_PS <- rowSums(cases_PS)
  
  return(list(cases_S=cases_S, cases_PS=cases_PS, 
              totalcases_S=totalcases_S, totalcases_PS=totalcases_PS, 
              prop_secondary=prop_secondary,
              susc=susc, mono=mono))
}

simcases2 <- function(nT, nA, lamH, lam, rho, gamma, pop, amin, amax){
  
  # with slight tweaks from o'driscoll et al
  # allowing a more dynamic historical lamH
  
  # immune profiles at beginning of time period
  age <- seq(0,99)
  age_idx <- 1:length(age)
  age_rev_idx <- length(age):1
  
  susc <- mono <- matrix(NA, nrow=nT+1, ncol=100)
  inc1 <- inc2 <- matrix(NA, nrow=nT, ncol=100)
  
  # initial susceptibility & monoinfection proportions
  for (i in 1:100){
    
    if (i > length(lamH)) {
      sum_lamH <- sum(lamH[1:length(lamH)])
    } else {
      sum_lamH <- sum(lamH[1:i])
    }
    
    susc[1,i] <- exp(-4*sum_lamH)
    mono[1,i] <- 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH))
    
  }
  
  # loop the catalytic model through time of observations
  for (t in 2:(nT+1)){
    
    susc[t,1] <- exp(-4*lam[t-1])
    mono[t,1] <- 4*exp(-3*lam[t-1])*(1-exp(-lam[t-1]))
    susc[t,2:100] <- susc[t-1,1:99] - 4*lam[t-1]*susc[t-1,1:99]
    mono[t,2:100] <- mono[t-1,1:99] + 4*lam[t-1]*susc[t-1,1:99] - 3*lam[t-1]*mono[t-1,1:99]
    inc1[t-1,1] <- 4*lam[t-1]*1
    inc2[t-1,1] <- 3*lam[t-1]*0
    inc1[t-1,2:100] <- 4*lam[t-1]*susc[t-1,1:99]
    inc2[t-1,2:100] <- 3*lam[t-1]*mono[t-1,1:99]
    
  }
  
  # reported cases
  cases_S <- matrix(NA, nrow=nT, ncol=nA)
  cases_PS <- matrix(NA, nrow=nT, ncol=nA)
  cases_S2 <- matrix(NA, nrow=nT, ncol=nA)
  cases_PS2 <- matrix(NA, nrow=nT, ncol=nA)
  susc_count <- matrix(NA, nrow=nT, ncol=nA)
  susc_prop <- matrix(NA, nrow=nT, ncol=nA)
  susc9_prop <- matrix(NA, nrow=nT, ncol=1)
  mono_count <- matrix(NA, nrow=nT, ncol=nA)
  mono_prop <- matrix(NA, nrow=nT, ncol=nA)
  prop_secondary <- matrix(NA, nrow=nT, ncol=nA)
  
  for(t in 1:nT) for(a in 1:nA){
    age_index_min <- which(age==amin[a])
    age_index_max <- which(age==amax[a])
    
    cases_S[t,a] = round(rho*(mean(inc2[t,age_index_min:age_index_max]))*pop[t,a],0)
    cases_PS[t,a] = round(rho*(mean(inc2[t,age_index_min:age_index_max]) + gamma*mean(inc1[t,age_index_min:age_index_max]))*pop[t,a],0)
    cases_S2[t,a] = round((mean(inc2[t,age_index_min:age_index_max]))*pop[t,a],0)
    cases_PS2[t,a] = round((mean(inc2[t,age_index_min:age_index_max]) + gamma*mean(inc1[t,age_index_min:age_index_max]))*pop[t,a],0)
    susc_count[t,a] = round((mean(susc[t,age_index_min:age_index_max]))*pop[t,a],0)
    susc_prop[t,a] = mean(susc[t,age_index_min:age_index_max])
    mono_count[t,a] = round((mean(mono[t,age_index_min:age_index_max]))*pop[t,a],0)
    mono_prop[t,a] = mean(mono[t,age_index_min:age_index_max])
    prop_secondary[t,a] =  round(rho*mean(inc2[t,age_index_min:age_index_max])*pop[t,a],0)/cases_PS[t,a]
  }
  
  for(t in 1:nT){
    susc9_prop[t,1] <- susc[t,10:10]
  }
  
  totalcases_S <- rowSums(cases_S)
  totalcases_PS <- rowSums(cases_PS)
  totalcases_S2 <- rowSums(cases_S2)
  totalcases_PS2 <- rowSums(cases_PS2)
  totalsusc <- rowSums(susc_count)
  totalmono <- rowSums(mono_count)
  
  return(list(cases_S=cases_S, cases_PS=cases_PS,
              cases_S2=cases_S2, cases_PS2=cases_PS2,
              totalcases_S=totalcases_S, totalcases_PS=totalcases_PS, 
              totalcases_S2=totalcases_S2, totalcases_PS2=totalcases_PS2, 
              susc_count=susc_count, mono_count=mono_count, 
              susc_prop=susc_prop, mono_prop=mono_prop, 
              susc9_prop=susc9_prop,
              totalsusc=totalsusc, totalmono=totalmono))
}

# add stochasticity using poisson distribution for case simulations
simcases_sto <- function(nT, nA, lamH, lam, rho, gamma, phi=NULL, dist="Poisson", 
                         pop, amin, amax, seed=12345){
  
  # with slight tweaks from o'driscoll et al
  # allowing a more dynamic historical lamH
  
  # immune profiles at beginning of time period
  age <- seq(0,99)
  age_idx <- 1:length(age)
  age_rev_idx <- length(age):1
  
  susc <- mono <- matrix(NA, nrow=nT+1, ncol=100)
  inc1 <- inc2 <- matrix(NA, nrow=nT, ncol=100)
  
  # initial susceptibility & monoinfection proportions
  for (i in 1:100){
    
    if (i > length(lamH)) {
      sum_lamH <- sum(lamH[1:length(lamH)])
    } else {
      sum_lamH <- sum(lamH[1:i])
    }
    
    susc[1,i] <- exp(-4*sum_lamH)
    mono[1,i] <- 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH))
    
  }
  
  # loop the catalytic model through time of observations
  for (t in 2:(nT+1)){
    
    susc[t,1] <- exp(-4*lam[t-1])
    mono[t,1] <- 4*exp(-3*lam[t-1])*(1-exp(-lam[t-1]))
    susc[t,2:100] <- susc[t-1,1:99] - 4*lam[t-1]*susc[t-1,1:99]
    mono[t,2:100] <- mono[t-1,1:99] + 4*lam[t-1]*susc[t-1,1:99] - 3*lam[t-1]*mono[t-1,1:99]
    inc1[t-1,1] <- 4*lam[t-1]*1
    inc2[t-1,1] <- 3*lam[t-1]*0
    inc1[t-1,2:100] <- 4*lam[t-1]*susc[t-1,1:99]
    inc2[t-1,2:100] <- 3*lam[t-1]*mono[t-1,1:99]
    
  }
  
  # reported cases
  cases_S <- matrix(NA, nrow=nT, ncol=nA)
  cases_PS <- matrix(NA, nrow=nT, ncol=nA)
  prop_secondary <- matrix(NA, nrow=nT, ncol=nA)
  
  for(t in 1:nT) for(a in 1:nA){
    age_index_min <- which(age==amin[a])
    age_index_max <- which(age==amax[a])
    
    if (dist=="Poisson"){
      set.seed(seed)
      cases_S[t,a] = rpois(1,round(rho*(mean(inc2[t,age_index_min:age_index_max]))*pop[t,a],0))
      set.seed(seed)
      cases_PS[t,a] = rpois(1,round(rho*(mean(inc2[t,age_index_min:age_index_max]) + gamma*mean(inc1[t,age_index_min:age_index_max]))*pop[t,a],0))
      prop_secondary[t,a] =  cases_S[t,a]/cases_PS[t,a]
    } else if(dist=="Negbin"){
      set.seed(seed)
      cases_S[t,a] = rnbinom(1,mu=round(rho*(mean(inc2[t,age_index_min:age_index_max]))*pop[t,a],0),size=phi)
      set.seed(seed)
      cases_PS[t,a] = rnbinom(1,mu=round(rho*(mean(inc2[t,age_index_min:age_index_max]) + gamma*mean(inc1[t,age_index_min:age_index_max]))*pop[t,a],0),size=phi)
      prop_secondary[t,a] =  cases_S[t,a]/cases_PS[t,a]
    } else{
      stop("dist should be either Poisson or Negbin")
    }
    
  }
  
  totalcases_S <- rowSums(cases_S)
  totalcases_PS <- rowSums(cases_PS)
  
  return(list(cases_S=cases_S, cases_PS=cases_PS, 
              totalcases_S=totalcases_S, totalcases_PS=totalcases_PS, 
              prop_secondary=prop_secondary,
              susc=susc, mono=mono))
}

#### customised function to generate seroprevalence in 2014
#### run with historical lam_H up to 2013
#### first column index means 1-year old
simseroprev_hist <- function(lamH){
  
  # with slight tweaks from o'driscoll et al
  # allowing a more dynamic historical lamH
  
  # immune profiles at beginning of time period
  age <- seq(0,99)
  age_idx <- 1:length(age)
  age_rev_idx <- length(age):1
  
  susc <- mono <- matrix(NA, nrow=1, ncol=100)
  
  # initial susceptibility & monoinfection proportions
  for (i in 1:100){
    
    if (i > length(lamH)) {
      sum_lamH <- sum(lamH[1:length(lamH)])
    } else {
      sum_lamH <- sum(lamH[1:i])
    }
    
    susc[1,i] <- exp(-4*sum_lamH)
    mono[1,i] <- 4*exp(-3*sum_lamH)*(1-exp(-sum_lamH))
    
  }
  
  return(list(susc=as.numeric(susc), mono=as.numeric(mono)))
  
}

#### function to get posterior samples from stan fit output
#### input:
#### stan_output
#### params
#### sample_size: sample size from posterior samples for further simulations
posterior_samples_output <- function(stan_output,params,sample_size){
  
  output <- rstan::extract(stan_output,pars=params)
  
  index_of_lam_H <- which(params=="lam_H")
  index_of_samples <- sample(1:length(output[[index_of_lam_H]]),sample_size)
  
  param_samples_list <- list()
  
  for (i in seq_len(length(params))){
    
    if (params[i] != "lam_t"){
      param_samples_list[[i]] <- output[[i]][index_of_samples] 
    } else {
      param_samples_list[[i]] <- matrix(NA,ncol=ncol(output[[i]]),nrow=sample_size)
      for (j in seq_len(ncol(output[[i]]))){
        param_samples_list[[i]][,j] <- output[[i]][index_of_samples,j]
      }
    }
  }
  
  names(param_samples_list) <- params
  
  return(param_samples_list)
}

#### function to summarise all params from stan output
#### input:
#### stan_output
#### params
#### original_values: named list should be in the same order as params
summarise_posterior <- function(stan_output,params,original_values){
  
  if(sum(names(original_values)==params) != length(params)){
    stop("params and original values name not in the same order")
  } else{
    output <- as_tibble(summary(stan_output,params)$summary) %>% 
      mutate(pars=rownames(summary(stan_output,params)$summary)) %>% 
      select(pars,everything())
  }
  
  plot_params <- list()
  
  for(i in seq_len(length(params))){
    
    output_filtered <- output %>% 
      filter(str_detect(pars,params[i]))
    
    if (params[i] == "lam_H"){
      
      if (sum(is.na(original_values[["lam_H"]]))){
        lam_H_ori <- tibble(pars="lam_H",t=-length(original_values[["lam_H"]]):-1)
        output_plot <- lam_H_ori %>% left_join(output_filtered)
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=t)) +
          geom_ribbon(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),fill="#fe0332",alpha=0.10) +
          geom_ribbon(aes(ymin=4*`25%`,ymax=4*`75%`),fill="#fe0332",alpha=0.50) +
          geom_line(aes(y=4*`50%`),col="#fe0332",linewidth=1.5) +
          theme_classic(base_size=12) +
          labs(x="Time",y="4*lam_H")
      } else {
        lam_H_ori <- tibble(pars="lam_H",t=-length(original_values[["lam_H"]]):-1,
                            original_values=original_values[["lam_H"]])
        output_plot <- lam_H_ori %>% left_join(output_filtered)
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=t)) +
          geom_ribbon(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),fill="#fe0332",alpha=0.10) +
          geom_ribbon(aes(ymin=4*`25%`,ymax=4*`75%`),fill="#fe0332",alpha=0.50) +
          geom_line(aes(y=4*`50%`),col="#fe0332",linewidth=1.5) +
          geom_point(aes(y=4*original_values)) +
          theme_classic(base_size=12) +
          labs(x="Time",y="4*lam_H")
      }
      
    } else {
      
      params_values <- output_filtered$pars
      
      if (sum(is.na(original_values[[params[i]]]))){
        output_plot <- output_filtered %>% 
          mutate(pars=factor(pars,levels=params_values))
        if (params[i]=="lam_t"){
          plot_params[[i]] <- output_plot %>% 
            ggplot(aes(x=pars)) +
            geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),col="#fe0332",linewidth=1) +
            geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),col="#fe0332",linewidth=2) +
            geom_point(aes(y=4*`50%`),col="#fe0332",size=2.5) +
            theme_classic(base_size=12) +
            labs(x="Params",y=paste0("4*",params[[i]]))
        } else{
          plot_params[[i]] <- output_plot %>% 
            ggplot(aes(x=pars)) +
            geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),col="#fe0332",linewidth=1) +
            geom_linerange(aes(ymin=`25%`,ymax=`75%`),col="#fe0332",linewidth=2) +
            geom_point(aes(y=`50%`),col="#fe0332",size=2.5) +
            theme_classic(base_size=12) +
            labs(x="Params",y=params[[i]])
        }
      } else {
        output_plot <- output_filtered %>% 
          mutate(pars=factor(pars,levels=params_values),
                 original_values=original_values[[params[i]]])
        
        if (params[i]=="lam_t"){
          plot_params[[i]] <- output_plot %>% 
            ggplot(aes(x=pars)) +
            geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),col="#fe0332",linewidth=1) +
            geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),col="#fe0332",linewidth=2) +
            geom_point(aes(y=4*`50%`),col="#fe0332",size=2.5) +
            geom_point(aes(y=4*original_values),size=2.5) +
            theme_classic(base_size=12) +
            labs(x="Params",y=paste0("4*",params[[i]]))
        } else{
          plot_params[[i]] <- output_plot %>% 
            ggplot(aes(x=pars)) +
            geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),col="#fe0332",linewidth=1) +
            geom_linerange(aes(ymin=`25%`,ymax=`75%`),col="#fe0332",linewidth=2) +
            geom_point(aes(y=`50%`),col="#fe0332",size=2.5) +
            geom_point(aes(y=original_values),size=2.5) +
            theme_classic(base_size=12) +
            labs(x="Params",y=params[[i]])
        }
      }
      
    }
    
  }
  
  return(list(summary=output,plot=plot_params))
}

#### function to simulate cases and compare to original data
#### input:
#### posterior_samples: from posterior_samples_output
#### cases: cases data fitted to the model
#### pop
#### amin
#### amax
#### hist_length
#### model_type: model type cases_S or cases_PS
#### dist: distribution base: Poisson will compare counts, multinomial will compare proportions
simulate_cases_from_posterior <- function(posterior_samples,
                                          cases,
                                          pop,
                                          amin,
                                          amax,
                                          hist_length,
                                          model_type,
                                          dist){
  
  if (dist=="Poisson"){
    
    nT <- ncol(posterior_samples[["lam_t"]])
    nA <- length(amin)
    age_groups <- paste(amin,amax,sep='-')
    
    simulated_cases_list <- list()
    for (i in seq_len(length(posterior_samples[["lam_H"]]))){
      
      lamH <- rep(posterior_samples[["lam_H"]][i],hist_length)
      lam <- posterior_samples[["lam_t"]][i,]
      rho <- posterior_samples[["rho"]][i]
      gamma <- posterior_samples[["gamma"]][i]
      simulated_output <- simcases(nT=nT, nA=nA, lamH=lamH, lam=lam, 
                                   rho=rho, gamma=gamma, pop=pop, amin=amin, amax=amax)
      simulated_cases_df <- simulated_output[[model_type]]
      colnames(simulated_cases_df) <- age_groups
      simulated_cases_list[[i]] <- simulated_cases_df %>% 
        as_tibble() %>% 
        mutate(t=1:nT,iter=i) %>% 
        pivot_longer(-c(t,iter),names_to="age_group",values_to="simulated_cases") %>% 
        mutate(age_group=factor(age_group,levels=age_groups))
      
    }
    
    simulated_cases <- bind_rows(simulated_cases_list)
    simulated_cases_summary <- simulated_cases %>% 
      group_by(t,age_group) %>% 
      summarise(median=median(simulated_cases,na.rm=TRUE),
                lo1=quantile(simulated_cases,probs=0.025,na.rm=TRUE),
                lo2=quantile(simulated_cases,probs=0.25,na.rm=TRUE),
                up1=quantile(simulated_cases,probs=0.975,na.rm=TRUE),
                up2=quantile(simulated_cases,probs=0.75,na.rm=TRUE)) %>% 
      ungroup()
    
    colnames(cases) <- age_groups
    reported_cases <- cases %>% 
      as_tibble() %>% 
      mutate(t=1:nT) %>% 
      pivot_longer(-t,names_to="age_group",values_to="reported_cases") %>% 
      mutate(age_group=factor(age_group,levels=age_groups))
    
    simulated_cases_summary <- simulated_cases_summary %>% 
      left_join(reported_cases)
    
    plot_compare <- simulated_cases_summary %>% 
      ggplot(aes(x=t)) +
      geom_linerange(aes(ymin=lo1,ymax=up1),col="#fe0332",linewidth=1) +
      geom_linerange(aes(ymin=lo2,ymax=up2),col="#fe0332",linewidth=2) +
      geom_point(aes(y=median),col="#fe0332",size=1.5) +
      geom_point(aes(y=reported_cases),size=1.5) +
      theme_classic(base_size=12) +
      facet_grid(age_group~.) +
      scale_x_continuous(breaks=1:nT) +
      labs(x="Time",y="Cases")
    
    return(list(simulated_cases=simulated_cases,
                simulated_cases_summary=simulated_cases_summary,
                plot=plot_compare))
    
  } else if (dist=="multinomial"){
    
    nT <- ncol(posterior_samples[["lam_t"]])
    nA <- length(amin)
    age_groups <- paste(amin,amax,sep='-')
    
    simulated_cases_list <- list()
    for (i in seq_len(length(posterior_samples[["lam_H"]]))){
      
      lamH <- rep(posterior_samples[["lam_H"]][i],hist_length)
      lam <- posterior_samples[["lam_t"]][i,]
      rho <- 1
      gamma <- posterior_samples[["gamma"]][i]
      simulated_output <- simcases(nT=nT, nA=nA, lamH=lamH, lam=lam, 
                                   rho=rho, gamma=gamma, pop=pop, amin=amin, amax=amax)
      simulated_cases_df <- simulated_output[[model_type]]
      colnames(simulated_cases_df) <- age_groups
      simulated_cases_list[[i]] <- simulated_cases_df %>% 
        as_tibble() %>% 
        mutate(t=1:nT,iter=i) %>% 
        pivot_longer(-c(t,iter),names_to="age_group",values_to="simulated_cases") %>% 
        mutate(age_group=factor(age_group,levels=age_groups)) %>% 
        group_by(iter,t) %>% 
        mutate(simulated_cases_prop=simulated_cases/sum(simulated_cases),
               simulated_cases_prop=ifelse(is.nan(simulated_cases_prop),NA,simulated_cases_prop)) %>% 
        ungroup()
    }
    
    simulated_cases <- bind_rows(simulated_cases_list)
    simulated_cases_summary <- simulated_cases %>% 
      group_by(t,age_group) %>% 
      summarise(median=median(simulated_cases_prop,na.rm=TRUE),
                lo1=quantile(simulated_cases_prop,probs=0.025,na.rm=TRUE),
                lo2=quantile(simulated_cases_prop,probs=0.25,na.rm=TRUE),
                up1=quantile(simulated_cases_prop,probs=0.975,na.rm=TRUE),
                up2=quantile(simulated_cases_prop,probs=0.75,na.rm=TRUE)) %>% 
      ungroup()
    
    colnames(cases) <- age_groups
    reported_cases <- cases %>% 
      as_tibble() %>% 
      mutate(t=1:nT) %>% 
      pivot_longer(-t,names_to="age_group",values_to="reported_cases") %>% 
      mutate(age_group=factor(age_group,levels=age_groups)) %>% 
      group_by(t) %>% 
      mutate(reported_cases_prop=reported_cases/sum(reported_cases)) %>% 
      ungroup()
    
    simulated_cases_summary <- simulated_cases_summary %>% 
      left_join(reported_cases)
    
    plot_compare <- simulated_cases_summary %>% 
      ggplot(aes(x=t)) +
      geom_linerange(aes(ymin=lo1,ymax=up1),col="#fe0332",linewidth=1) +
      geom_linerange(aes(ymin=lo2,ymax=up2),col="#fe0332",linewidth=2) +
      geom_point(aes(y=median),col="#fe0332",size=1.5) +
      geom_point(aes(y=reported_cases_prop),size=1.5) +
      theme_classic(base_size=12) +
      facet_grid(age_group~.) +
      scale_x_continuous(breaks=1:nT) +
      labs(x="Time",y="Cases")
    
    return(list(simulated_cases=simulated_cases,
                simulated_cases_summary=simulated_cases_summary,
                plot=plot_compare))
  
  } else{
    stop("set dist as either Poisson or multinomial")
  }
  
}

#### function to combine all plots for model checks
#### input:
#### posterior_summary: from summarise_posterior
#### simulated_cases: from simulate_cases_from_posterior
#### dist: distribution base: Poisson will compare counts, multinomial will compare proportions
plot_model_summary <- function(posterior_summary,simulated_cases,dist){
  
  if (dist=="Poisson"){
    
    plot1 <- plot_grid(posterior_summary$plot[[3]],posterior_summary$plot[[4]],nrow=1)
    plot2 <- plot_grid(posterior_summary$plot[[1]],posterior_summary$plot[[2]],
                       plot1,ncol=1)
    plot3 <- plot_grid(plot2,simulated_cases$plot,nrow=1)
    
    return(plot3)
    
  } else if(dist=="multinomial"){
    
    plot1 <- plot_grid(posterior_summary$plot[[3]],nrow=1)
    plot2 <- plot_grid(posterior_summary$plot[[1]],posterior_summary$plot[[2]],
                       plot1,ncol=1)
    plot3 <- plot_grid(plot2,simulated_cases$plot,nrow=1)
    
    return(plot3)
    
  } else{
    stop("set dist as either Poisson or multinomial")
  }
  
}

#### function to compare estimates from Poisson and multinomial model
#### input:
#### summary_poisson
#### summary_multinomial
#### params
#### hist_length: for lam_H
compare_estimates <- function(summary_poisson,summary_multinomial,params,hist_length=50){
  
  summary_poisson$dist <- "Poisson"
  summary_multinomial$dist <- "multinomial"
  
  summary_all <- bind_rows(summary_poisson,summary_multinomial)
  
  plot_params <- list()
  for (i in seq_len(length(params))){
    
    summary_filtered <- summary_all %>% 
      filter(str_detect(pars,params[i]))
    
    if (params[i] == "lam_H"){
      
      lam_H_ori <- tibble(pars="lam_H",t=rep(-hist_length:-1,2),
                          dist=rep(c("Poisson","multinomial"),each=hist_length))
      summary_plot <- lam_H_ori %>% left_join(summary_filtered)
      plot_params[[i]] <- summary_plot %>% 
        ggplot(aes(x=t,group=dist,col=dist,fill=dist)) +
        geom_ribbon(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),alpha=0.10) +
        geom_ribbon(aes(ymin=4*`25%`,ymax=4*`75%`),alpha=0.50) +
        geom_line(aes(y=4*`50%`),linewidth=1.5) +
        scale_colour_manual(breaks=c("Poisson","multinomial"),values=c("#D81B60","#004D40")) +
        scale_fill_manual(breaks=c("Poisson","multinomial"),values=c("#D81B60","#004D40")) +
        theme_classic(base_size=12) +
        labs(x="Time",y="4*lam_H")
      
    } else {
      
      params_values <- unique(summary_filtered$pars)
      
      summary_plot <- summary_filtered %>% 
        mutate(pars=factor(pars,levels=params_values))
      
      if (params[i] == "lam_t"){
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=4*`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=c("Poisson","multinomial"),values=c("#D81B60","#004D40")) +
          theme_classic(base_size=12) +
          labs(x="Params",y=paste0("4*",params[[i]]))
      } else{
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=c("Poisson","multinomial"),values=c("#D81B60","#004D40")) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]]) 
      }
    }
  
  }
  
  plot1 <- plot_grid(plot_params[[3]],plot_params[[4]],nrow=1)
  plot2 <- plot_grid(plot_params[[1]],plot_params[[2]],
                     plot1,ncol=1)
  
  return(list(plot_each=plot_params,plot_all=plot2))
  
}

#### function to summarise all params from stan output - this ver with yearly historical lam_H
#### input:
#### stan_output
#### params
#### original_values: named list should be in the same order as params
summarise_posterior_new <- function(stan_output,params,original_values){
  
  if(sum(names(original_values)==params) != length(params)){
    stop("params and original values name not in the same order")
  } else{
    output <- as_tibble(summary(stan_output,params)$summary) %>% 
      mutate(pars=rownames(summary(stan_output,params)$summary)) %>% 
      select(pars,everything())
  }
  
  plot_params <- list()
  
  for(i in seq_len(length(params))){
    
    output_filtered <- output %>% 
      filter(str_detect(pars,params[i]))
    
    params_values <- output_filtered$pars
    
    if (sum(is.na(original_values[[params[i]]]))){
      output_plot <- output_filtered %>% 
        mutate(pars=factor(pars,levels=params_values))
      if (params[i] %in% c("lam_H","lam_t")){
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=pars)) +
          geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),col="#fe0332",linewidth=1) +
          geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),col="#fe0332",linewidth=2) +
          geom_point(aes(y=4*`50%`),col="#fe0332",size=2.5) +
          theme_classic(base_size=12) +
          labs(x="Params",y=paste0("4*",params[[i]]))
      } else{
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=pars)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),col="#fe0332",linewidth=1) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),col="#fe0332",linewidth=2) +
          geom_point(aes(y=`50%`),col="#fe0332",size=2.5) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]])
      }
    } else {
      output_plot <- output_filtered %>% 
        mutate(pars=factor(pars,levels=params_values),
               original_values=original_values[[params[i]]])
      
      if (params[i] %in% c("lam_H","lam_t")){
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=pars)) +
          geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),col="#fe0332",linewidth=1) +
          geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),col="#fe0332",linewidth=2) +
          geom_point(aes(y=4*`50%`),col="#fe0332",size=2.5) +
          geom_point(aes(y=4*original_values),size=2.5) +
          theme_classic(base_size=12) +
          labs(x="Params",y=paste0("4*",params[[i]]))
      } else{
        plot_params[[i]] <- output_plot %>% 
          ggplot(aes(x=pars)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),col="#fe0332",linewidth=1) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),col="#fe0332",linewidth=2) +
          geom_point(aes(y=`50%`),col="#fe0332",size=2.5) +
          geom_point(aes(y=original_values),size=2.5) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]])
      }
    }
  }
  
  return(list(summary=output,plot=plot_params))
}

#### function to get posterior samples from stan fit output - this ver with yearly historical lam_H
#### input:
#### stan_output
#### params
#### sample_size: sample size from posterior samples for further simulations
posterior_samples_output_new <- function(stan_output,params,sample_size){
  
  output <- rstan::extract(stan_output,pars=params)
  
  index_of_gamma <- which(params=="gamma")
  index_of_samples <- sample(1:length(output[[index_of_gamma]]),sample_size)
  
  param_samples_list <- list()
  
  for (i in seq_len(length(params))){
    
    if (!(params[i] %in% c("lam_H","lam_t"))){
      param_samples_list[[i]] <- output[[i]][index_of_samples] 
    } else {
      param_samples_list[[i]] <- matrix(NA,ncol=ncol(output[[i]]),nrow=sample_size)
      for (j in seq_len(ncol(output[[i]]))){
        param_samples_list[[i]][,j] <- output[[i]][index_of_samples,j]
      }
    }
  }
  
  names(param_samples_list) <- params
  
  return(param_samples_list)
}

#### function to simulate cases and compare to original data - this ver with yearly historical lam_H
#### input:
#### posterior_samples: from posterior_samples_output
#### cases: cases data fitted to the model
#### pop
#### amin
#### amax
#### hist_length
#### model_type: model type cases_S or cases_PS
#### dist: distribution base: Poisson will compare counts, multinomial will compare proportions
simulate_cases_from_posterior_new <- function(posterior_samples,
                                              cases,
                                              pop,
                                              amin,
                                              amax,
                                              hist_length,
                                              model_type,
                                              dist){
  
  if (dist=="Poisson"){
    
    nT <- ncol(posterior_samples[["lam_t"]])
    nA <- length(amin)
    age_groups <- paste(amin,amax,sep='-')
    
    simulated_cases_list <- list()
    for (i in seq_len(length(posterior_samples[["gamma"]]))){
      
      lamH <- posterior_samples[["lam_H"]][i,]
      lam <- posterior_samples[["lam_t"]][i,]
      rho <- posterior_samples[["rho"]][i]
      gamma <- posterior_samples[["gamma"]][i]
      simulated_output <- simcases(nT=nT, nA=nA, lamH=lamH, lam=lam, 
                                   rho=rho, gamma=gamma, pop=pop, amin=amin, amax=amax)
      simulated_cases_df <- simulated_output[[model_type]]
      colnames(simulated_cases_df) <- age_groups
      simulated_cases_list[[i]] <- simulated_cases_df %>% 
        as_tibble() %>% 
        mutate(t=1:nT,iter=i) %>% 
        pivot_longer(-c(t,iter),names_to="age_group",values_to="simulated_cases") %>% 
        mutate(age_group=factor(age_group,levels=age_groups))
      
    }
    
    simulated_cases <- bind_rows(simulated_cases_list)
    simulated_cases_summary <- simulated_cases %>% 
      group_by(t,age_group) %>% 
      summarise(median=median(simulated_cases,na.rm=TRUE),
                lo1=quantile(simulated_cases,probs=0.025,na.rm=TRUE),
                lo2=quantile(simulated_cases,probs=0.25,na.rm=TRUE),
                up1=quantile(simulated_cases,probs=0.975,na.rm=TRUE),
                up2=quantile(simulated_cases,probs=0.75,na.rm=TRUE)) %>% 
      ungroup()
    
    colnames(cases) <- age_groups
    reported_cases <- cases %>% 
      as_tibble() %>% 
      mutate(t=1:nT) %>% 
      pivot_longer(-t,names_to="age_group",values_to="reported_cases") %>% 
      mutate(age_group=factor(age_group,levels=age_groups))
    
    simulated_cases_summary <- simulated_cases_summary %>% 
      left_join(reported_cases)
    
    plot_compare <- simulated_cases_summary %>% 
      ggplot(aes(x=t)) +
      geom_linerange(aes(ymin=lo1,ymax=up1),col="#fe0332",linewidth=1) +
      geom_linerange(aes(ymin=lo2,ymax=up2),col="#fe0332",linewidth=2) +
      geom_point(aes(y=median),col="#fe0332",size=1.5) +
      geom_point(aes(y=reported_cases),size=1.5) +
      theme_classic(base_size=12) +
      facet_grid(age_group~.) +
      scale_x_continuous(breaks=1:nT) +
      labs(x="Time",y="Cases")
    
    return(list(simulated_cases=simulated_cases,
                simulated_cases_summary=simulated_cases_summary,
                plot=plot_compare))
    
  } else if (dist=="multinomial"){
    
    nT <- ncol(posterior_samples[["lam_t"]])
    nA <- length(amin)
    age_groups <- paste(amin,amax,sep='-')
    
    simulated_cases_list <- list()
    for (i in seq_len(length(posterior_samples[["gamma"]]))){
      
      lamH <- posterior_samples[["lam_H"]][i,]
      lam <- posterior_samples[["lam_t"]][i,]
      rho <- 1
      gamma <- posterior_samples[["gamma"]][i]
      simulated_output <- simcases(nT=nT, nA=nA, lamH=lamH, lam=lam, 
                                   rho=rho, gamma=gamma, pop=pop, amin=amin, amax=amax)
      simulated_cases_df <- simulated_output[[model_type]]
      colnames(simulated_cases_df) <- age_groups
      simulated_cases_list[[i]] <- simulated_cases_df %>% 
        as_tibble() %>% 
        mutate(t=1:nT,iter=i) %>% 
        pivot_longer(-c(t,iter),names_to="age_group",values_to="simulated_cases") %>% 
        mutate(age_group=factor(age_group,levels=age_groups)) %>% 
        group_by(iter,t) %>% 
        mutate(simulated_cases_prop=simulated_cases/sum(simulated_cases),
               simulated_cases_prop=ifelse(is.nan(simulated_cases_prop),NA,simulated_cases_prop)) %>% 
        ungroup()
    }
    
    simulated_cases <- bind_rows(simulated_cases_list)
    simulated_cases_summary <- simulated_cases %>% 
      group_by(t,age_group) %>% 
      summarise(median=median(simulated_cases_prop,na.rm=TRUE),
                lo1=quantile(simulated_cases_prop,probs=0.025,na.rm=TRUE),
                lo2=quantile(simulated_cases_prop,probs=0.25,na.rm=TRUE),
                up1=quantile(simulated_cases_prop,probs=0.975,na.rm=TRUE),
                up2=quantile(simulated_cases_prop,probs=0.75,na.rm=TRUE)) %>% 
      ungroup()
    
    colnames(cases) <- age_groups
    reported_cases <- cases %>% 
      as_tibble() %>% 
      mutate(t=1:nT) %>% 
      pivot_longer(-t,names_to="age_group",values_to="reported_cases") %>% 
      mutate(age_group=factor(age_group,levels=age_groups)) %>% 
      group_by(t) %>% 
      mutate(reported_cases_prop=reported_cases/sum(reported_cases)) %>% 
      ungroup()
    
    simulated_cases_summary <- simulated_cases_summary %>% 
      left_join(reported_cases)
    
    plot_compare <- simulated_cases_summary %>% 
      ggplot(aes(x=t)) +
      geom_linerange(aes(ymin=lo1,ymax=up1),col="#fe0332",linewidth=1) +
      geom_linerange(aes(ymin=lo2,ymax=up2),col="#fe0332",linewidth=2) +
      geom_point(aes(y=median),col="#fe0332",size=1.5) +
      geom_point(aes(y=reported_cases_prop),size=1.5) +
      theme_classic(base_size=12) +
      facet_grid(age_group~.) +
      scale_x_continuous(breaks=1:nT) +
      labs(x="Time",y="Cases")
    
    return(list(simulated_cases=simulated_cases,
                simulated_cases_summary=simulated_cases_summary,
                plot=plot_compare))
    
  } else{
    stop("set dist as either Poisson or multinomial")
  }
  
}

#### function to compare estimates from Poisson and multinomial model - multiple models
#### input:
#### summary_poisson
#### summary_multinomial
#### params
#### hist_length: for lam_H
compare_estimates_new <- function(summary_list,model_name,params,original_values,hist_length=50){
  
  col_palette <- c("#D81B60","#1E88E5","#FFC107","#004D40")
  
  if (length(summary_list) != length(model_name)){
    stop("the length of the summary_list has to be the same as model_name")
  }
  
  for (i in seq_len(length(summary_list))){
    summary_list[[i]]$dist <- model_name[i]
  }
  
  summary_all <- bind_rows(summary_list)
  
  plot_params <- list()
  for (i in seq_len(length(params))){
    
    if (params[i] == "lam_H"){
      
      summary_filtered1 <- summary_all %>% 
        filter(str_detect(pars,params[i])) %>% 
        filter(pars == "lam_H")
      dist1 <- unique(summary_filtered1$dist)
      lam_H_ori <- tibble(pars="lam_H",t=rep(-hist_length:-1,2),
                          dist=rep(dist1,each=hist_length))
      summary_plot1 <- lam_H_ori %>% left_join(summary_filtered1)
      
      summary_plot2 <- summary_all %>% 
        filter(str_detect(pars,params[i])) %>% 
        filter(pars != "lam_H") %>% 
        mutate(t=parse_number(pars)-(hist_length+1))
      
      original_values_df <- tibble(pars=params[i],
                                   t=rep(-hist_length:-1),
                                   dist="Original",
                                   `50%`=original_values[[params[i]]],
                                   `2.5%`=original_values[[params[i]]],
                                   `25%`=original_values[[params[i]]],
                                   `75%`=original_values[[params[i]]],
                                   `97.5%`=original_values[[params[i]]])
      
      summary_plot <- bind_rows(summary_plot1,summary_plot2,original_values_df)
      length_model <- length(unique(summary_plot$dist))-1
      
      plot_params[[i]] <- summary_plot %>% 
        ggplot(aes(x=t,group=dist,col=dist,fill=dist)) +
        geom_ribbon(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),alpha=0.10) +
        geom_ribbon(aes(ymin=4*`25%`,ymax=4*`75%`),alpha=0.50) +
        geom_line(aes(y=4*`50%`),linewidth=1.5) +
        scale_colour_manual(breaks=c(model_name,"Original"),values=c(col_palette[1:length_model],"black")) +
        scale_fill_manual(breaks=c(model_name,"Original"),values=c(col_palette[1:length_model],"black")) +
        theme_classic(base_size=12) +
        labs(x="Time",y="4*lam_H")
      
    } else {
      
      summary_filtered <- summary_all %>% 
        filter(str_detect(pars,params[i])) 
      
      params_values <- unique(summary_filtered$pars)
      
      summary_plot <- summary_filtered %>% 
        mutate(pars=factor(pars,levels=params_values))
      
      length_model <- length(unique(summary_plot$dist))
      
      if (params[i] == "lam_t"){
        original_values_df <- tibble(pars=params_values,
                                     dist="Original",
                                     `50%`=original_values[[params[i]]],
                                     `2.5%`=original_values[[params[i]]],
                                     `25%`=original_values[[params[i]]],
                                     `75%`=original_values[[params[i]]],
                                     `97.5%`=original_values[[params[i]]])
        
        summary_plot <- bind_rows(summary_plot,original_values_df)
        
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=4*`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=c(model_name,"Original"),values=c(col_palette[1:length_model],"black")) +
          theme_classic(base_size=12) +
          labs(x="Params",y=paste0("4*",params[[i]]))
      } else if (params[i] == "rho"){
        original_values_df <- tibble(pars=params[i],
                                     dist="Original",
                                     `50%`=original_values[[params[i]]],
                                     `2.5%`=original_values[[params[i]]],
                                     `25%`=original_values[[params[i]]],
                                     `75%`=original_values[[params[i]]],
                                     `97.5%`=original_values[[params[i]]])
        
        summary_plot <- bind_rows(summary_plot,original_values_df)
        
        model_name2 <- unique(summary_plot$dist)
        col_index <- which(model_name %in% model_name2)
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=c(model_name2,"Original"),values=c(col_palette[col_index],"black")) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]])
      } else {
        original_values_df <- tibble(pars=params[i],
                                     dist="Original",
                                     `50%`=original_values[[params[i]]],
                                     `2.5%`=original_values[[params[i]]],
                                     `25%`=original_values[[params[i]]],
                                     `75%`=original_values[[params[i]]],
                                     `97.5%`=original_values[[params[i]]])
        
        summary_plot <- bind_rows(summary_plot,original_values_df)
        
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=c(model_name,"Original"),values=c(col_palette[1:length_model],"black")) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]]) 
      }
    }
    
  }
  
  plot1 <- plot_grid(plot_params[[3]],plot_params[[4]],nrow=1)
  plot2 <- plot_grid(plot_params[[1]],plot_params[[2]],
                     plot1,ncol=1)
  
  return(list(plot_each=plot_params,plot_all=plot2))
  
}

#### function to compare estimates from different districts - use multiple historical lambda model
#### input:
#### summary_poisson
#### summary_multinomial
#### params
#### hist_length: for lam_H
compare_estimates_district <- function(summary_list,district_name,params,lam_t_year,
                                       hist_length=50){
  
  col_palette <- c("#332288","#117733","#44AA99","#88CCEE","#DDCC77","#CC6677","#AA4499","#882255")
  
  if (length(summary_list) != length(district_name)){
    stop("the length of the summary_list has to be the same as district_name")
  }
  
  for (i in seq_len(length(summary_list))){
    summary_list[[i]]$dist <- district_name[i]
  }
  
  summary_all <- bind_rows(summary_list)
  
  plot_params <- list()
  for (i in seq_len(length(params))){
    
    if (params[i] == "lam_H"){
      
      summary_plot <- summary_all %>% 
        filter(str_detect(pars,params[i])) %>% 
        mutate(t=parse_number(pars)-(hist_length+1))
    
      length_model <- length(unique(summary_plot$dist))
      
      plot_params[[i]] <- summary_plot %>% 
        ggplot(aes(x=t,group=dist,col=dist,fill=dist)) +
        geom_ribbon(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),alpha=0.10) +
        geom_ribbon(aes(ymin=4*`25%`,ymax=4*`75%`),alpha=0.50) +
        geom_line(aes(y=4*`50%`),linewidth=1.5) +
        scale_colour_manual(breaks=district_name,values=col_palette[1:length_model]) +
        scale_fill_manual(breaks=district_name,values=col_palette[1:length_model]) +
        theme_classic(base_size=12) +
        labs(x="Time",y="4*lam_H")
      
    } else {
      
      summary_filtered <- summary_all %>% 
        filter(str_detect(pars,params[i])) 
      
      params_values <- unique(summary_filtered$pars)
      
      summary_plot <- summary_filtered %>% 
        mutate(pars=factor(pars,levels=params_values))
      
      length_model <- length(unique(summary_plot$dist))
      
      if (params[i] == "lam_t"){
        summary_plot <- summary_plot %>% left_join(lam_t_year)
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=year,group=dist,col=dist)) +
          geom_linerange(aes(ymin=4*`2.5%`,ymax=4*`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=4*`25%`,ymax=4*`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=4*`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=district_name,values=col_palette[1:length_model]) +
          theme_classic(base_size=12) +
          labs(x="Params",y=paste0("4*",params[[i]]))
      } else {
        plot_params[[i]] <- summary_plot %>% 
          ggplot(aes(x=pars,group=dist,col=dist)) +
          geom_linerange(aes(ymin=`2.5%`,ymax=`97.5%`),linewidth=1,position=position_dodge(width=0.5)) +
          geom_linerange(aes(ymin=`25%`,ymax=`75%`),linewidth=2,position=position_dodge(width=0.5)) +
          geom_point(aes(y=`50%`),size=2.5,position=position_dodge(width=0.5)) +
          scale_colour_manual(breaks=district_name,values=col_palette[1:length_model]) +
          theme_classic(base_size=12) +
          labs(x="Params",y=params[[i]]) 
      }
    }
    
  }
  
  plot1 <- plot_grid(plot_params[[3]],plot_params[[4]],nrow=1)
  plot2 <- plot_grid(plot_params[[1]],plot_params[[2]],
                     plot1,ncol=1)
  
  return(list(plot_each=plot_params,plot_all=plot2))
  
}

theme_bare <- theme(
  axis.line = element_blank(), 
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(),
  legend.text=element_text(size=9),
  legend.title=element_text(size=10),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "gray", fill=NA, linewidth=0.5)
)

theme_bare2 <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(), 
  axis.text.y = element_blank(),
  axis.ticks = element_blank(), 
  axis.title.x = element_blank(), 
  axis.title.y = element_blank(),
  legend.text=element_text(size=9),
  legend.title=element_text(size=10),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "white", fill=NA, linewidth=0.5)
)

# https://rpubs.com/koundy/71792

theme_Publication <- function(base_size=12, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold", hjust = 0.5),
                                      # size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.5, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(1.5,1.0,1.0,1.0),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

#### functions to extract stan output
extract_lam_H <- function(stan_fit,cases_all,province,
                          hist_length,n_samples,prob=c(0.025,0.5,0.975)){
  
  extract_stanfit <- extract(stan_fit,pars="lam_H") 
  
  cases_province <- cases_all %>% filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  lam_H_prov_list <- list()
  
  year_province <- cases_province %>% pull(year) %>% min()
  
  for (j in seq_len(length(admin2_province))){
    lam_H_prov_list[[j]] <- as_tibble(extract_stanfit$lam_H[,j,]) 
    colnames(lam_H_prov_list[[j]]) <- year_province + -1:-hist_length
    lam_H_prov_list[[j]] <- lam_H_prov_list[[j]] %>% 
      mutate(admin2=admin2_province[j],
             rep=1:n_samples) %>% 
      select(admin2,rep,everything()) %>% 
      pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_H")
  }
  
  lam_H_prov <- bind_rows(lam_H_prov_list)
  lam_H_output <- lam_H_prov %>% 
    group_by(admin2,t) %>% 
    summarise(as_tibble_row(quantile(lam_H,probs=prob))) %>% 
    dplyr::select(admin2,t,median_foi=`50%`,lo_foi=`2.5%`,up_foi=`97.5%`)
  
  return(lam_H_output)
}

extract_lam_t <- function(stan_fit,cases_all,province,n_samples,prob=c(0.025,0.5,0.975)){
  
  extract_stanfit <- extract(stan_fit,pars="lam_t") 
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  lam_t_prov_list <- list()
  
  year_first <- cases_province %>% pull(year) %>% min()
  
  for (j in seq_len(length(admin2_province))){
    lam_t_prov_list[[j]] <- as_tibble(extract_stanfit$lam_t[,j,]) 
    colnames(lam_t_prov_list[[j]]) <- year_first:2023
    lam_t_prov_list[[j]] <- lam_t_prov_list[[j]] %>% 
      mutate(admin2=admin2_province[j],
             rep=1:n_samples) %>% 
      select(admin2,rep,everything()) %>% 
      pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_t")
  }
  
  lam_t_prov <- bind_rows(lam_t_prov_list)
  lam_t_output <- lam_t_prov %>% 
    group_by(admin2,t) %>% 
    summarise(as_tibble_row(quantile(lam_t,probs=prob))) %>% 
    dplyr::select(admin2,t,median_foi=`50%`,lo_foi=`2.5%`,up_foi=`97.5%`)
  
  return(lam_t_output)
}

extract_rho <- function(stan_fit,cases_all,province,
                        n_samples,prob=c(0.025,0.5,0.975)){
  
  extract_stanfit <- extract(stan_fit,pars="rho") 
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  rho_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    rho_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                 rep=1:n_samples,
                                 rho=extract_stanfit$rho[,j])
  }
  
  rho_prov <- bind_rows(rho_prov_list)
  rho_output <- rho_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(rho,probs=prob))) %>% 
    dplyr::select(admin2,median_rho=`50%`,lo_rho=`2.5%`,up_rho=`97.5%`)
  
  return(rho_output)
}

extract_gamma <- function(stan_fit,cases_all,province,
                          n_samples,prob=c(0.025,0.5,0.975)){
  
  extract_stanfit <- extract(stan_fit,pars="gamma") 
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  gamma_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    gamma_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                   rep=1:n_samples,
                                   gamma=extract_stanfit$gamma[,j])
  }
  
  gamma_prov <- bind_rows(gamma_prov_list)
  gamma_output <- gamma_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(gamma,probs=prob))) %>% 
    dplyr::select(admin2,median_gamma=`50%`,lo_gamma=`2.5%`,up_gamma=`97.5%`)
  
  return(gamma_output)
}

extract_phi <- function(stan_fit,cases_all=cases_all,province,
                        n_samples,prob=c(0.025,0.5,0.975)){
  
  extract_stanfit <- extract(stan_fit,pars="phi") 
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  phi_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    phi_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                 rep=1:n_samples,
                                 phi=extract_stanfit$phi[,j])
  }
  
  phi_prov <- bind_rows(phi_prov_list)
  phi_output <- phi_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(phi,probs=prob))) %>% 
    dplyr::select(admin2,median_phi=`50%`,lo_phi=`2.5%`,up_phi=`97.5%`)
  
  return(phi_output)
}

#### function to extract province-level log-likelihoods
extract_province_loglik <- function(log_lik_hier, province_districts, nT, nA) {
  n_iter <- dim(log_lik_hier)[1]
  province_loglik <- array(NA, dim = c(n_iter, length(province_districts) * nT * nA))
  
  curr_idx <- 1
  for(d in province_districts) {
    start_idx <- ((d-1) * nT * nA) + 1
    end_idx <- d * nT * nA
    province_loglik[, curr_idx:(curr_idx + nT*nA - 1)] <- log_lik_hier[, start_idx:end_idx]
    curr_idx <- curr_idx + nT * nA
  }
  return(province_loglik)
}

extract_lam_H_indiv <- function(stan_fit,cases_all,province,
                                hist_length,n_samples,prob=c(0.025,0.5,0.975)){
  
  cases_province <- cases_all %>% filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  lam_H_prov_list <- list()
  
  year_province <- cases_province %>% pull(year) %>% min()
  
  for (j in seq_len(length(admin2_province))){
    extract_stanfit <- extract(stan_fit[[j]],pars="lam_H")
    lam_H_prov_list[[j]] <- as_tibble(extract_stanfit$lam_H) 
    colnames(lam_H_prov_list[[j]]) <- year_province + -1:-hist_length
    lam_H_prov_list[[j]] <- lam_H_prov_list[[j]] %>% 
      mutate(admin2=admin2_province[j],
             rep=1:n_samples) %>% 
      select(admin2,rep,everything()) %>% 
      pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_H")
  }
  
  lam_H_prov <- bind_rows(lam_H_prov_list)
  lam_H_output <- lam_H_prov %>% 
    group_by(admin2,t) %>% 
    summarise(as_tibble_row(quantile(lam_H,probs=prob))) %>% 
    dplyr::select(admin2,t,median_foi=`50%`,lo_foi=`2.5%`,up_foi=`97.5%`)
  
  return(lam_H_output)
}

extract_lam_t_indiv <- function(stan_fit,cases_all,province,n_samples,prob=c(0.025,0.5,0.975)){
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  lam_t_prov_list <- list()
  
  year_first <- cases_province %>% pull(year) %>% min()
  
  for (j in seq_len(length(admin2_province))){
    extract_stanfit <- extract(stan_fit[[j]],pars="lam_t") 
    lam_t_prov_list[[j]] <- as_tibble(extract_stanfit$lam_t) 
    colnames(lam_t_prov_list[[j]]) <- year_first:2023
    lam_t_prov_list[[j]] <- lam_t_prov_list[[j]] %>% 
      mutate(admin2=admin2_province[j],
             rep=1:n_samples) %>% 
      select(admin2,rep,everything()) %>% 
      pivot_longer(-c(admin2,rep),names_to="t",values_to="lam_t")
  }
  
  lam_t_prov <- bind_rows(lam_t_prov_list)
  lam_t_output <- lam_t_prov %>% 
    group_by(admin2,t) %>% 
    summarise(as_tibble_row(quantile(lam_t,probs=prob))) %>% 
    dplyr::select(admin2,t,median_foi=`50%`,lo_foi=`2.5%`,up_foi=`97.5%`)
  
  return(lam_t_output)
}

extract_rho_indiv <- function(stan_fit,cases_all,province,
                              n_samples,prob=c(0.025,0.5,0.975)){
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  rho_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    extract_stanfit <- extract(stan_fit[[j]],pars="rho") 
    rho_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                 rep=1:n_samples,
                                 rho=extract_stanfit$rho)
  }
  
  rho_prov <- bind_rows(rho_prov_list)
  rho_output <- rho_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(rho,probs=prob))) %>% 
    dplyr::select(admin2,median_rho=`50%`,lo_rho=`2.5%`,up_rho=`97.5%`)
  
  return(rho_output)
}

extract_gamma_indiv <- function(stan_fit,cases_all,province,
                                n_samples,prob=c(0.025,0.5,0.975)){
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  gamma_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    extract_stanfit <- extract(stan_fit[[j]],pars="gamma") 
    gamma_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                   rep=1:n_samples,
                                   gamma=extract_stanfit$gamma)
  }
  
  gamma_prov <- bind_rows(gamma_prov_list)
  gamma_output <- gamma_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(gamma,probs=prob))) %>% 
    dplyr::select(admin2,median_gamma=`50%`,lo_gamma=`2.5%`,up_gamma=`97.5%`)
  
  return(gamma_output)
}

extract_phi_indiv <- function(stan_fit,cases_all=cases_all,province,
                              n_samples,prob=c(0.025,0.5,0.975)){
  
  cases_province <- cases_all %>% 
    filter(admin1 == province)
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  phi_prov_list <- list()
  
  for (j in seq_len(length(admin2_province))){
    extract_stanfit <- extract(stan_fit[[j]],pars="phi") 
    phi_prov_list[[j]] <- tibble(admin2=admin2_province[j],
                                 rep=1:n_samples,
                                 phi=extract_stanfit$phi)
  }
  
  phi_prov <- bind_rows(phi_prov_list)
  phi_output <- phi_prov %>% 
    group_by(admin2) %>% 
    summarise(as_tibble_row(quantile(phi,probs=prob))) %>% 
    dplyr::select(admin2,median_phi=`50%`,lo_phi=`2.5%`,up_phi=`97.5%`)
  
  return(phi_output)
}



