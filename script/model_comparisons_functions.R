#### compare with seroprevalence
compare_seroprev_2014 <- function(model_names,catalytic_fit_data,province_idx,year_first){
  
  if (province_idx == 1){
    admin2_index <- c(2,3,5) 
  } else if (province_idx == 2){
    admin2_index <- c(1,4,6,9,16,21,23)  
  } else{
    stop("provide province_idx as either 1 or 2")
  }
  
  cases_province <- cases_all %>% 
    filter(admin1 == provinces[province_idx])
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  
  #### read catalytic fit
  catalytic_fit <- readRDS(paste0("output/stan_fit/",catalytic_fit_data))
  
  #### read model fit
  lam_H_samples <- list()
  if (grepl("non-hierarchical", model_names, fixed = TRUE)){
    for (i in admin2_index){
      lam_H_samples[[i]] <- extract(catalytic_fit[[i]],pars="lam_H")$lam_H
    }
  } else {
    for (i in admin2_index){
      lam_H_samples[[i]] <- extract(catalytic_fit,pars="lam_H")$lam_H[,i,]
    }
  }
  
  #### calculate seroprev 2014 from model results
  #### simulate seroprev up to 2013 historical lam_H to get 2014 seroprevalence
  seroprev_2014 <- list()
  for (i in admin2_index){
    
    seroprev_estim <- list()
    
    for (j in seq_len(n_samples)){
      
      len_lam_H <- length(lam_H_samples[[i]][j,])
      diff_to_2013 <- year_first - 1 - 2013
      idx_2013 <- len_lam_H - diff_to_2013
      lam_H_sim <- lam_H_samples[[i]][j,1:idx_2013] # up to 2013
      
      seroprev_estim_vec <- (1 - simseroprev_hist(lam_H_sim)$susc)[1:18]
      seroprev_estim[[j]] <- tibble(rep=j,age=1:18,seroprev=seroprev_estim_vec)
      
    }
    
    seroprev_2014[[i]] <- bind_rows(seroprev_estim)
    
  }
  
  seroprev_summary_2014 <- list()
  seroprev_overall_summary_2014 <- list()
  for (i in admin2_index){
    
    seroprev_summary_2014[[i]] <- seroprev_2014[[i]] %>% 
      group_by(age) %>% 
      summarise(median_seroprev=median(seroprev),
                lo_seroprev=quantile(seroprev,probs=0.025),
                up_seroprev=quantile(seroprev,probs=0.975))
    
    seroprev_overall_summary_2014[[i]] <- seroprev_2014[[i]] %>% 
      group_by(rep) %>%
      summarise(seroprev=mean(seroprev)) %>% 
      ungroup() %>% 
      summarise(median_seroprev=median(seroprev),
                lo_seroprev=quantile(seroprev,probs=0.025),
                up_seroprev=quantile(seroprev,probs=0.975))
    
  }
  
  seroprev_data_select <- seroprev_data %>% 
    filter(admin2 %in% admin2_province[admin2_index])
  
  if (province_idx == 1){
    seroprev_summary_admin2 <- bind_rows(seroprev_summary_2014[[2]] %>%
                                           mutate(admin2=admin2_province[2]),
                                         seroprev_summary_2014[[3]] %>%
                                           mutate(admin2=admin2_province[3]),
                                         seroprev_summary_2014[[5]] %>%
                                           mutate(admin2=admin2_province[5])) %>%
      mutate(type=model_names) %>%
      bind_rows(seroprev_data_select)
    
  } else{
    seroprev_summary_admin2 <- bind_rows(seroprev_summary_2014[[1]] %>%
                                           mutate(admin2=admin2_province[1]),
                                         seroprev_summary_2014[[4]] %>%
                                           mutate(admin2=admin2_province[4]),
                                         seroprev_summary_2014[[6]] %>%
                                           mutate(admin2=admin2_province[6]),
                                         seroprev_summary_2014[[9]] %>%
                                           mutate(admin2=admin2_province[9]),
                                         seroprev_summary_2014[[16]] %>%
                                           mutate(admin2=admin2_province[16]),
                                         seroprev_summary_2014[[21]] %>%
                                           mutate(admin2=admin2_province[21]),
                                         seroprev_summary_2014[[23]] %>%
                                           mutate(admin2=admin2_province[23])) %>%
      mutate(type=model_names) %>%
      bind_rows(seroprev_data_select)
    
  }
  
  plot_output <- seroprev_summary_admin2 %>% 
    ggplot(aes(x=age,colour=type)) +
    geom_point(aes(y=median_seroprev),position=position_dodge(width=1)) +
    geom_linerange(aes(ymin=lo_seroprev,ymax=up_seroprev),position=position_dodge(width=1)) +
    theme_bw() +
    facet_wrap(.~admin2) +
    theme(legend.position="bottom") +
    ylim(0,1) +
    scale_x_continuous(breaks=1:18)
  
  if (province_idx == 1){
    seroprev_summary_admin2 <- bind_rows(seroprev_summary_2014[[2]] %>%
                                           mutate(admin2=admin2_province[2]),
                                         seroprev_summary_2014[[3]] %>%
                                           mutate(admin2=admin2_province[3]),
                                         seroprev_summary_2014[[5]] %>%
                                           mutate(admin2=admin2_province[5])) %>%
      mutate(type=model_names)
    
    seroprev_overall_summary_admin2 <- bind_rows(seroprev_overall_summary_2014[[2]] %>%
                                                   mutate(admin2=admin2_province[2]),
                                                 seroprev_overall_summary_2014[[3]] %>%
                                                   mutate(admin2=admin2_province[3]),
                                                 seroprev_overall_summary_2014[[5]] %>%
                                                   mutate(admin2=admin2_province[5])) %>%
      mutate(type=model_names)
  } else{
    seroprev_summary_admin2 <- bind_rows(seroprev_summary_2014[[1]] %>%
                                           mutate(admin2=admin2_province[1]),
                                         seroprev_summary_2014[[4]] %>%
                                           mutate(admin2=admin2_province[4]),
                                         seroprev_summary_2014[[6]] %>%
                                           mutate(admin2=admin2_province[6]),
                                         seroprev_summary_2014[[9]] %>%
                                           mutate(admin2=admin2_province[9]),
                                         seroprev_summary_2014[[16]] %>%
                                           mutate(admin2=admin2_province[16]),
                                         seroprev_summary_2014[[21]] %>%
                                           mutate(admin2=admin2_province[21]),
                                         seroprev_summary_2014[[23]] %>%
                                           mutate(admin2=admin2_province[23])) %>%
      mutate(type=model_names)
    
    seroprev_overall_summary_admin2 <- bind_rows(seroprev_overall_summary_2014[[1]] %>%
                                                   mutate(admin2=admin2_province[1]),
                                                 seroprev_overall_summary_2014[[4]] %>%
                                                   mutate(admin2=admin2_province[4]),
                                                 seroprev_overall_summary_2014[[6]] %>%
                                                   mutate(admin2=admin2_province[6]),
                                                 seroprev_overall_summary_2014[[9]] %>%
                                                   mutate(admin2=admin2_province[9]),
                                                 seroprev_overall_summary_2014[[16]] %>%
                                                   mutate(admin2=admin2_province[16]),
                                                 seroprev_overall_summary_2014[[21]] %>%
                                                   mutate(admin2=admin2_province[21]),
                                                 seroprev_overall_summary_2014[[23]] %>%
                                                   mutate(admin2=admin2_province[23])) %>%
      mutate(type=model_names)
  }
  
  return(list(plot=plot_output,
              table=seroprev_summary_admin2,
              table2=seroprev_overall_summary_admin2)
  )
  
}

#### compare parameter estimates
#### lam_H
compare_lam_H <- function(model_names,province_idx,sero_only=1){
  
  summary_lam_H_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_lam_H_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_H_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  if (sero_only==1){
    
    if (province_idx == 1){
      admin2_index <- c(2,3,5) 
    } else if (province_idx == 2){
      admin2_index <- c(1,4,6,9,16,21,23)  
    } else{
      stop("provide province_idx as either 1 or 2")
    }
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    admin2_province <- admin2_province[admin2_index]
    
  } else {
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
  }
  
  lambda_estim <- readRDS("output/model_output/lambda_estim.rds")
  
  lambda_estim_select <- lambda_estim %>% 
    filter(admin2 %in% admin2_province[admin2_index]) %>% 
    mutate(t=2014)
  
  summary_lam_H <- bind_rows(summary_lam_H_model) %>% 
    filter(admin2 %in% admin2_province) %>% 
    mutate(t=as.numeric(t)) %>% 
    filter(t>=1990)
  
  return(summary_lam_H %>% 
           ggplot() +
           geom_line(aes(x=t,y=median_foi*4,colour=model,group=model),linewidth=1) +
           geom_ribbon(aes(x=t,ymin=lo_foi*4,ymax=up_foi*4,fill=model,group=model),alpha=0.25) +
           geom_point(data=lambda_estim_select,aes(x=t,y=lambda_median),colour="blue",size=2.5) +
           geom_linerange(data=lambda_estim_select,aes(x=t,ymin=lambda_lo,ymax=lambda_up),colour="blue") +
           facet_wrap(.~admin2) +
           theme_Publication(base_size=11) +
           scale_x_continuous(breaks=seq(1990,2020,5))
  )
  
  
}

#### lam_t
compare_lam_t <- function(model_names,province_idx,sero_only=1){
  
  summary_lam_t_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_lam_t_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_t_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  if (sero_only==1){
    
    if (province_idx == 1){
      admin2_index <- c(2,3,5) 
    } else if (province_idx == 2){
      admin2_index <- c(1,4,6,9,16,21,23)  
    } else{
      stop("provide province_idx as either 1 or 2")
    }
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    admin2_province <- admin2_province[admin2_index]
    
  } else {
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
  }
  
  lambda_estim <- readRDS("output/model_output/lambda_estim.rds")
  
  lambda_estim_select <- lambda_estim %>% 
    filter(admin2 %in% admin2_province[admin2_index]) %>% 
    mutate(t=2014)
  
  summary_lam_t <- bind_rows(summary_lam_t_model) %>% 
    filter(admin2 %in% admin2_province) %>% 
    mutate(t=as.numeric(t)) %>% 
    filter(t>=1990)
  
  return(summary_lam_t %>% 
           ggplot() +
           geom_point(aes(x=t,y=median_foi*4,colour=model,group=model),linewidth=1,
                      position=position_dodge(width=1)) +
           geom_linerange(aes(x=t,ymin=lo_foi*4,ymax=up_foi*4,colour=model,group=model),
                          position=position_dodge(width=1)) +
           geom_point(data=lambda_estim_select,aes(x=t,y=lambda_median),colour="blue") +
           geom_linerange(data=lambda_estim_select,aes(x=t,ymin=lambda_lo,ymax=lambda_up),colour="blue") +
           facet_wrap(.~admin2) +
           theme_Publication(base_size=11) +
           scale_x_continuous(breaks=seq(1990,2020,5))
  )
  
  
}

#### rho
compare_rho <- function(model_names,province_idx){
  
  summary_rho_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_rho_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_rho_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  summary_rho <- bind_rows(summary_rho_model)
  
  return(list(plot=summary_rho %>% 
           ggplot() +
           geom_point(aes(x=admin2,y=median_rho,colour=model,group=model),
                      position=position_dodge(width=1)) +
           geom_linerange(aes(x=admin2,ymin=lo_rho,ymax=up_rho,colour=model,group=model),
                          position=position_dodge(width=1)) +
           theme_Publication(base_size=11) +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
           labs(x="Districts",y=expression(rho)) +
           scale_colour_colorblind(),
           table=summary_rho)
  )
  
  
}

#### gamma
compare_gamma <- function(model_names,province_idx){
  
  summary_gamma_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_gamma_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_gamma_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  summary_gamma <- bind_rows(summary_gamma_model)
  
  return(list(plot=summary_gamma %>% 
           ggplot() +
           geom_point(aes(x=admin2,y=median_gamma,colour=model,group=model),
                      position=position_dodge(width=1)) +
           geom_linerange(aes(x=admin2,ymin=lo_gamma,ymax=up_gamma,colour=model,group=model),
                          position=position_dodge(width=1)) +
           theme_Publication(base_size=11) +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
           labs(x="Districts",y=expression(gamma)) +
           scale_colour_colorblind(),
           table=summary_gamma)
  )
  
  
}

#### all lambda
compare_all_lambda <- function(model_names,province_idx,sero_only=1){
  
  summary_lam_H_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_lam_H_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_H_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  summary_lam_t_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_lam_t_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_t_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  if (sero_only==1){
    
    if (province_idx == 1){
      admin2_index <- c(2,3,5) 
    } else if (province_idx == 2){
      admin2_index <- c(1,4,6,9,16,21,23)  
    } else{
      stop("provide province_idx as either 1 or 2")
    }
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    admin2_province <- admin2_province[admin2_index]
    
  } else {
    
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[province_idx])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
  }
  
  lambda_estim <- readRDS("output/model_output/lambda_estim.rds")
  
  lambda_estim_select <- lambda_estim %>% 
    filter(admin2 %in% admin2_province) %>% 
    mutate(t=2014) %>% 
    mutate(model="2014 survey")
  
  summary_lam_H <- bind_rows(summary_lam_H_model) %>% 
    filter(admin2 %in% admin2_province) %>% 
    mutate(t=as.numeric(t)) %>% 
    filter(t>=2000)
  
  summary_lam_t <- bind_rows(summary_lam_t_model) %>% 
    filter(admin2 %in% admin2_province) %>% 
    mutate(t=as.numeric(t))
  
  return(list(plot=ggplot() +
                geom_line(data=summary_lam_H,aes(x=t,y=median_foi*4,colour=model,group=model),linewidth=1) +
                geom_ribbon(data=summary_lam_H,aes(x=t,ymin=lo_foi*4,ymax=up_foi*4,fill=model,group=model),alpha=0.25) +
                geom_line(data=summary_lam_t,aes(x=t,y=median_foi*4,colour=model,group=model),linewidth=1) +
                geom_ribbon(data=summary_lam_t,aes(x=t,ymin=lo_foi*4,ymax=up_foi*4,fill=model,group=model),
                            alpha=0.25) +
                geom_point(data=lambda_estim_select,aes(x=t,y=lambda_median,colour=model,fill=model)) +
                geom_linerange(data=lambda_estim_select,aes(x=t,ymin=lambda_lo,ymax=lambda_up,colour=model,fill=model)) +
                facet_wrap(.~admin2,ncol=4) +
                theme_Publication(base_size=11) +
                scale_colour_colorblind() +
                scale_fill_colorblind() +
                scale_x_continuous(breaks=seq(1990,2020,5)) +
                labs(x="Year",y=expression("4"*lambda),colour=NULL,fill=NULL),
              table=summary_lam_t,
              table2=summary_lam_H)
  )
  
}

#### compare elpd
#### province
compare_elpd_province <- function(model_names,province_idx){
  
  loo_models <- list()
  
  for (i in seq_len(length(model_names))){
    loo_models[[i]] <- 
      readRDS(paste0("output/model_output/log_likelihood/loo_province_",province_idx,"_",model_names[i],".rds"))
  }
  
  elpd_models <- list()
  for (i in seq_len(length(model_names))){
    elpd_models[[i]] <- tibble(model=model_names[i],
                               elpd_mean=loo_models[[i]]$estimates[1,1],
                               elpd_lo=loo_models[[i]]$estimates[1,1]-loo_models[[i]]$estimates[1,2],
                               elpd_up=loo_models[[i]]$estimates[1,1]+loo_models[[i]]$estimates[1,2])
    
  }
  elpd_models <- bind_rows(elpd_models) %>% 
    mutate(model=factor(model,levels=model_name_all))
  
  if (province_idx==1){
    n_admin2 = 6
  } else if (province_idx==2){
    n_admin2 = 27
  }
  
  return(
    elpd_models %>% 
      ggplot(aes(x=model,col=model)) +
      geom_point(aes(y=elpd_mean/n_admin2)) +
      geom_linerange(aes(ymin=elpd_lo/n_admin2,ymax=elpd_up/n_admin2)) +
      scale_colour_colorblind() +
      theme_Publication(base_size=11) +
      theme(legend.position = "none") +
      labs(x="Model",y="Province-level ELPD")
  )
  
}

#### district
compare_elpd_district <- function(model_names,province_idx){
  
  loo_models <- list()
  
  for (i in seq_len(length(model_names))){
    loo_models[[i]] <- 
      readRDS(paste0("output/model_output/log_likelihood/loo_district_",province_idx,"_",model_names[i],".rds"))
  }
  
  elpd_models <- list()
  for (i in seq_len(length(model_names))){
    elpd_district <- list()
    for (j in seq_len(length(loo_models[[i]]))){
      elpd_district[[j]] <- tibble(model=paste0("Model ",i), # model_names[i],
                                   district=j,
                                   elpd_mean=loo_models[[i]][[j]]$estimates[1,1])
    }
    
    elpd_models[[i]] <- bind_rows(elpd_district)
  }
  elpd_models <- bind_rows(elpd_models) %>% 
    # mutate(model=factor(model,levels=model_name_all))
    mutate(model=factor(model))
  
  return(
    list(plot = elpd_models %>% 
      ggplot(aes(x=model,fill=model)) +
      geom_boxplot(aes(y=elpd_mean)) +
      scale_fill_colorblind(labels=model_names) +
      theme_Publication(base_size=11) +
      theme(legend.position = "none") +
      labs(x="Model",y="District-level ELPD"),
    data = elpd_models)
  )
  
}

#### correlations with incidence rate
correlation_incidence_rate <- function(model_names,province_idx,params=c("lam_t","rho","gamma")){
  
  if (params=="lam_t"){
    
    summary_lam_t_model <- list()
    
    for (i in seq_len(length(model_names))){
      summary_lam_t_model[[i]] <- 
        readRDS(paste0("output/model_output/parameters/fit_lam_t_",province_idx,"_",model_names[i],".rds"))
      
    }
    
    incidence_rate <- cases_all %>% 
      group_by(idadmin1,idadmin2,admin1,admin2,year) %>% 
      summarise(reported_cases=sum(reported_cases),
                pop=mean(pop)) %>% 
      ungroup() %>% 
      mutate(incidence_rate=reported_cases/pop*100000,
             t=as.character(year))
    
    summary_lam_t <- bind_rows(summary_lam_t_model) %>% 
      left_join(incidence_rate)
    
    corr_lam_t <- summary_lam_t %>% 
      group_by(model,admin2) %>% 
      summarise(correlation=cor(median_foi,incidence_rate))
    
    return(
      list(
        plot1=summary_lam_t %>% 
          ggplot(aes(x=incidence_rate,y=median_foi*4,colour=admin2)) +
          geom_point() +
          geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
          # geom_linerange(aes(ymin=lo_foi*4,ymax=up_foi*4)) +
          facet_wrap(.~model) +
          theme_Publication(base_size=11) +
          labs(x="Incidence rate per 100,000",y=expression("4"*lambda)),
        
        plot2=corr_lam_t %>% 
          ggplot(aes(y=correlation,x=model,fill=model)) +
          geom_boxplot() +
          theme_Publication(base_size=11) +
          scale_fill_colorblind() +
          theme(legend.position = "none") +
          ylim(0.5,1) +
          labs(x="Model",y="Incidence rate & FOI\ncorrelations (district-level)")
      )
    )
    
    
  } else if(params=="rho"){
    
    rho_model <- list()
    
    for (i in seq_len(length(model_names))){
      rho_model[[i]] <- 
        readRDS(paste0("output/model_output/parameters/fit_rho_",province_idx,"_",model_names[i],".rds"))
      
    }
    
    incidence_rate <- cases_all %>% 
      group_by(idadmin1,idadmin2,admin1,admin2) %>% 
      summarise(reported_cases=mean(reported_cases),
                pop=mean(pop)) %>% 
      ungroup() %>% 
      mutate(incidence_rate=reported_cases/pop*100000)
    
    summary_rho <- bind_rows(rho_model) %>% 
      left_join(incidence_rate)
    
    return(
      summary_rho %>% 
        ggplot(aes(x=incidence_rate,y=median_rho)) +
        geom_point() +
        # geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
        geom_linerange(aes(ymin=lo_rho,ymax=up_rho)) +
        facet_wrap(.~model) +
        theme_Publication(base_size=11) +
        labs(x="Incidence rate per 100,000",y=expression(rho)) +
        ylim(0,0.2) +
        xlim(0,180)
      
    )
    
  } else if(params=="gamma"){
    
    gamma_model <- list()
    
    for (i in seq_len(length(model_names))){
      gamma_model[[i]] <- 
        readRDS(paste0("output/model_output/parameters/fit_gamma_",province_idx,"_",model_names[i],".rds"))
      
    }
    
    incidence_rate <- cases_all %>% 
      group_by(idadmin1,idadmin2,admin1,admin2) %>% 
      summarise(reported_cases=mean(reported_cases),
                pop=mean(pop)) %>% 
      ungroup() %>% 
      mutate(incidence_rate=reported_cases/pop*100000)
    
    summary_gamma <- bind_rows(gamma_model) %>% 
      left_join(incidence_rate)
    
    return(
      summary_gamma %>% 
        ggplot(aes(x=incidence_rate,y=median_gamma)) +
        geom_point() +
        # geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
        geom_linerange(aes(ymin=lo_gamma,ymax=up_gamma)) +
        facet_wrap(.~model) +
        theme_Publication(base_size=11) +
        labs(x="Incidence rate per 100,000",y=expression(gamma)) +
        ylim(0,1) +
        xlim(0,180)
      
    )
    
  } else{
    stop("choose params as either 'lam_t', 'rho' or 'gamma'")
  }
  
}

#### combine provinces
correlation_incidence_rate_combine <- function(model_names,province_idx=c(1,2),params=c("rho","gamma")){
  
  if(params=="rho"){
    
    rho_model <- list()
    
    for (i in seq_len(length(model_names))){
      rho_province <- list()
      for (j in province_idx){
        rho_province[[j]] <- 
          readRDS(paste0("output/model_output/parameters/fit_rho_",j,"_",model_names[i],".rds"))
      }
      rho_model[[i]] <- bind_rows(rho_province)
    }
    
    incidence_rate <- cases_all %>% 
      filter(year > 2017) %>% 
      group_by(idadmin1,idadmin2,admin1,admin2) %>% 
      summarise(reported_cases=mean(reported_cases),
                pop=mean(pop)) %>% 
      ungroup() %>% 
      mutate(incidence_rate=reported_cases/pop*100000)
    
    summary_rho <- bind_rows(rho_model) %>% 
      left_join(incidence_rate)
    
    return(
      summary_rho %>% 
        ggplot(aes(x=incidence_rate,y=median_rho,colour=admin1)) +
        geom_point() +
        # geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
        geom_linerange(aes(ymin=lo_rho,ymax=up_rho)) +
        facet_wrap(.~model) +
        theme_Publication(base_size=11) +
        labs(x="Incidence rate per 100,000",y=expression(rho)) +
        scale_colour_colorblind() +
        ylim(0,0.2) +
        xlim(0,180)
      
    )
    
  } else if(params=="gamma"){
    
    gamma_model <- list()
    
    for (i in seq_len(length(model_names))){
      gamma_province <- list()
      for (j in province_idx){
        gamma_province[[j]] <- 
          readRDS(paste0("output/model_output/parameters/fit_gamma_",j,"_",model_names[i],".rds"))
      }
      gamma_model[[i]] <- bind_rows(gamma_province)
    }
    
    incidence_rate <- cases_all %>% 
      filter(year > 2017) %>% 
      group_by(idadmin1,idadmin2,admin1,admin2) %>% 
      summarise(reported_cases=mean(reported_cases),
                pop=mean(pop)) %>% 
      ungroup() %>% 
      mutate(incidence_rate=reported_cases/pop*100000)
    
    summary_gamma <- bind_rows(gamma_model) %>% 
      left_join(incidence_rate)
    
    return(
      summary_gamma %>% 
        ggplot(aes(x=incidence_rate,y=median_gamma,colour=admin1)) +
        geom_point() +
        # geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
        geom_linerange(aes(ymin=lo_gamma,ymax=up_gamma)) +
        facet_wrap(.~model) +
        theme_Publication(base_size=11) +
        labs(x="Incidence rate per 100,000",y=expression(gamma)) +
        scale_colour_colorblind() +
        ylim(0,1) +
        xlim(0,180)
      
    )
    
  } else{
    stop("choose params as either 'lam_t', 'rho' or 'gamma'")
  }
  
}

#### functions for model comparisons
#### posterior predictive
compare_posterior_predictive <- function(model_names,province_idx){
  
  summary_post_pred_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_post_pred_model[[i]] <- 
      readRDS(paste0("output/model_output/posterior_predictive/postpred_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  summary_post_pred <- bind_rows(summary_post_pred_model) %>% 
    left_join(cases_all)
  
  return(
    summary_post_pred %>% 
      ggplot(aes(x=reported_cases,y=median_cases,colour=model)) +
      geom_abline(slope=1, intercept = 0, linetype=2, colour="gray") +
      geom_point() +
      geom_linerange(aes(ymin=lo_cases,ymax=up_cases)) +
      facet_grid(age_group~model,scales="free_y") +
      theme_Publication(base_size=11) +
      scale_colour_colorblind() +
      theme(legend.position = "none") +
      labs(x="Reported cases",y="Modelled cases")
  )
  
}

correlation_foi_incidence_rate_adj <- function(model_names,province_idx,cases_prov){
    
  summary_lam_t_model <- list()
  
  for (i in seq_len(length(model_names))){
    summary_lam_t_model[[i]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_t_",province_idx,"_",model_names[i],".rds"))
    
  }
  
  incidence_rate <- cases_prov
  
  summary_lam_t <- bind_rows(summary_lam_t_model) %>% 
    left_join(incidence_rate)
  
  corr_lam_t <- summary_lam_t %>% 
    group_by(model,admin2) %>% 
    summarise(correlation=cor(median_foi,incidence_rate))
  
  return(
    list(
      plot1=summary_lam_t %>% 
        ggplot(aes(x=incidence_rate,y=median_foi*4,colour=admin2)) +
        geom_point() +
        geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
        # geom_linerange(aes(ymin=lo_foi*4,ymax=up_foi*4)) +
        facet_wrap(.~model) +
        theme_Publication(base_size=11) +
        labs(x="Incidence rate per 100,000",y=expression("4"*lambda)),
      
      plot2=corr_lam_t %>% 
        ggplot(aes(y=correlation,x=model,fill=model)) +
        geom_boxplot() +
        theme_Publication(base_size=11) +
        scale_fill_colorblind() +
        theme(legend.position = "none") +
        ylim(0.5,1) +
        labs(x="Model",y="Incidence rate & FOI\ncorrelations (district-level)")
    )
  )  
    
}
