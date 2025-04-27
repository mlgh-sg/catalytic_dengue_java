#### library
library(tidyverse)
library(janitor)
library(rstan)
library(cowplot)
library(bayesplot)
library(readxl)
library(sf)
library(patchwork)
library(loo)
library(ggh4x)
source("script/functions.R")

#### read data
cases_all <- read_csv("data/cases_jakarta.csv")
cases_jakarta_4age <- read_csv("data/cases_jakarta_4age.csv")

#### params & variables
model_name <- c("Poisson\nhierarchical","Negative binomial\nhierarchical",
                "Poisson\nnon-hierarchical","Negative binomial\nnon-hierarchical")

model_name_jkt2 <- c("Poisson\nhierarchical (4 age)","Negative binomial\nhierarchical (4 age)",
                     "Poisson\nnon-hierarchical (4 age)","Negative binomial\nnon-hierarchical (4 age)")

model_name_2year <- c("Poisson\nhierarchical (2 years)","Negative binomial\nhierarchical (2 years)",
                      "Poisson\nnon-hierarchical (2 years)","Negative binomial\nnon-hierarchical (2 years)")

model_name_4year <- c("Poisson\nhierarchical (4 years)","Negative binomial\nhierarchical (4 years)",
                      "Poisson\nnon-hierarchical (4 years)","Negative binomial\nnon-hierarchical (4 years)")

provinces <- c("DKI JAKARTA","JAWA BARAT")

#### extracting and processing the output
#### output summary:
#### parameters: lam_H, lam_t, rho, gamma
#### posterior predictions
#### log-likelihood: province, district

#### the general model:
#### hierarchical model
for (i in 1:2){
  
  for (j in 1:1){
    
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,".rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(80,79)
    
    fit_lam_H <- extract_lam_H(catalytic_fit,cases_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name[i],".rds"))
    fit_lam_t <- extract_lam_t(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name[i],".rds"))
    fit_rho <- extract_rho(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name[i],".rds"))
    fit_gamma <- extract_gamma(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name[i],".rds"))
    if (i %in% c(10)){
      fit_phi <- extract_phi(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    
    post_samples_model <- extract(catalytic_fit)$y_rep
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_year <- list()
      for (l in seq_len(length(year_province))){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,k,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name[i],".rds"))
    
    #### log-likelihood
    #### province
    #### extract hierarchical model log-lik for this province
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name[i],".rds"))
    
    #### district
    loo_district <- list()
    
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_list[[j]]$nA
    nT_prov <- data_list[[j]]$nT
    
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    
    for (d in province_districts){
      # extract hierarchical model log-lik for this province
      log_lik_district <- extract_province_loglik(log_lik, d, nT_prov, nA_prov)
      loo_district[[d]] <- loo(log_lik_district)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name[i],".rds"))
    
  }
  
}

#### non-hierarchical model
for (i in 3:4){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,".rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(80,79)
    
    fit_lam_H <- extract_lam_H_indiv(catalytic_fit,cases_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name[i],".rds"))
    fit_lam_t <- extract_lam_t_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name[i],".rds"))
    fit_rho <- extract_rho_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name[i],".rds"))
    fit_gamma <- extract_gamma_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name[i],".rds"))
    if (i %in% c(14)){
      fit_phi <- extract_phi_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_model <- extract(catalytic_fit[[k]])$y_rep
      post_samples_year <- list()
      for (l in seq_len(length(year_province))){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    rm(post_samples_model)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name[i],".rds"))
    
    #### log-likelihood
    #### province
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- data_indiv_list[[j]][[1]]$nT
    
    # combine non-hierarchical log-liks for all districts in province
    log_lik <- do.call(cbind, lapply(province_districts, function(d) {
      extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
    }))
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name[i],".rds"))
    
    #### district
    loo_district <- list()
    
    # get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- data_indiv_list[[j]][[1]]$nT
    
    for (d in province_districts){
      # combine non-hierarchical log-liks for all districts in province
      log_lik <- extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
      loo_district[[d]] <- loo(log_lik)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name[i],".rds"))
    
  }
  
}

#### jakarta with four age groups instead of nine
#### hierarchical model
for (i in 1:2){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_jkt2.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(80,79)
    
    fit_lam_H <- extract_lam_H(catalytic_fit,cases_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_jkt2[i],".rds"))
    fit_lam_t <- extract_lam_t(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_jkt2[i],".rds"))
    fit_rho <- extract_rho(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_jkt2[i],".rds"))
    fit_gamma <- extract_gamma(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_jkt2[i],".rds"))
    if (i %in% c(10)){
      fit_phi <- extract_phi(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_jkt2[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_jkt2[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    
    post_samples_model <- extract(catalytic_fit)$y_rep
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_year <- list()
      for (l in seq_len(length(year_province))){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,k,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    
    age_group_prov <- c("0-4","5-14","15-44","45+")
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_jkt2[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_jkt2[i],".rds"))
    
    #### log-likelihood
    #### province
    #### extract hierarchical model log-lik for this province
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_jkt2[i],".rds"))
    
    #### district
    loo_district <- list()
    
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- 4
    nT_prov <- data_list[[j]]$nT
    
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    
    for (d in province_districts){
      # extract hierarchical model log-lik for this province
      log_lik_district <- extract_province_loglik(log_lik, d, nT_prov, nA_prov)
      loo_district[[d]] <- loo(log_lik_district)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_jkt2[i],".rds"))
    
  }
  
}

#### non-hierarchical model
for (i in 3:4){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_jkt2.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(80,79)
    
    fit_lam_H <- extract_lam_H_indiv(catalytic_fit,cases_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_jkt2[i],".rds"))
    fit_lam_t <- extract_lam_t_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_jkt2[i],".rds"))
    fit_rho <- extract_rho_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_jkt2[i],".rds"))
    fit_gamma <- extract_gamma_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_jkt2[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_jkt2[i],".rds"))
    if (i %in% c(14)){
      fit_phi <- extract_phi_indiv(catalytic_fit,cases_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_jkt2[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_jkt2[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_model <- extract(catalytic_fit[[k]])$y_rep
      post_samples_year <- list()
      for (l in seq_len(length(year_province))){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    rm(post_samples_model)
    
    age_group_prov <- c("0-4","5-14","15-44","45+")
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_jkt2[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_jkt2[i],".rds"))
    
    #### log-likelihood
    #### province
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- 4
    nT_prov <- data_indiv_list[[j]][[1]]$nT
    
    # combine non-hierarchical log-liks for all districts in province
    log_lik <- do.call(cbind, lapply(province_districts, function(d) {
      extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
    }))
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_jkt2[i],".rds"))
    
    #### district
    loo_district <- list()
    
    # get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- 4
    nT_prov <- data_indiv_list[[j]][[1]]$nT
    
    for (d in province_districts){
      # combine non-hierarchical log-liks for all districts in province
      log_lik <- extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
      loo_district[[d]] <- loo(log_lik)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_jkt2[i],".rds"))
    
  }
  
}

#### 4year
#### hierarchical model
for (i in 1:2){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_4year.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(83,83)
    
    cases_new_all <- cases_all %>% 
      filter(year>=2020)
    
    fit_lam_H <- extract_lam_H(catalytic_fit,cases_new_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_4year[i],".rds"))
    fit_lam_t <- extract_lam_t(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_4year[i],".rds"))
    fit_rho <- extract_rho(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_4year[i],".rds"))
    fit_gamma <- extract_gamma(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_4year[i],".rds"))
    if (i %in% c(10)){
      fit_phi <- extract_phi(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_4year[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_4year[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    
    post_samples_model <- extract(catalytic_fit)$y_rep
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_year <- list()
      for (l in seq_len(4)){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,k,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_4year[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_4year[i],".rds"))
    
    #### log-likelihood
    #### province
    #### extract hierarchical model log-lik for this province
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_4year[i],".rds"))
    
    #### district
    loo_district <- list()
    
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_list[[j]]$nA
    nT_prov <- 4
    
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    
    for (d in province_districts){
      # extract hierarchical model log-lik for this province
      log_lik_district <- extract_province_loglik(log_lik, d, nT_prov, nA_prov)
      loo_district[[d]] <- loo(log_lik_district)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_4year[i],".rds"))
    
  }
  
}

#### non-hierarchical model
for (i in 3:4){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_4year.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(83,83)
    
    cases_new_all <- cases_all %>% 
      filter(year>=2020)
    
    fit_lam_H <- extract_lam_H_indiv(catalytic_fit,cases_new_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_4year[i],".rds"))
    fit_lam_t <- extract_lam_t_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_4year[i],".rds"))
    fit_rho <- extract_rho_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_4year[i],".rds"))
    fit_gamma <- extract_gamma_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_4year[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_4year[i],".rds"))
    if (i %in% c(14)){
      fit_phi <- extract_phi_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_4year[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_4year[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_model <- extract(catalytic_fit[[k]])$y_rep
      post_samples_year <- list()
      for (l in seq_len(4)){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    rm(post_samples_model)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_4year[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_4year[i],".rds"))
    
    #### log-likelihood
    #### province
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- 4
    
    # combine non-hierarchical log-liks for all districts in province
    log_lik <- do.call(cbind, lapply(province_districts, function(d) {
      extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
    }))
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_4year[i],".rds"))
    
    #### district
    loo_district <- list()
    
    # get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- 4
    
    for (d in province_districts){
      # combine non-hierarchical log-liks for all districts in province
      log_lik <- extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
      loo_district[[d]] <- loo(log_lik)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_4year[i],".rds"))
    
  }
  
}

#### 2year
#### hierarchical model
for (i in 1:2){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_2year.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(85,85)
    
    cases_new_all <- cases_all %>% 
      filter(year>=2022)
    
    fit_lam_H <- extract_lam_H(catalytic_fit,cases_new_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_2year[i],".rds"))
    fit_lam_t <- extract_lam_t(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_2year[i],".rds"))
    fit_rho <- extract_rho(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_2year[i],".rds"))
    fit_gamma <- extract_gamma(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_2year[i],".rds"))
    if (i %in% c(10)){
      fit_phi <- extract_phi(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_2year[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_2year[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    
    post_samples_model <- extract(catalytic_fit)$y_rep
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_year <- list()
      for (l in seq_len(2)){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,k,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_2year[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_2year[i],".rds"))
    
    #### log-likelihood
    #### province
    #### extract hierarchical model log-lik for this province
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_2year[i],".rds"))
    
    #### district
    loo_district <- list()
    
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_list[[j]]$nA
    nT_prov <- 2
    
    log_lik <- extract_log_lik(catalytic_fit, merge_chains = TRUE)
    
    for (d in province_districts){
      # extract hierarchical model log-lik for this province
      log_lik_district <- extract_province_loglik(log_lik, d, nT_prov, nA_prov)
      loo_district[[d]] <- loo(log_lik_district)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_2year[i],".rds"))
    
  }
  
}

#### non-hierarchical model
for (i in 3:4){
  
  for (j in 1:1){
    
    #### read model fitting output
    catalytic_fit <- 
      readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_2year.rds"))
    
    #### summary statistics of parameter estimates
    n_samples <- 12000
    hist_length_province <- c(85,85)
    
    cases_new_all <- cases_all %>% 
      filter(year>=2022)
    
    fit_lam_H <- extract_lam_H_indiv(catalytic_fit,cases_new_all,provinces[j],hist_length_province[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_lam_H,paste0("output/model_output/parameters/fit_lam_H_",j,"_",model_name_2year[i],".rds"))
    fit_lam_t <- extract_lam_t_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_lam_t,paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_name_2year[i],".rds"))
    fit_rho <- extract_rho_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_rho,paste0("output/model_output/parameters/fit_rho_",j,"_",model_name_2year[i],".rds"))
    fit_gamma <- extract_gamma_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
      mutate(model=model_name_2year[i])
    saveRDS(fit_gamma,paste0("output/model_output/parameters/fit_gamma_",j,"_",model_name_2year[i],".rds"))
    if (i %in% c(14)){
      fit_phi <- extract_phi_indiv(catalytic_fit,cases_new_all,provinces[j],n_samples=n_samples) %>% 
        mutate(model=model_name_2year[i])
      saveRDS(fit_phi,paste0("output/model_output/parameters/fit_phi_",j,"_",model_name_2year[i],".rds"))
    }
    
    #### posterior summary
    cases_province <- cases_all %>% 
      filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    year_province <- cases_province %>% pull(year) %>% unique()
    age_group_province <- cases_province %>% pull(age_group) %>% unique()
    post_samples_admin2 <- list()
    
    for (k in seq_len(length(admin2_province))){
      post_samples_model <- extract(catalytic_fit[[k]])$y_rep
      post_samples_year <- list()
      for (l in seq_len(2)){
        post_samples_year[[l]] <- as_tibble(post_samples_model[,l,])
        colnames(post_samples_year[[l]]) <- age_group_province
        post_samples_year[[l]] <- post_samples_year[[l]] %>% 
          mutate(rep=1:n_samples) %>% 
          pivot_longer(-rep,names_to="age_group",values_to="y_rep") %>% 
          mutate(year=year_province[l])
      }
      post_samples_admin2[[k]] <- bind_rows(post_samples_year) %>% 
        mutate(admin2=admin2_province[k])
      rm(post_samples_year)
    }
    post_samples <- bind_rows(post_samples_admin2) %>% 
      mutate(admin1=provinces[j]) %>% 
      dplyr::select(rep,admin1,admin2,year,age_group,y_rep)
    rm(post_samples_admin2)
    rm(post_samples_model)
    
    if (j == 1){
      age_group_prov <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
    } else{
      age_group_prov <- c("0-4","5-14","15-44","45+")
    }
    
    post_summary <- post_samples %>% 
      group_by(admin1,admin2,year,age_group) %>% 
      summarise(median_cases=median(y_rep),
                lo_cases=quantile(y_rep,probs=c(0.025)),
                up_cases=quantile(y_rep,probs=c(0.975))) %>% 
      ungroup() %>% 
      left_join(cases_all) %>% 
      dplyr::select(idadmin1,idadmin2,admin1,admin2,year,age_group,pop,everything()) %>% 
      mutate(age_group=factor(age_group,levels=age_group_prov)) %>% 
      arrange(idadmin2,year,age_group) %>% 
      mutate(model=model_name_2year[i])
    rm(post_samples)
    saveRDS(post_summary,paste0("output/model_output/posterior_predictive/postpred_",j,"_",model_name_2year[i],".rds"))
    
    #### log-likelihood
    #### province
    #### get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- 2
    
    # combine non-hierarchical log-liks for all districts in province
    log_lik <- do.call(cbind, lapply(province_districts, function(d) {
      extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
    }))
    loo_province <- loo(log_lik)
    saveRDS(loo_province,paste0("output/model_output/log_likelihood/loo_province_",j,"_",model_name_2year[i],".rds"))
    
    #### district
    loo_district <- list()
    
    # get districts in this province
    cases_province <- cases_all %>% filter(admin1==provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique() 
    province_districts <- seq_len(length(admin2_province))
    
    nA_prov <- data_indiv_list[[j]][[1]]$nA
    nT_prov <- 2
    
    for (d in province_districts){
      # combine non-hierarchical log-liks for all districts in province
      log_lik <- extract_log_lik(catalytic_fit[[d]], merge_chains = TRUE)
      loo_district[[d]] <- loo(log_lik)
    }
    saveRDS(loo_district,paste0("output/model_output/log_likelihood/loo_district_",j,"_",model_name_2year[i],".rds"))
    
  }
  
}

}