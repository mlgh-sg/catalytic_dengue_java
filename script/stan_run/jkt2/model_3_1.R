library(rstan)

options(mc.cores = 4)

# default params
# model_name <- c("Poisson-PS-Hierarchical","Negbin-PS-Hierarchical",
#                 "Poisson-PS-NonHierarchical","Negbin-PS-NonHierarchical")

# province_name <- c("jakarta","westjava")
no_district <- c(6,27)

# change these params for model runs
model_no <- 3
province <- 1
iter_sampling_stan <- 3000
iter_warmup_stan <- 1000

# read data
if (model_no %in% c(1:2)){
  data_list <- readRDS("data/data_jkt2.rds") 
} else if (model_no %in% c(3:4)){
  data_list <- readRDS("data/data_indiv_jkt2.rds") 
}
params_list <- readRDS("data/params_list.rds")
seed_rand <- readRDS("data/seed_rand.rds")

# read model
if (model_no==1){
  model_stan <- stan_model("model/dengue_catalytic_poisson_PS_hierarchical.stan")  
} else if(model_no==2){
  model_stan <- stan_model("model/dengue_catalytic_negbin_PS_hierarchical.stan")  
} else if(model_no==3){
  model_stan <- stan_model("model/dengue_catalytic_poisson_PS_nonhierarchical.stan")  
} else if(model_no==4){
  model_stan <- stan_model("model/dengue_catalytic_negbin_PS_nonhierarchical.stan")  
} 

# run model
if (model_no %in% c(1:2)){
  
  stanfit <- sampling(model_stan, 
                      data = data_list[[province]], 
                      warmup = iter_warmup_stan, 
                      iter = iter_warmup_stan+iter_sampling_stan, 
                      chains = 4,
                      pars = params_list[[model_no]],
                      seed = seed_rand[province])
  
} else if (model_no %in% c(3:4)){
  
  stanfit <- list()
  set.seed(seed_rand[province])
  seed_indiv_rand <- round(runif(no_district[province],1,1000000),0)
  
  for (i in seq_len(no_district[province])){
    
    stanfit[[i]] <- sampling(model_stan, 
                             data = data_list[[province]][[i]], 
                             warmup = iter_warmup_stan, 
                             iter = iter_warmup_stan+iter_sampling_stan, 
                             chains = 4,
                             pars = params_list[[model_no]],
                             seed = seed_indiv_rand[i])
    
  }
}

# save output
saveRDS(stanfit,paste0("output/stan_fit/fit_",model_no,"_",province,"_jkt2.rds"))