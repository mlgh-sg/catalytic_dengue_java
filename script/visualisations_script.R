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
library(binom)
source("script/functions.R")
source("script/modeL_comparisons_functions.R")

#### figure 1 epidemiology in jakarta and west java
#### yearly incidence rate
#### province-level
cases_all <- read_csv("data/cases_jakarta.csv")
cases_jakarta_4age <- read_csv("data/cases_jakarta_4age.csv")

incidence_rate_admin1 <- cases_all %>% 
  filter(admin1 %in% c("DKI JAKARTA")) %>% 
  group_by(admin1,year) %>% 
  summarise(reported_cases=sum(reported_cases),
            pop=sum(pop)) %>% 
  ungroup() %>% 
  mutate(incidence_rate=reported_cases/pop*100000)

#### district-level
incidence_rate_admin2 <- cases_all %>% 
  filter(admin1 %in% c("DKI JAKARTA")) %>% 
  group_by(admin1,admin2,year) %>% 
  summarise(reported_cases=sum(reported_cases),
            pop=sum(pop)) %>% 
  ungroup() %>% 
  mutate(incidence_rate=reported_cases/pop*100000)

#### visualise
fig1_A <- incidence_rate_admin1 %>% 
  ggplot(aes(x=year,y=incidence_rate)) +
  geom_line(data=incidence_rate_admin2,aes(group=admin2),colour="darkgray",alpha=0.25) +
  geom_point(colour="#de2d26",size=2.5,linewidth=2.5) +
  geom_line(colour="#de2d26",linewidth=1) +
  facet_wrap(.~admin1) +
  theme_Publication(base_size=11) + 
  # theme(axis.title=element_text(size=8)) +
  labs(x="Year",y="Reported incidence rate\nper 100,000",tag="A")

#### geographical distributions
incidence_rate_overall_admin2 <- cases_all %>% 
  filter(admin1 %in% c("DKI JAKARTA")) %>%
  filter(year > 2016) %>% 
  group_by(idadmin1,idadmin2,admin1,admin2,year) %>% 
  summarise(reported_cases=sum(reported_cases),
            pop=sum(pop)) %>% 
  ungroup() %>% 
  group_by(idadmin1,idadmin2,admin1,admin2) %>% 
  summarise(reported_cases=sum(reported_cases),
            pop=mean(pop)) %>% 
  ungroup() %>% 
  mutate(incidence_rate=reported_cases/pop*100000)

admin2_jkt_shp <- read_sf(dsn="data/shapefiles/admin2",
                              layer="admin2", stringsAsFactors = FALSE) %>% 
  filter(A1N_BPS %in% c("DKI Jakarta"))

admin2_kep_seribu_shp <- read_sf(dsn="data/shapefiles/admin2",
                                 layer="admin2", stringsAsFactors = FALSE) %>% 
  filter(A2N_BPS %in% c("Kepulauan Seribu"))

# box: xmin: 106.3703 ymin: -7.820979 xmax: 108.8469 ymax: -5.184322

admin1_java_shp <- read_sf(dsn="data/shapefiles/admin1",
                           layer="admin1", stringsAsFactors = FALSE) %>% 
  filter(A1N_2013 %in% c("DKI JAKARTA","JAWA BARAT","BANTEN","JAWA TENGAH"))

admin2_jkt_ir_shp <- admin2_jkt_shp %>% 
  left_join(incidence_rate_overall_admin2 %>% 
              rename(A2C_BPS=idadmin2))

admin2_kep_seribu_ir_shp <- admin2_kep_seribu_shp %>% 
  left_join(incidence_rate_overall_admin2 %>% 
              rename(A2C_BPS=idadmin2))

fig1_C <- admin2_jkt_ir_shp %>% 
  ggplot() +
  geom_sf(data=admin1_java_shp,fill="gray") +
  geom_sf(aes(fill=incidence_rate),color = 'gray30', linewidth = 0.2) +
  geom_sf(data=admin1_java_shp,fill=NA,col="black", linewidth=0.5) +
  annotate(geom="rect",
           xmin=106.3831, xmax=106.8479, ymin=-6.040159, ymax=-5.084322, 
           fill=NA,col="black",linewidth=0.5,linetype=2) +
  annotate(geom="text",
           x=106.5831, y=-5.184322, 
           label="428.84",col="black",size=3) +
  theme_Publication(base_size=11) +
  scale_fill_viridis_c() +
  coord_sf(xlim=c(106.3703, 108.8469),
           ylim=c(-7.820979, -5.184322)) +
  labs(fill="Incidence rate per 100,000\n(2017-2023)",x=NULL,y=NULL,tag="C")

#### visualise age distribution
prop_age_admin1 <- cases_jakarta_4age %>% 
  group_by(admin1,age_group) %>% 
  summarise(reported_cases=sum(reported_cases)) %>% 
  ungroup() %>% 
  group_by(admin1) %>% 
  mutate(prop_cases=reported_cases/sum(reported_cases)*100) %>% 
  ungroup()

prop_age_admin2 <- cases_jakarta_4age %>% 
  group_by(admin1,admin2,age_group) %>% 
  summarise(reported_cases=sum(reported_cases)) %>% 
  ungroup() %>% 
  group_by(admin1,admin2) %>% 
  mutate(prop_cases=reported_cases/sum(reported_cases)*100) %>% 
  ungroup()

fig1_B <- prop_age_admin1 %>% 
  ggplot(aes(x=age_group,y=prop_cases)) +
  geom_line(data=prop_age_admin2,aes(group=admin2),colour="darkgray",alpha=0.25) +
  geom_point(colour="#de2d26",size=2.5) +
  geom_line(aes(group=admin1),colour="#de2d26",linewidth=1) +
  facet_wrap(.~admin1) +
  theme_Publication(base_size=11) + 
  # theme(axis.title=element_text(size=11)) +
  ylim(0,60) +
  labs(x="Age group",y="Proportion of\nreported cases (%)",tag="B")

fig1 <- (fig1_A / fig1_B) | fig1_C
# ggsave("output/figures/fig1.jpg",fig1,height=11,width=17,unit="cm")

#### figure 2 elpd
plot_elpd_district_jkt <- compare_elpd_district(model_name[1:4],1)$plot + 
  labs(tag="A") + ylim(-800,0)

plot_versus_elpd_district_jkt <- plot_elpd_district_jkt$data %>% 
  pivot_wider(values_from=elpd_mean,names_from=model) %>% 
  ggplot() +
  geom_point(aes(x=`Model 2`,y=`Model 4`),
             colour="#de2d26") +
  theme_Publication(base_size=11) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  ylim(-400,-50) + xlim(-400,-50) +
  labs(tag="C",x="Negative binomial\nhierarchical",y="Negative binomial\nnon-hierarchical")

fig2 <- (plot_elpd_district_jkt) / (plot_versus_elpd_district_jkt)
# ggsave("output/figures/fig2.jpg",fig2,height=11,width=17,unit="cm")

#### seroprevalence comparisons
#### read seroprevalence data
seroprev_data <- read_csv("data/jakarta_seroprev_data.csv") %>% 
  mutate(admin2=toupper(location_admin2),admin3=toupper(location_admin3)) %>% 
  dplyr::select(admin2,admin3,n,n_pos,age=mid_age)

median_seroprev_data <- vector()
lo_seroprev_data <- vector()
up_seroprev_data <- vector()

for (i in seq_len(nrow(seroprev_data))){
  
  x <- seroprev_data[i,]$n_pos
  N <- seroprev_data[i,]$n
  
  median_seroprev_data[i] <- binom.bayes(x,N)$mean
  lo_seroprev_data[i] <- binom.bayes(x,N)$lower
  up_seroprev_data[i] <- binom.bayes(x,N)$upper
  
  
}

seroprev_data <- seroprev_data %>% 
  mutate(median_seroprev=median_seroprev_data,
         lo_seroprev=lo_seroprev_data,
         up_seroprev=up_seroprev_data,
         type="2014 survey") %>% 
  dplyr::select(-admin3,-n,-n_pos)

#### filter seroprev data by province
#### jakarta
seroprev_data1 <- seroprev_data %>% 
  filter(admin2 %in% c("KOTA JAKARTA BARAT","KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR"))

#### overall seroprevalence 1-18 years old
seroprev_overall_data <- read_csv("data/jakarta_seroprev_data.csv") %>% 
  mutate(admin2=toupper(location_admin2),admin3=toupper(location_admin3)) %>% 
  dplyr::select(admin2,admin3,n,n_pos,age=mid_age) %>% 
  group_by(admin2,admin3) %>% 
  summarise(n=sum(n),
            n_pos=sum(n_pos)) %>% ungroup()

median_seroprev_data <- vector()
lo_seroprev_data <- vector()
up_seroprev_data <- vector()

for (i in seq_len(nrow(seroprev_overall_data))){
  
  x <- seroprev_overall_data[i,]$n_pos
  N <- seroprev_overall_data[i,]$n
  
  median_seroprev_data[i] <- binom.bayes(x,N)$mean
  lo_seroprev_data[i] <- binom.bayes(x,N)$lower
  up_seroprev_data[i] <- binom.bayes(x,N)$upper
  
  
}

seroprev_overall_data <- seroprev_overall_data %>% 
  mutate(median_seroprev=median_seroprev_data,
         lo_seroprev=lo_seroprev_data,
         up_seroprev=up_seroprev_data,
         type="2014 survey") %>% 
  dplyr::select(-admin3,-n,-n_pos)

seroprev_overall_data1 <- seroprev_overall_data %>% 
  filter(admin2 %in% c("KOTA JAKARTA BARAT","KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR"))

#### run functions for comparisons
seroprev_nb_hier_jkt <- compare_seroprev_2014("Negative binomial\nhierarchical","fit_2_1.rds",1,2017)
seroprev_nb_nonhier_jkt <- compare_seroprev_2014("Negative binomial\nnon-hierarchical","fit_4_1.rds",1,2017)

seroprev_jkt <- bind_rows(seroprev_nb_hier_jkt$table,seroprev_nb_nonhier_jkt$table)

seroprev_hier_jkt <- bind_rows(seroprev_nb_hier_jkt$table %>% mutate(type="model"),
                               seroprev_data1 %>% mutate(type="data")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=type,values_from=c(median_seroprev,lo_seroprev,up_seroprev))

seroprev_nonhier_jkt <- bind_rows(seroprev_nb_nonhier_jkt$table %>% mutate(type="model"),
                                  seroprev_data1 %>% mutate(type="data")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=type,values_from=c(median_seroprev,lo_seroprev,up_seroprev))

seroprev_all_jkt <- bind_rows(seroprev_hier_jkt,seroprev_nonhier_jkt)

seroprev_overall_jkt <- bind_rows(seroprev_nb_hier_jkt$table2,seroprev_nb_nonhier_jkt$table2,
                                  seroprev_overall_data1)

plot_seroprev_jkt <- seroprev_jkt %>% 
  ggplot(aes(x=age)) +
  geom_line(aes(y=median_seroprev*100,colour=type),linewidth=1.25) +
  geom_ribbon(aes(ymin=lo_seroprev*100,ymax=up_seroprev*100,fill=type),alpha=0.25,colour=NA) +
  geom_point(data=seroprev_data1,aes(y=median_seroprev*100,colour=type,fill=type)) +
  geom_linerange(data=seroprev_data1,
                 aes(ymin=lo_seroprev*100,ymax=up_seroprev*100,colour=type,fill=type)) +
  theme_bw() +
  facet_wrap(.~admin2,ncol=4) +
  labs(x="Age",y="Seroprevalence (%)",colour=NULL,tag="A") +
  theme_Publication(base_size = 10) +
  theme(legend.position="none") +
  scale_colour_colorblind() +
  scale_fill_colorblind() +
  ylim(0,100) +
  scale_x_continuous(breaks=seq(1,18,2))

fig3_AB <- (plot_seroprev_jkt) + plot_layout(heights = c(0.5,1))

seroprev_overall_all <- bind_rows(seroprev_overall_jkt %>% mutate(admin1="DKI JAKARTA"))

fig3_C <- seroprev_overall_all %>% 
  ggplot(aes(y=admin2)) +
  geom_point(aes(x=median_seroprev*100,colour=type),
             position=position_dodge(width=0.75)) +
  geom_linerange(aes(xmin=lo_seroprev*100,xmax=up_seroprev*100,colour=type),
                 position=position_dodge(width=0.75)) +
  labs(x="Age",y="Seroprevalence overall 1-18 years (%)",colour=NULL,tag="C") +
  facet_wrap(.~admin1,ncol=1,scales="free_y") +
  force_panelsizes(rows = c(0.5, 1)) +
  theme_Publication(base_size = 10) +
  theme(legend.position="none",axis.text.y=element_text(size=6)) +
  scale_colour_colorblind() +
  xlim(0,100)

plot_versus_seroprev_jkt <- seroprev_all_jkt %>% 
  ggplot() +
  geom_point(aes(x=median_seroprev_data*100,y=median_seroprev_model*100,colour=model)) +
  geom_linerange(aes(x=median_seroprev_data*100,
                     ymin=lo_seroprev_model*100,ymax=up_seroprev_model*100,colour=model)) +
  geom_linerange(aes(y=median_seroprev_model*100,
                     xmin=lo_seroprev_data*100,xmax=up_seroprev_data*100,colour=model)) +
  scale_color_manual(values=c("#E69F00","#56B4E9")) +
  facet_wrap(.~model,ncol=1) +
  theme_Publication(base_size=10) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  xlim(0,100) +
  ylim(0,100) +
  labs(x="Seroprevalence (data)",
       y="Seroprevalence (model)",
       tag="D") +
  theme(legend.position = "none")

fig3_DE <- (plot_versus_seroprev_jkt)

fig3 <- (((free(plot_seroprev_jkt)) + plot_layout(heights = c(0.5,1),nrow=1)) +
           ((fig3_C | plot_versus_seroprev_jkt) +
              plot_layout(ncol=2,widths=c(0.70,1,1)))) + plot_layout(nrow=2,heights=c(0.4,0.8,1))

# ggsave("output/figures/fig3.jpg",fig3,height=20,width=17,unit="cm")

#### FoI
plot_foi_jkt <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nnon-hierarchical"),
                                   1,1)$plot +
  labs(tag="A") + theme(legend.position="none")

estim_foi_jkt <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nnon-hierarchical"),
                                    1,2)$table %>% 
  pivot_wider(names_from=model,values_from=c(median_foi,lo_foi,up_foi)) %>% 
  mutate(admin1="DKI JAKARTA")

plot_versus_foi <- estim_foi_jkt %>%
  ggplot() +
  geom_point(aes(x=`median_foi_Negative binomial\nhierarchical`*4,
                 y=`median_foi_Negative binomial\nnon-hierarchical`*4),
             colour="#de2d26") +
  geom_linerange(aes(x=`median_foi_Negative binomial\nhierarchical`*4,
                     ymin=`lo_foi_Negative binomial\nnon-hierarchical`*4,
                     ymax=`up_foi_Negative binomial\nnon-hierarchical`*4),
                 colour="#de2d26") +
  geom_linerange(aes(y=`median_foi_Negative binomial\nnon-hierarchical`*4,
                     xmin=`lo_foi_Negative binomial\nhierarchical`*4,
                     xmax=`up_foi_Negative binomial\nhierarchical`*4),
                 colour="#de2d26") +
  theme_Publication(base_size=11) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_wrap(.~admin1) +
  labs(tag="C",x="Negative binomial\nhierarchical",y="Negative binomial\nnon-hierarchical") +
  ylim(0,0.6) +
  xlim(0,0.6)

estim_foi_jkt2 <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nnon-hierarchical"),
                                     1,2)$table %>% 
  mutate(admin1="DKI JAKARTA")

estim_foi_ir_all <- estim_foi_jkt2 %>% 
  mutate(year=as.numeric(t)) %>% 
  left_join(incidence_rate_admin2)

plot_foi_ir <- estim_foi_ir_all %>% 
  ggplot(aes(x=incidence_rate,y=median_foi*4,colour=admin2)) +
  geom_point() +
  geom_smooth(se = FALSE,formula = y ~ x,method='lm') +
  theme_Publication(base_size=11) +
  facet_grid(admin1~model) +
  theme(legend.position = "none") +
  labs(x="Model",y="Incidence rate & FOI\ncorrelations (district-level)",tag="D")

plot_corr_foi_ir <- estim_foi_ir_all %>% 
  group_by(model,admin2,admin1) %>% 
  summarise(correlation=cor(median_foi,incidence_rate)) %>% 
  ggplot(aes(y=correlation,x=model,fill=model)) +
  geom_boxplot() +
  theme_Publication(base_size=11) +
  facet_wrap(.~admin1,ncol=1) +
  scale_fill_manual(values=c("#E69F00","#56B4E9")) +
  theme(legend.position = "none") +
  ylim(0.5,1) +
  labs(x="Model",y="Incidence rate & FOI\ncorrelations (district-level)",tag="D")

# Create the plot
plot_example_for_ir <- estim_foi_ir_all %>% 
  filter(admin2 %in% c("KEPULAUAN SERIBU","KOTA JAKARTA BARAT")) %>% 
  ggplot(aes(x = t)) +
  geom_ribbon(aes(ymin = lo_foi * 4, ymax = up_foi * 4, fill = model), alpha = 0.2) +
  geom_line(aes(y = median_foi * 4, color = model)) +
  geom_line(aes(y = incidence_rate/333), # Adjusted scaling factor for new FoI range
            color = "#de2d26", linetype = "dashed") +
  geom_point(aes(y = incidence_rate/333),
             color = "#de2d26", size = 2) +
  scale_y_continuous(
    name = expression("4"*lambda),
    limits = c(0, 0.6), # Adjusted for 4x multiplication
    sec.axis = sec_axis(~ . * 333, 
                        name = "Incidence rate per 100,000",
                        breaks = seq(0, 200, 50))
  ) +
  scale_x_continuous(breaks = 2017:2023) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9")) +
  facet_wrap(~admin2, scales = "free_y",ncol=2) +
  theme_Publication(base_size = 11) +
  theme(legend.position="bottom") +
  labs(x = "Year",colour=NULL,fill=NULL,tag="E")

fig4 <- (((free(plot_foi_jkt))) +
           ((plot_versus_foi | plot_corr_foi_ir) + 
              plot_layout(ncol=2,widths=c(1,1))) + 
           free(plot_example_for_ir)) + plot_layout(nrow=3,heights=c(0.5,0.8,1))

# ggsave("output/figures/fig4.jpg",fig4,height=22,width=17,unit="cm")

#### rho dan gamma
#### A scatter plot versus for rho and gamma
#### B correlation between rho and gamma with IR
#### C correlation between rho and gamma and FoI
#### supplementary forest plot for all districts

model_names <- c("Negative binomial\nhierarchical","Negative binomial\nnon-hierarchical")
province_idx <- c(1)

lam_t_model <- list()
gamma_model <- list()
rho_model <- list()

for (i in seq_len(length(model_names))){
  lam_t_province <- list()
  gamma_province <- list()
  rho_province <- list()
  for (j in province_idx){
    lam_t_province[[j]] <- 
      readRDS(paste0("output/model_output/parameters/fit_lam_t_",j,"_",model_names[i],".rds"))
    gamma_province[[j]] <- 
      readRDS(paste0("output/model_output/parameters/fit_gamma_",j,"_",model_names[i],".rds"))
    rho_province[[j]] <- 
      readRDS(paste0("output/model_output/parameters/fit_rho_",j,"_",model_names[i],".rds"))
  }
  lam_t_model[[i]] <- bind_rows(lam_t_province)
  gamma_model[[i]] <- bind_rows(gamma_province)
  rho_model[[i]] <- bind_rows(rho_province)
}

gamma_model1 <- list()
gamma_model2 <- list()

rho_model1 <- list()
rho_model2 <- list()

lam_t_model1 <- list()
lam_t_model2 <- list()

#### compare models
gamma_model1[[1]] <- gamma_model[[1]] %>% 
  dplyr::select(admin2,
                median_gamma_1=median_gamma,lo_gamma_1=lo_gamma,up_gamma_1=up_gamma) %>% 
  ungroup()
gamma_model1[[2]] <- gamma_model[[2]] %>% 
  dplyr::select(admin2,
                median_gamma_2=median_gamma,lo_gamma_2=lo_gamma,up_gamma_2=up_gamma) %>% 
  ungroup()

rho_model1[[1]] <- rho_model[[1]] %>% 
  dplyr::select(admin2,
                median_rho_1=median_rho,lo_rho_1=lo_rho,up_rho_1=up_rho) %>% 
  ungroup()
rho_model1[[2]] <- rho_model[[2]] %>% 
  dplyr::select(admin2,
                median_rho_2=median_rho,lo_rho_2=lo_rho,up_rho_2=up_rho) %>% 
  ungroup()

gamma_compare_model <- gamma_model1[[1]] %>% left_join(gamma_model1[[2]]) %>%
  left_join(incidence_rate_admin2)
rho_compare_model <- rho_model1[[1]] %>% left_join(rho_model1[[2]]) %>%
  left_join(incidence_rate_admin2)

plot_rho_compare_model <- rho_compare_model %>% 
  ggplot() +
  geom_point(aes(x=median_rho_1,
                 y=median_rho_2),
             colour="#de2d26") +
  geom_linerange(aes(x=median_rho_1,
                     ymin=lo_rho_2,
                     ymax=up_rho_2),
                 colour="#de2d26") +
  geom_linerange(aes(y=median_rho_2,
                     xmin=lo_rho_1,
                     xmax=up_rho_1),
                 colour="#de2d26") +
  theme_Publication(base_size=11) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_wrap(.~admin1) +
  labs(tag="A",x="Negative binomial\nhierarchical",y="Negative binomial\nnon-hierarchical",
       colour=NULL) +
  scale_x_continuous(breaks=c(0,0.1,0.2),limits=c(0,0.2)) +
  scale_y_continuous(breaks=c(0,0.1,0.2),limits=c(0,0.2))

plot_gamma_compare_model <- gamma_compare_model %>% 
  ggplot() +
  geom_point(aes(x=median_gamma_1,
                 y=median_gamma_2),
             colour="#de2d26") +
  geom_linerange(aes(x=median_gamma_1,
                     ymin=lo_gamma_2,
                     ymax=up_gamma_2),
                 colour="#de2d26") +
  geom_linerange(aes(y=median_gamma_2,
                     xmin=lo_gamma_1,
                     xmax=up_gamma_1),
                 colour="#de2d26") +
  theme_Publication(base_size=11) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_wrap(.~admin1) +
  labs(tag="B",x="Negative binomial\nhierarchical",y="Negative binomial\nnon-hierarchical",
       colour=NULL) +
  theme(axis.text=element_text(size=7.5)) +
  scale_x_continuous(breaks=c(0,0.5,1)) +
  scale_y_continuous(breaks=c(0,0.5,1))

#### compare with ir
gamma_model2[[1]] <- gamma_model[[1]] %>% 
  mutate(model=model_names[1]) %>% ungroup()
gamma_model2[[2]] <- gamma_model[[2]] %>% 
  mutate(model=model_names[2]) %>% ungroup()

rho_model2[[1]] <- rho_model[[1]] %>% 
  mutate(model=model_names[1]) %>% ungroup()
rho_model2[[2]] <- rho_model[[2]] %>% 
  mutate(model=model_names[2]) %>% ungroup()

lam_t_model2[[1]] <- lam_t_model[[1]] %>% 
  mutate(model=model_names[1]) %>% ungroup()
lam_t_model2[[2]] <- lam_t_model[[2]] %>% 
  mutate(model=model_names[2]) %>% ungroup()

gamma_compare_ir <- bind_rows(gamma_model2) %>% 
  left_join(incidence_rate_overall_admin2)
rho_compare_ir <- bind_rows(rho_model2) %>%
  left_join(incidence_rate_overall_admin2)
lam_t_compare_ir <- bind_rows(lam_t_model2) %>%
  left_join(incidence_rate_overall_admin2)

plot_rho_ir <- rho_compare_ir %>% 
  ggplot(aes(x=incidence_rate,y=median_rho,colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=lo_rho,ymax=up_rho)) +
  geom_smooth(se = FALSE,formula = y ~ x,method='lm',linetype=2) +
  facet_wrap(.~model) +
  theme_Publication(base_size=11) +
  labs(y=expression(rho),x="Incidence rate per 100,000",tag="C",colour=NULL) +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

plot_gamma_ir <- gamma_compare_ir %>% 
  ggplot(aes(x=incidence_rate,y=median_gamma,colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=lo_gamma,ymax=up_gamma)) +
  geom_smooth(se = FALSE,formula = y ~ x,method='lm',linetype=2) +
  facet_wrap(.~model) +
  theme_Publication(base_size=11) +
  labs(y=expression(gamma),x="Incidence rate per 100,000",tag="D",colour=NULL) +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

lam_t_max <- lam_t_compare_ir %>%
  group_by(admin2,model) %>%
  slice_max(order_by = median_foi, n = 1) %>% 
  arrange(idadmin2) %>% 
  dplyr::select(model,admin2,median_foi,lo_foi,up_foi)

rho_compare_foi <- lam_t_compare_ir %>% left_join(rho_compare_ir)
gamma_compare_foi <- lam_t_compare_ir %>% left_join(gamma_compare_ir)

plot_rho_foi <- rho_compare_foi %>% 
  ggplot(aes(x=median_foi*4,y=median_rho,colour=t)) +
  geom_point() +
  geom_linerange(aes(ymin=lo_rho,ymax=up_rho)) +
  geom_linerange(aes(xmin=lo_foi*4,xmax=up_foi*4)) +
  # geom_smooth(se = FALSE,formula = y ~ x,method='lm',linetype=2) +
  facet_grid(admin1~model) +
  theme_Publication(base_size=11) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 1)) +
  labs(y=expression(rho),x=expression("4"*lambda),tag="E",colour=NULL) +
  scale_colour_colorblind()

plot_gamma_foi <- gamma_compare_foi %>% 
  ggplot(aes(x=median_foi*4,y=median_gamma,colour=t)) +
  geom_point() +
  geom_linerange(aes(ymin=lo_gamma,ymax=up_gamma)) +
  geom_linerange(aes(xmin=lo_foi*4,xmax=up_foi*4)) +
  # geom_smooth(se = FALSE,formula = y ~ x,method='lm',linetype=2) +
  facet_grid(admin1~model) +
  theme_Publication(base_size=11) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 1)) +
  labs(y=expression(gamma),x=expression("4"*lambda),tag="E",colour=NULL) +
  scale_colour_colorblind()

fig5 <- (((free(plot_rho_compare_model) | free(plot_gamma_compare_model)) + plot_layout(ncol=2)) /
           ((plot_rho_ir | plot_gamma_ir) + plot_layout(ncol=2)) / 
           free(plot_rho_foi)) + plot_layout(nrow=3,heights=c(0.8,0.5,1))

# ggsave("output/figures/fig5.jpg",fig5,height=22,width=17,unit="cm")

#### different age grouping
foi_H_jkt_age_hier <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 age)"),1,1)$table2 %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nhierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nhierarchical")
foi_H_jkt_age_nonhier <- compare_all_lambda(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 age)"),1,1)$table2 %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nnon-hierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical")

foi_H_jkt_age <- bind_rows(foi_H_jkt_age_hier,foi_H_jkt_age_nonhier)

lambda_estim <- readRDS("output/model_output/lambda_estim.rds")

lambda_estim_select <- lambda_estim %>% 
  filter(admin2 %in% c("KOTA JAKARTA BARAT","KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR")) %>% 
  mutate(t=2014)

fig6_A <- foi_H_jkt_age %>% 
  ggplot() +
  geom_line(aes(x=t,y=median_foi*4,colour=age_groups,group=age_groups),linewidth=1) +
  geom_ribbon(aes(x=t,ymin=lo_foi*4,ymax=up_foi*4,fill=age_groups,group=age_groups),alpha=0.25) +
  geom_point(data=lambda_estim_select,aes(x=t,y=lambda_median)) +
  geom_linerange(data=lambda_estim_select,aes(x=t,ymin=lambda_lo,ymax=lambda_up)) +
  facet_grid(model~admin2) +
  theme_Publication(base_size=11) +
  scale_colour_manual(values=c("#0072B2","#CC79A7")) +
  scale_fill_manual(values=c("#0072B2","#CC79A7")) +
  scale_x_continuous(breaks=seq(1990,2016,5)) +
  theme(legend.position="right") +
  guides(color = guide_legend(ncol = 1)) +
  labs(x="Year",y=expression("4"*lambda),colour=NULL,fill=NULL,tag="A")

seroprev_nb_hier_jkt_9 <- compare_seroprev_2014("Negative binomial\nhierarchical","fit_2_1.rds",1,2017)
seroprev_nb_nonhier_jkt_9 <- compare_seroprev_2014("Negative binomial\nnon-hierarchical","fit_4_1.rds",1,2017)
seroprev_nb_hier_jkt_4 <- compare_seroprev_2014("Negative binomial\nhierarchical (4 age)","fit_2_1_jkt2.rds",1,2017)
seroprev_nb_nonhier_jkt_4 <- compare_seroprev_2014("Negative binomial\nnon-hierarchical (4 age)","fit_4_1_jkt2.rds",1,2017)

seroprev_jkt_age_hier <- bind_rows(seroprev_nb_hier_jkt_9$table,seroprev_nb_hier_jkt_4$table) %>% 
  pivot_wider(names_from=type,values_from=c(median_seroprev,lo_seroprev,up_seroprev))
seroprev_jkt_age_nonhier <- bind_rows(seroprev_nb_nonhier_jkt_9$table,seroprev_nb_nonhier_jkt_4$table) %>% 
  pivot_wider(names_from=type,values_from=c(median_seroprev,lo_seroprev,up_seroprev))

seroprev_jkt_age_hier <- bind_rows(seroprev_nb_hier_jkt_9$table,seroprev_nb_hier_jkt_4$table) %>% 
  mutate(age_groups=ifelse(type=="Negative binomial\nhierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(type="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_seroprev,lo_seroprev,up_seroprev))
seroprev_jkt_age_nonhier <- bind_rows(seroprev_nb_nonhier_jkt_9$table,seroprev_nb_nonhier_jkt_4$table) %>% 
  mutate(age_groups=ifelse(type=="Negative binomial\nnon-hierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(type="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_seroprev,lo_seroprev,up_seroprev))

seroprev_jkt_age <- bind_rows(seroprev_jkt_age_hier,seroprev_jkt_age_nonhier)

fig6_B <- seroprev_jkt_age %>% 
  ggplot(aes(x=`median_seroprev_9 age groups`,
             y=`median_seroprev_4 age groups`)) +
  geom_point(colour="#de2d26") +
  geom_linerange(aes(ymin=`lo_seroprev_4 age groups`,
                     ymax=`up_seroprev_4 age groups`),
                 colour="#de2d26") +
  geom_linerange(aes(xmin=`lo_seroprev_9 age groups`,
                     xmax=`up_seroprev_9 age groups`),
                 colour="#de2d26") +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(type~.) +
  ylim(0,1) +
  xlim(0,1) +
  theme_Publication(base_size=11) +
  labs(x="Seroprevalence (9 age groups)",
       y="Seroprevalence (4 age groups)",tag="B")

foi_jkt_age_hier <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 age)"),1,2)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nhierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_foi,lo_foi,up_foi))
foi_jkt_age_nonhier <- compare_all_lambda(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 age)"),1,2)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nnon-hierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_foi,lo_foi,up_foi))

foi_jkt_age <- bind_rows(foi_jkt_age_hier,foi_jkt_age_nonhier)

fig6_C <- foi_jkt_age %>% 
  ggplot(aes(x=`median_foi_9 age groups`*4,
             y=`median_foi_4 age groups`*4)) +
  geom_point(colour="#de2d26") +
  geom_linerange(aes(ymin=`lo_foi_4 age groups`*4,
                     ymax=`up_foi_4 age groups`*4),
                 colour="#de2d26") +
  geom_linerange(aes(xmin=`lo_foi_9 age groups`*4,
                     xmax=`up_foi_9 age groups`*4),
                 colour="#de2d26") +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.6) +
  xlim(0,0.6) +
  theme_Publication(base_size=11) +
  labs(x=expression("4"*lambda*" (9 age groups)"),
       y=expression("4"*lambda*" (4 age groups)"),tag="C")

rho_jkt_age_hier <- compare_rho(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 age)"),1)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nhierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_rho,lo_rho,up_rho))
rho_jkt_age_nonhier <- compare_rho(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 age)"),1)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nnon-hierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_rho,lo_rho,up_rho))
rho_jkt_age <- bind_rows(rho_jkt_age_hier,rho_jkt_age_nonhier)

fig6_D <- rho_jkt_age %>% 
  ggplot(aes(x=`median_rho_9 age groups`,
             y=`median_rho_4 age groups`)) +
  geom_point(colour="#de2d26") +
  geom_linerange(aes(ymin=`lo_rho_4 age groups`,
                     ymax=`up_rho_4 age groups`),
                 colour="#de2d26") +
  geom_linerange(aes(xmin=`lo_rho_9 age groups`,
                     xmax=`up_rho_9 age groups`),
                 colour="#de2d26") +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.25) +
  xlim(0,0.25) +
  theme_Publication(base_size=11) +
  labs(x=expression(rho*" (9 age groups)"),
       y=expression(rho*" (4 age groups)"),tag="D")

gamma_jkt_age_hier <- compare_gamma(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 age)"),1)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nhierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_gamma,lo_gamma,up_gamma))
gamma_jkt_age_nonhier <- compare_gamma(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 age)"),1)$table %>% 
  mutate(age_groups=ifelse(model=="Negative binomial\nnon-hierarchical (4 age)","4 age groups","9 age groups")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=age_groups,values_from=c(median_gamma,lo_gamma,up_gamma))
gamma_jkt_age <- bind_rows(gamma_jkt_age_hier,gamma_jkt_age_nonhier)

fig6_E <- gamma_jkt_age %>% 
  ggplot(aes(x=`median_gamma_9 age groups`,
             y=`median_gamma_4 age groups`)) +
  geom_point(colour="#de2d26") +
  geom_linerange(aes(ymin=`lo_gamma_4 age groups`,
                     ymax=`up_gamma_4 age groups`),
                 colour="#de2d26") +
  geom_linerange(aes(xmin=`lo_gamma_9 age groups`,
                     xmax=`up_gamma_9 age groups`),
                 colour="#de2d26") +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,1) +
  xlim(0,1) +
  theme_Publication(base_size=11) +
  labs(x=expression(gamma*" (9 age groups)"),
       y=expression(gamma*" (4 age groups)"),tag="E")

fig6 <- (free(fig6_A) / (free(fig6_B) | free(fig6_C)) / (free(fig6_D) | free(fig6_E)))
# ggsave("output/figures/fig6.jpg",fig6,height=24,width=17,unit="cm")

#### supplementary materials
#### foi 4 years
#### jkt
foi_jkt_4year_hier <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 years)"),1,2)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_foi,lo_foi,up_foi)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
foi_jkt_4year_nonhier <- compare_all_lambda(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 years)"),1,2)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_foi,lo_foi,up_foi)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

foi_4year <- bind_rows(foi_jkt_4year_hier,foi_jkt_4year_nonhier)

s_fig3_A <- foi_4year %>% 
  ggplot(aes(x=`median_foi_Complete`*4,
             y=`median_foi_4 years`*4,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_foi_4 years`*4,
                     ymax=`up_foi_4 years`*4)) +
  geom_linerange(aes(xmin=`lo_foi_Complete`*4,
                     xmax=`up_foi_Complete`*4)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.6) +
  xlim(0,0.6) +
  theme_Publication(base_size=11) +
  labs(x=expression("4"*lambda*" (complete)"),
       y=expression("4"*lambda*" (4 years)"),tag="A") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

#### foi 2 years
#### jkt
foi_jkt_2year_hier <- compare_all_lambda(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (2 years)"),1,2)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_foi,lo_foi,up_foi)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
foi_jkt_2year_nonhier <- compare_all_lambda(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (2 years)"),1,2)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_foi,lo_foi,up_foi)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

foi_2year <- bind_rows(foi_jkt_2year_hier,foi_jkt_2year_nonhier)

s_fig3_B <- foi_2year %>% 
  ggplot(aes(x=`median_foi_Complete`*4,
             y=`median_foi_2 years`*4,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_foi_2 years`*4,
                     ymax=`up_foi_2 years`*4)) +
  geom_linerange(aes(xmin=`lo_foi_Complete`*4,
                     xmax=`up_foi_Complete`*4)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.6) +
  xlim(0,0.6) +
  theme_Publication(base_size=11) +
  labs(x=expression("4"*lambda*" (complete)"),
       y=expression("4"*lambda*" (2 years)"),tag="B") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

#### rho 4 years
#### jkt
rho_jkt_4year_hier <- compare_rho(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_rho,lo_rho,up_rho)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
rho_jkt_4year_nonhier <- compare_rho(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_rho,lo_rho,up_rho)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

rho_4year <- bind_rows(rho_jkt_4year_hier,rho_jkt_4year_nonhier)

s_fig3_C <- rho_4year %>% 
  ggplot(aes(x=`median_rho_Complete`,
             y=`median_rho_4 years`,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_rho_4 years`,
                     ymax=`up_rho_4 years`)) +
  geom_linerange(aes(xmin=`lo_rho_Complete`,
                     xmax=`up_rho_Complete`)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.4) +
  xlim(0,0.4) +
  theme_Publication(base_size=11) +
  labs(x=expression(rho*" (complete)"),
       y=expression(rho*" (4 years)"),tag="C") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

#### rho 2 years
#### jkt
rho_jkt_2year_hier <- compare_rho(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (2 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_rho,lo_rho,up_rho)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
rho_jkt_2year_nonhier <- compare_rho(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (2 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_rho,lo_rho,up_rho)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

rho_2year <- bind_rows(rho_jkt_2year_hier,rho_jkt_2year_nonhier)

s_fig3_D <- rho_2year %>% 
  ggplot(aes(x=`median_rho_Complete`,
             y=`median_rho_2 years`,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_rho_2 years`,
                     ymax=`up_rho_2 years`)) +
  geom_linerange(aes(xmin=`lo_rho_Complete`,
                     xmax=`up_rho_Complete`)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,0.4) +
  xlim(0,0.4) +
  theme_Publication(base_size=11) +
  labs(x=expression(rho*" (complete)"),
       y=expression(rho*" (2 years)"),tag="D") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

#### gamma 4 years
#### jkt
gamma_jkt_4year_hier <- compare_gamma(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (4 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_gamma,lo_gamma,up_gamma)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
gamma_jkt_4year_nonhier <- compare_gamma(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (4 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (4 years)","4 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_gamma,lo_gamma,up_gamma)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

gamma_4year <- bind_rows(gamma_jkt_4year_hier,gamma_jkt_4year_nonhier)

s_fig3_E <- gamma_4year %>% 
  ggplot(aes(x=`median_gamma_Complete`,
             y=`median_gamma_4 years`,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_gamma_4 years`,
                     ymax=`up_gamma_4 years`)) +
  geom_linerange(aes(xmin=`lo_gamma_Complete`,
                     xmax=`up_gamma_Complete`)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,1) +
  xlim(0,1) +
  theme_Publication(base_size=11) +
  labs(x=expression(gamma*" (complete)"),
       y=expression(gamma*" (4 years)"),tag="E") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

#### gamma 2 years
#### jkt
gamma_jkt_2year_hier <- compare_gamma(c("Negative binomial\nhierarchical","Negative binomial\nhierarchical (2 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nhierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nhierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_gamma,lo_gamma,up_gamma)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()
gamma_jkt_2year_nonhier <- compare_gamma(c("Negative binomial\nnon-hierarchical","Negative binomial\nnon-hierarchical (2 years)"),1)$table %>% 
  mutate(no_years=ifelse(model=="Negative binomial\nnon-hierarchical (2 years)","2 years","Complete")) %>% 
  mutate(model="Negative binomial\nnon-hierarchical") %>% 
  pivot_wider(names_from=no_years,values_from=c(median_gamma,lo_gamma,up_gamma)) %>% 
  mutate(admin1="DKI JAKARTA") %>% 
  drop_na()

gamma_2year <- bind_rows(gamma_jkt_2year_hier,gamma_jkt_2year_nonhier)

s_fig3_F <- gamma_2year %>% 
  ggplot(aes(x=`median_gamma_Complete`,
             y=`median_gamma_2 years`,
             colour=admin1)) +
  geom_point() +
  geom_linerange(aes(ymin=`lo_gamma_2 years`,
                     ymax=`up_gamma_2 years`)) +
  geom_linerange(aes(xmin=`lo_gamma_Complete`,
                     xmax=`up_gamma_Complete`)) +
  geom_abline(intercept=0,slope=1,linetype=2) +
  facet_grid(model~.) +
  ylim(0,1) +
  xlim(0,1) +
  theme_Publication(base_size=11) +
  labs(x=expression(gamma*" (complete)"),
       y=expression(gamma*" (2 years)"),tag="F") +
  scale_colour_manual(values=c("#009E73","#CC79A7"))

s_fig3 <- (free(s_fig3_A) | free(s_fig3_B)) / (free(s_fig3_C) | free(s_fig3_D)) / (free(s_fig3_E) | free(s_fig3_F))
# ggsave("output/figures/s_fig3.jpg",s_fig3,height=24,width=17,unit="cm")

#### generated susceptible 9 years old, 1 - susceptible = seroprevalence
#### high transmisson = seroprevalence > 60%
#### the probability of a district is high transmission

#### models
#### all hierarchicals listed first
model_names_policy <- c("Negative binomial\nhierarchical","Negative binomial\nnon-hierarchical")
model_names_jkt2_policy <- c("Negative binomial\nhierarchical (4 age)","Negative binomial\nnon-hierarchical (4 age)")
model_names_4year_policy <- c("Negative binomial\nhierarchical (4 years)","Negative binomial\nnon-hierarchical (4 years)")
model_names_2year_policy <- c("Negative binomial\nhierarchical (2 years)","Negative binomial\nnon-hierarchical (2 years)")

#### sample the posterior for simulations
posterior_seed <- list()
n_samples_sim <- 1000
for (i in c(2,4)){
  
  set.seed(i+1000)
  posterior_seed[[i]] <- sample(1:n_samples,n_samples_sim,replace=FALSE)
  
}

#### extract params posterior samples for simulations
#### create a list for each cases and population data:
#### list index for each admin2, row is age group, and column is year
provinces <- c("DKI JAKARTA")
provinces_filename <- c("jakarta")
cases_list <- list()
pop_list <- list()
for (i in seq_len(length(provinces))){
  cases_province <- cases_all %>% 
    filter(admin1 == provinces[i])
  admin2_province <- cases_province %>% pull(admin2) %>% unique()
  cases_list[[i]] <- list()
  pop_list[[i]] <- list()
  
  for (j in seq_len(length(admin2_province))){
    
    cases_list[[i]][[j]] <- cases_province %>% 
      filter(admin2==admin2_province[j]) %>% 
      dplyr::select(year,age_group,reported_cases) %>% 
      pivot_wider(names_from=year,values_from=reported_cases) %>% 
      dplyr::select(-age_group) %>% 
      as.matrix()
    pop_list[[i]][[j]] <- cases_province %>% 
      filter(admin2==admin2_province[j]) %>% 
      dplyr::select(year,age_group,pop) %>% 
      pivot_wider(names_from=year,values_from=pop) %>% 
      dplyr::select(-age_group) %>% 
      as.matrix()
    
  }
}

cases_jkt2_list <- list()
pop_jkt2_list <- list()
admin2_jkt <- cases_jakarta_4age %>% pull(admin2) %>% unique()

for (i in seq_len(length(admin2_jkt))){
  
  cases_jkt2_list[[i]] <- cases_jakarta_4age %>% 
    filter(admin2==admin2_jkt[i]) %>% 
    dplyr::select(year,age_group,reported_cases) %>% 
    pivot_wider(names_from=year,values_from=reported_cases) %>% 
    dplyr::select(-age_group) %>% 
    as.matrix()
  pop_jkt2_list[[i]] <- cases_jakarta_4age %>% 
    filter(admin2==admin2_jkt[i]) %>% 
    dplyr::select(year,age_group,pop) %>% 
    pivot_wider(names_from=year,values_from=pop) %>% 
    dplyr::select(-age_group) %>% 
    as.matrix()
  
}

#### standard model
params_simulation <- list()

for (i in c(2,4)){
  
  params_simulation[[i]] <- list()
  
  for (j in 1){
    
    if (i == 2){
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,".rds"))
      
      extract_stanfit <- extract(catalytic_fit,
                                 pars=c("lam_H","lam_t","rho","gamma")) 
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation[[i]][[j]] <- list()
      
      for (k in seq_len(length(admin2_province))){
        
        params_simulation[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- t(pop_list[[j]][[k]])
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],k,]
        params_simulation[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],k,]
        params_simulation[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]],k]
        params_simulation[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]],k]
        params_simulation[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    } else{
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation[[i]][[j]] <- list()
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,".rds"))
      
      for (k in seq_len(length(admin2_province))){
        
        extract_stanfit <- extract(catalytic_fit[[k]],
                                   pars=c("lam_H","lam_t","rho","gamma")) 
        
        params_simulation[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- t(pop_list[[j]][[k]])
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],]
        params_simulation[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],]
        params_simulation[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]]]
        params_simulation[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]]]
        params_simulation[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    }
    
  }
  
}

#### jakarta 4 age groups
params_simulation_jkt2 <- list()

for (i in c(2,4)){
  
  params_simulation_jkt2[[i]] <- list()
  
  for (j in c(1)){
    
    if (i == 2){
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_jkt2.rds"))
      
      extract_stanfit <- extract(catalytic_fit,
                                 pars=c("lam_H","lam_t","rho","gamma")) 
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_jkt2[[i]][[j]] <- list()
      
      for (k in seq_len(length(admin2_province))){
        
        params_simulation_jkt2[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- t(pop_jkt2_list[[k]])
        
        if (j == 1) {
          amin_sim <- c(0,5,15,45)
          amax_sim <- c(4,14,44,99)
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_jkt2[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_jkt2[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_jkt2[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],k,]
        params_simulation_jkt2[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],k,]
        params_simulation_jkt2[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]],k]
        params_simulation_jkt2[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]],k]
        params_simulation_jkt2[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_jkt2[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_jkt2[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    } else{
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_jkt2[[i]][[j]] <- list()
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_jkt2.rds"))
      
      for (k in seq_len(length(admin2_province))){
        
        extract_stanfit <- extract(catalytic_fit[[k]],
                                   pars=c("lam_H","lam_t","rho","gamma")) 
        
        params_simulation_jkt2[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- t(pop_jkt2_list[[k]])
        
        if (j == 1) {
          amin_sim <- c(0,5,15,45)
          amax_sim <- c(4,14,44,99)
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_jkt2[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_jkt2[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_jkt2[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],]
        params_simulation_jkt2[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],]
        params_simulation_jkt2[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]]]
        params_simulation_jkt2[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]]]
        params_simulation_jkt2[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_jkt2[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_jkt2[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    }
    
  }
  
}

#### 4 years of data
params_simulation_4year <- list()

for (i in c(2,4)){
  
  params_simulation_4year[[i]] <- list()
  
  for (j in 1){
    
    if (i == 2){
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_4year.rds"))
      
      extract_stanfit <- extract(catalytic_fit,
                                 pars=c("lam_H","lam_t","rho","gamma")) 
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_4year[[i]][[j]] <- list()
      
      for (k in seq_len(length(admin2_province))){
        
        params_simulation_4year[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- tail(t(pop_list[[j]][[k]]),4)
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_4year[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_4year[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_4year[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],k,]
        params_simulation_4year[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],k,]
        params_simulation_4year[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]],k]
        params_simulation_4year[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]],k]
        params_simulation_4year[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_4year[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_4year[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    } else{
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_4year[[i]][[j]] <- list()
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_4year.rds"))
      
      for (k in seq_len(length(admin2_province))){
        
        extract_stanfit <- extract(catalytic_fit[[k]],
                                   pars=c("lam_H","lam_t","rho","gamma")) 
        
        params_simulation_4year[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- tail(t(pop_list[[j]][[k]]),4)
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_4year[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_4year[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_4year[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],]
        params_simulation_4year[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],]
        params_simulation_4year[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]]]
        params_simulation_4year[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]]]
        params_simulation_4year[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_4year[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_4year[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    }
    
  }
  
}

#### 2 years of data
params_simulation_2year <- list()

for (i in c(2,4)){
  
  params_simulation_2year[[i]] <- list()
  
  for (j in 1){
    
    if (i == 2){
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_2year.rds"))
      
      extract_stanfit <- extract(catalytic_fit,
                                 pars=c("lam_H","lam_t","rho","gamma")) 
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_2year[[i]][[j]] <- list()
      
      for (k in seq_len(length(admin2_province))){
        
        params_simulation_2year[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- tail(t(pop_list[[j]][[k]]),2)
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_2year[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_2year[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_2year[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],k,]
        params_simulation_2year[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],k,]
        params_simulation_2year[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]],k]
        params_simulation_2year[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]],k]
        params_simulation_2year[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_2year[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_2year[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    } else{
      
      cases_province <- cases_all %>% filter(admin1 == provinces[j])
      admin2_province <- cases_province %>% pull(admin2) %>% unique()
      
      params_simulation_2year[[i]][[j]] <- list()
      
      catalytic_fit <- 
        readRDS(paste0("output/stan_fit/fit_",i,"_",j,"_2year.rds"))
      
      for (k in seq_len(length(admin2_province))){
        
        extract_stanfit <- extract(catalytic_fit[[k]],
                                   pars=c("lam_H","lam_t","rho","gamma")) 
        
        params_simulation_2year[[i]][[j]][[k]] <- list()
        
        pop_admin2_sim <- tail(t(pop_list[[j]][[k]]),2)
        
        if (j == 1) {
          amin_sim <- data_list[[1]]$ageLims[1,]
          amax_sim <- data_list[[1]]$ageLims[2,]
        } else{
          amin_sim <- data_list[[2]]$ageLims[1,]
          amax_sim <- data_list[[2]]$ageLims[2,]
        }
        
        params_simulation_2year[[i]][[j]][[k]]$nT <- nrow(pop_admin2_sim)
        params_simulation_2year[[i]][[j]][[k]]$nA <- ncol(pop_admin2_sim)
        params_simulation_2year[[i]][[j]][[k]]$lam_H <- extract_stanfit$lam_H[posterior_seed[[i]],]
        params_simulation_2year[[i]][[j]][[k]]$lam_t <- extract_stanfit$lam_t[posterior_seed[[i]],]
        params_simulation_2year[[i]][[j]][[k]]$rho <- extract_stanfit$rho[posterior_seed[[i]]]
        params_simulation_2year[[i]][[j]][[k]]$gamma <- extract_stanfit$gamma[posterior_seed[[i]]]
        params_simulation_2year[[i]][[j]][[k]]$pop <- pop_admin2_sim
        params_simulation_2year[[i]][[j]][[k]]$amin <- amin_sim
        params_simulation_2year[[i]][[j]][[k]]$amax <- amax_sim
        
      }
      
    }
    
  }
  
}

#### simulate seroprevalence at 9 years old, susceptible and monotypic
#### posterior simulations
#### standard models
init_year <- c(2017,2016)
age_groups <- list()
age_groups[[1]] <- c("0-4","5-9","10-14","15-19","20-44","45-54","55-64","65-74","75+")
age_groups[[2]] <- c("0-4","5-14","15-44","45+")
age_groups_jkt2 <- list()
age_groups_jkt2[[1]] <- c("0-4","5-14","15-44","45+")

posterior_simulation <- list()

for (i in c(2,4)){
  
  posterior_simulation[[i]] <- list()
  
  for (j in 1){
    
    posterior_simulation[[i]][[j]] <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      
      posterior_simulation[[i]][[j]][[k]] <- list()
      post_sim_summary <- list()
      post_sim_sero_summary <- list()
      
      for (l in seq_len(n_samples_sim)){
        
        nT_sim <- params_simulation[[i]][[j]][[k]]$nT
        
        post_sim <- simcases2(params_simulation[[i]][[j]][[k]]$nT,
                              params_simulation[[i]][[j]][[k]]$nA,
                              params_simulation[[i]][[j]][[k]]$lam_H[l,],
                              params_simulation[[i]][[j]][[k]]$lam_t[l,],
                              params_simulation[[i]][[j]][[k]]$rho[l],
                              params_simulation[[i]][[j]][[k]]$gamma[l],
                              params_simulation[[i]][[j]][[k]]$pop,
                              params_simulation[[i]][[j]][[k]]$amin,
                              params_simulation[[i]][[j]][[k]]$amax)
        
        sero9 <- 1 - post_sim$susc9_prop[1:nT_sim,]
        high_transmission <- ifelse(sero9>0.6,1,0)
        susc <- post_sim$susc_prop[1:nT_sim,] 
        mono <- post_sim$mono_prop[1:nT_sim,]
        colnames(susc) <- age_groups[[j]]
        colnames(mono) <- age_groups[[j]]
        
        susc_df <- as_tibble(susc) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="susc") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        mono_df <- as_tibble(mono) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="mono") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        post_sim_summary[[l]] <- susc_df %>% left_join(mono_df, 
                                                       by=c("rep","admin1","admin2","year","age_group"))
        
        post_sim_sero_summary[[l]] <- tibble(rep=l,
                                             admin1=provinces[j],
                                             admin2=admin2_province[k],
                                             year=init_year[j]-1+(1:nT_sim),
                                             sero9=sero9,
                                             high_transmission=high_transmission)
        
      }
      
      posterior_simulation[[i]][[j]][[k]] <- list()
      posterior_simulation[[i]][[j]][[k]]$susc_mono <- bind_rows(post_sim_summary)
      posterior_simulation[[i]][[j]][[k]]$sero <- bind_rows(post_sim_sero_summary)
      
    }
    
  }
  
}
saveRDS(posterior_simulation,"output/model_output/posterior_simulation.rds")

#### jakarta four age group
posterior_simulation_jkt2 <- list()

for (i in c(2,4)){
  
  posterior_simulation_jkt2[[i]] <- list()
  
  for (j in c(1)){
    
    posterior_simulation_jkt2[[i]][[j]] <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      
      posterior_simulation_jkt2[[i]][[j]][[k]] <- list()
      post_sim_summary <- list()
      post_sim_sero_summary <- list()
      
      for (l in seq_len(n_samples_sim)){
        
        nT_sim <- params_simulation_jkt2[[i]][[j]][[k]]$nT
        
        post_sim <- simcases2(params_simulation_jkt2[[i]][[j]][[k]]$nT,
                              params_simulation_jkt2[[i]][[j]][[k]]$nA,
                              params_simulation_jkt2[[i]][[j]][[k]]$lam_H[l,],
                              params_simulation_jkt2[[i]][[j]][[k]]$lam_t[l,],
                              params_simulation_jkt2[[i]][[j]][[k]]$rho[l],
                              params_simulation_jkt2[[i]][[j]][[k]]$gamma[l],
                              params_simulation_jkt2[[i]][[j]][[k]]$pop,
                              params_simulation_jkt2[[i]][[j]][[k]]$amin,
                              params_simulation_jkt2[[i]][[j]][[k]]$amax)
        
        sero9 <- 1 - post_sim$susc9_prop[1:nT_sim,]
        high_transmission <- ifelse(sero9>0.6,1,0)
        susc <- post_sim$susc_prop[1:nT_sim,] 
        mono <- post_sim$mono_prop[1:nT_sim,]
        colnames(susc) <- age_groups_jkt2[[j]]
        colnames(mono) <- age_groups_jkt2[[j]]
        
        susc_df <- as_tibble(susc) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="susc") %>% 
          mutate(age_group=factor(age_group,levels=age_groups_jkt2[[j]]))
        
        mono_df <- as_tibble(mono) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="mono") %>% 
          mutate(age_group=factor(age_group,levels=age_groups_jkt2[[j]]))
        
        post_sim_summary[[l]] <- susc_df %>% left_join(mono_df, 
                                                       by=c("rep","admin1","admin2","year","age_group"))
        
        post_sim_sero_summary[[l]] <- tibble(rep=l,
                                             admin1=provinces[j],
                                             admin2=admin2_province[k],
                                             year=init_year[j]-1+(1:nT_sim),
                                             sero9=sero9,
                                             high_transmission=high_transmission)
        
      }
      
      posterior_simulation_jkt2[[i]][[j]][[k]] <- list()
      posterior_simulation_jkt2[[i]][[j]][[k]]$susc_mono <- bind_rows(post_sim_summary)
      posterior_simulation_jkt2[[i]][[j]][[k]]$sero <- bind_rows(post_sim_sero_summary)
      
    }
    
  }
  
}
saveRDS(posterior_simulation_jkt2,"output/model_output/posterior_simulation_jkt2.rds")

#### 4 years of data
posterior_simulation_4year <- list()
init_year_4year <- c(2020,2020)

for (i in c(2,4)){
  
  posterior_simulation_4year[[i]] <- list()
  
  for (j in 1){
    
    posterior_simulation_4year[[i]][[j]] <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      
      posterior_simulation_4year[[i]][[j]][[k]] <- list()
      post_sim_summary <- list()
      post_sim_sero_summary <- list()
      
      for (l in seq_len(n_samples_sim)){
        
        nT_sim <- params_simulation_4year[[i]][[j]][[k]]$nT
        
        post_sim <- simcases2(params_simulation_4year[[i]][[j]][[k]]$nT,
                              params_simulation_4year[[i]][[j]][[k]]$nA,
                              params_simulation_4year[[i]][[j]][[k]]$lam_H[l,],
                              params_simulation_4year[[i]][[j]][[k]]$lam_t[l,],
                              params_simulation_4year[[i]][[j]][[k]]$rho[l],
                              params_simulation_4year[[i]][[j]][[k]]$gamma[l],
                              params_simulation_4year[[i]][[j]][[k]]$pop,
                              params_simulation_4year[[i]][[j]][[k]]$amin,
                              params_simulation_4year[[i]][[j]][[k]]$amax)
        
        sero9 <- 1 - post_sim$susc9_prop[1:nT_sim,]
        high_transmission <- ifelse(sero9>0.6,1,0)
        susc <- post_sim$susc_prop[1:nT_sim,] 
        mono <- post_sim$mono_prop[1:nT_sim,]
        colnames(susc) <- age_groups[[j]]
        colnames(mono) <- age_groups[[j]]
        
        susc_df <- as_tibble(susc) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year_4year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="susc") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        mono_df <- as_tibble(mono) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year_4year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="mono") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        post_sim_summary[[l]] <- susc_df %>% left_join(mono_df, 
                                                       by=c("rep","admin1","admin2","year","age_group"))
        
        post_sim_sero_summary[[l]] <- tibble(rep=l,
                                             admin1=provinces[j],
                                             admin2=admin2_province[k],
                                             year=init_year_4year[j]-1+(1:nT_sim),
                                             sero9=sero9,
                                             high_transmission=high_transmission)
        
      }
      
      posterior_simulation_4year[[i]][[j]][[k]] <- list()
      posterior_simulation_4year[[i]][[j]][[k]]$susc_mono <- bind_rows(post_sim_summary)
      posterior_simulation_4year[[i]][[j]][[k]]$sero <- bind_rows(post_sim_sero_summary)
      
    }
    
  }
  
}
saveRDS(posterior_simulation_4year,"output/model_output/posterior_simulation_4year.rds")

#### 2 years of data
posterior_simulation_2year <- list()
init_year_2year <- c(2022,2022)

for (i in c(2,4)){
  
  posterior_simulation_2year[[i]] <- list()
  
  for (j in 1){
    
    posterior_simulation_2year[[i]][[j]] <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      
      posterior_simulation_2year[[i]][[j]][[k]] <- list()
      post_sim_summary <- list()
      post_sim_sero_summary <- list()
      
      for (l in seq_len(n_samples_sim)){
        
        nT_sim <- params_simulation_2year[[i]][[j]][[k]]$nT
        
        post_sim <- simcases2(params_simulation_2year[[i]][[j]][[k]]$nT,
                              params_simulation_2year[[i]][[j]][[k]]$nA,
                              params_simulation_2year[[i]][[j]][[k]]$lam_H[l,],
                              params_simulation_2year[[i]][[j]][[k]]$lam_t[l,],
                              params_simulation_2year[[i]][[j]][[k]]$rho[l],
                              params_simulation_2year[[i]][[j]][[k]]$gamma[l],
                              params_simulation_2year[[i]][[j]][[k]]$pop,
                              params_simulation_2year[[i]][[j]][[k]]$amin,
                              params_simulation_2year[[i]][[j]][[k]]$amax)
        
        sero9 <- 1 - post_sim$susc9_prop[1:nT_sim,]
        high_transmission <- ifelse(sero9>0.6,1,0)
        susc <- post_sim$susc_prop[1:nT_sim,] 
        mono <- post_sim$mono_prop[1:nT_sim,]
        colnames(susc) <- age_groups[[j]]
        colnames(mono) <- age_groups[[j]]
        
        susc_df <- as_tibble(susc) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year_2year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="susc") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        mono_df <- as_tibble(mono) %>% 
          mutate(rep=l,
                 admin1=provinces[j],
                 admin2=admin2_province[k],
                 year=init_year_2year[j]-1+(1:nT_sim)) %>% 
          pivot_longer(-c(rep,admin1,admin2,year),names_to="age_group",values_to="mono") %>% 
          mutate(age_group=factor(age_group,levels=age_groups[[j]]))
        
        post_sim_summary[[l]] <- susc_df %>% left_join(mono_df, 
                                                       by=c("rep","admin1","admin2","year","age_group"))
        
        post_sim_sero_summary[[l]] <- tibble(rep=l,
                                             admin1=provinces[j],
                                             admin2=admin2_province[k],
                                             year=init_year_2year[j]-1+(1:nT_sim),
                                             sero9=sero9,
                                             high_transmission=high_transmission)
        
      }
      
      posterior_simulation_2year[[i]][[j]][[k]] <- list()
      posterior_simulation_2year[[i]][[j]][[k]]$susc_mono <- bind_rows(post_sim_summary)
      posterior_simulation_2year[[i]][[j]][[k]]$sero <- bind_rows(post_sim_sero_summary)
      
    }
    
  }
  
}
saveRDS(posterior_simulation_2year,"output/model_output/posterior_simulation_2year.rds")

#### read posterior simulations
posterior_simulation <- readRDS("output/model_output/posterior_simulation.rds")
posterior_simulation_jkt2 <- readRDS("output/model_output/posterior_simulation_jkt2.rds")
posterior_simulation_4year <- readRDS("output/model_output/posterior_simulation_4year.rds")
posterior_simulation_2year <- readRDS("output/model_output/posterior_simulation_2year.rds")

#### summarise posterior simulation
#### standard models
posterior_simulation_summary <- list()
posterior_simulation_sero_summary <- list()

for (i in c(2,4)){
  
  posterior_simulation_summary[[i]] <- list()
  posterior_simulation_sero_summary[[i]] <- list()
  
  for (j in 1){
    
    post_sim <- list()
    post_sim_sero <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      post_sim[[k]] <- posterior_simulation[[i]][[j]][[k]]$susc_mono
      post_sim_sero[[k]] <- posterior_simulation[[i]][[j]][[k]]$sero
    }
    
    simulation_summary <- bind_rows(post_sim) %>% 
      mutate(model=model_name[i]) %>% 
      group_by(model,admin1,admin2,year,age_group) %>% 
      summarise(susc_median=median(susc),
                susc_lo=quantile(susc,probs=c(0.025)),
                susc_up=quantile(susc,probs=c(0.975)),
                mono_median=median(mono),
                mono_lo=quantile(mono,probs=c(0.025)),
                mono_up=quantile(mono,probs=c(0.975))) %>% 
      ungroup()
    simulation_sero_summary <- bind_rows(post_sim_sero) %>% 
      mutate(model=model_name[i]) %>% 
      group_by(model,admin1,admin2,year) %>% 
      summarise(sero9_median=median(sero9),
                sero9_lo=quantile(sero9,probs=c(0.025)),
                sero9_up=quantile(sero9,probs=c(0.975)),
                high_transmission=sum(high_transmission)/n_samples_sim) %>% 
      ungroup()
    posterior_simulation_summary[[i]][[j]] <- simulation_summary %>% 
      mutate(model=factor(model,levels=model_name[c(2,4)]))
    posterior_simulation_sero_summary[[i]][[j]] <- simulation_sero_summary %>% 
      mutate(model=factor(model,levels=model_name[c(2,4)]))
  }
}

#### jakarta 4 age groups
posterior_simulation_summary_jkt2 <- list()
posterior_simulation_sero_summary_jkt2 <- list()

for (i in c(2,4)){
  
  posterior_simulation_summary_jkt2[[i]] <- list()
  posterior_simulation_sero_summary_jkt2[[i]] <- list()
  
  for (j in c(1)){
    
    post_sim <- list()
    post_sim_sero <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      post_sim[[k]] <- posterior_simulation_jkt2[[i]][[j]][[k]]$susc_mono
      post_sim_sero[[k]] <- posterior_simulation_jkt2[[i]][[j]][[k]]$sero
    }
    
    simulation_summary <- bind_rows(post_sim) %>% 
      mutate(model=model_name_jkt2[i]) %>% 
      group_by(model,admin1,admin2,year,age_group) %>% 
      summarise(susc_median=median(susc),
                susc_lo=quantile(susc,probs=c(0.025)),
                susc_up=quantile(susc,probs=c(0.975)),
                mono_median=median(mono),
                mono_lo=quantile(mono,probs=c(0.025)),
                mono_up=quantile(mono,probs=c(0.975))) %>% 
      ungroup()
    simulation_sero_summary <- bind_rows(post_sim_sero) %>% 
      mutate(model=model_name_jkt2[i]) %>% 
      group_by(model,admin1,admin2,year) %>% 
      summarise(sero9_median=median(sero9),
                sero9_lo=quantile(sero9,probs=c(0.025)),
                sero9_up=quantile(sero9,probs=c(0.975)),
                high_transmission=sum(high_transmission)/n_samples_sim) %>% 
      ungroup()
    posterior_simulation_summary_jkt2[[i]][[j]] <- simulation_summary %>% 
      mutate(model=factor(model,levels=model_name_jkt2[c(2,4)]))
    posterior_simulation_sero_summary_jkt2[[i]][[j]] <- simulation_sero_summary %>% 
      mutate(model=factor(model,levels=model_name_jkt2[c(2,4)]))
  }
}

#### 4 years
posterior_simulation_summary_4year <- list()
posterior_simulation_sero_summary_4year <- list()

for (i in c(2,4)){
  
  posterior_simulation_summary_4year[[i]] <- list()
  posterior_simulation_sero_summary_4year[[i]] <- list()
  
  for (j in 1){
    
    post_sim <- list()
    post_sim_sero <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      post_sim[[k]] <- posterior_simulation_4year[[i]][[j]][[k]]$susc_mono
      post_sim_sero[[k]] <- posterior_simulation_4year[[i]][[j]][[k]]$sero
    }
    
    simulation_summary <- bind_rows(post_sim) %>% 
      mutate(model=model_name_4year[i]) %>% 
      group_by(model,admin1,admin2,year,age_group) %>% 
      summarise(susc_median=median(susc),
                susc_lo=quantile(susc,probs=c(0.025)),
                susc_up=quantile(susc,probs=c(0.975)),
                mono_median=median(mono),
                mono_lo=quantile(mono,probs=c(0.025)),
                mono_up=quantile(mono,probs=c(0.975))) %>% 
      ungroup()
    simulation_sero_summary <- bind_rows(post_sim_sero) %>% 
      mutate(model=model_name_4year[i]) %>% 
      group_by(model,admin1,admin2,year) %>% 
      summarise(sero9_median=median(sero9),
                sero9_lo=quantile(sero9,probs=c(0.025)),
                sero9_up=quantile(sero9,probs=c(0.975)),
                high_transmission=sum(high_transmission)/n_samples_sim) %>% 
      ungroup()
    posterior_simulation_summary_4year[[i]][[j]] <- simulation_summary %>% 
      mutate(model=factor(model,levels=model_name_4year[c(2,4)]))
    posterior_simulation_sero_summary_4year[[i]][[j]] <- simulation_sero_summary %>% 
      mutate(model=factor(model,levels=model_name_4year[c(2,4)]))
  }
}

#### 2 years
posterior_simulation_summary_2year <- list()
posterior_simulation_sero_summary_2year <- list()

for (i in c(2,4)){
  
  posterior_simulation_summary_2year[[i]] <- list()
  posterior_simulation_sero_summary_2year[[i]] <- list()
  
  for (j in 1){
    
    post_sim <- list()
    post_sim_sero <- list()
    
    cases_province <- cases_all %>% filter(admin1 == provinces[j])
    admin2_province <- cases_province %>% pull(admin2) %>% unique()
    
    for (k in seq_len(length(admin2_province))){
      post_sim[[k]] <- posterior_simulation_2year[[i]][[j]][[k]]$susc_mono
      post_sim_sero[[k]] <- posterior_simulation_2year[[i]][[j]][[k]]$sero
    }
    
    simulation_summary <- bind_rows(post_sim) %>% 
      mutate(model=model_name_2year[i]) %>% 
      group_by(model,admin1,admin2,year,age_group) %>% 
      summarise(susc_median=median(susc),
                susc_lo=quantile(susc,probs=c(0.025)),
                susc_up=quantile(susc,probs=c(0.975)),
                mono_median=median(mono),
                mono_lo=quantile(mono,probs=c(0.025)),
                mono_up=quantile(mono,probs=c(0.975))) %>% 
      ungroup()
    simulation_sero_summary <- bind_rows(post_sim_sero) %>% 
      mutate(model=model_name_2year[i]) %>% 
      group_by(model,admin1,admin2,year) %>% 
      summarise(sero9_median=median(sero9),
                sero9_lo=quantile(sero9,probs=c(0.025)),
                sero9_up=quantile(sero9,probs=c(0.975)),
                high_transmission=sum(high_transmission)/n_samples_sim) %>% 
      ungroup()
    posterior_simulation_summary_2year[[i]][[j]] <- simulation_summary %>% 
      mutate(model=factor(model,levels=model_name_2year[c(2,4)]))
    posterior_simulation_sero_summary_2year[[i]][[j]] <- simulation_sero_summary %>% 
      mutate(model=factor(model,levels=model_name_2year[c(2,4)]))
  }
}

#### combine all posterior simulation summary
posterior_simulation_summary <- bind_rows(posterior_simulation_summary)
posterior_simulation_summary_jkt2 <- bind_rows(posterior_simulation_summary_jkt2)
posterior_simulation_summary_4year <- bind_rows(posterior_simulation_summary_4year)
posterior_simulation_summary_2year <- bind_rows(posterior_simulation_summary_2year)

posterior_simulation_summary_all <- bind_rows(posterior_simulation_summary,
                                              posterior_simulation_summary_jkt2,
                                              posterior_simulation_summary_4year,
                                              posterior_simulation_summary_2year)

posterior_simulation_sero_summary <- bind_rows(posterior_simulation_sero_summary)
posterior_simulation_sero_summary_jkt2 <- bind_rows(posterior_simulation_sero_summary_jkt2)
posterior_simulation_sero_summary_4year <- bind_rows(posterior_simulation_sero_summary_4year)
posterior_simulation_sero_summary_2year <- bind_rows(posterior_simulation_sero_summary_2year)

posterior_simulation_sero_summary_all <- bind_rows(posterior_simulation_sero_summary,
                                                   posterior_simulation_sero_summary_jkt2,
                                                   posterior_simulation_sero_summary_4year,
                                                   posterior_simulation_sero_summary_2year)

#### visualise
jkt_admin2_sero <- c("KOTA JAKARTA BARAT","KOTA JAKARTA SELATAN","KOTA JAKARTA TIMUR")

#### seroprevalence age 9
#### hierarchical
s_fig4_A <- posterior_simulation_sero_summary_all %>% 
  filter(admin2 %in% jkt_admin2_sero) %>%
  filter(!str_detect(model, 'non')) %>% 
  ggplot(aes(x=year,y=sero9_median,colour=model,fill=model)) +
  geom_hline(yintercept=0.6,linetype=2) +
  geom_line() +
  geom_ribbon(aes(ymin=sero9_lo,ymax=sero9_up),alpha=0.25) +
  facet_wrap(.~admin2) +
  theme(legend.position="bottom") +
  ylim(0,1) +
  theme_Publication(base_size=11) +
  labs(y="Seroprevalence at age 9",x="Year",tag="A",colour=NULL,fill=NULL) +
  scale_colour_manual(values=c("#000000","#E69F00","#56B4E9","#009E73")) +
  scale_fill_manual(values=c("#000000","#E69F00","#56B4E9","#009E73"))

#### non-hierarchical
s_fig5_A <- posterior_simulation_sero_summary_all %>% 
  filter(admin2 %in% jkt_admin2_sero) %>%
  filter(str_detect(model, 'non')) %>% 
  ggplot(aes(x=year,y=sero9_median,colour=model,fill=model)) +
  geom_hline(yintercept=0.6,linetype=2) +
  geom_line() +
  geom_ribbon(aes(ymin=sero9_lo,ymax=sero9_up),alpha=0.25) +
  facet_wrap(.~admin2) +
  theme(legend.position="bottom") +
  ylim(0,1) +
  theme_Publication(base_size=11) +
  labs(y="Seroprevalence at age 9",x="Year",tag="A",colour=NULL,fill=NULL) +
  scale_colour_manual(values=c("#000000","#E69F00","#56B4E9","#009E73")) +
  scale_fill_manual(values=c("#000000","#E69F00","#56B4E9","#009E73"))

s_fig4 <- s_fig4_A
# ggsave("output/figures/s_fig4.jpg",s_fig4,height=15,width=17,unit="cm")

s_fig5 <- s_fig5_A
# ggsave("output/figures/s_fig5.jpg",s_fig5,height=15,width=17,unit="cm")

#### map of probability of high transmission
#### neg binomial non-hierarchical
#### 2017
admin2_jkt_sero9_2017_shp <- admin2_jkt_shp %>% 
  left_join(posterior_simulation_sero_summary %>% 
              filter(model=="Negative binomial\nnon-hierarchical") %>% 
              filter(year==2017) %>% 
              left_join(cases_all) %>% 
              rename(A2C_BPS=idadmin2))

s_fig6_A <- admin2_jkt_sero9_2017_shp %>% 
  ggplot() +
  geom_sf(data=admin1_java_shp,fill="gray") +
  geom_sf(aes(fill=high_transmission),color = 'gray30', linewidth = 0.2) +
  geom_sf(data=admin1_java_shp,fill=NA,col="black", linewidth=0.5) +
  annotate(geom="rect",
           xmin=106.3831, xmax=106.8479, ymin=-6.040159, ymax=-5.084322, 
           fill=NA,col="black",linewidth=0.5,linetype=2) +
  annotate(geom="text",
           x=106.5831, y=-5.184322, 
           label="1.000",col="black",size=3) +
  theme_Publication(base_size=11) +
  scale_fill_viridis_c(limits = c(0.5, 1), oob = scales::squish) +
  coord_sf(xlim=c(106.3703, 108.8469),
           ylim=c(-7.820979, -5.184322)) +
  labs(fill="Prob. of high transmisson (2017)",x=NULL,y=NULL,tag="A")

#### 2023
admin2_jkt_sero9_2023_shp <- admin2_jkt_shp %>% 
  left_join(posterior_simulation_sero_summary %>% 
              filter(model=="Negative binomial\nnon-hierarchical") %>% 
              filter(year==2023) %>% 
              left_join(cases_all) %>% 
              rename(A2C_BPS=idadmin2))

s_fig6_B <- admin2_jkt_sero9_2023_shp %>% 
  ggplot() +
  geom_sf(data=admin1_java_shp,fill="gray") +
  geom_sf(aes(fill=high_transmission),color = 'gray30', linewidth = 0.2) +
  geom_sf(data=admin1_java_shp,fill=NA,col="black", linewidth=0.5) +
  annotate(geom="rect",
           xmin=106.3831, xmax=106.8479, ymin=-6.040159, ymax=-5.084322, 
           fill=NA,col="black",linewidth=0.5,linetype=2) +
  annotate(geom="text",
           x=106.5831, y=-5.184322, 
           label="0.996",col="black",size=3) +
  theme_Publication(base_size=11) +
  scale_fill_viridis_c(limits = c(0.5, 1), oob = scales::squish) +
  coord_sf(xlim=c(106.3703, 108.8469),
           ylim=c(-7.820979, -5.184322)) +
  labs(fill="Prob. of high transmisson (2023)",x=NULL,y=NULL,tag="B")

s_fig6 <- s_fig6_A / s_fig6_B
# ggsave("output/figures/s_fig6.jpg",s_fig6,height=20,width=17,unit="cm")

#### map of probability of high transmission
#### neg binomial hierarchical
#### 2017
admin2_jkt_sero9_hier_2017_shp <- admin2_jkt_shp %>% 
  left_join(posterior_simulation_sero_summary %>% 
              filter(model=="Negative binomial\nhierarchical") %>% 
              filter(year==2017) %>% 
              left_join(cases_all) %>% 
              rename(A2C_BPS=idadmin2))

s_fig7_A <- admin2_jkt_sero9_hier_2017_shp %>% 
  ggplot() +
  geom_sf(data=admin1_java_shp,fill="gray") +
  geom_sf(aes(fill=high_transmission),color = 'gray30', linewidth = 0.2) +
  geom_sf(data=admin1_java_shp,fill=NA,col="black", linewidth=0.5) +
  annotate(geom="rect",
           xmin=106.3831, xmax=106.8479, ymin=-6.040159, ymax=-5.084322, 
           fill=NA,col="black",linewidth=0.5,linetype=2) +
  annotate(geom="text",
           x=106.5831, y=-5.184322, 
           label="1.000",col="black",size=3) +
  theme_Publication(base_size=11) +
  scale_fill_viridis_c(limits = c(0.5, 1), oob = scales::squish) +
  coord_sf(xlim=c(106.3703, 108.8469),
           ylim=c(-7.820979, -5.184322)) +
  labs(fill="Prob. of high transmisson (2017)",x=NULL,y=NULL,tag="A")

#### 2023
admin2_jkt_sero9_hier_2023_shp <- admin2_jkt_shp %>% 
  left_join(posterior_simulation_sero_summary %>% 
              filter(model=="Negative binomial\nhierarchical") %>% 
              filter(year==2023) %>% 
              left_join(cases_all) %>% 
              rename(A2C_BPS=idadmin2))

s_fig7_B <- admin2_jkt_sero9_hier_2023_shp %>% 
  ggplot() +
  geom_sf(data=admin1_java_shp,fill="gray") +
  geom_sf(aes(fill=high_transmission),color = 'gray30', linewidth = 0.2) +
  geom_sf(data=admin1_java_shp,fill=NA,col="black", linewidth=0.5) +
  annotate(geom="rect",
           xmin=106.3831, xmax=106.8479, ymin=-6.040159, ymax=-5.084322, 
           fill=NA,col="black",linewidth=0.5,linetype=2) +
  annotate(geom="text",
           x=106.5831, y=-5.184322, 
           label="0.951",col="black",size=3) +
  theme_Publication(base_size=11) +
  scale_fill_viridis_c(limits = c(0.5, 1), oob = scales::squish) +
  coord_sf(xlim=c(106.3703, 108.8469),
           ylim=c(-7.820979, -5.184322)) +
  labs(fill="Prob. of high transmisson (2023)",x=NULL,y=NULL,tag="B")

s_fig7 <- s_fig7_A / s_fig7_B
# ggsave("output/figures/s_fig7.jpg",s_fig7,height=20,width=17,unit="cm")
