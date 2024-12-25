

simple <- function(threshold=30,
                   mu=4,
                   sigma2=1*1){
  
dOne <- generateBiomarkerData(
  mu = c("logRES" = mu),
  sigma2 = sigma2,
  rho = 0.7,
  seed = floor(.Machine$integer.max * runif(1))
)

dOneWide <- dOne %>% 
  pivot_wider(
    values_from = Value,
    names_from = TimePoint
  ) 

dOut <- dOneWide %>% analyseTrial(threshold = threshold)

dOut
}

#simple()

dOut1 <- plyr::rdply(1000, simple()) %>% add_column(scenario="mu=4, sd=1")
dOut2 <- plyr::rdply(1000, simple(mu=1.25, sigma2=0.8*0.8))%>% add_column(scenario="mu=1.25, sd=0.8")
dOut3 <- plyr::rdply(1000, simple(mu=10, sigma2=1*1))%>% add_column(scenario="mu=10, sd=1")
dOut4 <- plyr::rdply(1000, simple(mu=4, sigma2=2*2))%>% add_column(scenario="mu=4, sd=2")
dOut5 <- plyr::rdply(1000, simple(mu=1.25, sigma2=1*1))%>% add_column(scenario="mu=1.25, sd=1")

dAll <- bind_rows(dOut1, dOut2, dOut3, dOut4, dOut5) %>% mutate(Effect=DeltaPct_Control - DeltaPct_Test)

dSummary <- dAll  %>% 
  group_by(scenario) %>% 
  summarise(across(
  c(Effect), 
    #Baseline_SD_Control, Baseline_SD_Test, Endpoint_SD_Control, Endpoint_SD_Test, Delta_SD_Control, Delta_SD_Test),
  list(Mean = \(x) mean(x, na.rm = TRUE),
       Median = \(x) quantile(x, 0.5, na.rm = TRUE),
       SD = \(x) sd(x, na.rm = TRUE)
))) 

dAll %>% group_by(scenario) %>% 
  ggplot(aes(x=Effect,  y=after_stat(density), color=scenario))  + geom_density(alpha=0.3) +
  xlim(-200, 200)
