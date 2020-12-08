#################
# COVID-19 infection date versus report date 
# Citation: Goldstein ND, Quick H, Burstyn I. The impact of adjusting for misclassification in case ascertainment and uncertainty about date of infection in COVID-19 time series surveillance data on the estimated effective reproduction number. Manuscript in preparation.
# 9/1/20 -- Neal Goldstein
#################


### FUNCTIONS ###

library("R2jags") #JAGS interface
library("MCMCvis") #visualizations

#function to add CIs to plot
add.ci <- function(x, est, se, lvl=.95, trnc.lo=0) {
  for (i in 1:length(x)) {
    points(rep(x[i],2),
           pmax(est[i]+qnorm(1-(1-lvl)/2)*c(-1,1)*se[i],rep(trnc.lo,2)),
           type="l")
  }
}

#Bayesian model specification for true counts based on reported positives; see: https://pubmed.ncbi.nlm.nih.gov/32505172/
genmod.JAGS = function()
{
  
  ### prior distribution
  
  ### how much weight on the linear component for sens
  ### (earlier version corresponds to no weight)
  
  ### sn.wt <- 0 no smoothness (as per earlier versions)
  ### sn.wt <- 1 perfectly linear
  sn.wt ~ dunif(0.5,0.9)  ### some unspecified amount of smoothing
  
  ### linear component endpoints
  sn.LL ~ dunif(sn.wt*sn.lo[1], sn.wt*sn.hi[1])
  sn.LR ~ dunif(sn.wt*sn.lo[num.kn], sn.wt*sn.hi[num.kn])
  
  ### prev, sn piecewise linear
  ### parameterized by value at knots
  for (i in 1:num.kn) {
    r.kn[i] ~ dunif(0, r.hi[i])
    
    ### linear component
    sn.L[i] <- ((knts[i]-knts[1])*sn.LR+(knts[num.kn]-knts[i])*sn.LL)/
      (knts[num.kn]-knts[1])
    
    ### jumpy component
    sn.J[i] ~ dunif((1-sn.wt)*sn.lo[i],(1-sn.wt)*sn.hi[i]) 
    
    ### two together
    sn.kn[i] ~ dsum(sn.L[i], sn.J[i])
  }
  
  sp ~ dunif(sp.lo,1)
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      r[knts[i]+j] <- ((spc.kn[i]-j)*r.kn[i]+j*r.kn[i+1])/(spc.kn[i])
    }  
  }    
  r[knts[num.kn]] <- r.kn[num.kn]
  
  ### these imply the daily values
  for (i in 1:(num.kn-1)) {
    for (j in 0:(spc.kn[i]-1)) {
      sn[knts[i]+j] <- ((spc.kn[i]-j)*sn.kn[i]+j*sn.kn[i+1])/(spc.kn[i])
    }  
  }    
  sn[knts[num.kn]] <- sn.kn[num.kn]
  
  for (i in 1:(knts[num.kn])) {
    y[i] ~ dbinom(r[i], n[i])            ### true positives
    ystr1[i] ~ dbinom(sn[i], y[i])       ### correct positives
    ystr2[i] ~ dbinom(1-sp, n[i]-y[i])   ### false positives
    ystr[i] ~ sum(ystr1[i], ystr2[i])
  }
  
}


### READ DATA ###

#Data may be downloaded from https://covid19.colorado.gov/data
#See link at bottom "Access the data files" which takes one to a Google Drive. Go to "Case Data". Each file corresponds to daily counts of # tested (cumulative), cases by report date, cases by onset date. The most recent file includes all data to date, EXCEPT number of tested which must be obtained individually.

#to create data, obtain most recent cases_demographics_tests and abstract:
#date: align onset to reported dates
#reported: Cases of COVID-19 in Colorado by Date Reported to the State
#onset: Cases of COVID-19 in Colorado by Date of Illness Onset 
colorado_data = read.csv("colorado_data.csv", as.is=T, stringsAsFactors=F)

#set date type
colorado_data$date = as.Date(colorado_data$date, format="%m/%d/%y")

#trim to date ranges that are compatible with misclassification code (remove 0 counts at beginning and large increase at end)
colorado_data = colorado_data[5:210, ]


### BAYESIAN MODELING for TRUE CASE COUNTS ###

#active data set
dta = colorado_data

#population size of CO
pop_size = 5758736

#number of time points
T.end = nrow(dta)

#total observed positive results
ystr_report = dta$report[1:T.end]
ystr_onset = dta$onset[1:T.end]

#susceptible population at each time point
#n = dta$tested[1:T.end]
n_report = c(pop_size - (cumsum(dta$report) - dta$report[1]))
n_onset = c(pop_size - (cumsum(dta$onset) - dta$onset[1]))

#observed prevalence per time point
q.hat_report = ystr_report/n_report
q.hat_onset = ystr_onset/n_onset

#error for observed test positivity
se_report = sqrt(q.hat_report*(1-q.hat_report)/n_report)
se_onset = sqrt(q.hat_onset*(1-q.hat_onset)/n_onset)

#select interior knots based on linear trends
#knts = c(1, 15, 26, 55, 89, 140, T.end)
knts = c(1, T.end)
num.kn = length(knts)
spc.kn = knts[-1] - knts[-num.kn]

#plot observed data and knots (report)
plot(1:T.end, q.hat_report, 
     xlab=paste("Dates (",dta$date[1]," - ",dta$date[T.end],")",sep=""),
     ylab="Proportion Positive", ylim=c(0, max(q.hat_report+2.1*se_report)))
add.ci((1:T.end), q.hat_report, se_report)
points(knts, rep(0,num.kn),pch=17,col="red")

#specify hyperparameters: reported cases
sn.lo.report = rep(0.6, num.kn)  # same bound at each knot
sn.hi.report = rep(0.9, num.kn)  # same upper-bound at each knot
r.hi.report = rep(0.05, num.kn)   # same upper-bound at each knot
sp.lo.report = 0.95
line_data_report = list(knts=knts, num.kn=num.kn, spc.kn=spc.kn, sp.lo=sp.lo.report, sn.lo=sn.lo.report, sn.hi=sn.hi.report, r.hi=r.hi.report, ystr=ystr_report, n=n_report)

#specify hyperparameters: observed cases
sn.lo.onset = rep(0.6, num.kn)  # same bound at each knot
sn.hi.onset = rep(0.9, num.kn)  # same upper-bound at each knot
r.hi.onset = rep(0.05, num.kn)   # same upper-bound at each knot
sp.lo.onset = 0.95
line_data_onset = list(knts=knts, num.kn=num.kn, spc.kn=spc.kn, sp.lo=sp.lo.onset, sn.lo=sn.lo.onset, sn.hi=sn.hi.onset, r.hi=r.hi.onset, ystr=ystr_onset, n=n_onset)

#initial values (assigned within JAGS)
line_inits = function() {list(y=round(1.2*ystr), ystr1=ystr, ystr2=rep(0,length(ystr)))}
#line_inits_onset = function() {list(y=round(1.2*ystr_onset), ystr1=ystr_onset, ystr2=rep(0,length(ystr_onset)))}

#run JAGS with benchmarking
time1 = Sys.time()
model_report_posterior = jags.parallel(data=line_data_report, inits=line_inits, parameters.to.save=c("sp","sn.wt","sn.kn","r.kn","y"), model.file=genmod.JAGS, n.chains=4, n.thin=200, n.iter=400000, n.burnin=1000)
#model_onset_posterior = jags.parallel(data=line_data_onset, inits=line_inits, parameters.to.save=c("sp","sn.wt","sn.kn","r.kn","y"), model.file=genmod.JAGS, n.chains=4, n.thin=200, n.iter=400000, n.burnin=1000)
time2 = Sys.time()
time2-time1
rm(time1, time2)

#diagnostics
options(max.print=9999)
print(model_report_posterior)
MCMCtrace(model_report_posterior, params=c("sp","sn.wt","sn.kn","r.kn","y"), wd="./", filename="Appendix.pdf")

#extract true counts and add to original (dated) data
y_post_report = data.frame(model_report_posterior$BUGSoutput$summary[grep("y",substr(row.names(model_report_posterior$BUGSoutput$summary),1,1)),], stringsAsFactors=F)
#y_post_onset = data.frame(model_onset_posterior$BUGSoutput$summary[grep("y",substr(row.names(model_onset_posterior$BUGSoutput$summary),1,1)),], stringsAsFactors=F)

colorado_data$true_report_median = round(y_post_report$X50.)
colorado_data$true_report_lo = round(y_post_report$X2.5.)
colorado_data$true_report_hi = round(y_post_report$X97.5.)
#colorado_data$true_onset_median = round(y_post_onset$X50.)
#colorado_data$true_onset_lo = round(y_post_onset$X2.5.)
#colorado_data$true_onset_hi = round(y_post_onset$X97.5.)

#extract posterior counts from all chains for reported date
count_posterior = function(covid_posterior) {
  #extract each chain as a dataframe
  #nrow = (n.iter - n.burnin) / n.thin
  chain1 = as.data.frame(as.mcmc(covid_posterior)[[1]])
  chain2 = as.data.frame(as.mcmc(covid_posterior)[[2]])
  chain3 = as.data.frame(as.mcmc(covid_posterior)[[3]])
  chain4 = as.data.frame(as.mcmc(covid_posterior)[[4]])
  
  #add zip code specific posteriors to dataframe
  count_posterior = data.frame(matrix(NA,ncol=nrow(colorado_data),nrow=nrow(chain1)*4))
  col_mapping = colnames(count_posterior)[order(names(count_posterior))]
  for (i in 1:ncol(count_posterior)) {
    
    #determine correct column (offset is based on first y.t column in mcmc output)
    col_index = which(col_mapping==colnames(count_posterior)[i]) + (which(colnames(chain1)=="y[1]") - 1)
    
    #add data
    count_posterior[,i] = c(chain1[,col_index], chain2[,col_index], chain3[,col_index], chain4[,col_index])
  }
  rm(i,col_mapping,col_index,chain1,chain2,chain3,chain4)
  
  return(count_posterior)
}

model_report_posterior_counts = count_posterior(model_report_posterior)

#clean up
rm(dta, line_data_onset, line_data_report, y_post_report, y_post_onset, knts, n, n_onset, n_report, num.kn, q.hat_onset, q.hat_report, r.hi, r.hi.onset, r.hi.report, se_onset, se_report, sn.hi, sn.hi.onset, sn.hi.report, sn.lo, sn.lo.onset, sn.lo.report, sp.lo, sp.lo.onset, sp.lo.report, spc.kn, T.end, ystr, ystr_onset, ystr_report, add.ci, genmod.JAGS, line_inits, count_posterior)

#save
save.image("bayes_posterior.RData")


### READ DATA ###

load("bayes_posterior.RData")


### FUNCTIONS ###

library("EpiNow2") #deconvolution; see: https://epiforecasts.io/EpiNow2/


### NAIVE CURVES ###

#plot naive curves
palette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(colorado_data$date, colorado_data$report, type="l", lwd=3, xaxt="n", xlab="Date", ylab="Cases", ylim=c(0,1500), col=palette[2])
axis(side=1, at=seq(min(colorado_data$date), max(colorado_data$date), by=7), labels=format(seq(min(colorado_data$date), max(colorado_data$date), by=7), "%m-%d"), las=2, cex.axis=0.7)
lines(colorado_data$date, colorado_data$onset, col=palette[3], lwd=3, lty=1)
legend("topright", legend=c("Report date", "Onset date"), lty=1, col=palette[2:3], lwd=3, cex=0.6)


### EPINOW2 APPROACH TO RECOVERING INFECTION DATE ###

#https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2020.25.17.2000257
generation_time = list(mean = covid_generation_times[1, ]$mean,
                       mean_sd = covid_generation_times[1, ]$mean_sd,
                       sd = covid_generation_times[1, ]$sd,
                       sd_sd = covid_generation_times[1, ]$sd_sd,
                       max = 30)

#https://annals.org/aim/fullarticle/2762808/incubation-period-coronavirus-disease-2019-covid-19-from-publicly-reported
incubation_period = list(mean = covid_incubation_period[1, ]$mean,
                         mean_sd = covid_incubation_period[1, ]$mean_sd,
                         sd = covid_incubation_period[1, ]$sd,
                         sd_sd = covid_incubation_period[1, ]$sd_sd,
                         max = 30)

#https://www.denverpost.com/2020/07/24/colorado-covid-testing-delays-july/
reporting_delay = bootstrapped_dist_fit(rlnorm(100, log(3.5), 1))
reporting_delay$max = 30

#random sample from posteriors
set.seed(777)
rand_samples = sample(1:nrow(model_report_posterior_counts), 100, replace=F)

#final posterior samples
epinow_posterior = data.frame("sample"=NA, "date"=as.Date("2020-01-01"), "value"=NA, "model"=NA, stringsAsFactors=F)

for (i in 1:length(rand_samples)) {
  
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #data frame for projections
  reported_cases = data.frame("date"=colorado_data$date, "confirm"=as.integer(model_report_posterior_counts[rand_samples[i], ]))
  
  #projections
  #time1 = Sys.time()
  model_report_posterior_epinow = epinow(reported_cases = reported_cases, generation_time = generation_time,
                                         delays = list(incubation_period, reporting_delay),
                                         rt_prior = list(mean = 1, sd = 1),
                                         samples = 2000, warmup = 200, cores = ifelse(interactive(), 4, 1), chains = 4,
                                         verbose = TRUE, return_fit = TRUE)
  #time2 = Sys.time()
  #time2-time1
  #rm(time1, time2)
  
  #extract posteriors
  model_report_posterior_epinow_counts = model_report_posterior_epinow$estimates$samples[model_report_posterior_epinow$estimates$samples$variable=="infections"]
  model_report_posterior_epinow_counts = model_report_posterior_epinow_counts[model_report_posterior_epinow_counts$sample %in% sample(1:2000, length(rand_samples), replace=F), c("sample","date","value")]
  model_report_posterior_epinow_counts$model = rand_samples[i]
  
  #add to dataframe
  epinow_posterior = rbind(epinow_posterior, model_report_posterior_epinow_counts)
  
  #clean up
  rm(model_report_posterior_epinow, reported_cases, model_report_posterior_epinow_counts)
  gc()
  
}
rm(i, incubation_period, generation_time, reporting_delay, model_report_posterior_counts)
epinow_posterior = epinow_posterior[-1, ]

#unique time-series indicator
epinow_posterior$sample_model = paste(epinow_posterior$sample, epinow_posterior$model, sep="-")

#save
save.image("bayes_posterior_EpiNow.RData")


### ADJUSTED CURVE ###

#compute 2.5, 50, 97.5 quantiles
plot_med = aggregate(epinow_posterior$value, list(epinow_posterior$date), median)
plot_lo = aggregate(epinow_posterior$value, list(epinow_posterior$date), quantile, probs=0.025)
plot_hi = aggregate(epinow_posterior$value, list(epinow_posterior$date), quantile, probs=0.975)
plot_dates = unique(epinow_posterior$date)

palette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(plot_dates, plot_med$x, type="l", lwd=3, xaxt="n", xlab="Date", ylab="Cases", ylim=c(0,1500), col=palette[4])
axis(side=1, at=seq(min(plot_dates), max(plot_dates), by=7), labels=format(seq(min(plot_dates), max(plot_dates), by=7), "%m-%d"), las=2, cex.axis=0.7)
polygon(x=c(plot_dates,rev(plot_dates)), y=c(plot_lo$x,rev(plot_hi$x)), col=rgb(t(col2rgb(palette[4])), alpha=100, maxColorValue=255), border=NA)
lines(colorado_data$date, colorado_data$report, col=palette[2], lwd=3, lty=1)
legend("topleft", legend=c("Reported", "Adjusted"), lty=1, col=palette[c(2,4)], lwd=3, cex=0.6)

rm(plot_med,plot_lo,plot_hi,plot_dates)


### READ DATA ###

load("bayes_posterior_EpiNow.RData")


### FUNCTIONS ###

library("R0") #estimate R; see: https://pubmed.ncbi.nlm.nih.gov/23249562/
library("DescTools") #confidence intervals for correlations


### COMPARISON of Rt ###

#estimate Rt from observed data
gen_time = generation.time("gamma", c(3.64, 3.08))
observed_report_r = est.R0.TD(colorado_data$report, GT=gen_time, begin=as.integer(1), end=nrow(colorado_data))
observed_onset_r = est.R0.TD(colorado_data$onset, GT=gen_time, begin=as.integer(1), end=nrow(colorado_data))

#estimate Rt from imputed data
infected_r = matrix(data=NA, nrow=1, ncol=length(unique(epinow_posterior$date)))
infected_r_lo = matrix(data=NA, nrow=1, ncol=length(unique(epinow_posterior$date)))
infected_r_hi = matrix(data=NA, nrow=1, ncol=length(unique(epinow_posterior$date)))
imputed_ts = unique(epinow_posterior$sample_model)

for (i in 1:length(imputed_ts)) {
  
  cat("\n\n************** ","Observation: ",i," **************\n",sep="")
  
  #some runs of EpiNow may have resulted in invalid data
  this_rt = tryCatch(est.R0.TD(round(epinow_posterior$value[epinow_posterior$sample_model==imputed_ts[i]]), GT=gen_time, begin=as.integer(1), end=length(epinow_posterior$value[epinow_posterior$sample_model==imputed_ts[i]])), error=function(e) NULL)
  
  if (length(this_rt)>0) {
    infected_r = rbind(infected_r, this_rt$R)
    infected_r_lo = rbind(infected_r_lo, this_rt$conf.int[,1])
    infected_r_hi = rbind(infected_r_hi, this_rt$conf.int[,2])
  }
}
rm(i, this_rt, imputed_ts)

#create Rt dataframe
compare_r = data.frame("date"=colorado_data$date, "observed_report_r"=NA, "observed_report_r_lo"=NA, "observed_report_r_hi"=NA, "observed_onset_r"=NA, "observed_onset_r_lo"=NA, "observed_onset_r_hi"=NA, "infect_report_r"=NA, "infect_report_r_lo"=NA, "infect_report_r_hi"=NA, stringsAsFactors=F)

#date offset to align observed data with imputed data
date_offset = (which(unique(epinow_posterior$date)==compare_r$date[1])-1)

for (i in 1:nrow(compare_r)) {
  
  #observed, report
  compare_r$observed_report_r[i] = observed_report_r$R[i]
  compare_r$observed_report_r_lo[i] = observed_report_r$conf.int[i,1]
  compare_r$observed_report_r_hi[i] = observed_report_r$conf.int[i,2]
  
  #observed, onset
  compare_r$observed_onset_r[i] = observed_onset_r$R[i]
  compare_r$observed_onset_r_lo[i] = observed_onset_r$conf.int[i,1]
  compare_r$observed_onset_r_hi[i] = observed_onset_r$conf.int[i,2]
  
  #infect, report
  compare_r$infect_report_r[i] = median(infected_r[, (i+date_offset)], na.rm=T)
  compare_r$infect_report_r_lo[i] = median(infected_r_lo[, (i+date_offset)], na.rm=T)
  compare_r$infect_report_r_hi[i] = median(infected_r_hi[, (i+date_offset)], na.rm=T)
}
rm(i, gen_time, date_offset, observed_report_r, observed_onset_r, infected_r, infected_r_lo, infected_r_hi)

#trim last day due to unstable estimates
compare_r = compare_r[-nrow(compare_r), ]

#save
save.image("bayes_posterior_EpiNow_R.RData")


### READ DATA ###

load("bayes_posterior_EpiNow_R.RData")


### Rt PLOTS ###

#full time series
palette = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
plot(compare_r$date, compare_r$observed_report_r, type="l", lwd=3, xaxt="n", xlab="Date", ylab="Rt", ylim=c(0.5,4), col=palette[2], main="A)")
axis(side=1, at=seq(min(compare_r$date), max(compare_r$date), by=7), labels=format(seq(min(compare_r$date), max(compare_r$date), by=7), "%m-%d"), las=2, cex.axis=0.7)
polygon(x=c(compare_r$date,rev(compare_r$date)), y=c(compare_r$observed_report_r_lo,rev(compare_r$observed_report_r_hi)), col=rgb(t(col2rgb(palette[2])), alpha=100, maxColorValue=255), border=NA)
lines(compare_r$date, compare_r$observed_onset_r, col=palette[3], lwd=3, lty=1)
polygon(x=c(compare_r$date,rev(compare_r$date)), y=c(compare_r$observed_onset_r_lo,rev(compare_r$observed_onset_r_hi)), col=rgb(t(col2rgb(palette[3])), alpha=100, maxColorValue=255), border=NA)
lines(compare_r$date, compare_r$infect_report_r, col=palette[4], lwd=3, lty=1)
polygon(x=c(compare_r$date,rev(compare_r$date)), y=c(compare_r$infect_report_r_lo,rev(compare_r$infect_report_r_hi)), col=rgb(t(col2rgb(palette[4])), alpha=100, maxColorValue=255), border=NA)
abline(h=1, lty=2, lwd=2)
legend("topright", legend=c("Report date", "Onset date", "Adjusted infection date"), lty=1, col=palette[2:4], lwd=3, cex=0.6)

#post March time series
plot(compare_r$date[compare_r$date>=as.Date("2020-04-01")], compare_r$observed_report_r[compare_r$date>=as.Date("2020-04-01")], type="l", lwd=3, xaxt="n", xlab="Date", ylab="Rt", ylim=c(0.5,1.8), col=palette[2], main="B)")
axis(side=1, at=seq(min(compare_r$date[compare_r$date>=as.Date("2020-04-01")]), max(compare_r$date[compare_r$date>=as.Date("2020-04-01")]), by=7), labels=format(seq(min(compare_r$date[compare_r$date>=as.Date("2020-04-01")]), max(compare_r$date[compare_r$date>=as.Date("2020-04-01")]), by=7), "%m-%d"), las=2, cex.axis=0.7)
polygon(x=c(compare_r$date[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$date[compare_r$date>=as.Date("2020-04-01")])), y=c(compare_r$observed_report_r_lo[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$observed_report_r_hi[compare_r$date>=as.Date("2020-04-01")])), col=rgb(t(col2rgb(palette[2])), alpha=100, maxColorValue=255), border=NA)
lines(compare_r$date[compare_r$date>=as.Date("2020-04-01")], compare_r$observed_onset_r[compare_r$date>=as.Date("2020-04-01")], col=palette[3], lwd=3, lty=1)
polygon(x=c(compare_r$date[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$date[compare_r$date>=as.Date("2020-04-01")])), y=c(compare_r$observed_onset_r_lo[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$observed_onset_r_hi[compare_r$date>=as.Date("2020-04-01")])), col=rgb(t(col2rgb(palette[3])), alpha=100, maxColorValue=255), border=NA)
lines(compare_r$date[compare_r$date>=as.Date("2020-04-01")], compare_r$infect_report_r[compare_r$date>=as.Date("2020-04-01")], col=palette[4], lwd=3, lty=1)
polygon(x=c(compare_r$date[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$date[compare_r$date>=as.Date("2020-04-01")])), y=c(compare_r$infect_report_r_lo[compare_r$date>=as.Date("2020-04-01")],rev(compare_r$infect_report_r_hi[compare_r$date>=as.Date("2020-04-01")])), col=rgb(t(col2rgb(palette[4])), alpha=100, maxColorValue=255), border=NA)
abline(h=1, lty=2, lwd=2)
legend("topright", legend=c("Report date", "Onset date", "Adjusted infection date"), lty=1, col=palette[2:4], lwd=3, cex=0.6)

summary(compare_r$observed_report_r)
summary(compare_r$observed_onset_r)
summary(compare_r$infect_report_r)

compare_r$observed_report_r[compare_r$date==max(compare_r$date)]
compare_r$observed_report_r_lo[compare_r$date==max(compare_r$date)]
compare_r$observed_report_r_hi[compare_r$date==max(compare_r$date)]
compare_r$observed_onset_r[compare_r$date==max(compare_r$date)]
compare_r$observed_onset_r_lo[compare_r$date==max(compare_r$date)]
compare_r$observed_onset_r_hi[compare_r$date==max(compare_r$date)]
compare_r$infect_report_r[compare_r$date==max(compare_r$date)]
compare_r$infect_report_r_lo[compare_r$date==max(compare_r$date)]
compare_r$infect_report_r_hi[compare_r$date==max(compare_r$date)]

#to compare: http://vassarstats.net/rdiff.html
SpearmanRho(compare_r$observed_report_r, compare_r$infect_report_r, conf.level=0.95)
SpearmanRho(compare_r$observed_onset_r, compare_r$infect_report_r, conf.level=0.95)


### DESCRIPTIVE STATS ###

#observed
sum(colorado_data$report)
sum(colorado_data$report) / pop_size
#sum(colorado_data$tested)
#sum(colorado_data$reported) / sum(colorado_data$tested)
#min(colorado_data$date[colorado_data$report>=1])
#min(colorado_data$date[colorado_data$onset>=1])

#posterior counts and missed cases
sum(colorado_data$true_report_median)
sum(colorado_data$true_report_lo)
sum(colorado_data$true_report_hi)
(sum(colorado_data$true_report_median) - sum(colorado_data$report))
(sum(colorado_data$true_report_hi) - sum(colorado_data$report))

#posterior prevalence
sum(colorado_data$true_report_median) / pop_size
sum(colorado_data$true_report_lo) / pop_size
sum(colorado_data$true_report_hi) / pop_size

#sensitivity and specificity
model_report_posterior$BUGSoutput$summary[grep("sp", row.names(model_report_posterior$BUGSoutput$summary)), ]
model_report_posterior$BUGSoutput$summary[grep("sn", row.names(model_report_posterior$BUGSoutput$summary)), ]

min(epinow_posterior$date)
