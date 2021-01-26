



load(file = '[set]GT_distribution.RData')
load(file = '[proc]CA_COVID_data.RData')






CA.data$daily.case
GT.lag = max(c(length(GT.dist), length(GT.dist.alt)))
#
CA.data$seed.candidate.0 = NA
CA.data$seed.candidate.1 = NA
CA.data$seed.candidate.all = NA
for(date.i in (GT.lag +1):nrow(CA.data)){#               date.i = 16
  temp.data = CA.data[(date.i -GT.lag) : (date.i -1),]
  CA.data$seed.candidate.0[date.i] = sum(rev(temp.data$daily.case *temp.data$prev.AA.0)[1:length(GT.dist)] *GT.dist)
  CA.data$seed.candidate.1[date.i] = sum(rev(temp.data$daily.case *temp.data$prev.AA.1)[1:length(GT.dist.alt)] *GT.dist.alt)
  CA.data$seed.candidate.all[date.i] = CA.data$seed.candidate.0[date.i] + CA.data$seed.candidate.1[date.i]
}



#                 save(CA.data, file = '[comp]CA_COVID_data.RData')




