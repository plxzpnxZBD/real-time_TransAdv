


# control.scale = 3.3
# eta = 1.3
# candidate.0 = 0
# candidate.1 = 0.3527102
# obs.case = 2
# num.AA.0 = 0
# num.AA.1 = 0


get.one.R.est.with.fixed.eta = function(
  num.AA.0, num.AA.1, 
  candidate.0, candidate.1, obs.case, 
  eta, control.scale = 3.3
){
  crude.R = obs.case / (candidate.0 + candidate.1 * eta)
  #
  this.R.min = crude.R /control.scale; this.R.min = ifelse(this.R.min < 0.0001, 0.0001, this.R.min)
  this.R.max = crude.R *control.scale#; this.R.max = ifelse(this.R.max > 10, 10, this.R.max)
  R.sample.array = c(seq(this.R.min, crude.R, length.out = 101),
                     seq(crude.R, this.R.max, length.out = 101))
  expect.case.array = (candidate.0 + candidate.1 * eta) *R.sample.array
  expect.case.0.array = (candidate.0) *R.sample.array
  expect.case.1.array = (candidate.1 * eta) *R.sample.array
  #
  expect.prev.AA.0 = (candidate.0) / (candidate.0 + candidate.1 * eta)
  expect.prev.AA.1 = (candidate.1 * eta) / (candidate.0 + candidate.1 * eta)
  expect.prev.array = c(expect.prev.AA.0, expect.prev.AA.1)
  #
  case.ll.array = dpois(x = obs.case, lambda = expect.case.array, log = T)
  prev.ll = log(expect.prev.AA.0)*num.AA.0 + log(expect.prev.AA.1)*num.AA.1; prev.ll = ifelse(is.infinite(prev.ll) | is.na(prev.ll), 0, prev.ll)
  ll.array = case.ll.array + prev.ll
  #   plot(R.sample.array, ll.array)
  cutoff.ll = max(ll.array, na.rm = T) - qchisq(p = 0.975, df = 1); sel.index.array = which(ll.array > cutoff.ll)
  R.est = median(R.sample.array[which.max(ll.array)]); R.ci.array = range(R.sample.array[sel.index.array])
  case.est = median(expect.case.array[which.max(ll.array)]); expect.case.ci.array = range(expect.case.array[sel.index.array])
  #
  expect.prev.AA.1.lwr = expect.prev.AA.1 +c(-1)*1.96*sqrt(expect.prev.AA.0 *expect.prev.AA.1 /case.est); expect.prev.AA.1.lwr = ifelse(expect.prev.AA.1.lwr <0, 0, expect.prev.AA.1.lwr)
  expect.prev.AA.1.upr = expect.prev.AA.1 +c(1)*1.96*sqrt(expect.prev.AA.0 *expect.prev.AA.1 /case.est); expect.prev.AA.1.upr = ifelse(expect.prev.AA.1.upr >1, 1, expect.prev.AA.1.upr)
  expect.prev.AA.1.array = c(expect.prev.AA.1, expect.prev.AA.1.lwr, expect.prev.AA.1.upr)
  #
  if((this.R.min == min(R.ci.array) & this.R.min > 0) | this.R.max == max(R.ci.array)){
    est.array = get.one.R.est.with.fixed.eta(
      num.AA.0 = num.AA.0, num.AA.1 = num.AA.1, 
      candidate.0 = candidate.0, candidate.1 = candidate.1, obs.case = obs.case, 
      eta = eta, control.scale = control.scale *2)
  } else if(crude.R == 0) {
    #, rep(0,3)
    est.array = c(rep(0,3), rep(NA,3), NA)
  } else {
    #case.est, expect.case.ci.array, 
    est.array = c(R.est, R.ci.array, expect.prev.AA.1.array, max(ll.array, na.rm = T))
  }
  names(est.array) = c('R.est', 'R.lwr', 'R.upr', 'prev.AA.1.est', 'prev.AA.1.lwr', 'prev.AA.1.upr', 'll')
  return(c(est.array))
}





# sel.data = na.omit(CA.data)
# num.AA.0.array = c(sel.data$num.AA.0)
# num.AA.1.array = c(sel.data$num.AA.1)
# candidate.0.array = c(sel.data$seed.candidate.0)
# candidate.1.array = c(sel.data$seed.candidate.1)
# obs.case.array = c(sel.data$daily.case)
# eta = 1.3
# control.scale = 3.3




get.all.R.est.with.fixed.eta = function(
  num.AA.0.array, num.AA.1.array, 
  candidate.0.array, candidate.1.array, obs.case.array, 
  eta, control.scale = 3.3
){
  #if(length(candidate.array) != length(obs.death.array)) stop('different lengths!')
  #
  est.mat = NULL
  for(case.j in 1:length(obs.case.array)){#      case.j = 13
    #
    temp.est.array = get.one.R.est.with.fixed.eta(
      num.AA.0 = num.AA.0.array[case.j], num.AA.1 = num.AA.1.array[case.j], 
      candidate.0 = candidate.0.array[case.j], candidate.1 = candidate.1.array[case.j], obs.case = obs.case.array[case.j], 
      eta = eta, control.scale = control.scale
    )
    est.mat = rbind(est.mat, c(temp.est.array))
  }
  est.mat = as.data.frame(est.mat)
  colnames(est.mat) = c('R.est', 'R.lwr', 'R.upr', 'prev.AA.1.est', 'prev.AA.1.lwr', 'prev.AA.1.upr', 'll')
  return(est.mat)
}#      sel.data[9,]









