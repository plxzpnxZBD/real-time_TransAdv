
source('functions_for_Rt.R')



#load(file = '[set]GT_distribution.RData')
load(file = '[comp]CA_COVID_data.RData')


#plot(CA.data$cal.date, CA.data$prev.AA.1)





eta.array = c(seq(0.1, 1, length.out = 101),
              seq(1, 10, length.out = 101))
sel.data = na.omit(CA.data)
#
ll.array = NULL
for(eta.i in 1:length(eta.array)){#            eta.i = 1
  temp.eta = eta.array[eta.i]
  temp.outcome.data = get.all.R.est.with.fixed.eta(
    num.AA.0.array = c(sel.data$num.AA.0), num.AA.1.array = c(sel.data$num.AA.1), 
    candidate.0.array = c(sel.data$seed.candidate.0), candidate.1.array = c(sel.data$seed.candidate.1), obs.case.array = c(sel.data$daily.case), 
    eta = temp.eta, control.scale = 3.3
  )
  temp.ll = sum(temp.outcome.data$ll, na.rm = T)
  ll.array = c(ll.array, temp.ll)
}#   par(las = 1); plot(eta.array, ll.array)
max(ll.array)
eta.array[which.max(ll.array)]
range(eta.array[ll.array >(max(ll.array) -5)])








sel.data = na.omit(CA.data)
#
outcome.data = get.all.R.est.with.fixed.eta(
  num.AA.0.array = c(sel.data$num.AA.0), num.AA.1.array = c(sel.data$num.AA.1), 
  candidate.0.array = c(sel.data$seed.candidate.0), candidate.1.array = c(sel.data$seed.candidate.1), obs.case.array = c(sel.data$daily.case), 
  eta = 1.54, control.scale = 3.3
)
outcome.data$cal.date = sel.data$cal.date
outcome.data$cal.time = sel.data$cal.time
sel.outcome.data = na.omit(subset(outcome.data, subset = R.est < 5))










best.eta.est = 1.52; lwr.eta.est = 1.36; upr.eta.est = 1.76
par(las = 1, mar = c(2.5,5,1.5,1), mfrow = c(3,1))

plot(1,1, type = 'n', xlim = c(2020.1,2020.8), ylim = c(0,14000), xaxs = 'i', yaxs = 'i', axe = F, ann = F, frame = T)
points(CA.data$cal.time, CA.data$daily.case, col = 'royalblue', type = 'S', lwd = 1, lend = 1)
axis(2); axis(1, at = 2020 +0:11/12 +1/24, tick = F, labels = month.abb)
axis(1, at = 2020 +0:12/12, tick = T, labels = F)
mtext(side = 2, line = 3.5, las = 0, text = 'daily # of cases')
mtext(side = 3, adj = 0, text = '(A) # of COVID-19 cases (by reporting date)', cex = 1.2)

#

plot(1,1, type = 'n', xlim = c(2020.1,2020.8), ylim = c(0,5), xaxs = 'i', axe = F, ann = F, frame = T)
for(date.i in 1:nrow(sel.outcome.data)){
  lines(rep(sel.outcome.data$cal.time[date.i], 2) -0.5/366, c(sel.outcome.data$R.lwr[date.i], sel.outcome.data$R.upr[date.i]), col = 'green3')
}
for(date.i in 1:nrow(sel.outcome.data)){
  lines(rep(sel.outcome.data$cal.time[date.i], 2) +0.5/366, c(sel.outcome.data$R.lwr[date.i], sel.outcome.data$R.upr[date.i]) *best.eta.est, col = 'darkorange')
}
points(sel.outcome.data$cal.time -0.5/366, sel.outcome.data$R.est, col = 'green3', pch = 15)
points(sel.outcome.data$cal.time +0.5/366, sel.outcome.data$R.est *best.eta.est, col = 'darkorange', pch = 15)
abline(h = 1, col = 'red', lty = 2)
axis(2); axis(1, at = 2020 +0:11/12 +1/24, tick = F, labels = month.abb)
axis(1, at = 2020 +0:12/12, tick = T, labels = F)
mtext(side = 2, line = 3, las = 0, text = 'reproduction number')
legend('topright', col = c('darkorange','green3'), lty = 1, pch = 15, legend = c('G', 'D'), bty = 'n')
mtext(side = 3, adj = 0, text = '(B) instantaneous reproduction number', cex = 1.2)

#

plot(1,1, type = 'n', xlim = c(2020.1,2020.8), ylim = c(0,1), xaxs = 'i', yaxs = 'r', axe = F, ann = F, frame = T)
polygon(
  x = c(sel.outcome.data$cal.time, rev(sel.outcome.data$cal.time)), 
  y = c(sel.outcome.data$prev.AA.1.lwr, rev(sel.outcome.data$prev.AA.1.upr)),
  border = NA, col = 'pink'
)
weekly.CA.data = CA.data[seq(1,nrow(CA.data), by = 7),]
points(weekly.CA.data$cal.time, weekly.CA.data$prev.AA.1, col = 'orange', pch = 20, cex = sqrt(weekly.CA.data$num.for.prev) /7)
#lines(sel.outcome.data$cal.time, sel.outcome.data$prev.AA.1, col = 'darkgreen')
lines(sel.outcome.data$cal.time, sel.outcome.data$prev.AA.1.est, col = 'red', lwd = 1)
axis(2); axis(1, at = 2020 +0:11/12 +1/24, tick = F, labels = month.abb)
axis(1, at = 2020 +0:12/12, tick = T, labels = F)
mtext(side = 2, line = 3, las = 0, text = 'proportion of G')
legend('bottom', col = c('orange','red'), lty = c(NA, 1), pch = c(19, NA), legend = c('observed (weekly)', 'fitted'), bty = 'n')
legend('bottomright', col = 'black', lwd = NA, pch = 1, pt.cex = sqrt(c(30,100,300)) /7, legend = c('n = 30', 'n = 100', 'n = 300'), bty = 'n')
mtext(side = 3, adj = 0, text = '(C) proportion of G at the 614th codon (by reporting date)', cex = 1.2)





