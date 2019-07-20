# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

btmle <- function(Y, A, W) {
  require(gtools)

  df = data.frame(W, A, Y)

  # Q estimation:: model Y | A,W, parameters: theta
  theta = coef(glm(Y ~ A + W1 + W2 + W3, family=binomial(link='logit'), data=df))
  #theta = d.glm(W = df[,c('A','W1','W2','W3')],Y = df[,'Y'], n.worker = 8)

  # compute the predicted Y hat for each patient
  Ya = apply(cbind(1,df[,c('A','W1','W2','W3')]), MARGIN = 1, FUN = function(x){inv.logit(x %*% theta)})
  df$Ya = Ya
  # compute Y1 hat and Y0 hat by setting A to 1 and 0 respectively
  Y1 = apply(cbind(1,1,df[,c('W1','W2','W3')]), MARGIN = 1, FUN = function(x){inv.logit(x %*% theta)})
  Y0 = apply(cbind(1,0,df[,c('W1','W2','W3')]), MARGIN = 1, FUN = function(x){inv.logit(x %*% theta)})
  df$Y1 = Y1
  df$Y0 = Y0
  # g estimation:: model A | W, parameters: beta
  beta = coef(glm(A ~ W1 + W2 + W3, family=binomial(link='logit'), data=df))
  #beta = d.glm(W = df[,c('W1','W2','W3')],Y = df[,'A'], n.worker = 8)

  # compute pi0 and pi1
  pi1 = apply(cbind(1,df[,c('W1','W2','W3')]), MARGIN = 1, FUN = function(x){inv.logit(x %*% beta)})
  pi0 = 1-pi1
  df$pi1 = pi1
  df$pi0 = pi0

  # compute H0 and H1
  H1 = 1/pi1
  H0 = -1/pi0
  df$h1 = H1
  df$h0 = H0

  # function to compute Ha based on the value of the treatment A
  ha = function(a, p1, p0) {
    if (a == 1){
      return(1/p1)
    } else {
      return(-1/p0)
    }
  }

  # compute Ha
  Ha = apply(df[,c('A','pi1','pi0')], MARGIN = 1, FUN = function(x){ha(x['A'],x['pi1'],x['pi0'])})
  df$Ha = Ha

  # model logit(Y*) = logit(Y) + sigma*Ha, parameter: sigma
  sigma = coef(glm(Y ~ Ha -1 + offset(logit(Ya)), family=binomial(link='logit'), data=df))
  #sigma = d.glm(W = df[,c('Ha')], Y = df[,'Y'], n.worker = 8, offset = logit(df[,'Ya']))

  # adjust Y0 and Y1 using the estimated fluctuation parameter sigma
  Y1.star = inv.logit(logit(Y1) + sigma*H1)
  Y0.star = inv.logit(logit(Y0) + sigma*H0)
  df$Y1.star = Y1.star
  df$Y0.star = Y0.star

  mu1.true = mean(Y1.true)
  mu0.true = mean(Y0.true)
  mu1 = mean(Y1.star)
  mu0 = mean(Y0.star)

  # compute Additive Treatment Effect
  ATE.true = mu1.true-mu0.true
  ATE.true

  ATE = mu1-mu0
  ATE
  IC.ATE = (A/pi1 - (1-A)/pi0)*(Y-Ya) + Y1 - Y0 - ATE
  ATE.var = var(IC.ATE)/n
  ATE.CI = c(ATE - (1.96 * sqrt(var(IC.ATE)/n)), ATE + (1.96 * sqrt(var(IC.ATE)/n)))
  ATE.CI
  ATE.p = 2*pnorm(-abs(ATE/sqrt(var(IC.ATE)/n)))
  ATE.p
  ATE.se = sqrt(var(IC.ATE)/n)
  ATE.se

  # compute Relative Risk
  RR.true = mu1.true/mu0.true
  RR.true
  RR = mu1/mu0
  RR
  log.RR = log(RR)
  IC.logRR = 1/mu1*((A/pi1)*(Y-Ya) + Y1 - mu1) - 1/mu0*((1-A)/pi0*(Y-Ya)+ Y0 - mu0)
  logRR.var = var(IC.logRR)/n
  RR.CI = c(exp(log(RR) -1.96 *sqrt(var(IC.logRR)/n)), exp(log(RR) +1.96 *sqrt(var(IC.logRR)/n)))
  logRR.CI = c(log(RR) -1.96 *sqrt(var(IC.logRR)/n), log(RR) +1.96 *sqrt(var(IC.logRR)/n))
  logRR.se = sqrt(var(IC.logRR)/n)

  # compute Odds Ratio
  OR.true = mu1.true/(1-mu1.true)/(mu0.true/(1-mu0.true))
  OR.true
  OR = mu1/(1-mu1)/(mu0/(1-mu0))
  OR
  log.OR = log(OR)
  IC.logOR = 1/(mu1*(1-mu1)) * (A/pi1*(Y-Ya) + Y1) - 1/(mu0*(1-mu0)) * ((1-A)/pi0*(Y-Ya) + Y0)
  OR.CI = c(exp(log(OR) -1.96 *sqrt(var(IC.logOR)/n)), exp(log(OR) +1.96 *sqrt(var(IC.logOR)/n)))
  logOR.CI = c(log(OR) -1.96 *sqrt(var(IC.logOR)/n), log(OR) +1.96 *sqrt(var(IC.logOR)/n))
  logOR.se = sqrt(var(IC.logOR)/n)

  print(ATE)
  print(ATE.CI)
  print(TMLE$ATE$psi)
  print(TMLE$ATE$CI)

  results = data.frame(psi=c("ATE","log(RR)","log(OR)"),
                       true=c(ATE.true,log(RR.true),log(OR.true)),
                       TMLE.estimate=c(TMLE$ATE$psi,TMLE$RR$log.psi,TMLE$OR$log.psi),
                       TMLE.se=c(sqrt(TMLE$ATE$var.psi),sqrt(TMLE$RR$var.log.psi),sqrt(TMLE$OR$var.log.psi)),
                       TMLE.l95=c(TMLE$ATE$CI[1],log(TMLE$RR$CI[1]),log(TMLE$OR$CI[1])),
                       TMLE.u95=c(TMLE$ATE$CI[2],log(TMLE$RR$CI[2]),log(TMLE$OR$CI[2])),
                       DTMLE.estimate=c(ATE,log.RR,log.OR),
                       DTMLE.se=c(ATE.se,logRR.se,logOR.se),
                       DTMLE.l95=c(ATE.CI[1],logRR.CI[1],logOR.CI[1]),
                       DTMLE.u95=c(ATE.CI[2],logRR.CI[2],logOR.CI[2])
  )

  results.tmle = results[,1:6]
  results.tmle$method = 'TMLE'
  results.dtmle = results[,c(1:2,7:10)]
  results.dtmle$method = 'DTMLE'
  names(results.tmle) = c('parameter','true','estimate','se','l95','u95','method')
  names(results.dtmle) = c('parameter','true','estimate','se','l95','u95','method')
  results = rbind(results.tmle,results.dtmle)

  return(results)
}
