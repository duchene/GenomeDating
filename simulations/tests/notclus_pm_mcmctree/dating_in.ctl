          seed = -1357
seqfile = notclus_pm.phy
treefile = notclus_pm.tre
outfile = notclus_pm_outfile
usedata = 2*
ndata = 4 *
         clock = 2   * 1: global clock; 2: independent rates; 3: correlated rates
       RootAge = '<1.1'
         model = 0    * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
         alpha = 0    * alpha for gamma rates at sites
         ncatG = 1    * No. categories in discrete gamma
     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?
       BDparas = 1 1 0  * lambda, mu, rho, psi for birth-death-sampling model
   kappa_gamma = 6 2   * gamma prior for kappa
   alpha_gamma = 1 1   * gamma prior for alpha
   rgene_gamma = 2 2   * gamma prior for rate for genes
  sigma2_gamma = 1 10    * gamma prior for sigma^2  (for clock=2)
      finetune = 1: .1 .1 .1 .1 .1 .1  * auto (0 or 1) : times, rates, mixing, paras, RateParas, FossilErr
         print = 1
        burnin = 5000
      sampfreq = 5
       nsample = 20000
