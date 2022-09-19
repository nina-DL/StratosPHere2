# Load Data: Stratosph part1 --------------------------------------------------------------------

load("StratosPHereI.RData")

# Primary endpoint: Panel of biomarkers
All_biomarkers_qPCR_dCT$All = apply(data.frame(All_biomarkers_qPCR_dCT), 1, mean)



# Allocation rule --------------------------------------------------------------------

# This function assumes a binary reward variable Y defined as > / <= a 30% increase

AP_binaryY <- function(Nsim = 5000, Alloc_rule, K = 3, Tmax = 20, b = c(6, 6, 8), stage, success = 0, failure = 0, prior_params = c(1,1)){
  # Nsim = simulations for the RAR algorithm
  # Alloc_rule = allocation rule to be used among:
  # 'TS': Thompson Sampling
  # 'Trippa': proceduree of Trippa et al
  # 'FR': fixed and equal randomization
  # K = number of arms
  # Tmax = max horizon (or sample size)
  # b = block size
  # prior_params = vector with prior success and failure 
  
  prior_s = array(prior_params[1],c(1,K))
  prior_f = array(prior_params[2],c(1,K))
  success = success + prior_s 
  failure = failure + prior_f
  n_patients <- (success+failure)-(prior_s+prior_f)  # total successes and failures (excluding priors)
  #stage = which.min(abs(cumsum(b)-sum(n_patients)))
  
  # posterior sampling: sample the success probability 
  postsample <- array(0,c(Nsim,K)) 
  for (k in 1:K) {
    postsample[,k] <- rbeta(Nsim,success[k],failure[k])
  }
  
  ap = array(0,c(1,K)) 
  
  # Thompson sampling
  # Alloc_rule = 'TS'
  if (Alloc_rule == 'TS') {
    best_arm = apply(postsample,1,which.max)
    c = (stage*b[stage])/(2*Tmax)
    for (k in 1:K){
      ap[,k] = mean(best_arm == k)^c                    
    }
    
    constant <- sum(ap)
    ap <- ap/constant
    allocation_prob <- t(ap)
  }
  
  # Trippa et al. Procedure
  # Alloc_rule = 'Trippa'
  else if (Alloc_rule == 'Trippa') {
    gamma = 13*((stage*b[stage])/Tmax)^2.75 # this governs the balance between the experimental arms
    eta = ((stage*b[stage])/Tmax)*0.25 # this governs the control arm
    
    # AP of treatments
    for (k in 2:K) {
      ap[k] = mean(postsample[,k]>postsample[,1])^gamma
    }
    constant_tr <- sum(ap[,2:3])
    
    for (k in 2:K){
      ap[k] <- ap[k]/constant_tr
    }
    
    # AP of control
    ap[1] <- 1/K*exp(max(n_patients[2:K])-n_patients[1])^eta
    
    constant <- sum(ap)
    ap <- ap/constant
    
    allocation_prob <- t(ap)
  } else {
    allocation_prob <- rep(1/K, K)
  }
  
  return(allocation_prob)
}

AP_binaryY(Alloc_rule = 'Trippa', stage = 1, success = rep(0,3), failure = rep(0,3))


# Inference function --------------------------------------------------------------------
library(DescTools)

do_inference = function(replicas = 5000, Nsim_alloc = 500, K = 3, b = c(6, 6, 8), burn_in_b1 = T, Nmax = 20, delta = 0.3, alpha = 0.05, my_data, H1 = 'Gain_score'){
  # details of the simulation studies
  # ss=20; #max sample size per subgroup
  # delta=0.3; #ESS: 30% increase
  # replica=5000;
  # b= 4; #block size
  # K=3; # arms = number of treatment arms in the trial (including the control)
  # Nsim=500; #simulations for RAR algorithm
  # my_data = real data of a biomarker to evaluate
  # H1 = 'Gain_score' or H1 = 'Follow-up'
  
  ap_arm0 = ap_arm1 = ap_arm2 = array(0,c(1,replicas))  
  pv_H1 = pv_H0arm1 = pv_H1arm2 = array(0, c(1, replicas))
  
  success = failure = array(0,c(replicas,K)) 
  outcome = array(0,c(replicas,Nmax))
  arm = array(0,c(replicas,Nmax))
  
  ### generate sample and trial replicas with our proposed RAR algorithm
  for (i in 1:replicas){
    
    if(burn_in_b1==T){
      # first interim = Fixed Restricted Randomization
      arm[i,1:b[1]] = sample(rep(c(1:3),b[1]/K), b[1], replace = F)
    } else {
      #get probabilities for block 1, based on equal uninformative priors. 
      alloc_prob = AP_binaryY(Nsim = Nsim_alloc, Alloc_rule = 'Trippa', K = K, b = b, Tmax = Nmax, stage = 1, success = success[i,], failure = failure[i,])
      arm[i,1:b[1]] = sample(1:3, b[1], replace = T, prob = alloc_prob)
    }
    
    ss = b[1]
    y_baseline = sample(my_data, Nmax, replace = TRUE) # we assume follow-up of arm0 is equal to baseline
    y_arm0 = y_arm1 = sample(my_data, Nmax, replace = TRUE) # we assume follow-up of arm0 and arm1 are equal to baseline
    y_arm2=sample((1+delta)*my_data, Nmax, replace = TRUE) # we assume follow-up of arm2 has a 30% increase
    y_obs = c()
    # assigning binary outcomes for the first block
    for(pts in 1:b[1]){
      outcome[i,pts] = 1*(arm[i,pts]==1)*((y_arm0[pts]-y_baseline[pts])/y_baseline[pts] >= delta) +
        1*(arm[i,pts]==2)*((y_arm1[pts]-y_baseline[pts])/y_baseline[pts] >= delta) +
        1*(arm[i,pts]==3)*((y_arm2[pts]-y_baseline[pts])/y_baseline[pts] >= delta)
      
      y_obs[pts] = 1*(arm[i,pts]==1)*(y_arm0[pts]) + 1*(arm[i,pts]==2)*(y_arm1[pts]) + 1*(arm[i,pts]==3)*(y_arm2[pts])
      
      for(k in 1:K){
        success[i,k] = success[i,k]+1*(arm[i,pts]==k)*outcome[i,pts]
        failure[i,k] = failure[i,k]+1*(arm[i,pts]==k)*(1-outcome[i,pts])
      }
    }
    
    for (j in 2:length(b)){
      ss = ss + b[j]
      alloc_prob = AP_binaryY(Nsim = Nsim_alloc, Alloc_rule = 'Trippa', K, b = b, Tmax = Nmax, stage = j, success = success[i,], failure = failure[i,])
      arm[i,(ss-b[j]+1):ss] = sample(1:3, b[j], replace = T, prob = alloc_prob)
      
      # assigning binary outcomes for the first block
      for(pts in (ss-b[j]+1):ss){
        outcome[i,pts] = 1*(arm[i,pts]==1)*((y_arm0[pts]-y_baseline[pts])/y_baseline[pts] >= delta) +
          1*(arm[i,pts]==2)*((y_arm1[pts]-y_baseline[pts])/y_baseline[pts] >= delta) +
          1*(arm[i,pts]==3)*((y_arm2[pts]-y_baseline[pts])/y_baseline[pts] >= delta)
        
        y_obs[pts] = 1*(arm[i,pts]==1)*(y_arm0[pts]) + 1*(arm[i,pts]==2)*(y_arm1[pts]) + 1*(arm[i,pts]==3)*(y_arm2[pts])
        
        for(k in 1:K){
          success[i,k] = success[i,k]+1*(arm[i,pts]==k)*outcome[i,pts]
          failure[i,k] = failure[i,k]+1*(arm[i,pts]==k)*(1-outcome[i,pts])
        }
      }
    }
    #success[i,];    failure[i,]; arm[i,]
    
    ap_arm0[i]=sum(arm[i,]==1)/Nmax
    ap_arm1[i]=sum(arm[i,]==2)/Nmax
    ap_arm2[i]=sum(arm[i,]==3)/Nmax
    
    mydf = data.frame(arm = arm[i,], baseline = y_baseline, followup = y_obs, diff = y_obs - y_baseline,
                      increase = (y_obs - y_baseline)/y_baseline)
    
    if (H1 == 'Gain_score'){
      pv_H0arm1[i]=try(wilcox.test(mydf$diff[mydf$arm==2], mydf$diff[mydf$arm==1], alternative = "greater")$p.value, silent = T)
      pv_H1arm2[i]=try(wilcox.test(mydf$diff[mydf$arm==3], mydf$diff[mydf$arm==1], alternative = "greater")$p.value, silent = T)
      
    } else {
      pv_H0arm1[i]=try(wilcox.test(mydf$followup[mydf$arm==2], mydf$followup[mydf$arm==1], alternative = "greater")$p.value, silent = T)
      pv_H1arm2[i]=try(wilcox.test(mydf$followup[mydf$arm==3], mydf$followup[mydf$arm==1], alternative = "greater")$p.value, silent = T)
    }
    
  }
  
  return(list(alloc_arm0 = mean(ap_arm0), alloc_arm1 = mean(ap_arm1), alloc_arm2 = mean(ap_arm2),
              power_arm2H1_wx = mean(pv_H1arm2<alpha),
              typeI_H0arm1 = mean(pv_H0arm1<alpha)))
}


qPCR_res_dCT = matrix(NA, nrow = 5, ncol = length(All_biomarkers_qPCR_dCT))
for (i in 1:length(All_biomarkers_qPCR_dCT)){
  qPCR_res_dCT[,i] = unlist(do_inference(Nmax = 20, my_data = All_biomarkers_qPCR_dCT[[i]], b = c(6, 6, 8), alpha = 0.067, H1 = 'Gain_score'))
}

colnames(qPCR_res_dCT) =c("ID3", "SMAD1", "SMAD5", "NOTCH1", "NOTCH2", "ID2", "ARL4C", "PTGS2", "All")
rownames(qPCR_res_dCT) =c("Alloc_arm0", "Alloc_arm1", "Alloc_arm2", "Power_arm2_wx",  "TypeIerror_wx")
round(qPCR_res_dCT, 3)
