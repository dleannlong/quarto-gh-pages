/*The following SAS NLMIXED code was used for the SafeTalk motivating example.*/

proc nlmixed data=work.for_analysis seed=31415 maxiter=500 qpoints=50 cov hess;

    /* null initial parameter values */
     parms g0 0 g1 0 g2 0 g3 0 g4 0
           a0 0 a1 0 a2 0 a3 0 a4 0;

    /* linear predictor for the zero-inflation probability  */
    /* logit(psi)=Z\gamma                                   */

     logit_psi = g0 + g1*arm + g2*site2 + g3*site3 + g4*baseline_uavi;

    /* Useful functions of psi */
     psi1  = exp(logit_psi)/(1+exp(logit_psi));  /*psi = exp(Z\gamma)/(1+exp(Z\gamma)) */
     psi2  = 1/(1+exp(logit_psi));               /*1-psi = (1+exp(Z\gamma))ˆ-1         */

    /* Overall mean \nu */
    /* log(nu) = X\alpha */

     log_nu   = a0 + a1*arm + a2*site2 + a3*site3 + a4*baseline_uavi;

     delta = log(psi2**(-1)) + log_nu;

    /* Build the mZIPlog likelihood */
 
     if outcome=0 then
           ll = log(psi1 + psi2*(exp(-exp(delta))));
      else ll = log(psi2) - exp(delta)  + outcome*(delta) - lgamma(outcome + 1);

     model outcome ˜  general(ll);
run;
