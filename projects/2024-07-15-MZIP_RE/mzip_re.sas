/*The following SAS NLMIXED code was used for the SafeTalk motivating example.*/

proc nlmixed data=safetalk seed=31415
	parms b0 0 b1 0 b2 0 b3 0 b4 0 b5 0 b6 0 b7 0 b8 0
	a0 0 a1 0 a2 0 a3 0 a4 0 a5 0 a6 0 a7 0 a8 0
	sigma1 1 sigma12 0 sigma2 1;

	/* linear predictor for the zeroinflation probability */
	logit__psi = a0 + a1*site2 + a2*site3 + a3*v2 + a4*v2*st + a5*v3 + a6*v3*st + a7*v4 + a8*v4*st + c1;

	*logit(\psi)=Z\gamma + c;

	/* useful functions of \psi */
	psi1 = exp(logit__psi)/(1+exp(logit__psi));

	*\psi = exp(Z\gamma+c)/(1+exp(Z\gamma+c));

	psi2 = 1/(1+exp(logit__psi));

	*1{\psi = (1+exp(Z\gamma+c));

	/* Overall mean \nu */
	log__nu = b0 + b1*site2 + b2*site3 + b3*v2 + b4*v2*st + b5*v3 + b6*v3*st + b7*v4 + b8*v4*st + d1;
	delta = log(psi2**(1)) + log__nu;

	/* Build the MZIP + RE loglikelihood */
	if outcome=0 then
	ll = log(psi1 + psi2*(exp(exp(delta))));
	else ll = log(psi2) exp(delta) + outcome*(delta) lgamma(outcome + 1);

	model outcome general(ll);

	random c1 d1 normal([0,0],[sigma1,sigma12,sigma2]) SUBJECT=urn;

	contrast "TX" b4, b6, b8;

run;
