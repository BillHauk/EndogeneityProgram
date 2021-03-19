clear
set more off

set matsize 200
*true value for alpha
global alpha=1/3
*true value for beta
global beta=1/3

global lambda=(1-$alpha-$beta)*(0.08)
global multiple=1-exp(-$lambda*5)

*set true value for gamma1 (rem first line and de-rem second line for structural parameters)
global gamma1=0.0579475
*global gamma1=$multiple*($alpha/(1-$alpha-$beta))
*set true value for gamma1 (rem first line and de-rem second line for structural parameters)
global gamma2=0.0407447
*global gamma2=$multiple*($beta/(1-$alpha-$beta))
*set true value for gamma3 (rem first line and de-rem second line for structural parameters)
global gamma3=-0.2141621
*global gamma3=-$multiple*(($alpha+$beta)/(1-$alpha-$beta))
*set true value for gamma4 (rem first line and de-rem second line for structural parameters)
global gamma4=0.7962023
*global gamma4=(1-$multiple)

*Set variance and correlation scaling for Fixed Effect term
global fevar=0.5
global fecorr=0.5

*number of countries in Monte Carlo sample draw
global N=2000

*Which residual correlations do you want to control?
* 1 = sk, sh and ndg
* 2 = sk and sh
* 3 = sk and ndg
* 4 = sh and ndg
* 5 = sk only
* 6 = sh only
* 7 = ndg only
global control=1

*correlation between residual term and sk
global skresidcorr=0.3
*correlation between residual term and sh
global shresidcorr=0.3
*correlation between residual term and ndg
global ndgresidcorr=-0.3

*Number of Monte Carlo draws
global drawnum=1000

use baselinedatabalnodelta

matrix Gamma=($gamma1, $gamma2, $gamma3, $gamma4, 1)'

sort countryindex period
iis countryindex
tis period
xtreg lyend lsk lsh lndg lyini i.period, fe
predict fes, u
reshape wide lyini lyend lsh lsk lndg fes, i(countryindex) j(period)
save baselinedatawide.dta, replace

sum
do residsgen$control

global t=1

*Each one of these loops generates coefficient values using fixed effects (i.e. the true model)
*shocked fixed effects, and shocked ols and saves them as data
*Note that the coefficients generated are of the form alpha/(1-alpha-beta), etc.

while $t==1{
	do MCgenerate
	* Define variables and tsset command for Arellano-Bond
	sort period country
	by period: egen meansk=mean(sk)
	by period: egen meansh=mean(sh)
	by period: egen meanndg=mean(ndg)
	by period: egen meanyini=mean(yini)
	by period: egen meanyend=mean(yend)
	gen skd=sk-meansk
	gen shd=sh-meansh
	gen ndgd=ndg-meanndg
	gen yinid=yini-meanyini
	gen yendd=yend-meanyend
	drop mean*
	sort country period
	tsset country period
	xtreg yend sk sh ndg yini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfe 	
	rename coefs2 shfe
	rename coefs3 ndgfe
	rename coefs4 yinife
	rename coefs5 consfe
	drop coefs*
	xtreg yend sk sh ndg yini i.period, be
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbe
	rename coefs2 shbe
	rename coefs3 ndgbe
	rename coefs4 yinibe
	rename coefs5 consbe
	drop coefs*
	xtreg yend sk sh ndg yini i.period, re
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skre
	rename coefs2 shre
	rename coefs3 ndgre
	rename coefs4 yinire
	rename coefs5 consre
	drop coefs*
	xtdpd yendd skd shd ndgd L.yendd, dgmmiv(skd shd ndgd L.yendd) noconstant
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skabond
	rename coefs2 shabond
	rename coefs3 ndgabond
	rename coefs4 yiniabond
	gen consabond=0
	xtdpd yendd skd shd ndgd L.yendd, dgmmiv(skd shd ndgd L.yendd) lgmmiv(skd shd ndgd L.yendd) noconstant
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbbond
	rename coefs2 shbbond
	rename coefs3 ndgbbond
	rename coefs4 yinibbond
	gen consbbond=0
	by country: replace sk=(sk+sk[_n+1]+sk[_n+2]+sk[_n+3]+sk[_n+4]+sk[_n+5]+sk[_n+6]+sk[_n+7]+sk[_n+8]+sk[_n+9])/10 if period==3
	by country: replace sh=(sh+sh[_n+1]+sh[_n+2]+sh[_n+3]+sh[_n+4]+sh[_n+5]+sh[_n+6]+sh[_n+7]+sh[_n+8]+sh[_n+9])/10 if period==3
	by country: replace ndg=(ndg+ndg[_n+1]+ndg[_n+2]+ndg[_n+3]+ndg[_n+4]+ndg[_n+5]+ndg[_n+6]+ndg[_n+7]+ndg[_n+8]+ndg[_n+9])/10 if period==3
	by country: replace yend=yend[_n+9] if period==3
	regress yend sk sh ndg yini if period==3
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skmrw
	rename coefs2 shmrw
	rename coefs3 ndgmrw
	rename coefs4 yinimrw
	rename coefs5 consmrw
	#delimit ;
	keep skfe shfe ndgfe yinife consfe 
	skre shre ndgre yinire consre
	skabond shabond ndgabond yiniabond consabond
	skbbond shbbond ndgbbond yinibbond consbbond
	skbe shbe ndgbe yinibe consbe
	skmrw shmrw ndgmrw yinimrw consmrw;
	#delimit cr
	keep if _n==1
	save montecarloresults, replace
	noisily disp $t
	global t=$t+1
}

quietly while $t<=$drawnum{
	do MCgenerate
	* Define variables and tsset command for Arellano-Bond
	sort period country
	by period: egen meansk=mean(sk)
	by period: egen meansh=mean(sh)
	by period: egen meanndg=mean(ndg)
	by period: egen meanyini=mean(yini)
	by period: egen meanyend=mean(yend)
	gen skd=sk-meansk
	gen shd=sh-meansh
	gen ndgd=ndg-meanndg
	gen yinid=yini-meanyini
	gen yendd=yend-meanyend
	drop mean*
	sort country period
	tsset country period
	xtreg yend sk sh ndg yini i.period, fe
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skfe 	
	rename coefs2 shfe
	rename coefs3 ndgfe
	rename coefs4 yinife
	rename coefs5 consfe
	drop coefs*
	xtreg yend sk sh ndg yini i.period, be
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbe
	rename coefs2 shbe
	rename coefs3 ndgbe
	rename coefs4 yinibe
	rename coefs5 consbe
	drop coefs*
	xtreg yend sk sh ndg yini i.period, re
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skre
	rename coefs2 shre
	rename coefs3 ndgre
	rename coefs4 yinire
	rename coefs5 consre
	drop coefs*
	xtdpd yendd skd shd ndgd L.yendd, dgmmiv(skd shd ndgd L.yendd) noconstant
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skabond
	rename coefs2 shabond
	rename coefs3 ndgabond
	rename coefs4 yiniabond
	gen consabond=0
	xtdpd yendd skd shd ndgd L.yendd, dgmmiv(skd shd ndgd L.yendd) lgmmiv(skd shd ndgd L.yendd) noconstant
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skbbond
	rename coefs2 shbbond
	rename coefs3 ndgbbond
	rename coefs4 yinibbond
	gen consbbond=0
	by country: replace sk=(sk+sk[_n+1]+sk[_n+2]+sk[_n+3]+sk[_n+4]+sk[_n+5]+sk[_n+6]+sk[_n+7]+sk[_n+8]+sk[_n+9])/10 if period==3
	by country: replace sh=(sh+sh[_n+1]+sh[_n+2]+sh[_n+3]+sh[_n+4]+sh[_n+5]+sh[_n+6]+sh[_n+7]+sh[_n+8]+sh[_n+9])/10 if period==3
	by country: replace ndg=(ndg+ndg[_n+1]+ndg[_n+2]+ndg[_n+3]+ndg[_n+4]+ndg[_n+5]+ndg[_n+6]+ndg[_n+7]+ndg[_n+8]+ndg[_n+9])/10 if period==3
	by country: replace yend=yend[_n+9] if period==3
	regress yend sk sh ndg yini if period==3
	matrix coefs=e(b)
	svmat coefs
	rename coefs1 skmrw
	rename coefs2 shmrw
	rename coefs3 ndgmrw
	rename coefs4 yinimrw
	rename coefs5 consmrw
	#delimit ;
	keep skfe shfe ndgfe yinife consfe 
	skre shre ndgre yinire consre
	skabond shabond ndgabond yiniabond consabond
	skbbond shbbond ndgbbond yinibbond consbbond
	skbe shbe ndgbe yinibe consbe
	skmrw shmrw ndgmrw yinimrw consmrw;
	#delimit cr
	keep if _n==1
	append using montecarloresults
	save montecarloresults, replace
	noisily disp $t
	global t=$t+1
}


gen lambdamrw=-log(yinimrw)/50
gen yinimrwadjust=exp(-lambdamrw*5)
replace yinimrwadjust=0 if yinimrwadjust==.
gen skmrwadjust=((1-yinimrwadjust)/(1-yinimrw))*skmrw
gen shmrwadjust=((1-yinimrwadjust)/(1-yinimrw))*shmrw
gen ndgmrwadjust=((1-yinimrwadjust)/(1-yinimrw))*ndgmrw
gen consmrwadjust=0


#delimit ;
order skfe shfe ndgfe yinife consfe skbe shbe ndgbe yinibe consbe
skre shre ndgre yinire consre skabond shabond ndgabond yiniabond consabond
skbbond shbbond ndgbbond yinibbond consbbond skmrw shmrw ndgmrw yinimrw consmrw
skmrwadjust shmrwadjust ndgmrwadjust yinimrwadjust consmrwadjust;
#delimit cr


