*This routine takes the modified raw data generated in datasetgenerate and
*calculates the mean, standard deviations, and correlations of those variables
*Monte Carlo draws are then generated using these moments

clear

*These loops get the means and standard deviations from the dataset generated
*in datasetgenerate

use baselinedatawide.dta
drop fes3 fes4 fes5 fes6 fes7 fes8 fes9 fes10 fes11
rename fes12 fes
gen lyini13=lyend12

global i 3
while $i==3{
	sum lyini$i	
	matrix N=(r(mean))
	matrix T=(r(sd))
	sum lsk$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	sum lsh$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	sum lndg$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	matrix M=N
	matrix S=T
	global i=$i+1
}

while $i<=12{
	sum lyini$i	
	matrix N=(r(mean))
	matrix T=(r(sd))
	sum lsk$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	sum lsh$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	sum lndg$i
	matrix N=(N, r(mean))
	matrix T=(T, r(sd))
	matrix M=(M, N)
	matrix S=(S, T)
	global i=$i+1
}

while $i==13{
	sum lyini$i	
	matrix N=(r(mean))
	matrix T=(r(sd))
	sum fes
	matrix N=(N, r(mean))
	matrix T=(T, ($fevar^0.5)*r(sd))
	matrix M=(M, N)
	matrix S=(S, T)
	global i=$i+1
}


* And then puts them into matricies to be used as parameters for a drawnorm command

#delimit ;
matrix accum X=lyini3 lsk3 lsh3 lndg3 lyini4 lsk4 lsh4 lndg4 lyini5 lsk5
lsh5 lndg5 lyini6 lsk6 lsh6 lndg6 lyini7 lsk7 lsh7 lndg7 lyini8 lsk8
lsh8 lndg8 lyini9 lsk9 lsh9 lndg9 lyini10 lsk10 lsh10 lndg10
lyini11 lsk11 lsh11 lndg11 lyini12 lsk12 lsh12 lndg12
lyini13 fes, nocons dev;
#delimit cr

matrix R=corr(X)
matrix Rmod=R[1..41,42]
matrix Rmod=$fecorr*Rmod
matrix R[1,42]=Rmod
matrix R[42,1]=Rmod'

#delimit ;
drawnorm yini3 sk3 sh3 ndg3 yini4 sk4 sh4 ndg4 yini5 sk5 sh5 ndg5 yini6 sk6 sh6
ndg6 yini7 sk7 sh7 ndg7 yini8 sk8 sh8 ndg8 yini9 sk9 sh9 ndg9 yini10 sk10 sh10
ndg10 yini11 sk11 sh11 ndg11 yini12 sk12 sh12 ndg12 yini13 fe, n($N) means(M) corr(R) sds(S) clear;
#delimit cr


generate country=_n
order country
sort country
keep country yini* sk* sh* ndg* fe

save generateddata, replace

* This section generates error terms that have the same means and variances as the data
* but make the error terms correlated with other RHS variables

matrix A=($meaneps3)
matrix B=($sigmaeps3)
drawnorm eps3, mean(A) sd(B) n($N)
gen nu3=$phi13*sk3+$phi23*sh3+$phi33*ndg3+$phi43*yini3+$phi53*fe+eps3
gen yend3=$gamma1*sk3+$gamma2*sh3+$gamma3*ndg3+$gamma4*yini3+fe+nu3
replace yini4=yend3

matrix A=($meaneps4)
matrix B=($sigmaeps4)
drawnorm eps4, mean(A) sd(B) n($N)
gen nu4=$phi14*sk4+$phi24*sh4+$phi34*ndg4+$phi44*yini4+$phi54*fe+eps4
gen yend4=$gamma1*sk4+$gamma2*sh4+$gamma3*ndg4+$gamma4*yini4+fe+nu4
replace yini5=yend4

matrix A=($meaneps5)
matrix B=($sigmaeps5)
drawnorm eps5, mean(A) sd(B) n($N)
gen nu5=$phi15*sk5+$phi25*sh5+$phi35*ndg5+$phi45*yini5+$phi55*fe+eps5
gen yend5=$gamma1*sk5+$gamma2*sh5+$gamma3*ndg5+$gamma4*yini5+fe+nu5
replace yini6=yend5

matrix A=($meaneps6)
matrix B=($sigmaeps6)
drawnorm eps6, mean(A) sd(B) n($N)
gen nu6=$phi16*sk6+$phi26*sh6+$phi36*ndg6+$phi46*yini6+$phi56*fe+eps6
gen yend6=$gamma1*sk6+$gamma2*sh6+$gamma3*ndg6+$gamma4*yini6+fe+nu6
replace yini7=yend6

matrix A=($meaneps7)
matrix B=($sigmaeps7)
drawnorm eps7, mean(A) sd(B) n($N)
gen nu7=$phi17*sk7+$phi27*sh7+$phi37*ndg7+$phi47*yini7+$phi57*fe+eps7
gen yend7=$gamma1*sk7+$gamma2*sh7+$gamma3*ndg7+$gamma4*yini7+fe+nu7
replace yini8=yend7

matrix A=($meaneps8)
matrix B=($sigmaeps8)
drawnorm eps8, mean(A) sd(B) n($N)
gen nu8=$phi18*sk8+$phi28*sh8+$phi38*ndg8+$phi48*yini8+$phi58*fe+eps8
gen yend8=$gamma1*sk8+$gamma2*sh8+$gamma3*ndg8+$gamma4*yini8+fe+nu8
replace yini9=yend8

matrix A=($meaneps9)
matrix B=($sigmaeps9)
drawnorm eps9, mean(A) sd(B) n($N)
gen nu9=$phi19*sk9+$phi29*sh9+$phi39*ndg9+$phi49*yini9+$phi59*fe+eps9
gen yend9=$gamma1*sk9+$gamma2*sh9+$gamma3*ndg9+$gamma4*yini9+fe+nu9
replace yini10=yend9

matrix A=($meaneps10)
matrix B=($sigmaeps10)
drawnorm eps10, mean(A) sd(B) n($N)
gen nu10=$phi110*sk10+$phi210*sh10+$phi310*ndg10+$phi410*yini10+$phi510*fe+eps10
gen yend10=$gamma1*sk10+$gamma2*sh10+$gamma3*ndg10+$gamma4*yini10+fe+nu10
replace yini11=yend10

matrix A=($meaneps11)
matrix B=($sigmaeps11)
drawnorm eps11, mean(A) sd(B) n($N)
gen nu11=$phi111*sk11+$phi211*sh11+$phi311*ndg11+$phi411*yini11+$phi511*fe+eps11
gen yend11=$gamma1*sk11+$gamma2*sh11+$gamma3*ndg11+$gamma4*yini11+fe+nu11
replace yini12=yend11

matrix A=($meaneps12)
matrix B=($sigmaeps12)
drawnorm eps12, mean(A) sd(B) n($N)
gen nu12=$phi112*sk12+$phi212*sh12+$phi312*ndg12+$phi412*yini12+$phi512*fe+eps12
gen yend12=$gamma1*sk12+$gamma2*sh12+$gamma3*ndg12+$gamma4*yini12+fe+nu12
replace yini13=yend12

gen yend2=yini3
keep country yend* sk* sh* ndg* yini*
reshape long yend sk sh ndg yini, i(country) j(period)

iis country
save generateddata, replace

