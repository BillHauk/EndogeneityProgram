*This is residsgen5

matrix accum XX=(lsk3 lsh3 lndg3 lyini3 fes3), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk3
global skmean=r(mean)
global sksd=r(sd)
sum lsh3
global shmean=r(mean)
global shsd=r(sd)
sum lndg3
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini3
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend3
global yendmean=r(mean)
global yendsd=r(sd)
sum fes3
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu3 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu3 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis3=invsym(XX)*cov

global phi13=phis3[1,1]
global phi23=phis3[2,1]
global phi33=phis3[3,1]
global phi43=phis3[4,1]
global phi53=phis3[5,1]


matrix FXXF=phis3'*XX*phis3
global sigma2exp2=FXXF[1,1]

global meaneps3=$yendmean-($gamma1+$phi13)*$skmean-($gamma2+$phi23)*$shmean-($gamma3+$phi33)*$ndgmean-($gamma4+$phi43)*$yinimean
global sigmaeps3=($sigmanu^2-$sigma2exp2)^.5


if $meaneps3==. {
	disp "meaneps3 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps3==. {
	disp "sigmaeps3 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk4 lsh4 lndg4 lyini4 fes4), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk4
global skmean=r(mean)
global sksd=r(sd)
sum lsh4
global shmean=r(mean)
global shsd=r(sd)
sum lndg4
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini4
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend4
global yendmean=r(mean)
global yendsd=r(sd)
sum fes4
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu4 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu4 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis4=invsym(XX)*cov

global phi14=phis4[1,1]
global phi24=phis4[2,1]
global phi34=phis4[3,1]
global phi44=phis4[4,1]
global phi54=phis4[5,1]

matrix FXXF=phis4'*XX*phis4
global sigma2exp2=FXXF[1,1]

global meaneps4=$yendmean-($gamma1+$phi14)*$skmean-($gamma2+$phi24)*$shmean-($gamma3+$phi34)*$ndgmean-($gamma4+$phi44)*$yinimean
global sigmaeps4=($sigmanu^2-$sigma2exp2)^.5


if $meaneps4==. {
	disp "meaneps4 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps4==. {
	disp "sigmaeps4 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk5 lsh5 lndg5 lyini5 fes5), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk5
global skmean=r(mean)
global sksd=r(sd)
sum lsh5
global shmean=r(mean)
global shsd=r(sd)
sum lndg5
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini5
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend5
global yendmean=r(mean)
global yendsd=r(sd)
sum fes5
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu5 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu5 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis5=invsym(XX)*cov

global phi15=phis5[1,1]
global phi25=phis5[2,1]
global phi35=phis5[3,1]
global phi45=phis5[4,1]
global phi55=phis5[5,1]

matrix FXXF=phis5'*XX*phis5
global sigma2exp2=FXXF[1,1]

global meaneps5=$yendmean-($gamma1+$phi15)*$skmean-($gamma2+$phi25)*$shmean-($gamma3+$phi35)*$ndgmean-($gamma4+$phi45)*$yinimean
global sigmaeps5=($sigmanu^2-$sigma2exp2)^.5


if $meaneps5==. {
	disp "meaneps5 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps5==. {
	disp "sigmaeps5 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk6 lsh6 lndg6 lyini6 fes6), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk6
global skmean=r(mean)
global sksd=r(sd)
sum lsh6
global shmean=r(mean)
global shsd=r(sd)
sum lndg6
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini6
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend6
global yendmean=r(mean)
global yendsd=r(sd)
sum fes6
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu6 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu6 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis6=invsym(XX)*cov

global phi16=phis6[1,1]
global phi26=phis6[2,1]
global phi36=phis6[3,1]
global phi46=phis6[4,1]
global phi56=phis6[5,1]

matrix FXXF=phis6'*XX*phis6
global sigma2exp2=FXXF[1,1]

global meaneps6=$yendmean-($gamma1+$phi16)*$skmean-($gamma2+$phi26)*$shmean-($gamma3+$phi36)*$ndgmean-($gamma4+$phi46)*$yinimean
global sigmaeps6=($sigmanu^2-$sigma2exp2)^.5


if $meaneps6==. {
	disp "meaneps6 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps6==. {
	disp "sigmaeps6 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk7 lsh7 lndg7 lyini7 fes7), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk7
global skmean=r(mean)
global sksd=r(sd)
sum lsh7
global shmean=r(mean)
global shsd=r(sd)
sum lndg7
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini7
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend7
global yendmean=r(mean)
global yendsd=r(sd)
sum fes7
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu7 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu7 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis7=invsym(XX)*cov

global phi17=phis7[1,1]
global phi27=phis7[2,1]
global phi37=phis7[3,1]
global phi47=phis7[4,1]
global phi57=phis7[5,1]

matrix FXXF=phis7'*XX*phis7
global sigma2exp2=FXXF[1,1]

global meaneps7=$yendmean-($gamma1+$phi17)*$skmean-($gamma2+$phi27)*$shmean-($gamma3+$phi37)*$ndgmean-($gamma4+$phi47)*$yinimean
global sigmaeps7=($sigmanu^2-$sigma2exp2)^.5

if $meaneps7==. {
	disp "meaneps7 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps7==. {
	disp "sigmaeps7 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk8 lsh8 lndg8 lyini8 fes8), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk8
global skmean=r(mean)
global sksd=r(sd)
sum lsh8
global shmean=r(mean)
global shsd=r(sd)
sum lndg8
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini8
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend8
global yendmean=r(mean)
global yendsd=r(sd)
sum fes8
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu8 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu8 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis8=invsym(XX)*cov

global phi18=phis8[1,1]
global phi28=phis8[2,1]
global phi38=phis8[3,1]
global phi48=phis8[4,1]
global phi58=phis8[5,1]

matrix FXXF=phis8'*XX*phis8
global sigma2exp2=FXXF[1,1]

global meaneps8=$yendmean-($gamma1+$phi18)*$skmean-($gamma2+$phi28)*$shmean-($gamma3+$phi38)*$ndgmean-($gamma4+$phi48)*$yinimean
global sigmaeps8=($sigmanu^2-$sigma2exp2)^.5

if $meaneps8==. {
	disp "meaneps8 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps8==. {
	disp "sigmaeps8 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk9 lsh9 lndg9 lyini9 fes9), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk9
global skmean=r(mean)
global sksd=r(sd)
sum lsh9
global shmean=r(mean)
global shsd=r(sd)
sum lndg9
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini9
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend9
global yendmean=r(mean)
global yendsd=r(sd)
sum fes9
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu9 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu9 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis9=invsym(XX)*cov

global phi19=phis9[1,1]
global phi29=phis9[2,1]
global phi39=phis9[3,1]
global phi49=phis9[4,1]
global phi59=phis9[5,1]

matrix FXXF=phis9'*XX*phis9
global sigma2exp2=FXXF[1,1]

global meaneps9=$yendmean-($gamma1+$phi19)*$skmean-($gamma2+$phi29)*$shmean-($gamma3+$phi39)*$ndgmean-($gamma4+$phi49)*$yinimean
global sigmaeps9=($sigmanu^2-$sigma2exp2)^.5

if $meaneps9==. {
	disp "meaneps9 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps9==. {
	disp "sigmaeps9 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk10 lsh10 lndg10 lyini10 fes10), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk10
global skmean=r(mean)
global sksd=r(sd)
sum lsh10
global shmean=r(mean)
global shsd=r(sd)
sum lndg10
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini10
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend10
global yendmean=r(mean)
global yendsd=r(sd)
sum fes10
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu10 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu10 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis10=invsym(XX)*cov

global phi110=phis10[1,1]
global phi210=phis10[2,1]
global phi310=phis10[3,1]
global phi410=phis10[4,1]
global phi510=phis10[5,1]

matrix FXXF=phis10'*XX*phis10
global sigma2exp2=FXXF[1,1]

global meaneps10=$yendmean-($gamma1+$phi110)*$skmean-($gamma2+$phi210)*$shmean-($gamma3+$phi310)*$ndgmean-($gamma4+$phi410)*$yinimean
global sigmaeps10=($sigmanu^2-$sigma2exp2)^.5

if $meaneps10==. {
	disp "meaneps10 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps10==. {
	disp "sigmaeps10 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

matrix accum XX=(lsk11 lsh11 lndg11 lyini11 fes11), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk11
global skmean=r(mean)
global sksd=r(sd)
sum lsh11
global shmean=r(mean)
global shsd=r(sd)
sum lndg11
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini11
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend11
global yendmean=r(mean)
global yendsd=r(sd)
sum fes11
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu11 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu11 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis11=invsym(XX)*cov

global phi111=phis11[1,1]
global phi211=phis11[2,1]
global phi311=phis11[3,1]
global phi411=phis11[4,1]
global phi511=phis11[5,1]

matrix FXXF=phis11'*XX*phis11
global sigma2exp2=FXXF[1,1]

global meaneps11=$yendmean-($gamma1+$phi111)*$skmean-($gamma2+$phi211)*$shmean-($gamma3+$phi311)*$ndgmean-($gamma4+$phi411)*$yinimean
global sigmaeps11=($sigmanu^2-$sigma2exp2)^.5

if $meaneps11==. {
	disp "meaneps11 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps11==. {
	disp "sigmaeps11 does not exist -- this is probably because the residual correlations are impossible"
	stop
}
matrix accum XX=(lsk12 lsh12 lndg12 lyini12 fes12), nocons dev
matrix XX=XX/(r(N)-1)
matrix A=($fevar^0.5)*$fecorr*XX[1..4,5]
matrix XX[1,5]=A
matrix XX[5,1]=A'
matrix B=$fevar*XX[5,5]
matrix XX[5,5]=B
matrix R=corr(XX)
matrix YXXY=Gamma'*XX*Gamma

sum lsk12
global skmean=r(mean)
global sksd=r(sd)
sum lsh12
global shmean=r(mean)
global shsd=r(sd)
sum lndg12
global ndgmean=r(mean)
global ndgsd=r(sd)
sum lyini12
global yinimean=r(mean)
global yinisd=r(sd)
sum lyend12
global yendmean=r(mean)
global yendsd=r(sd)
sum fes12
global fesmean=r(mean)
global fessd=r(sd)*($fevar^0.5)

global sigma2exp=YXXY[1,1]
global shresidcorr=(R[2,1]*$skresidcorr)
global ndgresidcorr=(R[3,1]*$skresidcorr)
global yiniresidcorr=(R[4,1]*$skresidcorr)
global feresidcorr=(R[5,1]*$skresidcorr)

global A=1
global B=2*($gamma1*$skresidcorr*$sksd+$gamma2*$shresidcorr*$shsd+$gamma3*$ndgresidcorr*$ndgsd+$gamma4*$yiniresidcorr*$yinisd+$feresidcorr*$fessd)
global C=$sigma2exp-(($yendsd)^2)
global sigmanu=(-$B+($B^2-4*$A*$C)^0.5)/(2*$A)

if $sigmanu==. {
	disp "sigmanu12 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmanu<=0 {
	disp "sigmanu12 is less than zero -- this is probably because the residual correlations are impossible"
	stop
}

matrix cov=($skresidcorr*$sigmanu*$sksd, $shresidcorr*$sigmanu*$shsd, $ndgresidcorr*$sigmanu*$ndgsd, $yiniresidcorr*$sigmanu*$yinisd, $feresidcorr*$sigmanu*$fessd)'
matrix phis12=invsym(XX)*cov

global phi112=phis12[1,1]
global phi212=phis12[2,1]
global phi312=phis12[3,1]
global phi412=phis12[4,1]
global phi512=phis12[5,1]

matrix FXXF=phis12'*XX*phis12
global sigma2exp2=FXXF[1,1]

global meaneps12=$yendmean-($gamma1+$phi112)*$skmean-($gamma2+$phi212)*$shmean-($gamma3+$phi312)*$ndgmean-($gamma4+$phi412)*$yinimean
global sigmaeps12=($sigmanu^2-$sigma2exp2)^.5

if $meaneps12==. {
	disp "meaneps12 does not exist -- this is probably because the residual correlations are impossible"
	stop
}

if $sigmaeps12==. {
	disp "sigmaeps12 does not exist -- this is probably because the residual correlations are impossible"
	stop
}
