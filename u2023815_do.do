
************************
* Univerisity of Warwick
* Department of Economics
* EC959 Msc Dissertation
"Does Bad Air Quality Contribute to Obesity? Evidence from China’s Central Heating System"
* ID u2023815 <Ma, Yuxuan(Sylvia)>
* Also see
"u2023815.txt" log file 
* Augest 2021
*************************


******************Data Processing*************************
{

*merge county datatset to individual dtataset 
use "D:\Home\Dissertation\obesity"
merge m:m countyid year using "D:\Home\Dissertation\pollution_modify"
drop _merge
save "C:\Users\maryh\Desktop\data\final proposal.dta"
*********************************************************
* + controls
destring countyid, replace
gen na = real(a)
encode b, gen(nb)
use "C:\Users\maryh\Desktop\data\final proposal.dta" 
merge m:m countyid using "C:\Users\maryh\Desktop\data\county cntrls_modify.dta" 
drop _merge
save, replace

}



*********************************************************
         ***** <obesity~air pollution> *****
*********************************************************

log using "C:\Users\maryh\Desktop\data\log\u2023815.smcl"
use "C:\Users\maryh\Desktop\data\final disser.dta" 
cd c:\Users\maryh\Desktop\results

*pre works
{
*set y & x
gen bmi = weight/(height*0.01)^2
rename mean pm

*log controls
gen lgdp=log(gdp)
gen lrev=log(generalrevenue)
gen ldep=log(savingsdeposits)
gen lexp=log(generalexpenditures)
gen lfa=log(invFA)
gen byte gender_dummy = 1 if gender=="Female"
replace gender_dummy=0 if gender=="Male"
//encode gender, gen(genderdummy)

*sum stat & select covariates
asdoc sum, save(disser.doc)
asdoc sum, detail
//drop gdp lgdp(cuz lack obs 3660/6000)

reg pm edu gender_dummy yob latitude longitude area lexp lrev employees population ldep lfa
estat vif
//drop generalrevenue lrev(cuz vif=11>10 multicollinearity)
pwcorr $Z $X
pwcorr $Z $X, sig
//correlation(|r|≥0.8,high;0.5≤|r|<0.8,medium;0.3≤|r|<0.5,low; |r|<0.3,irrelated)

*generate covariates
global x "edu gender_dummy yob"
global z "latitude longitude area lexp employees population ldep lfa"
global zz "longitude area lexp employees population ldep lfa"
//if global drop
//macro drop x
//macro drop_all
//raw controls: edu genderdummy yob longitude -gdp- generalexpenditures -generalrevenue- employees population area savingsdeposits invFA

}


**********************************************************
                    ** 1 OLS **
{
*validity check
**multicollinearity
reg bmi pm $x $z
asdoc hettest
asdoc estat vif
**correlation
asdoc pwcorr $Z $X
asdoc pwcorr $Z $X, sig

*separate
asdoc reg bmi pm
asdoc reg bmi pm, vce(cluster countyid)
asdoc reg bmi pm $x $z
asdoc reg bmi pm $x $z, vce(cluster countyid)
*combine
eststo x1: reg bmi pm
eststo x2: reg bmi pm, vce(cluster countyid)
eststo x3: reg bmi pm $x $z
eststo x4: reg bmi pm $x $z, vce(cluster countyid)
esttab x1 x2 x3 x4 using 1.ols.rtf, star(* .1 ** .05 *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)  

}


*********************************************************
                  ** 2 RD validity**
{

*treatment
gen d = 0
replace d = 1 if latitude > 33
*assignment standarlised
gen v = latitude - 33
gen vd = v*d
*loop for polynomial fit order
local i=1
forvalues i=2/4 {
gen v`i'=v^`i'
gen v`i'd=v`i'*d
}

}
	*** 2.1 obesity RD
{
*1 visualising the discontinuity
twoway (scatter bmi v, msymbol(+) msize(*0.4) mcolor(black*0.3)), tit(scatter)
graph save 2.1.1.b1.gph,  replace
cmogram bmi v,tit(Linear fit) cut(0) scatter lineat(0) lfitci
graph save 2.1.1.b2.gph,  replace
cmogram bmi v,tit(Polynomial fit) cut(0) scatter lineat(0) qfitci
graph save 2.1.1.b3.gph,  replace

graph combine 2.1.1.b1.gph 2.1.1.b2.gph 2.1.1.b3.gph
graph save 2.1.1.bcompare.gph,  replace

*2 Covariates are not discontinuous at the threshold
** non-parametric (placebo outcomes)<$x>
eststo x1: rdrobust edu v,c(0) kernel(uni) bwselect(mserd) all
eststo x2: rdrobust gender_dummy v,c(0) kernel(uni) bwselect(mserd) all
eststo x3: rdrobust yob v,c(0) kernel(uni) bwselect(mserd) all
esttab x1 x2 x3 using "C:\Users\maryh\Desktop\results\2.1.2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
//not significant= not reject H0 = satisfy

}
	*** 2.2 air pollution RD
{
*1 visualising the discontinuity
twoway (scatter pm v, msymbol(+) msize(*0.4) mcolor(black*0.3)), tit(scatter)
graph save 2.2.1.p1.gph,  replace
cmogram pm v,tit(Linear fit) cut(0) scatter lineat(0) lfitci
graph save 2.2.1.p2.gph,  replace
cmogram pm v,tit(Polynomial fit) cut(0) scatter lineat(0) qfitci
graph save 2.2.1.p3.gph,  replace

graph combine 2.2.1.p1.gph 2.2.1.p2.gph 2.2.1.p3.gph
graph save 2.2.1.pcompare.gph,  replace

*2 Covariates are not discontinuous at the threshold
** non-parametric (placebo outcomes)<$zz>
preserve
bys countyid v(pm): g n=_n
bys countyid v(pm): g N=_N
keep if n==N
drop n N

longitude area lexp employees population ldep lfa"
eststo x1: rdrobust longitude v,c(0) kernel(uni) bwselect(mserd) all
eststo x2: rdrobust area v,c(0) kernel(uni) bwselect(mserd) all
eststo x3: rdrobust lexp v,c(0) kernel(uni) bwselect(mserd) all
eststo x4: rdrobust employees v,c(0) kernel(uni) bwselect(mserd) all
eststo x5: rdrobust population v,c(0) kernel(uni) bwselect(mserd) all
eststo x6: rdrobust ldep v,c(0) kernel(uni) bwselect(mserd) all
eststo x7: rdrobust lfa v,c(0) kernel(uni) bwselect(mserd) all

esttab x1 x2 x3 x4 x5 x6 x7 using "C:\Users\maryh\Desktop\results\2.2.2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
//not significant= not reject H0 = satisfy
restore
}
	*** 2.3 latitude
//Imprecise Control over the Forcing Variable
{
preserve
bys countyid v: g n=_n
bys countyid v: g N=_N
keep if n==N
drop n N
*1 histogram
histogram v, w(1) lcolor(brown) fcolor(gs16) title("Latitude") xtitle("Degree to boundary")
graph save 2.3.1.gph, replace
*2 McCrary test
DCdensity v, breakpoint(0) gen(Xj Yj r0 fhat se_fhat)
graph save 2.3.2(1).gph, replace
return list
//scalars:
//          r(bandwidth) =  6.532940382404641
//            r(binsize) =  1.309181984318706
//                 r(se) =  .4804624939814367
//              r(theta) =  -.0831945215128846
gen t= -.0831945215128846/.4804624939814367
//t=-.1731551
display 2*ttail(99, t)
//p=0.86288306

gen upfhat=fhat+1.645*se_fhat
gen lowfhat=fhat-1.645*se_fhat
twoway (rarea upfhat lowfhat r0 if r0<0, sort fcolor(gs12) lcolor(gs12)) (rarea upfhat lowfhat r0 if r0>0, sort fcolor(gs12) lcolor(gs12)) (line fhat r0 if r0<0, lcolor(red)) (line fhat r0 if r0>0, lcolor(blue)) (scatter Yj Xj if Yj>0, mcolor(gs4) msymbol(circle_hollow)), ytitle("Density") xtitle("") xline(0) legend(off) 
graph save 2.3.2(2).gph, replace

*3 rddensity
asdoc rddensity v, c(0) 
//p=0.9727, reject H0, no manipulation

restore
}


**********************************************************
              ** 3 RD Estimation Results **

	*** 3.1 BMI
{
*1.non-parametric:LLR(local linear regression)
**Uniform/Triangle/Epanechnikov <x1*>
{
eststo x1: rdrobust bmi v,c(0) kernel(uni) bwselect(mserd) all
eststo x2: rdrobust bmi v,c(0) covs($x $zz) kernel(uni) bwselect(mserd) all
eststo x3: rdrobust bmi v, c(0) kernel(tri) bwselect(mserd) all
eststo x4: rdrobust bmi v,c(0) covs($x $zz) kernel(tri) bwselect(mserd) all
eststo x5: rdrobust bmi v, c(0) kernel(epa) bwselect(mserd) all 
eststo x6: rdrobust bmi v,c(0) covs($x $zz) kernel(epa) bwselect(mserd) all
}
{
esttab x1 x2 x3 x4 x5 x6 using "C:\Users\maryh\Desktop\results\3.1.1obllr.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

*2.parametric: nonlinear(local polynominal regression)
**full bandwidth <p=1** 2***>
{
eststo x11:reg bmi d v vd
eststo x12:reg bmi d v vd, vce(cluster countyid)
eststo x13:reg bmi d v vd $x $zz
eststo x14:reg bmi d v vd $x $zz, vce(cluster countyid)

eststo x21:reg bmi d v v2 vd v2d
eststo x22:reg bmi d v v2 vd v2d, vce(cluster countyid)
eststo x23:reg bmi d v v2 vd v2d $x $zz
eststo x24:reg bmi d v v2 vd v2d $x $zz, vce(cluster countyid)

eststo x31:reg bmi d v v2 v3 vd v2d v3d
eststo x32:reg bmi d v v2 v3 vd v2d v3d, vce(cluster countyid)
eststo x33:reg bmi d v v2 v3 vd v2d v3d $x $zz
eststo x34:reg bmi d v v2 v3 vd v2d v3d $x $zz, vce(cluster countyid)

eststo x41:reg bmi d v v2 v3 v4 vd v2d v3d v4d
eststo x42:reg bmi d v v2 v3 v4 vd v2d v3d v4d, vce(cluster countyid)
eststo x43:reg bmi d v v2 v3 v4 vd v2d v3d v4d $x $zz
eststo x44:reg bmi d v v2 v3 v4 vd v2d v3d v4d $x $zz, vce(cluster countyid)
}
{
esttab x11 x12 x13 x14 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_full1.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x21 x22 x23 x24 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_full2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x31 x32 x33 x34 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_full3.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x41 x42 x43 x44 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_full4.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}
**optimal bandwidth
rdbwselect bmi v, c(0) kernel(uni) bwselect(mserd) // h=3.067,b=6.103
{
preserve
keep if v>= -3.067 & v<= 3.067
eststo x11:reg bmi d v vd
eststo x12:reg bmi d v vd, vce(cluster countyid)

eststo x21:reg bmi d v v2 vd v2d
eststo x22:reg bmi d v v2 vd v2d, vce(cluster countyid)

eststo x31:reg bmi d v v2 v3 vd v2d v3d
eststo x32:reg bmi d v v2 v3 vd v2d v3d, vce(cluster countyid)

eststo x41:reg bmi d v v2 v3 v4 vd v2d v3d v4d
eststo x42:reg bmi d v v2 v3 v4 vd v2d v3d v4d, vce(cluster countyid)

esttab x11 x12 x21 x22 x31 x32 x41 x42 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_opt1.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
restore
}
**optimal bandwidth + controls
rdbwselect bmi v, c(0) covs($x $zz) kernel(uni) bwselect(mserd) // h=2.406,b=4.571
{
preserve
keep if v>= -2.406 & v<= 2.406
eststo x13:reg bmi d v vd $x $zz
eststo x14:reg bmi d v vd $x $zz, vce(cluster countyid)

eststo x23:reg bmi d v v2 vd v2d $x $zz
eststo x24:reg bmi d v v2 vd v2d $x $zz, vce(cluster countyid)

eststo x33:reg bmi d v v2 v3 vd v2d v3d $x $zz
eststo x34:reg bmi d v v2 v3 vd v2d v3d $x $zz, vce(cluster countyid)

eststo x43:reg bmi d v v2 v3 v4 vd v2d v3d v4d $x $zz
eststo x44:reg bmi d v v2 v3 v4 vd v2d v3d v4d $x $zz, vce(cluster countyid)

esttab x13 x14 x23 x24 x33 x34 x43 x44 using "C:\Users\maryh\Desktop\results\3.1.2obpoly_opt2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
restore
}
}
	*** 3.2 PM2.5
{
*1.county level <99obs>
{
preserve
bys countyid v(pm): g n=_n
bys countyid v(pm): g N=_N
keep if n==N
drop n N

*non-parametric:LLR(local linear regression)
**Uniform/Triangle/Epanechnikov <x1 x3***>
{
eststo x1: rdrobust pm v,c(0) kernel(uni) bwselect(mserd) all
eststo x2: rdrobust pm v,c(0) covs($zz) kernel(uni) bwselect(mserd) all
x
eststo x3: rdrobust pm v, c(0) kernel(tri) bwselect(mserd) all
eststo x4: rdrobust pm v,c(0) covs($zz) kernel(tri) bwselect(mserd) all
eststo x5: rdrobust pm v, c(0) kernel(epa) bwselect(mserd) all 
eststo x6: rdrobust pm v,c(0) covs($zz) kernel(epa) bwselect(mserd) all
}
{
esttab x1 x2 x3 x4 x5 x6 using "C:\Users\maryh\Desktop\results\3.2.1pmllr.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

*parametric: nonlinear(local polynominal regression)
**full bandwidth <x2**>
{
eststo x11:reg pm d v vd
eststo x12:reg pm d v vd, vce(cluster countyid)
eststo x13:reg pm d v vd $zz
eststo x14:reg pm d v vd $zz, vce(cluster countyid)

eststo x21:reg pm d v v2 vd v2d
eststo x22:reg pm d v v2 vd v2d, vce(cluster countyid)
eststo x23:reg pm d v v2 vd v2d $zz
eststo x24:reg pm d v v2 vd v2d $zz, vce(cluster countyid)

eststo x31:reg pm d v v2 v3 vd v2d v3d
eststo x32:reg pm d v v2 v3 vd v2d v3d, vce(cluster countyid)
eststo x33:reg pm d v v2 v3 vd v2d v3d $zz
eststo x34:reg pm d v v2 v3 vd v2d v3d $zz, vce(cluster countyid)

eststo x41:reg pm d v v2 v3 v4 vd v2d v3d v4d
eststo x42:reg pm d v v2 v3 v4 vd v2d v3d v4d, vce(cluster countyid)
eststo x43:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz
eststo x44:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz, vce(cluster countyid)
}
{
esttab x11 x12 x13 x14 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_full1.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x21 x22 x23 x24 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_full2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x31 x32 x33 x34 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_full3.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x41 x42 x43 x44 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_full4.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

**optimal bandwidth
rdbwselect pm v, c(0) kernel(uni) bwselect(mserd) // h=2.337,b=4.474
{
eststo x11:reg pm d v vd if v>= -2.337 & v<= 2.337
eststo x12:reg pm d v vd if v>= -2.337 & v<= 2.337, vce(cluster countyid)

eststo x21:reg pm d v v2 vd v2d if v>= -2.337 & v<= 2.337
eststo x22:reg pm d v v2 vd v2d if v>= -2.337 & v<= 2.337, vce(cluster countyid)

eststo x31:reg pm d v v2 v3 vd v2d v3d if v>= -2.337 & v<= 2.337
eststo x32:reg pm d v v2 v3 vd v2d v3d if v>= -2.337 & v<= 2.337, vce(cluster countyid)

eststo x41:reg pm d v v2 v3 v4 vd v2d v3d v4d if v>= -2.337 & v<= 2.337
eststo x42:reg pm d v v2 v3 v4 vd v2d v3d v4d if v>= -2.337 & v<= 2.337, vce(cluster countyid)

esttab x11 x12 x21 x22 x31 x32 x41 x42 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_opt1.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

**optimal bandwidth + controls
rdbwselect pm v, c(0) covs($zz) kernel(uni) bwselect(mserd) // h=2.112,b=3.211
{
eststo x13:reg pm d v vd $zz if v>= -2.112 & v<= 2.112
eststo x14:reg pm d v vd $zz if v>= -2.112 & v<= 2.112, vce(cluster countyid)

eststo x23:reg pm d v v2 vd v2d $zz if v>= -2.112 & v<= 2.112
eststo x24:reg pm d v v2 vd v2d $zz if v>= -2.112 & v<= 2.112, vce(cluster countyid)

eststo x33:reg pm d v v2 v3 vd v2d v3d $zz if v>= -2.112 & v<= 2.112
eststo x34:reg pm d v v2 v3 vd v2d v3d $zz if v>= -2.112 & v<= 2.112, vce(cluster countyid)

eststo x43:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz if v>= -2.112 & v<= 2.112
eststo x44:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz if v>= -2.112 & v<= 2.112, vce(cluster countyid)

esttab x13 x14 x23 x24 x33 x34 x43 x44 using "C:\Users\maryh\Desktop\results\3.2.1pmpoly_opt2.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

restore
}

*2.individual level <5940obs>
//no npnpara & optl bandwidth(h), cuz not enough variability to compute the loc. poly. h
*parametric: nonlinear(local polynominal regression)
**full bandwidth <***>
{
eststo x11:reg pm d v vd
eststo x12:reg pm d v vd, vce(cluster countyid)
eststo x13:reg pm d v vd $zz
eststo x14:reg pm d v vd $zz, vce(cluster countyid)

eststo x21:reg pm d v v2 vd v2d
eststo x22:reg pm d v v2 vd v2d, vce(cluster countyid)
eststo x23:reg pm d v v2 vd v2d $zz
eststo x24:reg pm d v v2 vd v2d $zz, vce(cluster countyid)

eststo x31:reg pm d v v2 v3 vd v2d v3d
eststo x32:reg pm d v v2 v3 vd v2d v3d, vce(cluster countyid)
eststo x33:reg pm d v v2 v3 vd v2d v3d $zz
eststo x34:reg pm d v v2 v3 vd v2d v3d $zz, vce(cluster countyid)

eststo x41:reg pm d v v2 v3 v4 vd v2d v3d v4d
eststo x42:reg pm d v v2 v3 v4 vd v2d v3d v4d, vce(cluster countyid)
eststo x43:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz
eststo x44:reg pm d v v2 v3 v4 vd v2d v3d v4d $zz, vce(cluster countyid)
}
{
esttab x11 x12 x13 x14 using "C:\Users\maryh\Desktop\results\3.2.2pmpoly_full11.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x21 x22 x23 x24 using "C:\Users\maryh\Desktop\results\3.2.2pmpoly_full22.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x31 x32 x33 x34 using "C:\Users\maryh\Desktop\results\3.2.2pmpoly_full33.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
esttab x41 x42 x43 x44 using "C:\Users\maryh\Desktop\results\3.2.2pmpoly_full44.rtf", star(* .1 ** .05  *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)
}

}


***********************************************************
                    ** 4 2SlS **
{
global s "v vd"
global s2 "v v2 vd v2d"
global s3 "v v2 v3 vd v2d v3d"
global s4 "v v2 v3 v4 vd v2d v3d v4d"

*separate
asdoc ivregress 2sls bmi $s (pm = d $s), first
asdoc ivregress 2sls bmi $s (pm = d $s), vce(cluster countyid) first
asdoc ivregress 2sls bmi $s $x $zz (pm = d $s $x $zz), first
asdoc ivregress 2sls bmi $s $x $zz (pm = d $s $x $zz), vce(cluster countyid) first


*combine
eststo x1: ivregress 2sls bmi $s (pm = d $s), first
eststo x2: ivregress 2sls bmi $s (pm = d $s), vce(cluster countyid) first
eststo x3: ivregress 2sls bmi $s $x $zz (pm = d $s $x $zz), first
eststo x4: ivregress 2sls bmi $s $x $zz (pm = d $s $x $zz), vce(cluster countyid) first

esttab x1 x2 x3 x4 using 4.2sls.rtf, star(* .1 ** .05 *** .01) nogap nonumber replace se(%5.4f) ar2 aic(%10.4f) bic(%10.4f)  

}


***********************************************************
                *** 5 RD Robust check***

	**** 5.1 Placebo (cutoff points)
{

**1 pm
//Quintiles of left and right (20%, 40%, 60%, 80%)
//8 placebo cutoffs not sig, no treatment effect
{
preserve
bys countyid v(pm): g n=_n
bys countyid v(pm): g N=_N
keep if n==N
drop n N

sum v
local vmax = r(max)
local vmin = r(min)

forvalues i=1(1)4{
local jr=`vmax' * (`i'/5)
local jl=`vmin' * (`i'/5)
rdrobust pm v,c(`jr') kernel(uni) bwselect(mserd)
estimates store jr`i'
rdrobust pm v,c(`jl') kernel(uni) bwselect(mserd)
estimates store jl`i'
}

rdrobust pm v,c(0) kernel(uni) bwselect(mserd)
estimates store jbaseline

local vlist "jl4 jl3 jl2 jl1 jbaseline jr1 jr2 jr3 jr4"

coefplot `vlist',  yline(0) drop(_cons) vertical
graph save 5.1.1.gph, replace

restore
}

**2 bmi <omit TAT>
{
preserve
forvalues i=1(1)4{
local jr=`vmax'*(`i'/5)
local jl=`vmin'*(`i'/5)
rdrobust bmi v,c(`jr') kernel(uni) bwselect(mserd)
estimates store jr`i'
rdrobust bmi v,c(`jl') kernel(uni) bwselect(mserd)
estimates store jl`i'
}

rdrobust bmi v,c(0) kernel(uni) bwselect(mserd)
estimates store jbaseline

local vlist "jl4 jl3 jl2 jl1 jbaseline jr1 jr2 jr3 jr4"
coefplot `vlist',  yline(0) drop(_cons) vertical
graph save 5.1.2.gph, replace

restore
}

}
	**** 5.2 bandwidths
**1 pm
{
preserve
bys countyid v(pm): g n=_n
bys countyid v(pm): g N=_N
keep if n==N
drop n N
rdbwselect pm v, kernel(uni) bwselect(mserd)
local h = 2.5
forvalues i=1/6{
local hrobust=`h'*0.25*`i'
rdrobust pm v,h(`hrobust') kernel(uni) bwselect(mserd)
estimates store hrob`i'
}

local vlist "hrob2 hrob3 hrob4 hrob5 hrob6"
coefplot `vlist',  yline(0) drop(_cons) vertical
graph save 5.2.1.gph, replace
restore
}

**2 bmi <omit TAT>
{
preserve
rdbwselect bmi v, kernel(uni) bwselect(mserd)
local h = 3
forvalues i=1/6{
local hrobust=`h'*0.1*`i'
rdrobust bmi v,h(`hrobust') kernel(uni) bwselect(mserd)
estimates store hrob`i'
}

local vlist "hrob2 hrob3 hrob4 hrob5 hrob6"
coefplot `vlist',  yline(0) drop(_cons) vertical
graph save 5.2.2.gph, replace
restore
}

log close
