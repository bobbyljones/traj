/* 
		traj	

		Copyright (C) 2023 Bobby L Jones <bobbyljones@gmail.com>
		
		This source code is subject to the terms of the 3-Clause
		BSD License (Open Source Initiative). The license can be
		obtained from https://opensource.org/licenses/BSD-3-Clause.
  	
		May 2023  replace model string with integer code on command line
				  for easier parsing
		Apr 2023  remove ocovbygroup, scoreci, reps options
		Feb 2023  rhobygroup option
		Aug 2022  AR1 option
		May 2020  add tol option
		Feb 2020  check model is present for all models in multmodel
		Nov 2019  altstart option
		Jun 2019  twostep option for joint distal outcome traj model.
		Nov 2018  support for exposure and tcov for models 4 - 6.
		Oct 2018  support for dropout2.
		May 2018  remove abbreviations.
		Feb 2018  add expos2-3, beta model, multinomial w/ baseoutcome
           		 (chg version 9 to 11, use _rmcoll)
*/

program traj, eclass

	version 11

	syntax [if] [in][, Model1(str) Var1(varlist numeric) Indep1(varlist numeric) Order1(numlist int >=-1 <=5) MIN1(real 0.0) MAX1(real 1e-38) Iorder1(numlist int >=-1 <=5) Rorder1(numlist int >=-1 <=1) Risk1(varlist numeric) Weight1(varlist numeric) Expos1(varlist numeric) Tcov1(varlist numeric) Plottcov1(string) Refgroup1(int 1) Dropout1(numlist int >=-1 <=2)	Dcov1(varlist numeric) OBsmar1(varlist numeric) Outcome1(varlist numeric) OModel1(str) Baseoutcome1(int -99998) Ocov1(varlist numeric) model2(str) var2(varlist numeric) indep2(varlist numeric) order2(numlist int >=-1 <=5) min2(real 0.0) max2(real 1e-38) iorder2(numlist int >=-1 <=5)  Dropout2(numlist int >=-1 <=2)	Dcov2(varlist numeric) expos2(varlist numeric) risk2(varlist numeric) refgroup2(int 1) model3(str) var3(varlist numeric) indep3(varlist numeric) order3(numlist int >=-1 <=5) min3(real 0.0) max3(real 1e-38) iorder3(numlist int >=-1 <=5) expos3(varlist numeric) model4(str) var4(varlist numeric) indep4(varlist numeric) order4(numlist int >=-1 <=5) min4(real 0.0) max4(real 1e-38) iorder4(numlist int >=-1 <=5) Expos4(varlist numeric)  model5(str) var5(varlist numeric) indep5(varlist numeric) order5(numlist int >=-1 <=5) min5(real 0.0) max5(real 1e-38) iorder5(numlist int >=-1 <=5) Expos5(varlist numeric) model6(str) var6(varlist numeric) indep6(varlist numeric) order6(numlist int >=-1 <=5) min6(real 0.0) max6(real 1e-38) iorder6(numlist int >=-1 <=5) Expos6(varlist numeric) MULTRisk(varlist numeric) Outofsample(varlist numeric) *]  
    
	local 0 `", `options'"'

 	marksample touse
	quietly count if `touse'

	syntax [ ,Tcov2(varlist numeric) Tcov3(varlist numeric) Tcov4(varlist numeric) Tcov5(varlist numeric) Tcov6(varlist numeric) Plottcov2(string) Plottcov3(string) Plottcov4(string) Plottcov5(string) Plottcov6(string) MULTGroups(integer 0) PROBUPDATES DETAIL NOVAR TRACE CI SIGMABYGROUP RHOBYGROUP AR1 NDERIV TWOSTEP ALTSTART Lower(string) Upper(string) Start(string) tol(real 0.0) ]
	
	local traj_start_len = 0
	local traj_lower_len = 0
	local traj_upper_len = 0
	local maxmodels = 6
  
	if `r(N)' == 0 {
		error 2000
	}

	local rN = `r(N)'

	cap qui drop _traj_*
	
	forval i = 3/`maxmodels' {
		local numdropout`i' = 0
	}
  
	forval i = 2/`maxmodels' {
		if ( length( "`model`i''" ) == 0 ) {
			local model`i' = "none"
			local mdl`i' = 0
		}
	}
		
	if ( length( "`omodel1'" ) == 0 ) {
			local omodel1 = "none"
			local omdl1 = 0
	}	
	
	local nummultrisk : word count `multrisk' 
	local numocov1 : word count `ocov1'
	local numdcov1 : word count `dcov1' 
	local numdcov2 : word count `dcov2' 
	local numoutcome1 : word count `outcome1' 
	local numobsmar1 : word count `obsmar1'
	local numdropout1 : word count `dropout1'
	local numdropout2 : word count `dropout2'
		
	local option_altstart = 0
	if "`altstart'" == "altstart" {
		local option_altstart = 1
	}

	local itdetail = 0
	if "`detail'" == "detail" {
		local itdetail = 1
	}

	local ittrace = 0
	if "`trace'" == "trace" {
		local ittrace = 1
	}

	local itnovar = 0
	if "`novar'" == "novar" {
		local itnovar = 1
	}

	local itprobupdates = 0
	if "`probupdates'" == "probupdates" {
		local itprobupdates = 1
	}
	
	local option_ci = 0
	if "`ci'" == "ci" {
		local option_ci = 1
	}

	local option_sigmabygroup = 0
	if "`sigmabygroup'" == "sigmabygroup" {
		local option_sigmabygroup = 1
	}
 
    local option_rhobygroup = 0
	if "`rhobygroup'" == "rhobygroup" {
		local option_rhobygroup = 1
	}
 
	local option_ar1 = 0
	if "`ar1'" == "ar1" {
		local option_ar1 = 1
	}

	local option_nderiv = 0
	if "`nderiv'" == "nderiv" {
		if "`twostep'" == "twostep" {
			di as err "deriv required with twostep. Option ignored."
		} 
		if "`twostep'" != "twostep" {
			local option_nderiv = 1
		}
	}
	
	local option_twostep = 0
	if "`twostep'" == "twostep" {
		local option_twostep = 1
	}
	
	if "`twostep'" == "twostep" & ! ("`omodel1'"=="logit" | "`omodel1'"=="mlogit" | ///
	 "`omodel1'"=="poisson" | "`omodel1'"=="normal" ) {  
		di as err "twostep is used with outcome modeling. Option ignored."
	}
	
	if "`start'" != "" {
	
		cap confirm matrix `start'
	
		if _rc {
			di as err "Start must be a matrix"
			exit 198
		}
	
		tempname traj_start
		local traj_start_len = colsof(`start')
		matrix traj_start = `start'
	}

	if "`lower'" != "" {
		
		cap confirm matrix `lower'
		
		if _rc {
			di as err "Lower must be a matrix"
			exit 198
		}
		
		tempname traj_lower
		local traj_lower_len = colsof(`lower')
		matrix traj_lower = `lower'
	}
	
	if "`upper'" != "" {
		
		cap confirm matrix `upper'
		
		if _rc {
			di as err "Upper must be a matrix"
			exit 198
		}
		
		tempname traj_upper
		local traj_upper_len = colsof(`upper')
		matrix traj_upper = `upper'
	}
	
	local numoutofsample : word count `outofsample'
	
	if ( `numoutofsample' > 1 ) {
		di as err "Only one variable is allowed for outofsample."
		exit 198
	}

	matrix traj_mat_opts = J( 1, 17, 0 )
	matrix traj_mat_opts[1,1] = `itprobupdates'
	matrix traj_mat_opts[1,3] = `itdetail'
	matrix traj_mat_opts[1,4] = `itnovar'
	matrix traj_mat_opts[1,5] = `ittrace'
	matrix traj_mat_opts[1,6] = `multgroups'
	matrix traj_mat_opts[1,7] = `traj_start_len'
	matrix traj_mat_opts[1,8] = `option_ci'
	matrix traj_mat_opts[1,9] = `option_sigmabygroup'
	matrix traj_mat_opts[1,10] = `traj_lower_len'
	matrix traj_mat_opts[1,11] = `traj_upper_len'
	matrix traj_mat_opts[1,13] = `option_nderiv'
	matrix traj_mat_opts[1,14] = `option_twostep'
	matrix traj_mat_opts[1,15] = `option_altstart'
	matrix traj_mat_opts[1,16] = `option_ar1'
	matrix traj_mat_opts[1,17] = `option_rhobygroup'
		
	matrix traj_mat_plottcov_len = J( 1, `maxmodels', 0 )
	
	matrix traj_mat_refgp = J( 1, 2, 0 )
	matrix traj_mat_refgp[1,1] = real("`refgroup1'")
	matrix traj_mat_refgp[1,2] = real("`refgroup2'")
	
	matrix traj_mat_min = J( 1, `maxmodels', 0 )
	matrix traj_mat_max = J( 1, `maxmodels', 0 )
	
	matrix traj_mat_tol = J( 1, 1, 0 )
	matrix traj_mat_tol[1,1] = real("`tol'")
		
	if "`omodel1'" == "normal" {
		local omdl1 = 1
	}
	
	if "`omodel1'" == "logit" {
		cap assert `outcome1' == 0 | `outcome1' == 1 | `outcome1' >= .
		if _rc { 
			di as err "`outcome1' is not 0/1 (logit model)."
			exit 198
		}
		local omdl1 = 2
	}			

	if "`omodel1'" == "mlogit" {
		if ( `baseoutcome1' != -99998 ) {	
	
			_rmcoll `outcome1' if `touse', mlogit baseoutcome(`baseoutcome1')		
			matrix _traj_Aout = r(out)
		}
		if ( `baseoutcome1' == -99998 )	{
			_rmcoll `outcome1' if `touse', mlogit			
			matrix _traj_Aout = r(out)
		}
		local noutc = r(k_out)
		if (`noutc' == 1) {
			exit 148	    
		}
		local base = r(baseout)
		local omdl1 = 5
	}
			
	if "`omodel1'" == "poisson" {
		cap assert `outcome1' >= 0
		if _rc { 
			di as err "Outcome has negative values (poisson model)."		
			exit 198
		}
		local omdl1 = 3
	}	

	if ( `numoutcome1' > 0 ) {
		if ( "`omodel1'" != "normal" & ///
				"`omodel1'" != "poisson" & "`omodel1'" != "mlogit" &  ///
				"`omodel1'" != "logit" ) {			
			di as err "Outcome model is not normal, poisson, mlogit or logit."
			exit 198
		} 
	}
	
	if ( `numobsmar1' > 1 ) {
		di as err "Only one variable is allowed for obsmar."
		exit 198
	}
		
	if "`omodel1'" == "mlogit" {
		
		local tmpxx = 	_traj_Aout[1,1]		
		matrix traj_mat_mlogitoutc = J( 1, 2, 0 )
	
	/* matrix traj_mat_mlogitoutc[1,1] = real("`base'") */
	
		matrix traj_mat_mlogitlevels = r(out)	
		matrix traj_mat_mlogitoutc[1,1] = real("`tmpxx'")		
		matrix traj_mat_mlogitoutc[1,2] = real("`noutc'")	
	}
	
	if ( `numocov1' > 0 & `numoutcome1' == 0 ) {
		di as err "Option ocov is only valid when an outcome is specified."
		exit 198	
	}
	
	local numweight1 : word count `weight1'

	if ( `numweight1' > 1 ) {
		di as err "Only one variable is allowed for weight."
		exit 198
	}

	if ( `nummultrisk' > 0 & `multgroups' < 2 ) {
		di as err "Multrisk is only valid when multgroups is > 1."
    	exit 198
	}

	local numrorder1 0
	tokenize "`rorder1'"
	
	while "`1'" != "" {
		local numrorder1 = `numrorder1' + 1
		local c`numrorder1' `1'
		local rordr1 "`rordr1' "`1'""
		mac shift
	}  
	
	if ( `numrorder1' != 0 ) {
		matrix traj_mat_rorder1 = J( 1, `numrorder1', 0 )
		tokenize "`rorder1'"
		forval i = 1/`numrorder1' {
			matrix traj_mat_rorder1[1,`i'] = real("`1'")
			mac shift
		}
	}

/*
	if ( `numrorder1' > 0 ) {
		di as err "rorder is not working at this time. Contact bjones@andrew.cmu.edu"
		exit 198
	}
*/	
	
	forval i = 1/`maxmodels' { 
	
		local numrisk`i' : word count `risk`i'' 
		local numexpos`i' : word count `expos`i'' 
		local numtcov`i' : word count `tcov`i''
		matrix traj_mat_min[1,`i'] = real("`min`i''")
		matrix traj_mat_max[1,`i'] = real("`max`i''")
		
		if "`model`i''" != "none" {
			if "`model`i''" != "cnorm" & "`model`i''" != "zip" & "`model`i''" != "beta" & "`model`i''" != "logit" {
				di as err "Model`i' is not cnorm, logit, zip, or beta."
				exit 198
			}

			if "`model`i''"=="beta" {
				tokenize "`var`i''"			
				while "`1'" != "" {
					local vr `1'
					cap assert ( `vr' >= 0 & `vr' <= 1 ) | `vr' >= .   
					if _rc { 
						di as err "`vr' is not within 0 and 1 (beta model)."
						exit 198
					}
					mac shift
				}
				local mdl`i' = 4
			}
			
			if "`model`i''" == "cnorm" {
				if `max`i'' == 1e-38 {
					di as err  "max`i' required for model`i' cnorm."
					exit 198
				}
				tokenize "`var`i''"	
				while "`1'" != "" {
					local vr `1'
					cap sum `vr'
					cap assert r(min) >= `min`i'' & r(max) <= `max`i''
					if _rc { 
						di as err "`vr' is not within min`i' = `min`i'' and max`i' = `max`i''."
						exit 198
					}
					mac shift
				}
				local mdl`i' = 1
			}
	
			if "`model`i''" == "logit" {
				tokenize "`var`i''"
				while "`1'" != "" {
					local vr `1'
					cap assert `vr' == 0 | `vr' == 1 | `vr' >= .   
					if _rc { 
						di as err "`vr' is not 0/1 (logit model)."
						exit 198
					}
					mac shift
				}
				local mdl`i' = 2
			}
   
			if "`model`i''" == "zip" {
				tokenize "`var`i''"
				while "`1'" != "" {
					local vr `1'
					cap assert `vr' >= 0  
					if _rc { 
						di as err "`vr' has negative values (zip model)."
						exit 198
					}
					mac shift
				}
				local mdl`i' = 3
			}
		}

		if "`ar1'" == "ar1" & !( "`model`i''" == "cnorm" | "`model`i''" == "none" ) {
			di as err "The ar1 option is only supported with the cnorm model."
			exit 198
		}
	
		local numindep`i' : word count `indep`i''
		local numvar`i' : word count `var`i'' 

		if ( `numindep`i'' != `numvar`i'' ) {
			di as err "The number of variables `numvar`i'' in var`i' and `numindep`i'' in indep`i' must match."
			exit 198
		}	
		
		local plottcov`i'_len = 0
 
		if "`plottcov`i''" != "" {			
			cap confirm matrix `plottcov`i''			
			if _rc {
				di as err "Plottcov`i' must specify a matrix"
				exit 198
			}    
			matrix traj_mat_plottcov`i' = `plottcov`i''
			matrix traj_mat_plottcov_len[1,`i'] = colsof(`plottcov`i'')
		}
  
 		if ( `numexpos`i'' > 1 ) {
			
			tokenize "`expos`i''"
			
			while "`1'" != "" {
				local vr `1'
				cap assert `vr' > 0  
				if _rc { 
					di as err "Exposure variable `vr' is not > 0."
					exit 198
				}
				mac shift
			}
			
			if ( `numexpos`i'' != `numvar`i'' ) {
				di as err "The number of variables in var`i' and expos`i' must match."
				exit 198
			}
	 		
			if ( "`model`i''" != "zip" ) {
				di as err "Exposure is only supported with the zip model."
				exit 198
			}	
		}
		
		if ( `numtcov`i'' != 0 & `numtcov`i'' / `numindep`i'' != int( `numtcov`i'' / `numindep`i'' ) ) {
			di as err "The number of variables in tcov`i' and var`i' must match or be a multiple."
			exit 198
		}
		
		local numorder`i' 0
		tokenize "`order`i''"
		
		while "`1'" != "" {
			local numorder`i' = `numorder`i'' + 1
			local c`numorder`i'' `1'
			local ordr`i' "`ordr`i'' "`1'""
			mac shift
		}
		
		if ( `multgroups' != 0 ) {
			if ( `numorder`i'' != 0 & `numorder`i'' != `multgroups' ) {
				di as err "The number of multgroups and orders for all models need to agree."
				exit 198
			}
		}
		
		if "`model`i''" != "none" & `numorder`i'' == 0 {
			di as err "model`i' has been specified but order`i' is missing."
			exit 198
		}
		
		if ( `numorder`i'' != 0 ) {
			matrix traj_mat_order`i' = J( 1, `numorder`i'', 0 )
			tokenize "`order`i''"
		
			forval j = 1/`numorder`i'' {
				matrix traj_mat_order`i'[1,`j'] = real("`1'")
				mac shift
			}
		}

		local numiorder`i' 0
		tokenize "`iorder`i''"
	
		while "`1'" != "" {
			local numiorder`i' = `numiorder`i'' + 1
			local c`numiorder`i'' `1'
			local iordr`i' "`iordr`i'' "`1'""
			mac shift
		}	  
		
		if ( `numiorder`i'' != 0 ) {
			matrix traj_mat_iorder`i' = J( 1, `numiorder`i'', 0 )
			tokenize "`iorder`i''"
			forval j = 1/`numiorder`i'' {
				matrix traj_mat_iorder`i'[1,`j'] = real("`1'")
				mac shift
			}
		}
				
 		if ( ( "`model`i''" != "zip" & "`model`i''" != "beta" ) & `numiorder`i'' != 0 ) {
			di as err "Option iorder is only allowed with the zip model."
			exit 198
		}
		
		if ( `numiorder`i'' != `numorder`i'' & `numiorder`i'' != 1 & `numiorder`i'' != 0 ) {
			di as err "There should be 1 or `numorder`i'' iorders for iorder`i'."
			exit 198
		}
	}
	
	if ( `numrisk1' > 0 & `multgroups' > 0 ) {
		di as err "Risk1 is not supported with multgroups. Try multrisk."
		exit 198
	}

	if ( `refgroup1' != 1 ) {
		if ( `numrisk1' == 0 & `nummultrisk' == 0 ) {
			di as err "Refgroup is only valid when risk or multrisk is specified."
			exit 198
		}
		if ( `refgroup1' > `numorder1' ) {
			di as err "Refgroup exceeds the number of groups, `numorder1'."
			exit 198
		}
		if ( `refgroup1' < 1 ) {
			di as err "Refgroup is less than 1."
			exit 198
		}
		if ( `nummultrisk' != 0 & `refgroup1' > `multgroups' ) {
			di as err "Refgroup exceeds the number of groups."
			exit 198
		}
	}
		
	if ( `refgroup2' != 1 ) {
		if ( `numrisk2' == 0 ) {
			di as err "Refgroup2 is only valid when risk2 is specified."
			exit 198
		}
		
		if ( `refgroup2' < 1 ) {
			di as err "Refgroup is less than 1."
			exit 198
		}

		if ( `refgroup2' > `numorder2' ) {
			di as err "Refgroup2 exceeds the number of groups, `numorder2'."
			exit 198
		}
	}

	if ( `numdcov1' != 0 & `numdcov1' / `numindep1' != int( `numdcov1' / `numindep1' ) ) {
   	
		di as err "The number of variables in dcov and var must match or be a multiple."
		exit 198
	}
	
	if ( `numdcov2' != 0 & `numdcov2' / `numindep2' != int( `numdcov2' / `numindep2' ) ) {
   	
		di as err "The number of variables in dcov2 and var2 must match or be a multiple."
		exit 198
	}
	
	if ( `numrisk1' > 0 & `numorder1' < 2 ) {
		di as err "Risk is only valid when the number of groups is > 1."
		exit 198
	}
	
	if ( `numdropout1' != 0 ) {
  	
		matrix traj_mat_dropoutorder1 = J( 1, `numdropout1', 0 )
		tokenize "`dropout1'"
		
		forval i = 1/`numdropout1' {
		
			matrix traj_mat_dropoutorder1[1,`i'] = real("`1'")
			
			if ( traj_mat_dropoutorder1[1,`i'] > 0 & `multgroups' > 0 ) {
				di as err "Previous response dropout rate dependence is not supported in multi-traj models."			
				exit 198
			}
			mac shift
		}
	}

	if ( `numdropout2' != 0 ) {
  	
		matrix traj_mat_dropoutorder2 = J( 1, `numdropout2', 0 )
		
		tokenize "`dropout2'"
		
		forval i = 1/`numdropout2' {		
			matrix traj_mat_dropoutorder2[1,`i'] = real("`1'")						
			mac shift
		}
	}

	if ( `plottcov1_len' != 0 & `plottcov1_len' / `numindep1' != int( `plottcov1_len' / `numindep1' ) ) {
		di as err "There should be `numtcov1' values for plottcov1."
		exit 198
	}

 	if ( `plottcov2_len' != 0 & `plottcov2_len' / `numindep2' != int( `plottcov2_len' / `numindep2' ) ) {
		di as err "There should be `numtcov2' values for plottcov2."
		exit 198
	}

 	if ( `plottcov3_len' != 0 & `plottcov3_len' / `numindep3' != int( `plottcov3_len' / `numindep3' ) ) {
		di as err "There should be `numtcov3' values for plottcov3."
		exit 198
	}

	if ( "`model1'" != "cnorm" & `numrorder1' != 0 ) {
		di as err "Option rorder is only supported for the cnorm model."
		exit 198
	}
	
	if ( ( "`model1'" == "beta" | "`model2'" == "beta" )  & `numorder2' != 0 & `multgroups' == 0 ) {
		di as err "The beta model is not supported with the joint traj model."
		exit 198
	}
	
	if ( `numorder2' != 0 & `numoutcome1' != 0 & `multgroups' == 0 ) {
		di as err "The outcome option is not supported with the joint traj model."
		exit 198
	}
	
	if ( "`ci'" == "ci" & `numorder2' != 0 & `multgroups' == 0 ) {
		di as err "Option ci is not supported with the joint traj model."
		exit 198
	}

	if ( "`probupdates'" == "probupdates" & `numorder2' != 0 & `multgroups' == 0 ) {
		di as err "Option probupdates is not supported for the joint traj model."
		exit 198
	}	  
	**********************************************************************************

	qui gen _traj_Group = .         
  
	local out1 "_traj_Group"     

	forval j = 1/`numorder1' {
		qui gen _traj_ProbG`j' = .        
		local out1 "`out1' _traj_ProbG`j'"
	}
	
	if ( `numoutcome1' > 0 & "`probupdates'" == "" ) {
	 
		qui gen _traj_Outcome = .		
		local out1 "`out1' _traj_Outcome" 	
		
		if "`omodel1'"=="mlogit" {
	  
			forval j = 1/`noutc' {	   
				local level = _traj_Aout[1,`j']     
				qui gen _traj_ProbO`level' = .     
				local out1 "`out1' _traj_ProbO`level'"
			}
		}
	}
	
	local outMat `out1'  
			
	if ( `numorder2' != 0 & `multgroups' == 0 ) {
	
		qui gen _traj_Model2_Group = .        	
		local out2 "_traj_Model2_Group"     
		
		forvalues j = 1/`numorder2' {		
			qui gen _traj_Model2_ProbG`j' = .        			
			local out2 "`out2' _traj_Model2_ProbG`j'"		
		}
		
		local outMat "`outMat' `out2'"
	}
	
	**********************************************************************************
	
	if "`probupdates'" == "probupdates" {
		
		local probupdatewaves = `numindep1' - 1
		
		forvalues j = 1/`probupdatewaves' {		
			qui gen _traj_Group_T`j' = .		
			local out1updates "`out1updates' _traj_Group_T`j'"    
			
			forvalues k = 1/`numorder1' {			
				qui gen _traj_ProbG`k'_T`j' = .        	
				local out1updates "`out1updates' _traj_ProbG`k'_T`j'"				
			}	
		}
		
		if ( `numoutcome1' > 0 ) {		
			qui gen _traj_Outcome = .		
			local out1updates "`out1updates' _traj_Outcome" 	
		
			if "`omodel1'"=="mlogit" {	  
				matrix out = r(out)
	  
				forval j = 1/`noutc' {	   
					local level = _traj_Aout[1,`j']     
					qui gen _traj_ProbO`level' = .     
					local out1updates "`out1updates' _traj_ProbO`level'"
				}
			}
			
			forvalues j = 1/`probupdatewaves' {	
				qui gen _traj_Outcome_T`j' = .				
				local out1updates "`out1updates' _traj_Outcome_T`j'"
				
				if "`omodel1'"=="mlogit" {  
					matrix out = r(out)
	  
					forval r = 1/`noutc' {
						local level = _traj_Aout[1,`r']
						qui gen _traj_ProbO`level'_T`j' = .     
						local out1updates "`out1updates' _traj_ProbO`level'_T`j'"						
					}
				}			
			}	
		}
		
		local outMat "`outMat' `out1updates'"
/*
		if ( `numorder2' != 0 & `multgroups' == 0 ) {		
			forvalues r = 1/`numindep2' {		
				qui gen _traj_Model2_Group_T`r' = .			
				local out2updates "`out2updates' _traj_Model2_Group_T`r'"		
				forvalues k = 1/`numorder2' {				
					qui gen _traj_Model2_ProbG`k'_T`r' = .        			
					local out2updates "`out2updates' _traj_Model2_ProbG`k'_T`r'"			
				}			
			}			
			local outMat "`outMat' `out2updates'"		
		}
*/
	}
	
	**********************************************************************************

	if ("`ci'" == "ci") {
		
		forvalues j = 1/`numorder1' {
			qui gen _traj_ProbG`j'_CI5 = .		
			qui gen _traj_ProbG`j'_CI95 = .						
			local prout1 "`prout1' _traj_ProbG`j'_CI5 _traj_ProbG`j'_CI95"			
		}
		
		if ( `numoutcome1' > 0 ) {
				
			if "`omodel1'"=="mlogit" {	  
				matrix out = r(out)
	  
				forval r = 1/`noutc' {	   
					local level = _traj_Aout[1,`r']     
					qui gen _traj_ProbO`level'_CI5 = .     
					qui gen _traj_ProbO`level'_CI95 = .
					local prout1 "`prout1' _traj_ProbO`level'_CI5 _traj_ProbO`level'_CI95"			
				}
			}
			else {
				qui gen _traj_Outcome_CI5 = .						
				qui gen _traj_Outcome_CI95 = .								
				local prout1 "`prout1' _traj_Outcome_CI5 _traj_Outcome_CI95"
			}
		}
		
		local outMat "`outMat' `prout1'"
			
		if "`probupdates'" == "probupdates" {
							
			forval j = 1/`probupdatewaves' {			
					
				forval k = 1/`numorder1' {			
					qui gen _traj_ProbG`k'_CI5_T`j' = .					
					qui gen _traj_ProbG`k'_CI95_T`j' = .        					
					local prout1 "`prout1' _traj_ProbG`k'_CI5_T`j' _traj_ProbG`k'_CI95_T`j'"				
				}								
							
				if ( `numoutcome1' > 0 ) {
				
					if "`omodel1'"=="mlogit" {	  
						matrix out = r(out)
	  
						forval r = 1/`noutc' {	   
							local level = _traj_Aout[1,`r']
							qui gen _traj_ProbO`level'_CI5_T`j' = .     
							qui gen _traj_ProbO`level'_CI95_T`j' = .	
							local prout1 "`prout1' _traj_ProbO`level'_CI5_T`j' _traj_ProbO`level'_CI95_T`j'"						
						}
					}						
					
					if "`omodel1'" != "mlogit" {						
						qui gen _traj_Outcome_CI5_T`j' = .						
						qui gen _traj_Outcome_CI95_T`j' = .								
						local prout1 "`prout1' _traj_Outcome_CI5_T`j' _traj_Outcome_CI95_T`j'"					
					}
				}					
			}
		}
	}
	
	local outMat "`outMat' `prout1'"
	
	**********************************************************************************

	forval i = 1/`maxmodels' {
		if `numorder`i'' != 0 {
			local plotVarLabels`i' "trajT"   
  
			forval j = 1/`numorder`i'' {
				local plotVarLabels`i' "`plotVarLabels`i'' Avg`j'"
			}
			forval j = 1/`numorder`i'' {
				local plotVarLabels`i' "`plotVarLabels`i'' Est`j'"
			}
			forval j = 1/`numorder`i'' {
				local plotVarLabels`i' "`plotVarLabels`i'' L95`j'"
				local plotVarLabels`i' "`plotVarLabels`i'' U95`j'"
			}
		}
		tokenize "`plotVarLabels`i''"
	} 

	if `numdropout1' != 0 {
		
		if ( `numdropout1' != 1 & `numdropout1' != `numorder1' ) {
			di as err "Error: there should be 1 or `numorder1' values for dropout."
			exit 198
		}
    
		forval j = 1/`numdropout1' {
			local plotVarLabels1 "`plotVarLabels1' Dropout`j'"
		} 
    
		tokenize "`plotVarLabels1'"
	} 

	if `numdropout2' != 0 {
		
		if ( `numdropout2' != 1 & `numdropout2' != `numorder2' ) {
			di as err "Error: there should be 1 or `numorder2' values for dropout."
			exit 198
		}
    
		forvalues j = 1/`numdropout2' {
		local plotVarLabels2 "`plotVarLabels2' Dropout`j'"
    } 
    
		tokenize "`plotVarLabels2'"
	}
	
	tempname b V parmData varData

	scalar BIC_n = 0
	scalar BIC_N = 0
	scalar aic = 0
	scalar loglike = 0 	
	
	matrix parmData = J( 1, 300, 0 )
	matrix varData = J( 300, 300, 0 )
	matrix groupPct1 = J( 1, `numorder1', . ) 
	
	if (`numorder2' != 0 & `multgroups' == 0 ) {
		matrix groupPct2 = J( 1, `numorder2', . ) 
	}
  
	forval i = 1/`maxmodels' {
		if `numorder`i'' != 0 {
			tempname outPlot`i' groupPct`i' 
			if `numindep`i'' == 4 * `numorder`i'' + 1 + `numdropout`i'' {
      	matrix outPlot`i' = J( 1 + `numindep`i'', 4 * `numorder`i'' + 1 + `numdropout`i'', . )
      }
      else {
      	matrix outPlot`i' = J( `numindep`i'', 4 * `numorder`i'' + 1 + `numdropout`i'', . )    
			}
		}
	}
  
  **********************************************************************************

	capture program _traj, plugin using("traj.plugin")	
	
	 plugin call _traj `var1' `indep1' `expos1'	`tcov1' `risk1' `weight1' ///
		`outcome1' `ocov1' `dcov1' `obsmar1' `outofsample' `multrisk'	///
    `var2' `indep2' `expos2' `tcov2' `risk2' 		///
    `var3' `indep3'	`expos3' `tcov3'						///
    `var4' `indep4'	`expos4' `tcov4'						///
    `var5' `indep5'	`expos5' `tcov5'						///
    `var6' `indep6'	`expos6' `tcov6'						///
    `outMat' if `touse' `in',					///									
    `mdl1'	 				///
    `numindep1' 			///
    `numrisk1' 				///
    `numweight1' 			///
    `numexpos1' 			///
    `numtcov1' 				///
    `numorder1' 			///
    `numiorder1'			///
    `numrorder1'			///
    `numdropout1' 			///
    `numdcov1'				///
    `numobsmar1'			///
    `omdl1'					///
	`numocov1'	  			///
	`numoutofsample'		///
    `mdl2'					///
    `numindep2' 			///
    `numrisk2' 				///
    `numexpos2' 			///
    `numtcov2' 				///
    `numorder2' 			///
    `numiorder2'			/// 
    `numdropout2' 			///
    `numdcov2'				///
    `mdl3'	 				///
    `numindep3' 			///
    `numorder3' 			///
    `numiorder3' 			///
    `numexpos3' 			///
    `numtcov3' 				///
    `mdl4' 					///
    `numindep4' 			///
    `numorder4' 			///
    `numiorder4' 			///		
    `numexpos4' 			///
    `numtcov4' 				///
    `mdl5'	 				///
    `numindep5' 			///
    `numorder5' 			///
    `numiorder5' 			///		
    `numexpos5' 			///
    `numtcov5' 				///
    `mdl6'	 				///
    `numindep6' 			///
    `numorder6' 			///
    `numiorder6' 			///		
	`numexpos6' 			///
    `numtcov6' 				///
    `nummultrisk'		

	ereturn clear

	local np : word count `outvars'

	matrix `b' = parmData[1,1..`np']
	matrix `V' = varData[1..`np', 1..`np']

	mat colnames `b' = `outvars' 
	mat rownames `V' = `outvars'
	mat colnames `V' = `outvars'

	ereturn post `b' `V', e(`touse')
	ereturn local cmd "traj"

	ereturn matrix groupSize1 = groupPct1
  
	if (`numorder2' != 0 & `multgroups' == 0 ) {
       ereturn matrix groupSize2 = groupPct2
	}
	
	forval i = 1/`maxmodels' {
		if `numorder`i'' != 0 {
			mat colnames outPlot`i' = `plotVarLabels`i''
			ereturn scalar numGroups`i' = `numorder`i'' 
			ereturn matrix plot`i' = outPlot`i'
		}
	}
	
	ereturn scalar BIC_N_data = BIC_N
	ereturn scalar BIC_n_subjects = BIC_n
	ereturn scalar AIC = aic
	ereturn scalar ll = loglike 
	ereturn scalar rc = rc
	ereturn local cmdline `"traj `0'"'
	trajentropy
end
/* end of traj.ado */
