/* 
	===========================================================================
	trajplot
	=========================================================================== 

	Copyright (C) 2022 Bobby L Jones <bobbyljones@gmail.com>
		
	This source code is subject to the terms of the 3-Clause
	BSD License (Open Source Initiative). The license can be
	obtained from https://opensource.org/licenses/BSD-3-Clause.
  	
	Aug 2018 summarize function to calculate group size
	Nov 2018 symxsize for traj symbol size in legend
	Feb 2019 e(groupSize1) and e(groupSize2) for percentages instead
			 of stata summarize function (problem with weights / outofsample)
*/

program trajplot
	version 9
	syntax ,[ model(integer 1) xlabel(str) ylabel(str) xtitle(str) ytitle(str) 	///
  	legendLabels(str) plotvars(varlist numeric) legendOpts(str)					///
  	CI NOLEGEND DROPOUT NOAVERAGES]   
  
	preserve

	tempname A plotVarNames    
  
	capture confirm matrix e(plot`model')         
  
	if _rc != 0 { 
	
		display as error "plot data not found" 
    
		exit _rc
	}

	if "`dropout'" == "dropout" {
  
		if "`ci'" == "ci" {
  	
			display "ci not supported for dropout probability"
  		
			local ci = "" 
		}
	}

	mat Piv = e(groupSize`model')
  
	local ng = e(numGroups`model')
  
	local numVar : word count `plotvars' 
	
	if ( `numVar' > 0 ) {
		mkmat trajT - U95`ng' if trajT != . , mat(A)
	}
	
	if ( `numVar' == 0 ) {
		
		mat A = e(plot`model')
  
		svmat A, names(col) 
	}

	local plotCmd = " " 
  
	local orderTxt = " "
  
	forvalues i = 1/`ng' {
	
		local orderTxt "`orderTxt' `i'" 
	
		local pi`i' : display %4.1f round( el( Piv, 1, `i' ), .1 )
	
		if "`dropout'" == "dropout" {
			local plotCmd "`plotCmd' (line Dropout`i' trajT, graphregion(color(white)) lwidth(thick))"
		}
		else {
			local plotCmd "`plotCmd' (line Est`i' trajT, graphregion(color(white)) lwidth(thick))"
		}
		
		local label`i' "`i'    `pi`i''%"	
	}

	if "`ci'" == "ci" {
		
		forvalues j = 1/`ng' {
			
			if "`noaverages'" != "noaverages" {
		
				local plotCmd "`plotCmd' (scatter Avg`j' trajT, mstyle(p`j'))"
			} 
			
			local plotCmd "`plotCmd' (line L95`j' trajT, lpattern(shortdash) lcolor(gs10))" 
			
			local plotCmd "`plotCmd' (line U95`j' trajT, lpattern(shortdash) lcolor(gs10))"
		}
	}
	else {
		
		forvalues j = 1/`ng' {
	
			if "`dropout'" == "dropout" {

				quietly capture confirm e Dropout`j'

				if !_rc {
			
					local plotCmd "`plotCmd' (scatter Dropout`j' trajT, mstyle(p`j'))" 
				} 
			}
			else {
			
				if "`noaverages'" != "noaverages" {
					
					local plotCmd "`plotCmd' (scatter Avg`j' trajT, mstyle(p`j'))"
				}	 
			}
		}
	}
	
	if "`nolegend'" == "nolegend" {
		local plotCmd "`plotCmd', ytitle(`ytitle') xtitle(`xtitle') yscale(titlegap(*5)) xscale(titlegap(*5)) legend(off) ylabel(`ylabel') xlabel(`xlabel') scale(1.2)"
	}
	else {
  
		if "`legendOpts'" != "" {
			local plotCmd "`plotCmd' ,legend(`legendOpts' region(lcolor(white)) symxsize(6) "
		}
		else {
			local plotCmd "`plotCmd' ,legend(region(lcolor(white)) style(zyx2) ring(1) cols(3) pos(6) symxsize(6) order(`"`orderTxt'"')"
		} 
	}

	if "`legendLabels'" != "" & "`nolegend'" == "" {	
		local plotCmd "`plotCmd' `legendLabels' )"
	}
	
	if "`legendLabels'" == "" & "`nolegend'" == "" {
	
		local plotCmd "`plotCmd' label(1 "`label1'") label(2 "`label2'") label(3 "`label3'") "
		local plotCmd "`plotCmd' label(4 "`label4'") label(5 "`label5'") label(6 "`label6'") "
		local plotCmd "`plotCmd' label(7 "`label7'") label(8 "`label8'") label(9 "`label9'") "
		local plotCmd "`plotCmd' label(10 "`label10'") label(11 "`label11'") label(12 "`label12'") "
		local plotCmd "`plotCmd' label(14 "`label14'") label(15 "`label15'") label(16 "`label16'") "
		local plotCmd "`plotCmd' label(17 "`label17'") label(18 "`label18'") label(19 "`label19'") )"
	}
	
	if "`nolegend'" != "nolegend" {
		local plotCmd "`plotCmd' ytitle(`ytitle') xtitle(`xtitle') ylabel(`ylabel') xlabel(`xlabel') yscale(titlegap(*5)) xscale(titlegap(*5)) scale(1.2)"
	}
    
	twoway `plotCmd'

	restore

end
/* end of trajplot.ado */
