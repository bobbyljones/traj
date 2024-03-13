*! multtrajplot Jan 2016  Bobby L Jones

program multtrajplot
  version 9
  syntax ,[ xlabel(str) xtitle(str) ///
	ylabel1(str) ylabel2(str) ylabel3(str) ///
	ylabel4(str) ylabel5(str) ylabel6(str) ///
	ytitle1(str) ytitle2(str) ytitle3(str) ///
	ytitle4(str) ytitle5(str) ytitle6(str) 	CI NOAVERAGES]   
  preserve

  tempname A1 A2 A3 A4 A5 A6
     
	forval j = 1/6 {
  	capture confirm scalar e(numGroups`j')         
  	if _rc == 0 { 
   	 	local nv "`j'" 
  	}   
	}
	
  local ng = e(numGroups1)
	
	forval j = 1/`nv' {
  	capture confirm matrix e(plot`j')         
  	if _rc != 0 { 
   	 display as error "plot data not found" 
   	 exit _rc 
  	}
		mat A`j' = e(plot`j')
  	svmat A`j', names(col)
		rename trajT	_`j'_T
		forval i = 1/`ng' {
			rename Est`i' _`j'_Est`i'
			rename Avg`i' _`j'_Avg`i'
			rename L95`i' _`j'_L95`i'
			rename U95`i' _`j'_U95`i'
		}
	}

  mat Piv = e(groupSize1)
	local titlel " "
	forval i = 1/`ng' {
		local pi`i' : display %4.1f round( el( Piv, 1, `i' ), .1 )
  	local titlel "`titlel' G`i' (`pi`i''%) "	
  }
  
  local graphs = ""
  
  forval j = 1/`nv' {		
		
		forval i = 1/`ng' {
			
			if `i' > 1 {
				local ytitle`j' : display ""
			}

			local graphs "`graphs' _plot_`j'_`i'" 
			local plotCmd "(line _`j'_Est`i' _`j'_T, lwidth(thick))"
		
			if "`ci'" == "ci" {
				if "`noaverages'" != "noaverages" {
					local plotCmd "`plotCmd' (scatter _`j'_Avg`i' _`j'_T, mstyle(p`i'))"
				} 
				local plotCmd "`plotCmd' (line _`j'_L95`i' _`j'_T, lpattern(shortdash) lcolor(gs10))" 
				local plotCmd "`plotCmd' (line _`j'_U95`i' _`j'_T, lpattern(shortdash) lcolor(gs10))"
			}
			else {
				if "`noaverages'" != "noaverages" {
					local plotCmd "`plotCmd' (scatter _`j'_Avg`i' _`j'_T, mstyle(p`i'))"
				} 
			}
			
			local plotCmd "`plotCmd' , ytitle(`ytitle`j'') xtitle(`xtitle') ylabel(`ylabel`j'')" 
			local plotCmd "`plotCmd' xlabel(`xlabel') yscale(titlegap(*5)) xscale(titlegap(*5))"
			local plotCmd "`plotCmd' legend(off) scale(1.2) name(_plot_`j'_`i', replace) nodraw"
			twoway `plotCmd'
		}
			
	}
	
  graph combine `graphs', cols(`ng') title( "`titlel'", span size(vsmall) )

  restore

end
/* end of multtrajplot.ado */
