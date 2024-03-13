/* 
		=================================================================================================
		Random Start Values for traj	03apr2018
		================================================================================================= 

		Copyright (C) 2018 Bobby L Jones <bobbyljones@gmail.com>
		
		This source code is subject to the terms of the 3-Clause BSD License (Open Source Initiative).
		The license can be obtained from https://opensource.org/licenses/BSD-3-Clause.
*/

program define trajstart, rclass

	version 9

	syntax [, Start(string) Sigma(real .1) Groups(integer 0) Risk]

	if "`start'" == "" {

		di as err "Start must specify a matrix"

		exit 198
	}

	if "`start'" != "" {

		cap confirm matrix `start'

		if _rc {
    
			di as err "Start must specify a matrix"
    
			exit 198
		}

		local traj_start_len = colsof(`start')
	}

	if (`groups' == 0 & "`risk'" != "risk") {

		di as err "Specify number of groups."

		exit 198
	}

	matrix traj_strt = J(1, `traj_start_len', 0)

	forval i = 1/`traj_start_len'{

		mat traj_strt[1,`i'] = `start'[1,`i'] + rnormal(`sigma') * `sigma'

	}

	local os = `traj_start_len' - `groups' + 1

	local denom = 0
 
	if "`risk'" == "risk" {
		/* nothing yet */
	} 
	else {
  
		forval i = `os'/`traj_start_len'{  
   
			local denom = `denom' + traj_strt[1,`i']
   
		} 

		forval i = `os'/`traj_start_len'{  
		
			if (el(traj_strt,1,`i') <= 0) {
				
				di as err "Group size percentages in start values must be greater than zero."

				exit 198
			}
			
			mat traj_strt[1,`i'] = traj_strt[1,`i']/`denom'*100
		}    
	}

	return matrix trajstart = traj_strt

end

/* end of trajstart.ado */
