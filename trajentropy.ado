/* 
		===========================================================================
		traj entropy	
		=========================================================================== 

		Copyright (C) 2022 Bobby L Jones <bobbyljones@gmail.com>
		
		This source code is subject to the terms of the 3-Clause
		BSD License (Open Source Initiative). The license can be
		obtained from https://opensource.org/licenses/BSD-3-Clause.
*/

program trajentropy

 quietly gen _ent_row = log(_traj_ProbG1) * _traj_ProbG1 if _traj_ProbG1 > 0

 forval j = 2 / `e(numGroups1)' {
 
  quietly replace _ent_row = _ent_row + log(_traj_ProbG`j') * _traj_ProbG`j' if _traj_ProbG`j' > 0
 
 }

 quietly sum _ent_row

 di as text _newline(1) " Entropy = " %5.3f 1 + `r(mean)' / log( `e(numGroups1)' )

 drop _ent_row

end
