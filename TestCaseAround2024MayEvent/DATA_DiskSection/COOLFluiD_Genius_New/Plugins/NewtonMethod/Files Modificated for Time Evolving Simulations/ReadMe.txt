In the "Plugins\NewtonMethod" folder, the following code files have been modified for time-evolving simulations and enhanced numerical stability.
1. BDF2.cxx
In Lines 148 and 149, Haopeng Wang changed "(m_data->freezeJacobian() && k > 1) ? m_data->setDoComputeJacobFlag(false) : m_data->setDoComputeJacobFlag(true);" 
to (m_data->freezeJacobian() && (k > 1) && (k % 4 != 0)) ? m_data->setDoComputeJacobFlag(false) : m_data->setDoComputeJacobFlag(true);

2. "StdUpdateSolPP.cxx" and "StdUpdateSolPP.cxx" are newly developed files, thererfore in "CMakeLists.txt" file, you should include the following lines:
"
StdUpdateSolPP.hh
StdUpdateSolPP.cxx
".

3. UpdateSolMHD.cxx
There are also some modification in "UpdateSolMHD.cxx", given that we didn't need to use it, don't pay attention to the modifications here.