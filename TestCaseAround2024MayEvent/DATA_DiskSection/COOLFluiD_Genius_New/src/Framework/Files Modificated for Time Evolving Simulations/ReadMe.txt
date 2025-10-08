In the Framework folder, the following code files have been modified for time-evolving simulations.

1. NodalStatesExtrapolator.ci
You need to modify this file in the sections marked with the the following comments, as described in the the manu:
"case-dependent_1", "case-dependent_2", "case-dependent_3", "case-dependent_4", and "case-dependent_5".

2. NodalStateExtrapolator.hh  2024.05.14
Between Lines 224 and 226, Haopeng Wang  added the following declaration, which also work for linear interpolation.
void readSurfaceDataextrapolator(std::vector<SurfaceData*>& surfaces,
		       const std::string& fileName, const std::string& fileName2,
			   double Coefi1, double Coefi2);	

3. ConvergenceMethod.cxx
Between Lines 447 and 455, Haopeng Wang transferd the time from code unite to physical unite in hours and print it to the output dialog.
L450 modified 'out << subSysStatus->getCurrentTimeDim();' to 'out << subSysStatus->getCurrentTimeDim()*1447.2/3600.0;'
L454 modified 'out << subSysStatus->getDTDim();' to 'out << subSysStatus->getDTDim()*1447.2/3600.0;'.

4. PathAppender.cxx 
Between Lines 113 and 115, Haopeng Wang transfered code unite to physical unite in hours.
L113     CFreal physicaltime = time*1447.2/3600.0;
L114     physicaltime=int(physicaltime*100.0)/100.0;
L115     std::string time_string = Common::StringOps::to_str(physicaltime);

5. SubSystemStatus.cxx
Between Lines 222 and 227, Haopeng Wang defined the start time in physical hour and then transfer it to code unite. It is related to "SubSystemStatusStack::getActive()->getCurrentTimeDim()"
L222 CFreal restarttime = 0.0
L227 return (PhysicalModelStack::getActive()->getImplementor()->getRefTime())*m_currentTime+restarttime*1447.2/3600.0;

With help from Andrea, this currently can be done in the CFcase, therefore no longer need to modify this file.
