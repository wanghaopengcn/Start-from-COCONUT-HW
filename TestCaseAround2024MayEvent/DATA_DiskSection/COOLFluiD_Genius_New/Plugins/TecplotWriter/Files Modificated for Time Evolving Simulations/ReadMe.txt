In the "Plugins\TecplotWriter" folder, the following code files have been modified for time-evolving simulations.
1. ParWriteSolutionBlock.cxx 
Line 546 Haopeng Wang transfer the time unit from code unit to physical hour
Line 546 CFreal timeDim = SubSystemStatusStack::getActive()->getCurrentTimeDim()*1447.2/3600.0;
In Line 503 Haopeng Wang tried to comment "cf_assert(tt.nodesOffset[iType].second == maxpos);" to avoid exit while saving data.

2. WriteSolutionBlockFV.cxx 
transfer code unites to physical unites in output files.
Line 248 CFreal solutiontime = subSysStatus->getCurrentTimeDim() > 0 ? subSysStatus->getCurrentTimeDim()*1447.2/3600.0 : subSysStatus->getNbIter();
Line 268 << ", AUXDATA PhysTime=\"" << subSysStatus->getCurrentTimeDim()*1447.2/3600.0 << "\""
Line 699 CFreal solutiontime = subSysStatus->getCurrentTimeDim() > 0 ? subSysStatus->getCurrentTimeDim()*1447.2/3600.0 : subSysStatus->getNbIter();