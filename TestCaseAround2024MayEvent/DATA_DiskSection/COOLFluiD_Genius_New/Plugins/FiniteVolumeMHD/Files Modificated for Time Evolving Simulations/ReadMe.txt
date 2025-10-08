In the "Plugins\FiniteVolumeMHD" folder, the following code files have been modified for time-evolving simulations.
1. SuperInletProjectionConstrained.hh
In Line 69, Haopeng Wang added "virtual void preProcess();"
In Line 83, Haopeng Wang comment "Common::CFMap<Framework::TopologicalRegionSet*, RealVector*> m_mapTrs2Twall;" 

2. SuperInletProjectionConstrained.cxx
It corresponding to the innerboundary conditions described in the time-evolving COCONUT paper.
Pay attention the difinition of  maxvA, currently it is defined as "CFreal maxvA = 3.0e6;", but it can also be changed to "CFreal maxvA = 3.0e6;" for solar minimum.

3. SuperInletProjection.hh
In line 65, Andrea added "virtual void preProcess();".

4. SuperInletProjection.cxx
Between Lines 343 and 354, Andrea defined the "virtual void preProcess();" function.

5. CMakelists.txt
Since "SuperInletProjectionConstrained.hh", "SuperInletProjectionConstrained.cxx", "SuperInletProjectionConstrainedST.hh", "SuperInletProjectionConstrainedST.cxx", "SuperInletProjectionConstrainedMichaela", and "SuperInletProjectionConstrainedMichaela" are newly developed files, you should add the following commands to "CMakelists.txt" file.
"
SuperInletProjectionConstrained.cxx
SuperInletProjectionConstrained.hh
SuperInletProjectionConstrainedST.cxx
SuperInletProjectionConstrainedST.hh
SuperInletProjectionConstrainedMichaela.cxx
SuperInletProjectionConstrainedMichaela.hh
".

6. Venktn3DStrict.cxx
Between Lines 185 and 186 Haopeng Wang comment out "psimin = min(psi, psimin);" and added "psimin = 1.0;".

7. MHDConsACASourceTerm.cxx 
Between Lines 329 and 494: Haopeng Wang changed the rediative cooling curve from Rosner1978 to chiant 9.
Lines 357-494 Haopeng Wang added the radiattive loss function derived from chiant 9. Correspondingly, Lines 329 and 353 are commented out.