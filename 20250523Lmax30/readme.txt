This folder containe the modifications made to source codes of the COCONUT version updated around July 2026 and a corresponding .CFcase.
For simulation containing transiation region, please 
1. Set Simulator.SubSystem.Flow.Data.MHDConsACAHWST.TransitionRegion = 1 in the .CFcase file. 
2. comment out "smooth_factor = 1.0;" in MHD3DProjectionDiffPrimHW.cxx.
3. Modify the corresponding boundary conditions and initial distribution of solutions in the .CFcase file.
4. Keep in mind simulation with transiation region require high resolution grid in bottom region of the grid.
