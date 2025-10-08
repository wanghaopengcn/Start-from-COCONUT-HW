This test case is part of the simulations described in the manuscript "H. P. Wang, S. Poedts, A. Lani, L. Linan, T. Baratashvili, F. Zhang, D. Sorokina, J. H.-J., Y. C. Li, N.-Z. Mahdi, and B. Schmieder. Time-evolving coronal modelling of solar maximum around the May 2024 storm by COCONUT. 2025.[arXiv:2505.11990]".

It reads 121 hourly-updated magnetograms at the start of the time-evolving coronal simulation, and is configured to run on the Level-5 mesh, the coarsest resolution, and simulates one day of coronal evolution in physical time. You may adjust the simulation duration by modifying the parameter "Simulator.SubSystem.MaxTime.maxTime" in the "map_unsteadyeclipse_fullMHD.CFcase" file, provided it remains less than 120.0*3600.0/1447.2.

Alternatively, Level-6 or Level-7 meshes can be selected by editing the CFcase files. We recommend using the Level-6 mesh, as it provides a well-balanced compromise between computational efficiency, accuracy, memory usage, and numerical stability. For Level-6, it is advisable to use 900 CPU cores instead of the 270 configured for Level-5 mesh.

For any questions or issues, please feel free to contact me at Haopeng.wang1@kuleuven.be.


