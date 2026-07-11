This folder containe the modifications made to source codes of the COCONUT version updated around July 2026 and a corresponding .CFcase.
For simulation containing transiation region, please 
1. Set Simulator.SubSystem.Flow.Data.MHDConsACAHWST.TransitionRegion = 1 in the .CFcase file. 
2. comment out "smooth_factor = 1.0;" in the following snippet in MHD3DProjectionDiffPrimHW.cxx.
"   CFreal qSound = -1.5*rho_temp * rho0*std::pow(C_sound, 3.0)*(qF_x*normal[XX] + qF_y*normal[YY] + qF_z*normal[ZZ]);
   CFreal smooth_factor = max(0.0, min(1.0 - 0.5*std::abs(qFlux_1 / qSound), 1.0));
   smooth_factor = 1.0;
   qFlux_1 = qFlux_1*smooth_factor + qSound*(1.0 - smooth_factor);
   // comment it out to cancel TRACE
   qFlux_1 = qFlux_1*std::pow((1 + std::pow(T / T_threshold, -10)), 0.25);
   // ------changed on 2026.06.06 << ----------
   qFlux_2=3.0/2.0*alpha_cond*state[0]*rho0/mu/mH*kB*T0_new*(Vx*normal[XX]+Vy*normal[YY]+Vz*normal[ZZ]);
   //------modified>> on 2026.07.11  ----------
   smooth_factor = max(0.0, min(1.0 - 0.5*std::abs(qFlux_2 / qSound), 1.0));
   smooth_factor = 1.0;
   qFlux_2 = qFlux_2*smooth_factor + qSound*(1.0 - smooth_factor);"
3. Modify the corresponding boundary conditions and initial distribution of solutions in the .CFcase file.
4. Keep in mind simulation with transiation region require high resolution grid in bottom region of the grid.
