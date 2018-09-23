# Makefile for source/physics/sourceTerms/Burn/nuclearBurn/XNet

Burn += bn_initNetwork.o bn_xnetData.o bn_xnetFinalize.o bn_xnetInit.o \
	xnet_abundances.o xnet_conditions.o xnet_constants.o xnet_controls.o \
	xnet_data.o xnet_eos.o xnet_evolve.o xnet_fd.o xnet_ffn.o xnet_flux.o \
	xnet_integrate.o xnet_integrate_bdf.o xnet_integrate_be.o xnet_jacobian_dense.o \
	xnet_match.o xnet_mpi.o xnet_nse.o xnet_output.o xnet_parallel.o xnet_preprocess.o \
	xnet_screening.o xnet_timers.o xnet_types.o xnet_util.o

Burn.o : bn_interface.o bn_xnetData.o Burn_data.o Driver_interface.o Eos_interface.o Grid_interface.o Hydro_interface.o Simulation_interface.o Timers_interface.o
bn_burner.o : Burn_dataEOS.o Burn_data.o bn_interface.o xnet_abundances.o xnet_conditions.o xnet_constants.o xnet_controls.o xnet_evolve.o xnet_timers.o
bn_initNetwork.o : Burn_data.o RuntimeParameters_interface.o bn_xnetData.o bn_interface.o
bn_finalizeNetwork.o : bn_interface.o
bn_xnetData.o : xnet_controls.o xnet_data.o xnet_eos.o
bn_xnetFinalize.o : Driver_interface.o bn_xnetData.o xnet_jacobian_dense.o xnet_timers.o
bn_xnetInit.o : Driver_data.o xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_eos.o xnet_flux.o xnet_integrate_bdf.o xnet_jacobian_dense.o xnet_match.o xnet_parallel.o xnet_preprocess.o xnet_screening.o xnet_util.o

xnet_controls.o: xnet_parallel.o xnet_types.o xnet_util.o
xnet_data.o: xnet_constants.o xnet_controls.o xnet_ffn.o xnet_parallel.o xnet_types.o xnet_util.o
xnet_ffn.o: xnet_constants.o xnet_controls.o xnet_types.o
xnet_flux.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_match.o xnet_timers.o xnet_types.o
xnet_match.o: xnet_controls.o xnet_data.o xnet_types.o
xnet_abundances.o: xnet_conditions.o xnet_controls.o xnet_data.o xnet_eos.o xnet_nse.o xnet_types.o xnet_util.o
xnet_conditions.o: xnet_controls.o xnet_types.o xnet_util.o
xnet_constants.o: xnet_types.o
xnet_evolve.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_integrate_bdf.o xnet_integrate_be.o xnet_output.o xnet_timers.o xnet_types.o xnet_util.o
xnet_integrate.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_eos.o xnet_ffn.o xnet_screening.o xnet_timers.o xnet_types.o
xnet_integrate_bdf.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_jacobian_dense.o xnet_timers.o xnet_types.o
xnet_integrate_be.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_integrate.o xnet_jacobian_dense.o xnet_timers.o xnet_types.o
xnet_nse.o: xnet_constants.o xnet_controls.o xnet_data.o xnet_eos.o xnet_types.o xnet_util.o
xnet_output.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_flux.o xnet_match.o xnet_timers.o xnet_types.o
xnet_preprocess.o: xnet_constants.o xnet_data.o xnet_types.o xnet_util.o
xnet_screening.o: xnet_abundances.o xnet_conditions.o xnet_constants.o xnet_controls.o xnet_data.o xnet_eos.o xnet_types.o
xnet_timers.o: xnet_types.o
xnet_util.o: Driver_interface.o xnet_types.o

xnet_fd.o: xnet_types.o
xnet_eos.o : Eos_interface.o Simulation_interface.o xnet_constants.o xnet_controls.o xnet_data.o xnet_fd.o xnet_types.o

xnet_jacobian_dense.o: xnet_abundances.o xnet_conditions.o xnet_controls.o xnet_data.o xnet_timers.o xnet_types.o

xnet_parallel.o: xnet_mpi.o xnet_types.o
