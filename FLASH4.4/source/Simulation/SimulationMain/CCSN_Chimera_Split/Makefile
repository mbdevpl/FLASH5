
Simulation += Simulation_data.o chimera_model.o model_interp.o quadrature.o hdf5_read_write.o

Simulation_init.o : Simulation_data.o chimera_model.o quadrature.o
Simulation_initBlock.o : Simulation_data.o chimera_model.o quadrature.o
chimera_model.o : Simulation_data.o model_interp.o quadrature.o hdf5_read_write.o
mesa_model.o : Simulation_data.o model_interp.o quadrature.o hdf5_read_write.o
