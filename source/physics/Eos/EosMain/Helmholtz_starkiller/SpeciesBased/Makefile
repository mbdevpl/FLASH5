#	Makefile for the Helmholtz equation of state, species based

Eos += bl_types.o bl_constants.o actual_network.o network.o eos_type.o actual_eos.o eos_initHelmholtz.o eos_helmholtz.o

Eos.o : eos_vecData.o  Eos_data.o
bl_constants.o : bl_types.o
actual_network.o : bl_types.o
network.o : bl_types.o bl_constants.o actual_network.o
eos_type.o : bl_types.o bl_constants.o network.o
actual_eos.o : bl_types.o bl_constants.o eos_type.o
eos_initHelmholtz.o : network.o actual_eos.o
eos_helmholtz.o : eos_type.o actual_eos.o
