IO += io_h5write_header.o\
	io_h5write_lists.o \
	io_h5read_header.o \
	io_h5read_lists.o \
	io_h5read_unknowns.o \
	io_h5write_unknowns.o \
	io_h5write_unknowns_sp.o \
	io_h5write_generic_int_arr.o \
	io_h5write_generic_real_arr.o \
	io_h5read_generic_int_arr.o \
	io_h5read_generic_real_arr.o \
	io_comp_decomp.o \
	io_h5_xfer_wrapper.o \
	io_h5_xfer.o \
	io_h5create_dataset.o \
	io_h5_type.o \
	io_h5_attribute.o \
	io_h5_type_matched_xfer.o \
	io_h5create_dataset_protons.o \
	io_h5create_raydset.o \
	io_h5write_protondata.o \
	io_h5write_raydata.o \
	io_h5_read_present_dims.o \
	io_h5_report_xfer_method.o \
	io_createDatasets.o

io_h5_attribute.o: constants.h Flash.h io_flash.h io_h5_type.h

# The following GNU make special prevents that apparent dependencies
# on the file hdf5.mod trigger unnecessary recompilation of files that
# refer to the HDF5 module. The file hdf5.mod does usually not
# actually exist in the object directory but refers to the HDF5 module
# that is part of the FORTRAN interface of HDF5 and often resides in a
# directory that should be pointed to with a '-I' FORTRAN compiler
# option.
.SECONDARY: hdf5.mod

