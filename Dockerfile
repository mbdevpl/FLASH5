FROM mbdevpl/usable-ubuntu:18.04

MAINTAINER Mateusz Bysiek <mateusz.bysiek.spam@gmail.com>

USER user
WORKDIR /home/user

RUN mkdir -p Projects && \
  git clone https://github.com/mbdevpl/spack.git --branch amrex_no_cmake_yes_configure --depth 1 Spack && \
  echo "source /home/user/Spack/share/spack/setup-env.sh" >> /home/user/.profile

RUN spack install mpich@3.2.1 && \
  spack install --no-checksum hdf5@1.8.20 +cxx +fortran +hl +mpi +szip +threadsafe ^mpich && \
  spack install --no-checksum hypre@2.14.0 +mpi ^mpich ^openblas@0.3.2 threads=openmp && \
  spack install amrex@develop dimensions=2 ~openmp +fortran +particles +mpi ^mpich && \
  spack install superlu@5.2.1 ^openblas threads=openmp && \
  echo "spack load -r mpich && spack load -r hdf5 && spack load -r hypre && spack load -r amrex && spack load -r superlu" >> /home/user/.bash_history && \
  echo "./setup Sod -auto -2d +Mode1 -site spack" >> /home/user/.bash_history && \
  echo "./setup Sod -auto -2d +Mode3 -site spack" >> /home/user/.bash_history && \
  echo "cd object" >> /home/user/.bash_history && \
  echo "make" >> /home/user/.bash_history && \
  echo "./flash4" >> /home/user/.bash_history && \
  echo "mpirun -np 1 ./flash4" >> /home/user/.bash_history && \
  echo "mpirun -np 2 ./flash4" >> /home/user/.bash_history

#spack install hpctoolkit@master +mpi ^mpich
