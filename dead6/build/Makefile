
SHELL=/bin/sh

FC=mpif90
CURDIR=.
NETCDF_INCLUDE=/soft/apps/packages/netcdf-3.5.1_pgi4.0/include
NETCDF_LIB=/soft/apps/packages/netcdf-3.5.1_pgi4.0/lib
MCT_MPEU_INCLUDE=/home/steder/software-pgi/include
MCT_MPEU_LIB=/home/steder/software-pgi/lib

DEFS=-DMODELNUMBER=1

INCLUDE_PATHS=-I${CURDIR} -I${MCT_MPEU_INCLUDE} -I{NETCDF_INCLUDE}
MODULES_PATHS=-p $(MCT_MPEU_INCLUDE) 
LIBRARY_PATHS=-L${MCT_MPEU_LIB} -L${NETCDF_LIB}
LIBRARIES=-lmct -lmpeu -lmph -lnetcdf
# SOURCE FILES

TARGET=dead
SRCS_F90=$(wildcard *.F90)

OBJS_ALL=$(SRCS_F90:.F90=.o)

# MACHINE AND COMPILER FLAGS
#include ../Makefile.conf

# TARGETS

all: ${TARGET}.bin

${TARGET}.bin: ${OBJS_ALL}
	${FC} -o $@ ${LIBRARY_PATHS} ${MODULE_PATHS} ${LIBRARIES} ${OBJS_ALL}

# ADDITIONAL FLAGS SPECIFIC FOR MCT COMPILATION

# RULES

.SUFFIXES:
.SUFFIXES: .F90 .o

.F90.o:
	$(FC) -c ${DEFS} ${INCLUDE_PATHS} ${MODULE_PATHS} $*.F90

clean:
	${RM} *.o *.mod ${TARGET}.bin

# COMPILATION ORDER

cpl_bundle_mod.o: cpl_mct_mod.o cpl_comm_mod.o cpl_domain_mod.o cpl_kind_mod.o cpl_control_mod.o shr_sys_mod.o shr_mpi_mod.o
cpl_comm_mod.o: cpl_fields_mod.o mph.o

mph.o: mph.F
	${FC} -c ${DEFS} ${INCLUDE_PATHS} ${MODULE_PATHS} mph.F

cpl_const_mod.o: cpl_kind_mod.o shr_const_mod.o
cpl_contract_mod.o: shr_timer_mod.o cpl_kind_mod.o cpl_mct_mod.o cpl_comm_mod.o cpl_fields_mod.o cpl_bundle_mod.o cpl_domain_mod.o cpl_infobuf_mod.o cpl_control_mod.o
cpl_control_mod.o: shr_sys_mod.o shr_cal_mod.o shr_date_mod.o shr_alarm_mod.o shr_orb_mod.o cpl_kind_mod.o
cpl_domain_mod.o: shr_sys_mod.o cpl_kind_mod.o cpl_mct_mod.o cpl_comm_mod.o cpl_control_mod.o
cpl_fields_mod.o: cpl_mct_mod.o cpl_kind_mod.o
cpl_infobuf_mod.o: cpl_kind_mod.o cpl_fields_mod.o shr_timer_mod.o shr_sys_mod.o shr_mpi_mod.o
cpl_interface_mod.o: cpl_mct_mod.o cpl_comm_mod.o cpl_fields_mod.o cpl_bundle_mod.o cpl_domain_mod.o cpl_infobuf_mod.o cpl_contract_mod.o cpl_kind_mod.o cpl_control_mod.o shr_sys_mod.o shr_timer_mod.o shr_mpi_mod.o
cpl_iobin_mod.o: cpl_mct_mod.o cpl_comm_mod.o cpl_fields_mod.o cpl_bundle_mod.o cpl_domain_mod.o cpl_kind_mod.o cpl_control_mod.o shr_sys_mod.o shr_date_mod.o shr_mpi_mod.o
cpl_iocdf_mod.o:cpl_mct_mod.o cpl_comm_mod.o cpl_fields_mod.o cpl_bundle_mod.o cpl_domain_mod.o cpl_kind_mod.o cpl_const_mod.o cpl_control_mod.o shr_sys_mod.o shr_date_mod.o
cpl_kind_mod.o: shr_kind_mod.o
cpl_map_mod.o: cpl_mct_mod.o cpl_domain_mod.o cpl_bundle_mod.o cpl_comm_mod.o cpl_kind_mod.o cpl_control_mod.o shr_sys_mod.o shr_mpi_mod.o

cpl_mct_mod.o: shr_sys_mod.o shr_mpi_mod.o cpl_kind_mod.o cpl_const_mod.o #MCT_MODULES

MCT_MODULES:
	cp ${MCT_MPEU_INCLUDE}/*.mod .

data_mod.o: shr_kind_mod.o cpl_fields_mod.o cpl_contract_mod.o
dead.o: cpl_interface_mod.o cpl_fields_mod.o cpl_contract_mod.o cpl_control_mod.o shr_const_mod.o shr_timer_mod.o shr_sys_mod.o shr_msg_mod.o shr_kind_mod.o data_mod.o

shr_alarm_mod.o: shr_cal_mod.o shr_date_mod.o shr_sys_mod.o shr_kind_mod.o
shr_cal_mod.o: shr_kind_mod.o
shr_const_mod.o: shr_kind_mod.o
shr_date_mod.o: shr_cal_mod.o shr_sys_mod.o shr_kind_mod.o
shr_file_mod.o: shr_kind_mod.o shr_sys_mod.o
shr_mpi_mod.o: shr_kind_mod.o
shr_msg_mod.o: shr_kind_mod.o shr_sys_mod.o
shr_orb_mod.o: shr_kind_mod.o shr_sys_mod.o shr_const_mod.o
shr_sys_mod.o: shr_kind_mod.o shr_mpi_mod.o
shr_timer_mod.o: shr_kind_mod.o
shr_vmath_mod.o: shr_kind_mod.o










