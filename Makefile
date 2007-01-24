
deadatmo: dead.F90 main.F90 cpl_comm_mod.F90
# Note, the cd command has to be on the same line as the command I want executed.
#cp dead.F90 ~/ccsm3_0/scripts/deadtest/SourceMods/src.xatm
#cp main.F90 ~/ccsm3_0/scripts/deadtest/SourceMods/src.cpl
#cp cpl_comm_mod.F90 ~/ccsm3_0/models/csm_share/cpl/
	cd ~/ccsm3_0/scripts/deadtest/; ./deadtest.joxaren.build

docs:
	doxygen Doxyfile