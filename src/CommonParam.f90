!######################################################################################
module CommonParam 
!######################################################################################
	implicit none
!!!!!!!!! descriptors for input files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_coord = 11
	character(5), public, parameter :: name_coord = 'COORD' ! for coordinates
	integer(4), public, parameter :: n_bonds = 12
	character(5), public, parameter :: name_bonds = 'BONDS' ! list of bonds
	integer(4), public, parameter :: n_veloc = 13
	character(5), public, parameter :: name_veloc = 'VELOC' ! initial velocities
	integer(4), public, parameter :: n_contr = 14
	character(5), public, parameter :: name_contr = 'CONTR' ! control parameters
	integer(4), public, parameter :: n_field = 15
	character(5), public, parameter :: name_field = 'FIELD' ! "force-field"
	integer(4), public, parameter :: n_fixed = 16
	character(5), public, parameter :: name_fixed = 'FIXED' ! list of the fixed beads
	integer(4), public, parameter :: n_angls = 17
	character(5), public, parameter :: name_angls = 'ANGLS' ! list of the angles	
!!!!!!!!! descriptors fot output files !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer(4), public, parameter :: n_ent = 20             ! descriptor for ent-files
	integer(4), public, parameter :: n_track = 21
	character(5), public, parameter :: name_track = 'TRACK'   ! print coordinates
	integer(4), public, parameter :: n_coordf = 22
	character(6), public, parameter :: name_coordf = 'COORDF' ! final configuration
	integer(4), public, parameter :: n_velocf = 23
	character(6), public, parameter :: name_velocf = 'VELOCF' ! final velocities
	integer(4), public, parameter :: n_infor = 24
	character(5), public, parameter :: name_infor = 'INFOR'   ! report about the code execution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	character(6), public, parameter :: palette = 'NOCSPF'
	real(4), public, parameter :: theta = 3.0
	real(4), public, parameter :: gamma = 4.5
	real(4), public, parameter :: r_cut = 1.0
	real(4), public, parameter :: lambda = 0.65
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	real(4) :: UH, UL, dH
	logical :: wall_on, solvent_on 
	character(len=128) :: path
	integer(4) :: solvent, num_solvent
	
Contains


end module CommonParam



