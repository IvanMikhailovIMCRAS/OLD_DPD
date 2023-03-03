!######################################################################################
module ForceEnergy
!######################################################################################
    use ErrorList
    use Lib
    use CommonParam
    use BoxMuller
    implicit none
		
Contains

!***************************************************************************************
subroutine calc_acceleration(N, NB, x, y, z, vx, vy, vz, BOX, FF, l_bond, k_bond,k_ang,&
&    b1,b2,N_dblist, ub1, ub2, ax, ay, az, typ, num_types, sqrdt,wall,num_ang, list_ang)
!***************************************************************************************
	integer(4), intent(in) :: N, NB, N_dblist, num_types, num_ang
	integer(4), intent(in) :: typ(1:N), list_ang(1:3,1:num_ang)
	real(4), intent(in) :: FF(1:num_types,1:num_types)
	integer(4), intent(in) :: b1(1:NB), b2(1:NB), ub1(0:N_dblist), ub2(1:N_dblist)
	real(4), intent(in) :: x(1:N), y(1:N), z(1:N), vx(1:N), vy(1:N), vz(1:N)
	real(4), intent(in) :: sqrdt, wall
	real(4), intent(in) :: l_bond, k_bond, k_ang, BOX(3)
	real(4), intent(out) :: ax(1:N), ay(1:N), az(1:N)
	integer(4) i
	real(4) dx, dy, dz, dvx, dvy, dvz, r, force, force_C, force_R, force_D, omega, ksi
	real(4) lower_border, upper_border
	real(4) ex, zt
	real(4) f_ang(1:3,1:3)

	ax(:) = 0.0; ay(:) = 0.0; az(:) = 0.0
	
          ! interactions for unbonded beads
	do i = 1, ub1(0)
		dx = x(ub1(i))-x(ub2(i)); dy = y(ub1(i))-y(ub2(i)); dz = z(ub1(i))-z(ub2(i))
		call periodic_test(dx,BOX(1));	call periodic_test(dy,BOX(2));	call periodic_test(dz,BOX(3))
		r = sqrt(dx**2 + dy**2 + dz**2)
		if (r.gt.r_cut) then
			force = 0.0
		else
			omega = 1.0 - r/r_cut
			call GaussRand(ksi)
			dvx = vx(ub1(i))-vx(ub2(i)); dvy = vy(ub1(i))-vy(ub2(i)); dvz = vz(ub1(i))-vz(ub2(i))	
											
			force_C = FF(typ(ub1(i)),typ(ub2(i))) * omega / r
			force_R = theta * omega * ksi / sqrdt / r
			force_D = - gamma * omega * omega * (dx*dvx + dy*dvy + dz*dvz) / r / r	
					
			force = force_C + force_R + force_D
			force = force 
		endif
		ax(ub1(i)) = ax(ub1(i)) + force * dx; ax(ub2(i)) = ax(ub2(i)) - force * dx 
		ay(ub1(i)) = ay(ub1(i)) + force * dy; ay(ub2(i)) = ay(ub2(i)) - force * dy
		az(ub1(i)) = az(ub1(i)) + force * dz; az(ub2(i)) = az(ub2(i)) - force * dz
	enddo

	! interactions for bonded beads
	
	do i = 1, NB
		dx = x(b1(i))-x(b2(i)); dy = y(b1(i))-y(b2(i)); dz = z(b1(i))-z(b2(i))
		call periodic_test(dx,BOX(1));	call periodic_test(dy,BOX(2));	call periodic_test(dz,BOX(3))
		r = sqrt(dx**2 + dy**2 + dz**2)
		force_C = k_bond * (l_bond/r - 1.0) 
				
		
		if (r.gt.r_cut) then
			force = force_C
		else
			omega = 1.0 - r/r_cut
			call GaussRand(ksi)
			dvx = vx(b1(i))-vx(b2(i)); dvy = vy(b1(i))-vy(b2(i)); dvz = vz(b1(i))-vz(b2(i))	
						
			force_R = theta * omega * ksi / sqrdt / r
			force_D = - gamma * omega * omega * (dx*dvx + dy*dvy + dz*dvz) / r / r	
					
			force = force_C + force_R + force_D
		endif
		ax(b1(i)) = ax(b1(i)) + force * dx; ax(b2(i)) = ax(b2(i)) - force * dx 
		ay(b1(i)) = ay(b1(i)) + force * dy; ay(b2(i)) = ay(b2(i)) - force * dy
		az(b1(i)) = az(b1(i)) + force * dz; az(b2(i)) = az(b2(i)) - force * dz
	enddo
	
	! Interaction between wall and beads
	if (wall_on) then
		lower_border = - BOX(3)*0.5 - l_bond*0.5 + r_cut
		upper_border = + BOX(3)*0.5 + l_bond*0.5 - r_cut
		do i = 1, N
				if (z(i).lt.lower_border) az(i) = az(i) + (lower_border-z(i))*wall
				if (z(i).gt.upper_border) az(i) = az(i) + (upper_border-z(i))*wall
				if (interact_on) then
					if (typ(i).ne.solvent) then
						zt = (z(i) + 0.5*BOX(3) - UL)/dH
						if (abs(zt).lt.10.0) then
						ex = exp(zt)
						az(i) = az(i) - UH * ex / dH / (1.0 + ex)**2
						endif
					endif
				endif
		enddo
	endif
	
	! Forces keeping angles in molecules
	do i = 1, num_ang
		call angular_projections(x(list_ang(1,i)), x(list_ang(2,i)), x(list_ang(3,i)), &
&                                y(list_ang(1,i)), y(list_ang(2,i)), y(list_ang(3,i)), &
&                                z(list_ang(1,i)), z(list_ang(2,i)), z(list_ang(3,i)), &
&                                f_ang(1,1),f_ang(1,2),f_ang(1,3), &
&                                f_ang(2,1),f_ang(2,2),f_ang(2,3), &
&                                f_ang(3,1),f_ang(3,2),f_ang(3,3) ) 
		ax(list_ang(1,i)) = ax(list_ang(1,i)) + k_ang * f_ang(1,1) 
		ax(list_ang(2,i)) = ax(list_ang(2,i)) + k_ang * f_ang(1,2)
		ax(list_ang(3,i)) = ax(list_ang(3,i)) + k_ang * f_ang(1,3)
		ay(list_ang(1,i)) = ay(list_ang(1,i)) + k_ang * f_ang(2,1)
		ay(list_ang(2,i)) = ay(list_ang(2,i)) + k_ang * f_ang(2,2)
		ay(list_ang(3,i)) = ay(list_ang(3,i)) + k_ang * f_ang(2,3)
		az(list_ang(1,i)) = az(list_ang(1,i)) + k_ang * f_ang(3,1)
		az(list_ang(2,i)) = az(list_ang(2,i)) + k_ang * f_ang(3,2)
		az(list_ang(3,i)) = az(list_ang(3,i)) + k_ang * f_ang(3,3)
	enddo
	
return
end subroutine calc_acceleration
!***************************************************************************************
    
!######################################################################################
end module ForceEnergy
!######################################################################################
