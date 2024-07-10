program main   
    implicit none
    ! Parameters and functions for main
    integer(8), parameter :: N_b = 2               !Number of bodies
    real(8), parameter :: pi = 3.1415926
    real(8), parameter :: G = 6.67259e-8        !gravitational constant (cm3/g/s2)
    real(8), parameter :: C = 2.998e10          !speed of light (cm/s)
    real(8), parameter :: C2 = C**2             !C^2
    real(8), parameter :: dC2 = 1/C**2          !divided by C^2
    real(8), parameter :: GdC = G/C             !G divided by C
    real(8), parameter :: Gdc2 = G/C2           !G divided by C^2
    character(len=128), parameter :: f_ilename = 'output_chks_.dat'
    !real(8), parameter :: M_shell = 1e-3*M0     !Mass of the shell

    ! Declaration of variables
    integer(8) :: ch_id,N_grav,Sch_type
    real(8), dimension(N_b,3) :: x,v,S
    real(8), dimension(N_b) :: M
    real(8) :: t_step,t_beg,t_end,t_wr,tol_err,T1,eps,t

    ! Some preparations
    eps = EPSILON(pi)
    t_beg = 0!0.          !simulation starts at
    t = t_beg    
    call CPU_TIME(T1)
        
    !ignition system starts!
    ch_id = 0.0      
    call init_read(Sch_type,N_grav,t_end,t_wr,t_step,tol_err,M,S,x,v)
    !delete existing output file
    open(unit=222, file = f_ilename,status = "unknown",position = "APPEND")        
    close(unit=222, status = "delete")
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Check if scheme type is correct
    if ((Sch_type.ne.1).and.(Sch_type.ne.2)) then
        write(*,*) "Wrong scheme type"
        ERROR STOP
    endif
     
    !Write initial stat
    call checkpoint_wr(x,v,t,t_step)
    call step(x,v,t,t_step)
    !Main cycle
    do while (t.lt.t_end)
        call step(x,v,t,t_step)
        call if_wr(x,v,t,t_step,ch_id,t_wr)
    enddo
    
    !Write final state
    call checkpoint_wr(x,v,t,t_step)
 
    contains  
    
    real(8) function norm_l2(v)
        implicit none
        !Euclidean norm of vector v
        !TEMP
        real(8), dimension(3) :: v
        !BODY
        norm_l2 = sqrt(sum(v**2))
        return
    endfunction norm_l2
    
    subroutine checkpoint_wr(x,v,t,t_step)  
        implicit none
        !Writes current state into checkpoint file
        
        !IN
        real(8), dimension(N_b,3), intent(in) :: x,v
        real(8), intent(in) :: t,t_step
        !OUT
        
        !INOUT
        
        !TEMP
        character(len=16) :: format_i ='(E15.7)'
        integer(8) :: ii,jj
        real(8) :: T_curr,T2,secs
        integer(8) :: hours,mins
        !BODY
        
        !checkpoint_ID
        open(unit=222, file = f_ilename,status = "unknown",position = "APPEND")        
        write (222,format_i,advance="no") t
        write (222,format_i,advance="no") t_step
        do ii = 1,N_b
            do jj = 1,3
                write (222,format_i,advance="no") x(ii,jj)    
            enddo
        enddo
        do ii = 1,N_b
            do jj = 1,3
                write (222,format_i,advance="no") v(ii,jj)    
            enddo
        enddo
        !write(222,*) 'current_time; ', 'current_time_step; ' ,'checkpoint_num'
        !write(222,*) t, t_step, ch_id
        !write(222,*) 'x_1(3);', ' v_1(3)'
        !write(222,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)') curr_state(1:6)
        !write(222,*) 'x_2(3);', ' v_2(3)'
        !write(222,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)') curr_state(7:12)
        !write(222,*) 'x_3(3);', ' v_3(3)'
        !write(222,'(E15.7,E15.7,E15.7,E15.7,E15.7,E15.7)') curr_state(13:18)
        if ((command_argument_count().eq.0).and.(t_end.ne.0.0)) then
            
            call CPU_TIME(T2)
            T_curr = T2-T1
            secs = mod(mod(T_curr,3600.0),60.0)
            mins = NINT((T_curr-secs)/60.0)
            hours = NINT((mins-mod(mins,60))/60.0)
            mins = mod(mins,60)
            
            !write(*,*) 100*REAL(ch_id,8)
            write(*,*) 'Progress:',100*REAL(ch_id,8)/REAL(CEILING(t_end/t_wr),8),'%',', current t_step = ',t_step
            write(*,'(A,I3,A4,I2,A5,F4.1,A3)') 'Elaped time:',hours,'(h);',mins,'(m); ',secs,'(s)'
            
        endif
        close(222)       

    endsubroutine checkpoint_wr
                       
    subroutine if_wr(x,v,t,t_step,ch_id,t_wr)
        implicit none
        !Check if timestep is a checkpoint
        !IN
        real(8), intent(in) :: t_wr,t,t_step
        real(8), dimension(N_b,3), intent(in) :: x,v
        !OUT
        
        !INOUT
        integer(8), intent(inout) ::  ch_id
        !TEMP
        integer(8) :: aux_ch_id
        !BODY
        
        aux_ch_id = int(t/t_wr,8)
        if (aux_ch_id .gt. ch_id) then
            ch_id = aux_ch_id
            call checkpoint_wr(x,v,t,t_step)
        endif      
        
    endsubroutine if_wr
    
    subroutine init_read(Sch_type,N_grav,t_end,t_wr,t_step,tol_err,M,S,x,v)
        implicit none
        !Read initial parameters from init.dat file
        !OUT
        integer(8), intent(out) :: N_grav,Sch_type
        real(8), intent(out) :: t_end,t_wr,t_step,tol_err
        real(8), dimension(N_b), intent(out) :: M
        real(8), dimension(N_b,3), intent(out) :: S,x,v
        !INOUT
    
        !TEMP
        integer(8) :: ipos,ios,i,j
        integer(8), parameter :: nlen = 100
        real(8) :: aux
        character (len=nlen) :: str
        !BODY
        open(unit=10,file='init.dat',status='old')
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) Sch_type
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) N_grav
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) t_end
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) t_wr
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) t_step
        
        read(10, '(A)', iostat=ios) str
        ipos = scan(str,"=",back=.true.)
        read (str(1+ipos:),*) tol_err
        
        do i = 1,N_b
            read(10, '(A)', iostat=ios) str
            ipos = scan(str,"=",back=.true.)
            read (str(1+ipos:),*) M(i)
        enddo 
        
        do i = 1,N_b
            do j = 1,3
                read(10, '(A)', iostat=ios) str
                ipos = scan(str,"=",back=.true.)
                read (str(1+ipos:),*) S(i,j)
            enddo
        enddo 
        
        do i = 1,N_b
            do j = 1,3
                read(10, '(A)', iostat=ios) str
                ipos = scan(str,"=",back=.true.)
                read (str(1+ipos:),*) x(i,j)
            enddo
            do j = 1,3
                read(10, '(A)', iostat=ios) str
                ipos = scan(str,"=",back=.true.)
                read (str(1+ipos:),*) v(i,j)
            enddo
        enddo 
        
        close(10)
    endsubroutine init_read
        
    subroutine step(x,v,t,h)
        implicit none
        !Step of the scheme
        !IN
        
        !OUT
        
        !INOUT
        real(8), dimension(N_b,3), intent(inout) :: x,v
        real(8), intent(inout) :: h,t
        !TEMP
        real(8), dimension(N_b,3) :: x1,v1,x2,v2,x3,v3,x4,v4,x5,v5,x6,v6
        real(8) :: k1_err,k2_err,k3_err,k4_err,k5_err,k6_err
        real(8) :: err,h_prev
        !BODY
        err = tol_err+1.0
        do while ((err.gt.tol_err) .and. (h.gt.eps))
            call r_h(x1,v1,t,x,v,k1_err)
            x1 = x1*h
            v1 = v1*h
            call r_h(x2,v2,t+h*2/9,x + x1*2/9,v + v1*2/9,k2_err)
            x2 = x2*h
            v2 = v2*h
            call r_h(x3,v3,t+h*1/3,x + x1*1/12 + x2/4,v + v1*1/12 + v2/4,k3_err)
            x3 = x3*h
            v3 = v3*h
            call r_h(x4,v4,t+h*3/4,x + x1*69/128 - x2*243/128 + x3*270/128,v + v1*69/128 - v2*243/128 + v3*270/128,k4_err)
            x4 = x4*h
            v4 = v4*h
            call r_h(x5,v5,t+h,x - x1*17/12 + x2*27/4 - x3*27/5 + x4*16/15,v - v1*17/12 + v2*27/4 - v3*27/5 + v4*16/15,k5_err)
            x5 = x5*h
            v5 = v5*h
            call r_h(x6,v6,t+h*5/6,x + x1*65/432 - x2*5/16 + x3*13/16 + x4*4/27 + x5*5/144,v + v1*65/432 - v2*5/16 + v3*13/16 &
				&+ v4*4/27 + v5*5/144,k6_err)
            x6 = x6*h
            v6 = v6*h
            !k6_err = h*(-k1_err*1/150 + 0 + k3_err*3/100 - k4_err*16/75 - k5_err*1/20 + k6_err*6/25)
            !err = k2_err(1)
            err = h*abs(-k1_err*1/150 + 0 + k3_err*3/100 - k4_err*16/75 - k5_err*1/20 + k6_err*6/25)+eps
            !write(*,*) err
            h_prev = h
            h = 0.9*h*(tol_err/err)**(0.2)
        enddo
        if (h_prev.le.eps) then
            !write(*,*) 'h is too small, increase tol_err!'
            stop 'h is too small, increase tol_err!' 
            h = 4*h_prev/3
        endif
        x = x + (x1*47/450 + 0 + x3*12/25 + x4*32/225 + x5*1/30 + x6*6/25)
        v = v + (v1*47/450 + 0 + v3*12/25 + v4*32/225 + v5*1/30 + v6*6/25)
        !write(*,*) state(4:6)
        t = t + h_prev
        !write(*,*) 'Progress:',100*t/t_end,'%',', current t_step = ',t_step
        !write(*,*) t
    endsubroutine step
    
    subroutine r_h(x_out,v_out,t,x,v,err)
        implicit none
        !Righ hand function of system of ODEs
        !rh_func_rez - rezult vector 
        !input vector
        !t - time
        !#d\vector{y}\dt = f(\vector{y},t)  
        !IN
        real(8), dimension(N_b,3), intent(in) :: x,v
        real(8), intent(in) :: t
        !OUT
        real(8), dimension(N_b,3), intent(out) :: x_out,v_out
        real(8), intent(out) :: err
        !INOUT
        
        !TEMP
        real(8) :: aux1,aux2,r_da,r_ak,r_bk,r_ij,u_aux
        real(8), dimension(3) :: n_da,n_ak,n_ka,n_bk,aux3,v_ak,v_ka
        real(8), dimension(3) :: a_eih_k,a_so_k,a_ss_k
        real(8), dimension(3) :: H_i,g_i,v_ij,n_ij,g_t
        real(8), dimension(N_b) :: v2,err_a
        real(8), dimension(3) :: R_i,U1n_i,U2n_i,Udn_i,Un_i
        real(8), dimension(N_b,3) :: an
        integer(8) :: i,j,k,a,b,d
        
        !BODY
        if (Sch_type.eq.1) then
              !set v**2 and newtone acceleration an
            do i = 1,N_b
                !v_squared
                v2(i) = sum(v(i,:)*v(i,:))
                an(i,:) = 0
                do j = 1,N_grav
                    if (i.ne.j) then
                        an(i,:) = an(i,:) - G*M(j)*(x(i,:)-x(j,:))/(norm_l2(x(i,:)-x(j,:))**3)
                    endif
                enddo
            enddo
        
            do i = 1,N_b
                !Gravimagnetic at i
                H_i = 0
            
                !Gravielectric at i
                g_i = 0
            
                !(dg at i)/dt
                g_t = 0
            
                !Newtone potential at i
                Un_i = 0
            
                !d(previous)/dt
                Udn_i = 0
            
                !Right hand side of EoM
                R_i = 0
            
                do j = 1,N_grav
                    if (i.ne.j) then
                        !r_ij and n_ji = -n_ij
                        v_ij = v(i,:) - v(j,:)
                        r_ij = norm_l2(x(i,:)-x(j,:))
                        n_ij = (x(i,:)-x(j,:))/r_ij
                    
                        !Newton potential
                        Un_i = Un_i + M(j)/r_ij
                        Udn_i = Udn_i - M(j)*sum(n_ij*v_ij)/(r_ij**2)
                    
                        !Gravimagnetic field
                        H_i = H_i + ( 4*M(j)*CP(n_ij,v(j,:))/(r_ij**2) + 2*( S(j,:) - 3*n_ij*sum(S(j,:)*n_ij) )/(r_ij**3) )*GdC
                    
                        !Gravielectric field
                        !u_aux non-convinient part of g: u_aux = sum_k.ne.j M(k)/r_kj
                        u_aux = 0
                        do k = 1,N_grav
                            if (j.ne.k) then
                                u_aux = u_aux + M(k)/(norm_l2(x(j,:)-x(k,:)))
                            endif
                        enddo
					
                        g_i = g_i - n_ij*M(j)*Gdc2/(r_ij**2)*( c2 + 2*v2(j) - 1.5*(sum(v(j,:)*n_ij))**2 -&
                           & u_aux*G - 0.5*r_ij*sum(an(j,:)*n_ij) + 6*sum(v(j,:)*CP(S(j,:),n_ij))/(M(j)*r_ij) ) +&
							    & ( (3*sum(n_ij*v(j,:))/(r_ij**2))*( M(j)*v(j,:) + 2*CP(S(j,:),n_ij)/(r_ij) ) +&
									& 3.5*M(j)*an(j,:)/r_ij - 4*CP(S(j,:),v(j,:))/(r_ij**3) )*Gdc2
                    
                        !dg/dt for part: -[S_i x dg_i/dt]
                        g_t = g_t + M(j)*( v_ij - 3*n_ij*sum(v_ij*n_ij) )/(r_ij**3)
                    
                        !spin-spin: grad of H_i's spin part: 0.5*grad(S_i*H^s_ji):
                        R_i = R_i + ( ( 15*n_ij*sum(n_ij*S(j,:))*sum(n_ij*S(i,:)) - 3*n_ij*sum(S(i,:)*S(j,:)) -&
						    3*S(i,:)*sum(n_ij*S(j,:)) - 3*S(j,:)*sum(n_ij*S(i,:)) )/(M(i)*r_ij) +&
							    !spin-orb: gradient of non-spin H_i's part: 0.5*grad(S_i*4*(M_j/r_ij**2)*[n_ij x v_j]):
							    ( -2*CP(S(i,:),v(j,:)) - 6*n_ij*sum(v(j,:)*CP(S(i,:),n_ij)) +&
                                     !spin-orb: grad(s.[v x g]):
								     CP(S(i,:),v(i,:))*2 + 6*n_ij*sum(v(i,:)*CP(S(i,:),n_ij)) )*M(j)/M(i) )*Gdc2/(r_ij**3)
                    endif
                enddo
                !Mult parts:
                U2n_i = 1 - Un_i*Gdc2 + 1.5*v2(i)*dc2
                U1n_i = 1 + 3*Un_i*Gdc2 + 0.5*v2(i)*dc2
            
                !Old parts of equation
                !-[S_i x dg_i/dt]
			    ![s_i x dg_i/dt]:
                !Final R_i:
                R_i = R_i + U2n_i*g_i + CP(v(i,:)/c,H_i) + CP(S(i,:)/M(i),g_t*Gdc2) - 3*Udn_i*v(i,:)*Gdc2
            
                !a_i*u1ni + (a_i.v_i)*v_i/c^2 = R_i ->
                !(a_i*v_i) = (R_i.v_i)*c^2 / (u1nic^2 + v_i^2)
                
                !out x = v
                x_out(i,:) = v(i,:)
                err_a(i) = norm_l2(x_out(i,:))
                !Acceleration of i-th mass with resepct to (a_i.v_i)
                !a(i,:) = ( R_i - v(i,:)*sum(R_i*v(i,:))/(U1n_i*c2+v2(i)) )/U1n_i
                v_out(i,:) = ( R_i - v(i,:)*sum(R_i*v(i,:))/(U1n_i*c2+v2(i)) )/U1n_i
                !v_out(i,:) = v_out(i,:)-an(i,:)
                !a(i,:) = ( R_i - v(i,:)*sum(R_i*v(i,:))/(U1n_i*c2+v2(i)) )/U1n_i               
            enddo    
            
        else
            do k = 1,N_b
                a_eih_k = 0
                a_so_k = 0
                a_ss_k = 0
                do a = 1,N_grav
                    if (a.ne.k) then                  
                        r_ak = norm_l2(x(a,:)-x(k,:))
                        n_ak = (x(a,:)-x(k,:))/r_ak
                        n_ka = -n_ak
                        v_ak = v(a,:) - v(k,:)
                        v_ka = -v_ak
                        
                        do b = 1,N_grav
                            aux1 = 0
                            if (b.ne.k) then
                                r_bk = norm_l2(x(b,:)-x(k,:))
                                aux1 = aux1 + M(b)/r_bk
                            endif
                        enddo
                        aux1 = aux1*Gdc2
                        
                        do d = 1,N_grav
                            aux2 = 0
                            aux3 = 0
                            if (d.ne.a) then
                                r_da = norm_l2(x(d,:)-x(a,:))
                                n_da = (x(d,:)-x(a,:))/r_da
                                aux2 = aux2 + (M(d)/r_da)*(1 - (r_ak/(2*r_da))*(sum(n_ak*n_da)))
                                aux3 = aux3 + 3.5*n_da*M(a)*M(d)/(r_ak*r_da**2)
                            endif
                        enddo
                        aux2 = aux2*Gdc2
                        aux3 = aux3*(G**2)*dc2
                        
                        a_eih_k = a_eih_k + (M(a)*G*n_ak/r_ak**2)*(1 - 4*aux1 - aux2 + (sum(v(k,:)*v(k,:)) + 2*sum(v(a,:)*v(a,:)) &
                            - 4*sum(v(a,:)*v(k,:)) - 1.5*sum(v(a,:)*n_ak)**2)*dc2) - (v_ak*M(a)/(r_ak**2)) &
                                *sum(n_ak*(3*v(a,:)-4*v(k,:)))*Gdc2 + 3.5*aux3
                        
                        a_so_k = a_so_k + (M(a)/r_ak**3)*(6*n_ka*(sum(CP(S(a,:),n_ka)*v_ka)) + 4*CP(S(a,:),v_ka) &
                            - 6*CP(S(a,:),n_ka)*sum(v_ka*n_ka) + 6*n_ka*(sum(CP(S(k,:),n_ka)*v_ka)) + 3*CP(S(k,:),v_ka) &
                                - 3*CP(S(k,:),n_ka)*sum(v_ka*n_ka))*Gdc2
                        
                        a_ss_k = a_ss_k + (M(a)/r_ak**4)*(15*n_ak*sum(n_ka*S(k,:))*sum(n_ka*S(a,:)) - 3*n_ka*sum(S(k,:)*S(a,:)) &
                            - 3*S(k,:)*sum(n_ka*S(a,:)) - 3*S(a,:)*sum(n_ka*S(k,:)))*Gdc2
                    endif        
                enddo
                x_out(k,:) = v(k,:)
                err_a(k) = norm_l2(x_out(k,:))
                v_out(k,:) = a_eih_k + a_so_k + a_ss_k
            enddo          
        endif
        err = MAXVAL(err_a)
    endsubroutine r_h
        
    subroutine check_clock(T)
        implicit none
        !Cross croduct of v_1 and v_2
        !IN
        real(8), intent(inout) :: T
        !OUT

        !INOUT
    
        !TEMP
        integer(8) :: h,m,s
        real(8) :: p
        !BODY
        p = 60.0
        h = int(T/p**2,8)
        T = T-real(h,8)*p**2
        m = int(T/p,8)
        s = int(T-real(m,8)*p,8)
        write(*,*) 'Elaped',h,'(h)',m,'(m)',s,'(s)'
    end subroutine check_clock
    
  	function CP(v_1, v_2) result (v_out)
		implicit none
        !Cross croduct of v_1 and v_2
		real(8), dimension(3), intent(in) :: v_1,v_2
		real(8), dimension(3) :: v_out
		!body
		v_out(1) = v_1(2)*v_2(3)-v_1(3)*v_2(2)
        v_out(2) = v_1(3)*v_2(1)-v_1(1)*v_2(3)
        v_out(3) = v_1(1)*v_2(2)-v_1(2)*v_2(1)
    end function

endprogram main
