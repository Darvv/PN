!Computational subroutine
module MySubs

integer(8), parameter :: N_b = 3               !Number of bodies
integer(8), parameter :: N_grav = 3               !Number of bodies
integer(8), parameter :: size = 3               !Size of the euclidean space
                                                    !only for matlab interlayer
real(8), parameter :: pi = 3.1415926
real(8), parameter :: G = 6.67259e-8        !gravitational constant (cm3/g/s2)
real(8), parameter :: C = 2.998e10          !speed of light (cm/s)
real(8), parameter :: C2 = C**2             !C^2
real(8), parameter :: dC2 = 1/C**2          !divided by C^2
real(8), parameter :: GdC = G/C             !G divided by C
real(8), parameter :: Gdc2 = G/C2           !G divided by C^2

contains

subroutine r_h(a_sm,x,v,M,S)
    implicit none
    !Righ hand function of system of ODEs
    !rh_func_rez - rezult vector 
    !input vector
    !t - time
    !#d\vector{y}\dt = f(\vector{y},t)  
    !IN
    real(8), dimension(N_b,3), intent(in) :: x,v,S
    real(8), dimension(N_b), intent(in) :: M
    !OUT
    real(8), dimension(N_b,3), intent(out) :: a_sm !small correction to a_n
    !INOUT

    !TEMP
    real(8), dimension(3) :: H_i,g_i,v_ij,r_ij,r_kj,n_ij
    real(8), dimension(N_b) :: v2,u_aux,g_t,err_a
    real(8), dimension(3) :: R_i,U1n_i,U2n_i,Udn_i,Un_i
    real(8), dimension(N_b,3) :: a,an
    integer(8) :: i,j,k

    !BODY
    !set v**2 and newtone acceleration an
    do i = 1,N_b
        !out x = v
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
                    u_aux*G - 0.5*r_ij*sum(an(j,:)*n_ij) + 6*sum(v(j,:)*CP(S(j,:),n_ij))/(M(j)*r_ij) ) +&
                        ( (3*sum(n_ij*v(j,:))/(r_ij**2))*( M(j)*v(j,:) + 2*CP(S(j,:),n_ij)/(r_ij) ) +&
                             3.5*M(j)*an(j,:)/r_ij - 4*CP(S(j,:),v(j,:))/(r_ij**3) )*Gdc2

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
        a(i,:) = ( R_i - v(i,:)*sum(R_i*v(i,:))/(U1n_i*c2+v2(i)) )/U1n_i
    a_sm(i,:) = a(i,:) - an(i,:)
    enddo    
endsubroutine r_h

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

real(8) function norm_l2(v)
    implicit none
    !Euclidean norm of vector v
    !TEMP
    real(8), dimension(3) :: v
    !BODY
    norm_l2 = sqrt(sum(v**2))
    return
endfunction norm_l2

end module MySubs

!The gateway routine
subroutine mexFunction(nlhs, plhs, nrhs, prhs)
    use MySubs

    integer :: mxGetM, mxGetN, mxGetPr
    integer :: mxIsNumeric, mxCreateDoubleMatrix
    integer :: plhs(*), prhs(*)
    !integer :: x_pr, y_pr
    integer :: nlhs, nrhs
    !real(8) :: 
    !real(8), dimension(3) :: x,y
    real(8), dimension(N_b,3) :: x_in,v_in,S,g_out
    real(8), dimension(N_b) :: M

    !Get the size of the input array.
    !m = mxGetM(prhs(1))
    !n = mxGetN(prhs(1))
    
    call mxCopyPtrToReal8(mxGetPr(prhs(1)), x_in, N_b*size)
    call mxCopyPtrToReal8(mxGetPr(prhs(2)), v_in, N_b*size)    
    call mxCopyPtrToReal8(mxGetPr(prhs(3)), M, N_b)    
    call mxCopyPtrToReal8(mxGetPr(prhs(4)), S, N_b*size)
    
    call r_h(g_out,x_in,v_in,M,S)

    !Create matrix for the return argument.
    
    plhs(1) = mxCreateDoubleMatrix(N_b, size, mxREAL)
    call mxCopyReal8ToPtr(g_out, mxGetPr(plhs(1)), N_b*size)     
    return
endsubroutine mexFunction
