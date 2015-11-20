!     Subroutines for basic 3D linear elastic elements


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
subroutine el_gurson_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
    n_state_variables, initial_state_variables, &                                                        ! Input variables
    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
!    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
               
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
          
    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element

    ! Local Variables
    integer      :: n_points,kint, i, j, k, l, a

    real (prec)  ::  strain(6), dstrain(3,3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, nu              ! Material properties

    real (prec)  ::  YY, edot0, mm, q1, q2, q3, fn, epsn, sn, fc, ff
    real (prec)  ::  dF(3,3), F_mid(3,3), F_midinv(3,3), dLbar(3,3), dw(3,3), dL(3,3)
    real (prec)  ::  ident(3,3), dR(3,3), temp(3,3), temp_inv(3,3)
    real (prec)  ::  eta, deta, volume, JJ, det_temp, stress1(6),sum1, strain_inc(3,3)
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0

    E = element_properties(1)
    nu = element_properties(2)
    YY = element_properties(3)
    edot0 = element_properties(4)
    mm = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    epsn = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)


! Initialize the variables
     dF = 0.d0
     F_mid = 0.d0
     volume = 0.d0
     eta = 0.d0
     deta = 0.d0
      dnbardy=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     --  Loop over integration points
do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

! calculate increment in deformation gradient delata F
df=0.d0
    do i = 1,3
        do j = 1,3
            do a = 1, n_nodes
            dF(i,j) = dF(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
            end do
        end do
    end do

! calculate the mid point deformation gradient F_mid
    F_mid=0.d0
     F_mid(1,1)=1.d0
     F_mid(2,2)=1.d0
     F_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

! calculate the volume averaged quantities
call invert_small(F_mid,F_midinv,JJ)
    eta = eta + JJ*w(kint)*determinant
    dL = matmul(dF,F_midinv)
    deta = deta + JJ*(dL(1,1)+dL(2,2)+dL(3,3))*w(kint)*determinant

! calculate the contribution from current integration point to the volume averaged spatial shape function derivatives
    dNdy = matmul(dNdx,F_midinv)
    dNbardy = dNbardy + JJ*dNdy*w(kint)*determinant

! compute the contribution from ITP to element volume
    volume = volume+w(kint)*determinant
end do

! Divide the volume averagted quantities by volume
    eta = eta/volume
    deta = deta/(volume*eta)
    dNbardy = dNbardy/(volume*eta)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do a second loop over the ITP
do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

! calculate increment in deformation gradient delata F
df=0.d0
    do i = 1,3
        do j = 1,3
            do a = 1,n_nodes
            dF(i,j) = dF(i,j)+dNdx(a,j)*(dof_increment(3*(a-1)+i))
            end do
        end do
    end do

! calculate the mid point deformation gradient F_mid
f_mid=0.d0
     F_mid(1,1)=1.d0
     F_mid(2,2)=1.d0
     F_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

call invert_small(F_mid,F_midinv,JJ)


! compute the corrected velocity gradient increment
dlbar=0.d0
    dLbar = matmul(dF,F_midinv)
!    write(6,*) 'dLbar', dLbar
    sum1=dLbar(1,1) + dLbar(2,2) + dLbar(3,3)
    do i = 1,3
       do j = 1,3
            if (i == j) then
               dLbar(i,j) = dLbar(i,j) + (deta-sum1)/3.d0
             end if
        end do
    end do

! compute the strain and spin increments

    strain_inc = (dLbar+transpose(dLbar))/2.d0
    dw = (dLbar-transpose(dLbar))/2.d0

! compute R increment
ident = 0.d0
ident(1,1) = 1.d0
ident(2,2) = 1.d0
ident(3,3) = 1.d0

    temp = ident-dw/2.d0

    call invert_small(temp,temp_inv,det_temp)
    dR = matmul(temp_inv,(ident+dw/2.d0))

call gurson(element_properties,n_properties,8,initial_state_variables(((kint-1)*8+1):((kint-1)*8+8)), &
updated_state_variables(((kint-1)*8+1):((kint-1)*8+8)),strain_inc,dR,stress)

!write(6,*) ' updated svars ',updated_state_variables((kint-1)*8+1:(kint-1)*8+8)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+(1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(1,3:3*n_nodes:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+(1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(2,3:3*n_nodes:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(3,1:3*n_nodes-2:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(3,2:3*n_nodes-1:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(3,3:3*n_nodes:3) = dNdy(1:n_nodes,3)+(1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

!
!        if (updated_state_variables((kint-1)*8+7) > ff) then
!            element_deleted = .true.
!        else
!            element_deleted = .false.
!       end if
!        write(6,*) stress
        stress1 = stress
        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress1)*w(kint)*determinant

    end do
  
    return
end subroutine el_gurson_dynamic

!=========================SUBTOUTINE gurson==============================================
subroutine gurson(element_properties,n_properties,n_state_variables,initial_state_variables,  &
updated_state_variables,dstrain,dR, stress1)
     use Types
     use ParamIO
     use Globals, only : TIME, dtime
     use Element_Utilities, only : invert_small
     implicit none

     integer, intent( in ) :: n_properties
     integer, intent( in ) :: n_state_variables


     real (prec), intent( in ) :: element_properties(n_properties)
     real (prec), intent( in ) :: initial_state_variables(n_state_variables)
     real (prec), intent( in ) :: dstrain(3,3)
     real (prec), intent( in ) :: dR(3,3)
     real (prec), intent( out ) :: stress1(6)
     real (prec), intent( out ) :: updated_state_variables(n_state_variables)

    integer      :: n_points,kint, i, j, k, l, a
    real (prec)  :: stress0(6), ematrix, Vf, tol
    real (prec)  :: taumat(3,3), taumatdiv(3,3), dstraindiv(3,3), pstar, Sstar(3,3), phi2, sigmastar, fstar
    real (prec)  :: ffbar, taumat1(3,3), tau1(6)
    real (prec)  :: E, nu, YY, edot0, mm, q1, q2, q3, fn, epsn, sn, fc, ff
    real (prec)  :: eps, error, dee,  dev, sigmae, pp, phi1
    real (prec)  :: dphi2_dsigmae, dphi2_dp, dsigma_ddee, dp_ddev, dphi2_ddee, dsigmae_ddee, dphi2_ddev
    real (prec)  :: XX, dXX_ddee, dXX_ddev, F1, F2, dF1_ddee, dF1_ddev, dF2_ddee, dF2_ddev
    real (prec)  :: matrix(2,2), matrix_inv(2,2), matrix_det, left(2), ddeeddev(2), debarmat, Vf_1, ematrix_1


    E = element_properties(1)
    nu = element_properties(2)
    YY = element_properties(3)
    edot0 = element_properties(4)
    mm = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    epsn = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)

    stress0 = initial_state_variables(1:6) ! Stress at start of increment
    ematrix = initial_state_variables(7)
    Vf = initial_state_variables(8)

! define 3*3 matrix taumat
    taumat = 0.d0
    taumat(1,1) = stress0(1)
    taumat(2,2) = stress0(2)
    taumat(3,3) = stress0(3)
    taumat(1,2) = stress0(4)
    taumat(2,1) = stress0(4)
    taumat(3,1) = stress0(5)
    taumat(1,3) = stress0(5)
    taumat(2,3) = stress0(6)
    taumat(3,2) = stress0(6)

! compute taumatdiv, dstraindiv, pstar and Sij star
   taumatdiv = 0.d0
    do i = 1,3
        do j = 1,3
            if (i == j) then
            taumatdiv(i,j) = taumat(i,j)-(taumat(1,1)+taumat(2,2)+taumat(3,3))/3.d0
            else
            taumatdiv(i,j) = taumat(i,j)
            end if
        end do
    end do


    dstraindiv = 0.d0
    do i = 1,3
        do j = 1,3
            if (i == j) then
            dstraindiv(i,j) = dstrain(i,j)-(dstrain(1,1)+dstrain(2,2)+dstrain(3,3))/3.d0
            else
            dstraindiv(i,j) = dstrain(i,j)
            end if
        end do
    end do


   ! write(6,*) dstrain
   !write(6,*) dstraindiv

    pstar = (taumat(1,1)+taumat(2,2)+taumat(3,3))/3.d0+E/(3.d0*(1.d0-2.d0*nu))*(dstrain(1,1)+ &
                             dstrain(2,2)+ dstrain(3,3))


    Sstar =  E/(1.d0+nu)*dstraindiv + matmul(dR, matmul(taumatdiv, transpose(dR)))


! define fstar
ffbar = (q1+dsqrt(q1**2.d0-q3))/q3

    if (Vf < fc) then
    fstar = Vf
    else if (fc < Vf .AND. Vf < ff) then
    fstar = fc + (ffbar-fc)/(ff-fc)*(Vf-fc)
    else if (Vf > ff) then
    fstar = ffbar
    end if

! define sigmastar and phi2
    sigmastar = dsqrt(3.d0/2.d0*(Sstar(1,1)**2.d0+Sstar(2,2)**2.d0+Sstar(3,3)**2.d0+ &
                       Sstar(1,2)**2.d0*2.d0+Sstar(2,3)**2.d0*2.d0+Sstar(1,3)**2.d0*2.d0))

    phi1 = sigmastar**2.d0/YY**2.d0+2.d0*q1*fstar*cosh(3.d0/2.d0*q2*pstar/YY)-(1.d0+q3*fstar**2.d0)

! Determine if phi2 is greater or smaller than eps
eps = 10.d0**(-8.d0)
taumat1=0.d0
    if (phi1 < eps) then
        do i = 1,3
            do j = 1,3
            if (i == j) then
            taumat1(i,j) = Sstar(i,j)+pstar
            else
            taumat1(i,j) = Sstar(i,j)
            end if
            end do
         end do
       debarmat=0.d0
       ematrix_1=ematrix

    else

    error = 1.d0
    dee=0.0d0
    dev=0.0d0
    tol = 10.d0**(-6.d0)
    do while (error > tol)
            sigmae = sigmastar - 1.5d0*E/(1.d0+nu)*dee
            pp = pstar - E/3.d0/(1.d0-2.d0*nu)*dev
            phi2 = sigmae**2.d0/YY**2.d0 + 2.d0*q1*fstar*cosh(1.5d0*q2*pp/YY)-(1.d0+q3*fstar**2.d0)
            if (phi2<0.d0) then
                dee=dee/10.d0
                dev=dev/10.d0
                cycle
            endif
            phi2 = dsqrt(phi2)
            dphi2_dsigmae = sigmae/phi2/YY**2.d0
            dphi2_dp = 1.5d0*q1*q2*fstar*sinh(1.5d0*q2*pp/YY)/phi2/YY
            dsigmae_ddee = -1.5d0*e/(1+nu)
            dp_ddev= -E/(3.d0*(1-2.d0*nu))
            dphi2_ddee = dphi2_dsigmae * dsigmae_ddee
            dphi2_ddev = dphi2_dp*dp_ddev

            XX = dsqrt(dphi2_dsigmae**2.d0 + 2.d0/9.d0*dphi2_dp**2.d0)
            dXX_ddee = (2.d0*sigmae*dsigmae_ddee/phi2/phi2/YY**4.d0 - sigmae**2.d0/YY**4.d0*2.d0/phi2**3.d0 &
               *dphi2_ddee + 2.d0/9.d0*dphi2_dp**2.d0*(-2.d0)/phi2*dphi2_ddee)/2.d0/XX
            dXX_ddev = (q1**2.d0*q2**2.d0*fstar**2.d0/phi2**2.d0/YY**2.d0*sinh(1.5d0*q2/YY*pp)*cosh(1.5d0*q2/YY*pp) &
                * 1.5d0*q2/YY*dp_ddev + 2.d0/9.d0*dphi2_dp**2.d0*(-2.d0)/phi2*dphi2_ddev-sigmae**2.d0/YY**4.d0*2.d0 &
                    /phi2**3.d0*dphi2_ddev)/2.d0/XX

            F1 = XX*dee/dtime/edot0 - dphi2_dsigmae*phi2**mm
            F2 = XX*dev/dtime/edot0 - dphi2_dp*phi2**mm

            dF1_ddee = XX/dtime/edot0 - dsigmae_ddee/YY**2.d0*phi2**(mm-1.d0) - &
                sigmae/YY**2.d0*(mm-1.d0)*phi2**(mm-2.d0)*dphi2_ddee + dee/dtime/edot0*dXX_ddee
            dF1_ddev = - sigmae/YY**2.d0*(mm-1.d0)*phi2**(mm-2.d0)*dphi2_ddev + dee/dtime/edot0*dXX_ddev
            dF2_ddee = dev/dtime/edot0 * dXX_ddee - dphi2_dp*(phi2**mm)*(mm-1.d0)/(phi2)*dphi2_ddee
            dF2_ddev = XX/dtime/edot0 - dphi2_dp*(phi2**mm)*(mm-1.d0)/(phi2)*dphi2_ddev - 1.5d0*q1*q2*fstar/YY*phi2**(mm-1.d0)* &
                cosh(1.5d0*q2*pp/YY)*dp_ddev*1.5d0*q2/YY + dev/dtime/edot0 * dXX_ddev

            matrix(1,1:2)=[dF1_ddee,df1_ddev]
            matrix(2,1:2)=[df2_ddee,df2_ddev]
    call invert_small(matrix,matrix_inv,matrix_det)
            left =[-F1,-F2]
            ddeeddev = matmul(matrix_inv, left)
            dee=dee+ddeeddev(1)
            dev=dev+ddeeddev(2)
            error=dsqrt(ddeeddev(1)**2.d0+ddeeddev(2)**2.d0)/sqrt(dee**2+dev**2.d0)
    end do

 taumat1 = 0.d0
        do i = 1,3
            do j = 1,3
                if (i == j) then
                if(sigmastar/=0.d0) then
                    taumat1(i,j) = Sstar(i,j)-dee*E/(1.d0+nu)*1.5d0*Sstar(i,j)/sigmastar &
                                   +(pstar-dev*E/(3.d0*(1.d0-2.d0*nu)))
                 else
                  taumat1(i,j) = Sstar(i,j) +(pstar-dev*E/(3.d0*(1.d0-2.d0*nu)))
                 end if
                else
                if(sigmastar/=0.d0) then
                    taumat1(i,j) = Sstar(i,j)-dee*E/(1.d0+nu)*1.5d0*Sstar(i,j)/sigmastar
                    else
                     taumat1(i,j) = Sstar(i,j)
                     end if
                end if
            end do
        end do

        debarmat =  edot0*dtime/(1.d0-Vf)*phi2**mm*(dphi2_dsigmae**2.d0+2.d0/9.d0*dphi2_dp**2.d0)**(-0.5d0) &
                   *(dphi2_dsigmae*sigmae+1.d0/3.d0*dphi2_dp*pp)

        ematrix_1 = ematrix +debarmat

    end if

        Vf_1 = 1.d0+(Vf-1.d0)*exp(-dev)+fn*debarmat/Sn/dsqrt(2.d0*3.1415d0)*exp(-0.5d0*((ematrix-epsn)/sn)**2.d0)

    stress1(1) = taumat1(1,1)
    stress1(2) = taumat1(2,2)
    stress1(3) = taumat1(3,3)
    stress1(4) = taumat1(1,2)
    stress1(5) = taumat1(1,3)
    stress1(6) = taumat1(2,3)

    updated_state_variables(1) = stress1(1)
    updated_state_variables(2) = stress1(2)
    updated_state_variables(3) = stress1(3)
    updated_state_variables(4) = stress1(4)
    updated_state_variables(5) = stress1(5)
    updated_state_variables(6) = stress1(6)
    updated_state_variables(7) = ematrix_1
    updated_state_variables(8) = Vf_1

    return
end subroutine gurson


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_gurson_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
    use Element_Utilities, only : dNbardy => vol_avg_shape_function_derivatives_3D
    implicit none

    integer, intent( in )         :: lmn                                                    ! Element number
    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
    integer, intent( in )         :: n_properties                                           ! # properties for the element
    integer, intent( in )         :: length_coord_array                                     ! Total # coords
    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
    integer, intent( in )         :: n_field_variables                                      ! # field variables

    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
    !  type node
    !      sequence
    !      integer :: flag                          ! Integer identifier
    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
    !      integer :: n_coords                      ! Total no. coordinates for the node
    !      integer :: dof_index                     ! Index of first DOF in dof array
    !      integer :: n_dof                         ! Total no. of DOF for node
    !   end type node
    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element

    character (len=100), intent(in) :: field_variable_names(n_field_variables)

    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,x2,(x3) for each node in turn
    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment

    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
    real( prec ), intent( inout )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step
             
    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
  
   integer      :: n_points,kint, i, j, k, l, a

    real (prec)  ::  strain(6), dstrain(3,3)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, nu              ! Material properties

    real (prec)  ::  YY, edot0, mm, q1, q2, q3, fn, epsn, sn, fc, ff, Vf
    real (prec)  ::  dF(3,3), F_mid(3,3), F_midinv(3,3), dLbar(3,3), dw(3,3), dL(3,3)
    real (prec)  ::  ident(3,3), dR(3,3), temp(3,3), temp_inv(3,3)
    real (prec)  ::  eta, deta, volume, JJ, det_temp, stress1(6),sum1, strain_inc(3,3)
    !     Subroutine to compute element force vector for a linear elastodynamic problem
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    nodal_fieldvariables = 0.d0

    E = element_properties(1)
    nu = element_properties(2)
    YY = element_properties(3)
    edot0 = element_properties(4)
    mm = element_properties(5)
    q1 = element_properties(6)
    q2 = element_properties(7)
    q3 = element_properties(8)
    fn = element_properties(9)
    epsn = element_properties(10)
    sn = element_properties(11)
    fc = element_properties(12)
    ff = element_properties(13)


! Initialize the variables
     dF = 0.d0
     F_mid = 0.d0
     volume = 0.d0
     eta = 0.d0
     deta = 0.d0
      dnbardy=0.d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     --  Loop over integration points
do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

! calculate increment in deformation gradient delata F
df=0.d0
    do i = 1,3
        do j = 1,3
            do a = 1, n_nodes
            dF(i,j) = dF(i,j)+dNdx(a,j)*dof_increment(3*(a-1)+i)
            end do
        end do
    end do

! calculate the mid point deformation gradient F_mid
    F_mid=0.d0
     F_mid(1,1)=1.d0
     F_mid(2,2)=1.d0
     F_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

! calculate the volume averaged quantities
call invert_small(F_mid,F_midinv,JJ)
    eta = eta + JJ*w(kint)*determinant
    dL = matmul(dF,F_midinv)
    deta = deta + JJ*(dL(1,1)+dL(2,2)+dL(3,3))*w(kint)*determinant

! calculate the contribution from current integration point to the volume averaged spatial shape function derivatives
    dNdy = matmul(dNdx,F_midinv)
    dNbardy = dNbardy + JJ*dNdy*w(kint)*determinant

! compute the contribution from ITP to element volume
    volume = volume+w(kint)*determinant
end do

! Divide the volume averagted quantities by volume
    eta = eta/volume
    deta = deta/(volume*eta)
    dNbardy = dNbardy/(volume*eta)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Do a second loop over the ITP
do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)

! calculate increment in deformation gradient delata F
df=0.d0
    do i = 1,3
        do j = 1,3
            do a = 1,n_nodes
            dF(i,j) = dF(i,j)+dNdx(a,j)*(dof_increment(3*(a-1)+i))
            end do
        end do
    end do

! calculate the mid point deformation gradient F_mid
f_mid=0.d0
     F_mid(1,1)=1.d0
     F_mid(2,2)=1.d0
     F_mid(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F_mid(i,j)=F_mid(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+0.5d0*dof_increment(3*(a-1)+i))
                end do
            end do
        end do

call invert_small(F_mid,F_midinv,JJ)


! compute the corrected velocity gradient increment
dlbar=0.d0
    dLbar = matmul(dF,F_midinv)
!    write(6,*) 'dLbar', dLbar
    sum1=dLbar(1,1) + dLbar(2,2) + dLbar(3,3)
    do i = 1,3
       do j = 1,3
            if (i == j) then
               dLbar(i,j) = dLbar(i,j) + (deta-sum1)/3.d0
             end if
        end do
    end do

! compute the strain and spin increments

    strain_inc = (dLbar+transpose(dLbar))/2.d0
    dw = (dLbar-transpose(dLbar))/2.d0

! compute R increment
ident = 0.d0
ident(1,1) = 1.d0
ident(2,2) = 1.d0
ident(3,3) = 1.d0

    temp = ident-dw/2.d0

    call invert_small(temp,temp_inv,det_temp)
    dR = matmul(temp_inv,(ident+dw/2.d0))

call gurson(element_properties,n_properties,8,initial_state_variables(((kint-1)*8+1):((kint-1)*8+8)), &
updated_state_variables(((kint-1)*8+1):((kint-1)*8+8)),strain_inc,dR,stress)

!write(6,*) ' updated svars ',updated_state_variables((kint-1)*8+1:(kint-1)*8+8)
        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)+(1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(1,2:3*n_nodes-1:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(1,3:3*n_nodes:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(2,1:3*n_nodes-2:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)+(1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(2,3:3*n_nodes:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(3,1:3*n_nodes-2:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,1)-dNdy(1:n_nodes,1))
        B(3,2:3*n_nodes-1:3) = (1.d0/3.d0)*(dNbardy(1:n_nodes,2)-dNdy(1:n_nodes,2))
        B(3,3:3*n_nodes:3) = dNdy(1:n_nodes,3)+(1.d0/3.d0)*(dNbardy(1:n_nodes,3)-dNdy(1:n_nodes,3))
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

       Vf = updated_state_variables(8*(kint-1)+8)
        stress1 = updated_state_variables(8*(kint-1)+1 : 8*kint-2)
!        write(6,*) ' stress1 in field proj '
 !       write(6,*) stress1(1:6)

        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress1(6)*N(1:n_nodes)*determinant*w(kint)
!            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
!                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'Vf',2) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + Vf*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
  end do
!  write(6,*) nodal_fieldvariables(1,1:n_nodes)
    return
end subroutine fieldvars_gurson_dynamic


