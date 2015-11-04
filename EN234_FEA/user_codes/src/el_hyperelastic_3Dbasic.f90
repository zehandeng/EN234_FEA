!==========================SUBROUTINE el_linelast_3dbasic ==============================
subroutine el_hyperelastic_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties, element_coords, length_coord_array, &                      ! Input variables
    dof_increment, dof_total, length_dof_array, &                                                ! Input variables
    n_state_variables, initial_state_variables, &                                                ! Input variables
    updated_state_variables,element_stiffness,element_residual, fail)      ! Output variables                          ! Output variables
    use Types
    use ParamIO
    !  use Globals, only: TIME,DTIME  For a time dependent problem uncomment this line to access the time increment and total time
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only : dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D

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
  
    logical, intent( out )        :: fail                                                   ! Set to .true. to force a timestep cutback
    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
    real( prec ), intent( out )   :: element_stiffness(length_dof_array,length_dof_array)   ! Element stiffness (ROW,COLUMN)
    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)

    ! Local Variables
    integer      :: n_points,kint, i, j,a

    real (prec)  ::  strain(6), dstrain(6), F(3,3), Finv(3,3), matrix(6,6)            ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6), Bcauchy1D(6), ident(6), Bcauchyinv1D(6),stressmat(3,3)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6), Bstar(9,3*n_nodes), G(6,9)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array), Bcauchy(3,3),Bcauchyinv(3,3), S(3,length_dof_array/3)            ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  volume, mu1, K1, detF,detB
    real (prec)  ::  ident_dyadic_Binv(6,6), ident_dyadic_ident(6,6), B_dyadic_Binv(6,6)
    real (prec)  ::  Pvec(3*n_nodes), Pmat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Svec(3*n_nodes), Smat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  doff(length_dof_array/3,3)
    !
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    fail = .false.
    
    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27

    call initialize_integration_points(n_points, n_nodes, xi, w)

    element_residual = 0.d0
    element_stiffness = 0.d0
!
    D = 0.d0
    mu1 = element_properties(1)
    K1 = element_properties(2)


    !     --  Loop over integration points


     do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
!        doff=transpose(reshape((dof_total+dof_increment),(/3,length_dof_array/3/)))

     F=0.d0
     F(1,1)=1.d0
     F(2,2)=1.d0
     F(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F(i,j)=F(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
                end do
            end do
        end do


    Bcauchy=matmul(F,TRANSPOSE(F))

   call invert_small(F,Finv,detF)
         dNdy = matmul(dNdx,Finv)

        B = 0.d0
        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)


     Bcauchy1D = 0.d0
     Bcauchy1D(1) = Bcauchy(1,1)
     Bcauchy1D(2) = Bcauchy(2,2)
     Bcauchy1D(3) = Bcauchy(3,3)
     Bcauchy1D(4) = Bcauchy(1,2)
     Bcauchy1D(5) = Bcauchy(1,3)
     Bcauchy1D(6) = Bcauchy(2,3)

    call invert_small(Bcauchy,Bcauchyinv,detB)
     Bcauchyinv1D(1) = Bcauchyinv(1,1)
     Bcauchyinv1D(2) = Bcauchyinv(2,2)
     Bcauchyinv1D(3) = Bcauchyinv(3,3)
     Bcauchyinv1D(4) = Bcauchyinv(1,2)
     Bcauchyinv1D(5) = Bcauchyinv(1,3)
     Bcauchyinv1D(6) = Bcauchyinv(2,3)


     ident = [1.d0, 1.d0, 1.d0, 0.d0, 0.d0, 0.d0]



     matrix = 0.d0
     do i = 1,3
        matrix(i,i) = 1.d0
        matrix(i+3,i+3) = 0.5d0
     end do

     ident_dyadic_Binv = spread(ident,dim=2,ncopies=6)*spread(Bcauchyinv1D,dim=1,ncopies=6)
     ident_dyadic_ident = spread(ident,dim=2,ncopies=6)*spread(ident,dim=1,ncopies=6)
     B_dyadic_Binv = spread(Bcauchy1D,dim=2,ncopies=6)*spread(Bcauchyinv1D,dim=1,ncopies=6)

     D = mu1/(detF**(2.d0/3.d0))*matrix+mu1/(3.d0*detF**(2.d0/3.d0))*((Bcauchy1D(1)+Bcauchy1D(2)+Bcauchy1D(3))/3.d0 &
     *ident_dyadic_Binv-ident_dyadic_ident-B_dyadic_Binv)+K1*detF*(detF-0.5d0)*ident_dyadic_Binv

     G = 0.d0
     G(1,1:9) = [2.d0*Bcauchy(1,1),0.d0,0.d0,2.d0*Bcauchy(1,2),0.d0,2.d0*Bcauchy(1,3),0.d0,0.d0,0.d0]
     G(2,1:9) = [0.d0,2.d0*Bcauchy(2,2),0.d0,0.d0,2.d0*Bcauchy(1,2),0.d0,0.d0,2.d0*Bcauchy(2,3),0.d0]
     G(3,1:9) = [0.d0,0.d0,2.d0*Bcauchy(3,3),0.d0,0.d0,0.d0,2.d0*Bcauchy(1,3),0.d0,2.d0*Bcauchy(1,3)]
     G(4,1:9) = [2.d0*Bcauchy(1,2),2.d0*Bcauchy(1,2),0.d0,2.d0*Bcauchy(2,2),2.d0*Bcauchy(1,1),2.d0*Bcauchy(2,3)&
                  ,0.d0,2.d0*Bcauchy(1,3),0.d0]
     G(5,1:9) = [2.d0*Bcauchy(1,3),0.d0,2.d0*Bcauchy(1,3),2.d0*Bcauchy(2,3),0.d0,2.d0*Bcauchy(3,3),2.d0*Bcauchy(1,1)&
                  ,0.d0,2.d0*Bcauchy(1,2)]
     G(6,1:9) = [0.d0,2.d0*Bcauchy(2,3),2.d0*Bcauchy(2,3),0.d0,2.d0*Bcauchy(1,3),0.d0,2.d0*Bcauchy(1,2),&
                  2.d0*Bcauchy(3,3),2.d0*Bcauchy(2,2)]


     Bstar = 0.d0
     Bstar(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
     Bstar(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
     Bstar(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
     Bstar(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
     Bstar(5,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
     Bstar(6,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
     Bstar(7,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
     Bstar(8,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
     Bstar(9,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)


     do i = 1,3
        do j = 1,3
        if (i ==j) then
        stressmat(i,j) = mu1/(detF**(5.d0/3.d0))*(Bcauchy(i,j)-1.d0/3.d0*(Bcauchy(1,1)+Bcauchy(2,2) &
        +Bcauchy(3,3)))+K1*(detF-1.d0)
        else
        stressmat(i,j) = mu1/(detF**(5.d0/3.d0))*(Bcauchy(i,j))
        end if
        end do
     end do

    stressmat = detF*stressmat

     stress=0.d0
     stress(1) = stressmat(1,1)
     stress(2) = stressmat(2,2)
     stress(3) = stressmat(3,3)
     stress(4) = stressmat(1,2)
     stress(5) = stressmat(1,3)
     stress(6) = stressmat(2,3)

     strain = matmul(B,dof_total)
     dstrain = matmul(B,dof_increment)

      S = reshape(matmul(transpose(B),stress),(/3,length_dof_array/3/))
          do i = 1,n_nodes
                Pvec = reshape(spread(transpose(dNdy(i:i,1:3)),dim=2,ncopies=n_nodes),(/3*n_nodes/))
                Pmat(3*i-2:3*i,1:3*n_nodes) = spread(Pvec,dim=1,ncopies=3)
                Svec = reshape(spread(S(1:3,i:i),dim=2,ncopies=n_nodes),(/3*n_nodes/))
                Smat(3*i-2:3*i,1:3*n_nodes) = spread(Svec,dim=1,ncopies=3)
          end do
     Sigma = Pmat*transpose(Smat)

        element_residual = element_residual - matmul(transpose(B),stress)*w(kint)*determinant
        element_stiffness = element_stiffness+ (matmul(transpose(B),matmul(D,matmul(G,Bstar))) &
        -sigma)*w(kint)*determinant
    end do
  
    return
end subroutine el_hyperelastic_3dbasic


!==========================SUBROUTINE el_linelast_3dbasic_dynamic ==============================
!subroutine el_linelast_3dbasic_dynamic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
!    n_properties, element_properties,element_coords, length_coord_array, &                               ! Input variables
!    dof_increment, dof_total, length_dof_array,  &                                                       ! Input variables
!    n_state_variables, initial_state_variables, &                                                        ! Input variables
!    updated_state_variables,element_residual,element_deleted)                                            ! Output variables
!    use Types
!    use ParamIO
!    use Mesh, only : node
!    use Element_Utilities, only : N => shape_functions_3D
!    use Element_Utilities, only:  dNdxi => shape_function_derivatives_3D
!    use Element_Utilities, only:  dNdx => shape_function_spatial_derivatives_3D
!    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
!    use Element_Utilities, only : dxdxi => jacobian_3D
!    use Element_Utilities, only : initialize_integration_points
!    use Element_Utilities, only : calculate_shapefunctions
!    use Element_Utilities, only : invert_small
!    implicit none
!
!    integer, intent( in )         :: lmn                                                    ! Element number
!    integer, intent( in )         :: element_identifier                                     ! Flag identifying element type (specified in .in file)
!    integer, intent( in )         :: n_nodes                                                ! # nodes on the element
!    integer, intent( in )         :: n_properties                                           ! # properties for the element
!    integer, intent( in )         :: length_coord_array                                     ! Total # coords
!    integer, intent( in )         :: length_dof_array                                       ! Total # DOF
!    integer, intent( in )         :: n_state_variables                                      ! # state variables for the element
!
!    type (node), intent( in )     :: node_property_list(n_nodes)                  ! Data structure describing storage for nodal variables - see below
!    !  type node
!    !      sequence
!    !      integer :: flag                          ! Integer identifier
!    !      integer :: coord_index                   ! Index of first coordinate in coordinate array
!    !      integer :: n_coords                      ! Total no. coordinates for the node
!    !      integer :: dof_index                     ! Index of first DOF in dof array
!    !      integer :: n_dof                         ! Total no. of DOF for node
!    !   end type node
!    !   Access these using node_property_list(k)%n_coords eg to find the number of coords for the kth node on the element
!
!    real( prec ), intent( in )    :: element_coords(length_coord_array)                     ! Coordinates, stored as x1,(x2),(x3) for each node in turn
!    real( prec ), intent( in )    :: dof_increment(length_dof_array)                        ! DOF increment, stored as du1,du2,du3,du4... for each node in turn
!    real( prec ), intent( in )    :: dof_total(length_dof_array)                            ! accumulated DOF, same storage as for increment
!
!    real( prec ), intent( in )    :: element_properties(n_properties)                       ! Element or material properties, stored in order listed in input file
!    real( prec ), intent( in )    :: initial_state_variables(n_state_variables)             ! Element state variables.  Defined in this routine
!
!    real( prec ), intent( inout ) :: updated_state_variables(n_state_variables)             ! State variables at end of time step
!    real( prec ), intent( out )   :: element_residual(length_dof_array)                     ! Element residual force (ROW)
!
!    logical, intent( inout )      :: element_deleted                                        ! Set to .true. to delete element
!
!    ! Local Variables
!    integer      :: n_points,kint
!
!    real (prec)  ::  strain(6), dstrain(6)             ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
!    real (prec)  ::  stress(6)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
!    real (prec)  ::  D(6,6)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
!    real (prec)  ::  B(6,length_dof_array)             ! strain = B*(dof_total+dof_increment)
!    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
!    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
!    real (prec)  :: E, xnu, D44, D11, D12              ! Material properties
!    !
!    !     Subroutine to compute element force vector for a linear elastodynamic problem
!    !     El props are:
!
!    !     element_properties(1)         Young's modulus
!    !     element_properties(2)         Poisson's ratio
!
!    x = reshape(element_coords,(/3,length_coord_array/3/))
!
!    if (n_nodes == 4) n_points = 1
!    if (n_nodes == 10) n_points = 4
!    if (n_nodes == 8) n_points = 8
!    if (n_nodes == 20) n_points = 27
!
!    call initialize_integration_points(n_points, n_nodes, xi, w)
!
!    element_residual = 0.d0
!
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44
!
!    !     --  Loop over integration points
!    do kint = 1, n_points
!        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdx(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdx(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdx(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdx(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdx(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdx(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdx(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdx(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdx(1:n_nodes,2)
!
!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)
!
!        stress = matmul(D,strain+dstrain)
!        element_residual(1:3*n_nodes) = element_residual(1:3*n_nodes) - matmul(transpose(B),stress)*w(kint)*determinant
!
!    end do
!
!    return
!end subroutine el_linelastc_3dbasic_dynamic


!==========================SUBROUTINE fieldvars_linelast_3dbasic ==============================
subroutine fieldvars_hyperelastic_3dbasic(lmn, element_identifier, n_nodes, node_property_list, &           ! Input variables
    n_properties, element_properties,element_coords,length_coord_array, &                                ! Input variables
    dof_increment, dof_total, length_dof_array,  &                                                      ! Input variables
    n_state_variables, initial_state_variables,updated_state_variables, &                               ! Input variables
    n_field_variables,field_variable_names, &                                                           ! Field variable definition
    nodal_fieldvariables)      ! Output variables
    use Types
    use ParamIO
    use Mesh, only : node
    use Element_Utilities, only : N => shape_functions_3D
    use Element_Utilities, only: dNdxi => shape_function_derivatives_3D
    use Element_Utilities, only: dNdx => shape_function_spatial_derivatives_3D
    use Element_Utilities, only:  dNdy => shape_function_spatial_derivatives_3D
    use Element_Utilities, only : xi => integrationpoints_3D, w => integrationweights_3D
    use Element_Utilities, only : dxdxi => jacobian_3D
    use Element_Utilities, only : initialize_integration_points
    use Element_Utilities, only : calculate_shapefunctions
    use Element_Utilities, only : invert_small
    use Element_Utilities, only : dNbardx => vol_avg_shape_function_derivatives_3D
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
    real( prec ), intent( in )    :: updated_state_variables(n_state_variables)             ! State variables at end of time step

    real( prec ), intent( out )   :: nodal_fieldvariables(n_field_variables,n_nodes)        ! Nodal field variables
  
    ! Local Variables
    logical      :: strcmp
    ! Local Variables
    integer      :: n_points,kint, i, j,a, k

    real (prec)  ::  strain(6), dstrain(6), F(3,3), Finv(3,3), matrix(6,6)            ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
    real (prec)  ::  stress(6), Bcauchy1D(6), ident(6), Bcauchyinv1D(6),stressmat(3,3)                        ! Stress vector contains [s11, s22, s33, s12, s13, s23]
    real (prec)  ::  D(6,6), Bstar(9,3*n_nodes), G(6,9)                            ! stress = D*(strain+dstrain)  (NOTE FACTOR OF 2 in shear strain)
    real (prec)  ::  B(6,length_dof_array), Bcauchy(3,3),Bcauchyinv(3,3),  S(3,length_dof_array/3)            ! strain = B*(dof_total+dof_increment)
    real (prec)  ::  dxidx(3,3), determinant           ! Jacobian inverse and determinant
    real (prec)  ::  x(3,length_coord_array/3)         ! Re-shaped coordinate array x(i,a) is ith coord of ath node
    real (prec)  ::  E, xnu, D44, D11, D12              ! Material properties
    real (prec)  ::  volume, mu1, K1, detF,detB
!    real (prec)  ::  ident_dyadic_Binv(6,6), ident_dyadic_ident(6,6), B_dyadic_Binv(6,6)
!    real (prec)  ::  Pvec(3*n_nodes), Pmat(3*n_nodes,3*n_nodes)
!    real (prec)  ::  Svec(3*n_nodes), Smat(3*n_nodes,3*n_nodes)
    real (prec)  ::  Sigma(3*n_nodes,3*n_nodes)
    real (prec)  ::  doff(length_dof_array/3,3)
    real (prec)  ::  p, sdev(6), smises
    !     Subroutine to compute element stiffness matrix and residual force vector for 3D linear elastic elements
    !     El props are:

    !     element_properties(1)         Young's modulus
    !     element_properties(2)         Poisson's ratio

    !fail = .false.

    x = reshape(element_coords,(/3,length_coord_array/3/))

    if (n_nodes == 4) n_points = 1
    if (n_nodes == 10) n_points = 4
    if (n_nodes == 8) n_points = 8
    if (n_nodes == 20) n_points = 27


    call initialize_integration_points(n_points, n_nodes, xi, w)
    nodal_fieldvariables = 0.d0
    D = 0.d0
    mu1 = element_properties(1)
    K1 = element_properties(2)
!     volume = 0.d0
!     dNbardx = 0.d0
!     do kint = 1, n_points
!        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
!        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
!        call invert_small(dxdxi,dxidx,determinant)
!        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
!        dNbardx(1:n_nodes,1) = dNbardx(1:n_nodes,1)+dNdx(1:n_nodes,1)*w(kint)*determinant
!        dNbardx(1:n_nodes,2) = dNbardx(1:n_nodes,2)+dNdx(1:n_nodes,2)*w(kint)*determinant
!        dNbardx(1:n_nodes,3) = dNbardx(1:n_nodes,3)+dNdx(1:n_nodes,3)*w(kint)*determinant
!        volume = volume+w(kint)*determinant
!        end do




!    element_residual = 0.d0
!    element_stiffness = 0.d0
!
!    D = 0.d0
!    E = element_properties(1)
!    xnu = element_properties(2)
!    d44 = 0.5D0*E/(1+xnu)
!    d11 = (1.D0-xnu)*E/( (1+xnu)*(1-2.D0*xnu) )
!    d12 = xnu*E/( (1+xnu)*(1-2.D0*xnu) )
!    D(1:3,1:3) = d12
!    D(1,1) = d11
!    D(2,2) = d11
!    D(3,3) = d11
!    D(4,4) = d44
!    D(5,5) = d44
!    D(6,6) = d44

    !     --  Loop over integration points


     do kint = 1, n_points
        call calculate_shapefunctions(xi(1:3,kint),n_nodes,N,dNdxi)
        dxdxi = matmul(x(1:3,1:n_nodes),dNdxi(1:n_nodes,1:3))
        call invert_small(dxdxi,dxidx,determinant)
        dNdx(1:n_nodes,1:3) = matmul(dNdxi(1:n_nodes,1:3),dxidx)
        B = 0.d0
        doff=transpose(reshape((dof_total+dof_increment),(/3,length_dof_array/3/)))
!    do i = 1,3
!       do  j = 1,3
!          do a = 1,n_nodes
!         if (i == j) then
!   F(i,j) = 1+dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
!         else
!   F(i,j) = dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
!         end if
!         end do
!       end do
!    end do
!
!    do i = 1,3
!       do j = 1,3
!    Bcauchy(i,j) = F(i,1)*F(j,1)+F(i,2)*F(j,2)+F(i,3)*F(j,3)
!       end do

    F=0.d0
     F(1,1)=1.d0
     F(2,2)=1.d0
     F(3,3)=1.d0
        do i=1,3
            do j=1,3
                do a=1,n_nodes
                    F(i,j)=F(i,j)+dNdx(a,j)*(dof_total(3*(a-1)+i)+dof_increment(3*(a-1)+i))
                end do
            end do
        end do

    call invert_small(F,Finv,detF)
!          dNdy(1:n_nodes,1:3) = matmul(dNdx(1:n_nodes,1:3),Finv)
!
!        B = 0.d0
!        B(1,1:3*n_nodes-2:3) = dNdy(1:n_nodes,1)
!        B(2,2:3*n_nodes-1:3) = dNdy(1:n_nodes,2)
!        B(3,3:3*n_nodes:3)   = dNdy(1:n_nodes,3)
!        B(4,1:3*n_nodes-2:3) = dNdy(1:n_nodes,2)
!        B(4,2:3*n_nodes-1:3) = dNdy(1:n_nodes,1)
!        B(5,1:3*n_nodes-2:3) = dNdy(1:n_nodes,3)
!        B(5,3:3*n_nodes:3)   = dNdy(1:n_nodes,1)
!        B(6,2:3*n_nodes-1:3) = dNdy(1:n_nodes,3)
!        B(6,3:3*n_nodes:3)   = dNdy(1:n_nodes,2)

        Bcauchy=matmul(F,TRANSPOSE(F))

     stressmat = 0.d0
     do i = 1,3
        do j = 1,3
        if (i ==j) then
        stressmat(i,j) = mu1/(detF**(5.d0/3.d0))*(Bcauchy(i,j)-1.d0/3.d0*(Bcauchy(1,1)+Bcauchy(2,2) &
        +Bcauchy(3,3)))+K1*(detF-1.d0)
        else
        stressmat(i,j) = mu1/(detF**(5.d0/3.d0))*(Bcauchy(i,j))
        end if
        end do
     end do

            stress(1)=stressmat(1,1)
            stress(2)=stressmat(2,2)
            stress(3)=stressmat(3,3)
            stress(4)=stressmat(1,2)
            stress(5)=stressmat(1,3)
            stress(6)=stressmat(2,3)

!
!        strain = matmul(B,dof_total)
!        dstrain = matmul(B,dof_increment)

        p = sum(stress(1:3))/3.d0
        sdev = stress
        sdev(1:3) = sdev(1:3)-p
        smises = dsqrt( dot_product(sdev(1:3),sdev(1:3)) + 2.d0*dot_product(sdev(4:6),sdev(4:6)) )*dsqrt(1.5d0)
        ! In the code below the strcmp( string1, string2, nchar) function returns true if the first nchar characters in strings match
        do k = 1,n_field_variables
            if (strcmp(field_variable_names(k),'S11',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(1)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S22',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(2)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S33',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(3)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S12',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(4)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S13',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(5)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'S23',3) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + stress(6)*N(1:n_nodes)*determinant*w(kint)
            else if (strcmp(field_variable_names(k),'SMISES',6) ) then
                nodal_fieldvariables(k,1:n_nodes) = nodal_fieldvariables(k,1:n_nodes) + smises*N(1:n_nodes)*determinant*w(kint)
            endif
        end do
 
    end do
  
    return
end subroutine fieldvars_hyperelastic_3dbasic

