!**********************************************************************************************************************************
! LICENSING
! Copyright (C) 2015-2016  National Renewable Energy Laboratory
!
!    This file is part of WakeDynamics.
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.
!
!**********************************************************************************************************************************
! File last committed: $Date$
! (File) Revision #: $Rev$
! URL: $HeadURL$
!**********************************************************************************************************************************
!> WakeDynamics is a time-domain aerodynamics module for horizontal-axis wind turbines.
module WakeDynamics
    
   use NWTC_Library
   use WakeDynamics_Types
   use WakeDynamics_IO
     
   implicit none

   private
         

   ! ..... Public Subroutines ...................................................................................................

   public :: WD_Init                           ! Initialization routine
   public :: WD_End                            ! Ending routine (includes clean up)
   public :: WD_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                               !   continuous states, and updating discrete states
   public :: WD_CalcOutput                     ! Routine for computing outputs
   public :: WD_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   
   
      ! Unit testing routines
   public :: WD_TEST_Init_BadData    
   public :: WD_TEST_Init_GoodData    
   public :: WD_TEST_UpdateStates
  
   contains  
   
real(ReKi) function EddyFilter(x_plane, D_rotor, C_Dmin, C_Dmax, C_Fmin, C_Exp)

   real(ReKi),          intent(in   ) :: x_plane
   real(ReKi),          intent(in   ) :: D_rotor
   real(ReKi),          intent(in   ) :: C_Dmin
   real(ReKi),          intent(in   ) :: C_Dmax
   real(ReKi),          intent(in   ) :: C_Fmin
   real(ReKi),          intent(in   ) :: C_Exp


   if ( x_plane <= C_Dmin*D_rotor ) then
      EddyFilter = C_Fmin
   else if (x_plane >= C_Dmax*D_rotor) then
      EddyFilter = 1
   else
      EddyFilter = C_Fmin + (1-C_Fmin)*( ( (x_plane/D_rotor) - C_DMin ) / (C_Dmax-C_Dmin) )**C_Exp
   end if


end function EddyFilter

real(ReKi) function WakeDiam( Mod_WakeDiam, nr, dr, rArr, Vx_wake, Vx_wind_disk, D_rotor, C_WakeDiam)

   integer(intKi),      intent(in   ) :: Mod_WakeDiam
   integer(intKi),      intent(in   ) :: nr
   real(ReKi),          intent(in   ) :: dr
   real(ReKi),          intent(in   ) :: rArr(0:)
   real(ReKi),          intent(in   ) :: Vx_wake(0:)
   real(ReKi),          intent(in   ) :: Vx_wind_disk
   real(ReKi),          intent(in   ) :: D_rotor
   real(ReKi),          intent(in   ) :: C_WakeDiam

   integer(IntKi) :: ILo
   real(ReKi) :: dArr(0:nr-1)
   real(ReKi) :: m(0:nr-1)
   integer(IntKi) :: i
   
   ILo = 0
   
   select case ( Mod_WakeDiam )
      case (WakeDiamMod_RotDiam) 
         WakeDiam = D_rotor
      case (WakeDiamMod_Velocity)  
               !if (C_WakeDiam <= 0 .or. C_WakeDiam >= 1) error
               ! Ensure the wake diameter is at least as large as the rotor diameter
            ! TODO: Add this to a parameter to speed up the calculation     
         dArr = rArr*2.0_ReKi 
         WakeDiam = max(D_rotor, InterpBin( (C_WakeDiam-1_ReKi)*Vx_wind_disk, darr, Vx_wake, ILo, nr )  )
      case (WakeDiamMod_MassFlux) 
         m(0) = 0.0
         do i = 0,nr-1
            m(i) = m(i-1) + pi*dr*Vx_wake(i)*darr(i)/2.0 + Vx_wake(i-1)*darr(i-1)/2.0
         end do
         WakeDiam = InterpBin( C_WakeDiam*m(nr-1), darr, m, ILo, nr )
      case (WakeDiamMod_MtmFlux)
      
      case default 
       !  ERROR: Invalid Mod_WakeDiam
   
   end select

      
end function WakeDiam


subroutine ThomasAlgorithm(nr, a, b, c, d, x)

   integer(IntKi),      intent(in   ) :: nr
   real(ReKi),          intent(inout) :: a(0:)
   real(ReKi),          intent(inout) :: b(0:)
   real(ReKi),          intent(inout) :: c(0:)
   real(ReKi),          intent(inout) :: d(0:)
   real(ReKi),          intent(inout) :: x(0:)

   real(ReKi)     :: m
   integer(IntKi) :: i
   
! Assumes all arrays are the same length
      ! Check that tridiagonal matrix is not diagonally dominant
  ! if ( abs(b(0)) <= abs(c(0)) ) then
      ! SET ERROR
  ! end if
  ! do i = 1,nr-2
   !   if ( abs(b(i)) <= abs(a(i)+c(i)) ) then
         ! Set error
   !   end if
   !end do
   !if ( abs(b(nr)) <= abs(a(nr)) ) then
      ! Set error
  ! end if
   
   do i = 1,nr-1 
      m = -a(i)/b(i-1)
      b(i) = b(i) +  m*c(i-1)
      d(i) = d(i) +  m*d(i-1)
   end do
   
   x(nr-1) = d(nr-1)/b(nr-1)
   do i = nr-2,0, -1
      x(i) = ( d(i) - c(i)*x(i+1) ) / b(i)
   end do
   

end subroutine ThomasAlgorithm

!----------------------------------------------------------------------------------------------------------------------------------   
!> This subroutine sets the initialization output data structure, which contains data to be returned to the calling program (e.g.,
!! FAST or WakeDynamics_Driver)   
subroutine WD_SetInitOut(p, InputInp, InitOut, errStat, errMsg)

   type(WD_InitOutputType),       intent(  out)  :: InitOut          ! output data
   type(WD_InitInputType),            intent(in   )  :: InputInp         ! initialization input  data 
   type(WD_ParameterType),        intent(in   )  :: p                ! Parameters
   integer(IntKi),                intent(  out)  :: errStat          ! Error status of the operation
   character(*),                  intent(  out)  :: errMsg           ! Error message if ErrStat /= ErrID_None


      ! Local variables
   integer(intKi)                               :: ErrStat2          ! temporary Error status
   character(ErrMsgLen)                         :: ErrMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'WD_SetInitOut'
   
   
   
   integer(IntKi)                               :: i, j, k, f
   integer(IntKi)                               :: NumCoords
#ifdef DBG_OUTS
   integer(IntKi)                               :: m
   character(5)                                 ::chanPrefix
#endif   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
      
                     
   
   InitOut%Ver = WD_Ver
   
   
   
   
end subroutine WD_SetInitOut
!----------------------------------------------------------------------------------------------------------------------------------   
!> This routine is called at the start of the simulation to perform initialization steps.
!! The parameters are set here and not changed during the simulation.
!! The initial states and initial guess for the input are defined.
subroutine WD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
!..................................................................................................................................

   type(WD_InitInputType),       intent(in   ) :: InitInp       !< Input data for initialization routine
   type(WD_InputType),           intent(  out) :: u             !< An initial guess for the input; input mesh must be defined
   type(WD_ParameterType),       intent(  out) :: p             !< Parameters
   type(WD_ContinuousStateType), intent(  out) :: x             !< Initial continuous states
   type(WD_DiscreteStateType),   intent(  out) :: xd            !< Initial discrete states
   type(WD_ConstraintStateType), intent(  out) :: z             !< Initial guess of the constraint states
   type(WD_OtherStateType),      intent(  out) :: OtherState    !< Initial other states
   type(WD_OutputType),          intent(  out) :: y             !< Initial system outputs (outputs are not calculated;
                                                                !!   only the output mesh is initialized)
   type(WD_MiscVarType),         intent(  out) :: m             !< Initial misc/optimization variables
   real(DbKi),                   intent(inout) :: interval      !< Coupling interval in seconds: the rate that
                                                                !!   (1) WD_UpdateStates() is called in loose coupling &
                                                                !!   (2) WD_UpdateDiscState() is called in tight coupling.
                                                                !!   Input is the suggested time from the glue code;
                                                                !!   Output is the actual coupling interval that will be used
                                                                !!   by the glue code.
   type(WD_InitOutputType),      intent(  out) :: InitOut       !< Output for initialization routine
   integer(IntKi),               intent(  out) :: errStat       !< Error status of the operation
   character(*),                 intent(  out) :: errMsg        !< Error message if ErrStat /= ErrID_None
   

      ! Local variables
   integer(IntKi)                              :: i             ! loop counter
   
   integer(IntKi)                              :: errStat2      ! temporary error status of the operation
   character(ErrMsgLen)                        :: errMsg2       ! temporary error message 
      
   !type(WD_InitInputType)                      :: InputFileData ! Data stored in the module's input file
   integer(IntKi)                              :: UnEcho        ! Unit number for the echo file
   
   character(*), parameter                     :: RoutineName = 'WD_Init'
   
   
      ! Initialize variables for this routine

   errStat = ErrID_None
   errMsg  = ""
   UnEcho  = -1

      ! Initialize the NWTC Subroutine Library

   call NWTC_Init( EchoLibVer=.FALSE. )

      ! Display the module information

   call DispNVD( WD_Ver )
   
   
   
   p%OutFileRoot  = TRIM(InitInp%RootName)//'.WD'
   
   
         
      
      ! Validate the initialization inputs
   call ValidateInitInputData( interval, InitInp, InitInp%InputFileData, ErrStat2, ErrMsg2 )
      call SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName ) 
      if (ErrStat >= AbortErrLev) then
         call Cleanup()
         return
      end if
      
      !............................................................................................
      ! Define parameters
      !............................................................................................
      
   
      
      ! set the rest of the parameters
   p%DT            = interval         
   p%NumPlanes     = InitInp%InputFileData%NumPlanes   
   p%NumRadii      = InitInp%InputFileData%NumRadii    
   p%dr            = InitInp%InputFileData%dr               
   p%C_NearWake    = InitInp%InputFileData%C_NearWake  
   p%C_vAmb_DMin   = InitInp%InputFileData%C_vAmb_DMin 
   p%C_vAmb_DMax   = InitInp%InputFileData%C_vAmb_DMax 
   p%C_vAmb_FMin   = InitInp%InputFileData%C_vAmb_FMin 
   p%C_vAmb_Exp    = InitInp%InputFileData%C_vAmb_Exp  
   p%C_vShr_DMin   = InitInp%InputFileData%C_vShr_DMin 
   p%C_vShr_DMax   = InitInp%InputFileData%C_vShr_DMax 
   p%C_vShr_FMin   = InitInp%InputFileData%C_vShr_FMin 
   p%C_vShr_Exp    = InitInp%InputFileData%C_vShr_Exp  
   p%k_vAmb        = InitInp%InputFileData%k_vAmb      
   p%k_vShr        = InitInp%InputFileData%k_vShr      
   p%Mod_WakeDiam  = InitInp%InputFileData%Mod_WakeDiam
   p%C_WakeDiam    = InitInp%InputFileData%C_WakeDiam  
   
   allocate( p%r(0:p%NumRadii-1),stat=errStat)
   
   do i = 0,p%NumRadii-1
      p%r(i)       = p%dr*i     
   end do
   
   p%filtParam     = exp(-2*pi*p%dt*InitInp%InputFileData%f_c)
   
      !............................................................................................
      ! Define and initialize inputs here 
      !............................................................................................
   
   allocate( u%V_plane       (3,0:p%NumPlanes-1),stat=errStat)
   allocate( u%Ct_azavg      (  0:p%NumRadii-1 ),stat=errStat)
   
   
   
   

         
      
      !............................................................................................
      ! Define outputs here
      !............................................................................................

   
   
      !............................................................................................
      ! Initialize states and misc vars : Note these are not the correct initializations because
      ! that would require valid input data, which we do not have here.  Instead we will check for
      ! an firstPass flag on the miscVars and if it is false we will properly initialize these state
      ! in CalcOutput or UpdateStates, as necessary.
      !............................................................................................
      
   allocate ( xd%xhat_plane       (3, 0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%p_plane          (3, 0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%Vx_wind_disk_filt(0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%x_plane          (0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%TI_amb_filt      (0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%D_rotor_filt     (0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%Ct_azavg_filt    (0:p%NumRadii-1) , STAT=ErrStat2 )
   allocate ( xd%Vx_wake     (0:p%NumRadii-1,0:p%NumPlanes-1) , STAT=ErrStat2 )
   allocate ( xd%Vr_wake     (0:p%NumRadii-1,0:p%NumPlanes-1) , STAT=ErrStat2 )
   

   xd%xhat_plane          = 0.0_ReKi
   xd%p_plane             = 0.0_ReKi
   xd%x_plane             = 0.0_ReKi
   xd%Vx_wake             = 0.0_ReKi
   xd%Vr_wake             = 0.0_ReKi
   xd%Vx_wind_disk_filt   = 0.0_ReKi
   xd%TI_amb_filt         = 0.0_ReKi
   xd%D_rotor_filt        = 0.0_ReKi
   xd%Vx_rel_disk_filt    = 0.0_ReKi
   xd%Ct_azavg_filt       = 0.0_ReKi
   m%firstPass            = .true.     
   
   
      !............................................................................................
      ! Define initialization output here
      !............................................................................................
   
   allocate ( y%xhat_plane(3,0:p%NumPlanes-1), STAT=ErrStat2 )
   allocate ( y%p_plane   (3,0:p%NumPlanes-1), STAT=ErrStat2 )
   allocate ( y%Vx_wake   (0:p%NumRadii-1,0:p%NumPlanes-1), STAT=ErrStat2 )
   allocate ( y%Vr_wake   (0:p%NumRadii-1,0:p%NumPlanes-1), STAT=ErrStat2 )
   allocate ( y%D_wake    (0:p%NumPlanes-1), STAT=ErrStat2 )
   
   y%xhat_plane = 0.0_Reki
   y%p_plane    = 0.0_Reki
   y%Vx_wake    = 0.0_Reki
   y%Vr_wake    = 0.0_Reki
   y%D_wake     = 0.0_Reki
   
   call Cleanup() 
      
contains
   subroutine Cleanup()

     ! CALL WD_DestroyInputFile( InputFileData, ErrStat2, ErrMsg2 )
      IF ( UnEcho > 0 ) CLOSE( UnEcho )
      
   end subroutine Cleanup

end subroutine WD_Init

!----------------------------------------------------------------------------------------------------------------------------------
!> This routine is called at the end of the simulation.
subroutine WD_End( u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
!..................................................................................................................................

      type(WD_InputType),           intent(inout)  :: u           !< System inputs
      type(WD_ParameterType),       intent(inout)  :: p           !< Parameters
      type(WD_ContinuousStateType), intent(inout)  :: x           !< Continuous states
      type(WD_DiscreteStateType),   intent(inout)  :: xd          !< Discrete states
      type(WD_ConstraintStateType), intent(inout)  :: z           !< Constraint states
      type(WD_OtherStateType),      intent(inout)  :: OtherState  !< Other states
      type(WD_OutputType),          intent(inout)  :: y           !< System outputs
      type(WD_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
      integer(IntKi),               intent(  out)  :: ErrStat     !< Error status of the operation
      character(*),                 intent(  out)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None



         ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = ""


         ! Place any last minute operations or calculations here:


         ! Close files here:



         ! Destroy the input data:

      call WD_DestroyInput( u, ErrStat, ErrMsg )


         ! Destroy the parameter data:

      call WD_DestroyParam( p, ErrStat, ErrMsg )


         ! Destroy the state data:

      call WD_DestroyContState(   x,           ErrStat, ErrMsg )
      call WD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
      call WD_DestroyConstrState( z,           ErrStat, ErrMsg )
      call WD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )
      call WD_DestroyMisc(        m,           ErrStat, ErrMsg ) 

         ! Destroy the output data:

      call WD_DestroyOutput( y, ErrStat, ErrMsg )




end subroutine WD_End
!----------------------------------------------------------------------------------------------------------------------------------
!> Loose coupling routine for solving for constraint states, integrating continuous states, and updating discrete and other states.
!! Continuous, constraint, discrete, and other states are updated for t + Interval
subroutine WD_UpdateStates( t, n, u, p, x, xd, z, OtherState, m, errStat, errMsg )
!..................................................................................................................................

   real(DbKi),                     intent(in   ) :: t          !< Current simulation time in seconds
   integer(IntKi),                 intent(in   ) :: n          !< Current simulation time step n = 0,1,...
   type(WD_InputType),             intent(inout) :: u       !< Inputs at utimes (out only for mesh record-keeping in ExtrapInterp routine)
  ! real(DbKi),                     intent(in   ) :: utimes   !< Times associated with u(:), in seconds
   type(WD_ParameterType),         intent(in   ) :: p          !< Parameters
   type(WD_ContinuousStateType),   intent(inout) :: x          !< Input: Continuous states at t;
                                                               !!   Output: Continuous states at t + Interval
   type(WD_DiscreteStateType),     intent(inout) :: xd         !< Input: Discrete states at t;
                                                               !!   Output: Discrete states at t  + Interval
   type(WD_ConstraintStateType),   intent(inout) :: z          !< Input: Constraint states at t;
                                                               !!   Output: Constraint states at t+dt
   type(WD_OtherStateType),        intent(inout) :: OtherState !< Input: Other states at t;
                                                               !!   Output: Other states at t+dt
   type(WD_MiscVarType),           intent(inout) :: m          !< Misc/optimization variables
   integer(IntKi),                 intent(  out) :: errStat    !< Error status of the operation
   character(*),                   intent(  out) :: errMsg     !< Error message if ErrStat /= ErrID_None

   ! local variables
   type(WD_InputType)                           :: uInterp     ! Interpolated/Extrapolated input
   integer(intKi)                               :: errStat2          ! temporary Error status
   character(ErrMsgLen)                         :: errMsg2           ! temporary Error message
   character(*), parameter                      :: RoutineName = 'WD_UpdateStates'
   real(ReKi)                                   :: lstar, dx, Vx_wake_min, r_wake, V_planeDT(3), a_interp
   real(ReKi), allocatable                      :: dvdr(:)  , dvtdr(:), vt(:), a(:), b(:), c(:), d(:)  
   integer(intKi)                               :: i,j, ILo
   ErrStat = ErrID_None
   ErrMsg  = ""
     
   ILo = 0
   
      ! Check if we are fully initialized
   if ( m%firstPass ) then
      call InitStatesWithInputs(p%NumPlanes, p%NumRadii, u , xd, errStat, errMsg)
      m%firstPass = .false.
      return
   end if
   
   
   !call WD_CopyInput( u, uInterp, MESH_NEWCOPY, errStat2, errMsg2)
   !   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
   !   if (ErrStat >= AbortErrLev) then
   !      call Cleanup()
   !      return
   !   end if

      ! set values of u(2) from inputs interpolated at t+dt:
  ! call WD_Input_ExtrapInterp(u,utimes,uInterp,t+p%DT, errStat2, errMsg2)
  !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !call SetInputs(p, uInterp, m, 2, errStat2, errMsg2)      
   !   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
      
      ! set values of m%BEMT_u from inputs (uInterp) interpolated at t:
      ! I'm doing this second in case we want the other misc vars at t as before, but I don't think it matters      
  ! call WD_Input_ExtrapInterp(u,utimes,uInterp, t, errStat2, errMsg2)
  !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)

   !call SetInputs(p, uInterp, m, 1, errStat2, errMsg2)      
   !   call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
  
      
   !======================================================================
   !
   ! Near-wake Calculations
   !
   !======================================================================
      ! a uses the previous value of xd%Ct_azavg_filt so this needs to be before the 
   
   
      ! we probably want these to be miscvars to avoid the allocation per timestep
   allocate ( dvdr(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate ( dvtdr(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate (   vt(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate (    a(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate (    b(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate (    c(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   allocate (    d(0:p%NumRadii-1 ) , STAT=ErrStat2 )
   
      !IF ( ErrStat2 /= 0 )  THEN
      !   CALL SetErrStat( ErrID_Fatal, 'Error allocating memory for the WaveTime array.', ErrStat, ErrMsg, 'WAMIT_Init')
      !   CALL Cleanup()
      !   RETURN
      !END IF
   
   !======================================================================
   !
   ! Wake Advection, Deflection, and Meandering Calculations
   !
   !======================================================================
   
      ! We are going to update Vx_Wake
      ! The quantities in these loops are all at time [n], so we need to compute prior to updating the states to [n+1]
   do i = p%NumPlanes-1, 1, -1  
      
      lstar = WakeDiam( p%Mod_WakeDiam, p%numRadii, p%dr, p%r, xd%Vx_wake(i-1,:), xd%Vx_wind_disk_filt(i-1), xd%D_rotor_filt(i-1), p%C_WakeDiam ) 
      ! The following two quantities need to by for the time increments:
      !           [n+1]             [n]
      !dx      = xd%x_plane(i) - xd%x_plane(i-1)
      ! This is equivalent to
      
      V_planeDT(1)            =  u%V_plane   (1,i-1)*p%DT
      V_planeDT(2)            =  u%V_plane   (2,i-1)*p%DT
      V_planeDT(3)            =  u%V_plane   (3,i-1)*p%DT
      dx = dot_product(xd%xhat_plane(:,i-1),V_planeDT)
      
      Vx_wake_min = 9e9_ReKi
      do j = 0,p%NumRadii-1
         Vx_wake_min = min(Vx_wake_min, xd%Vx_wake(i-1,j))
      end do
         
      do j = 0,p%NumRadii-1      
         if ( j == 0 ) then
          dvdr(j) =   0.0_ReKi
         elseif (j <= p%NumRadii-2) then
            dvdr(j) = ( xd%Vx_wake(i-1,j+1) - xd%Vx_wake(i-1,j-1) ) / (2*p%dr)
         else
            dvdr(j) = - xd%Vx_wake(i-1,j-1)  / (2*p%dr)
         end if
            ! All of the following states are at [n] 
         vt(j) =  EddyFilter(xd%x_plane(i-1),xd%D_rotor_filt(i-1), p%C_vAmb_DMin, p%C_vAmb_DMax, p%C_vAmb_FMin, p%C_vAmb_Exp) * p%k_vAmb * xd%TI_amb_filt(i-1) * xd%Vx_wind_disk_filt(i-1) * xd%D_rotor_filt(i-1)/2.0_ReKi &
                - EddyFilter(xd%x_plane(i-1),xd%D_rotor_filt(i-1), p%C_vShr_DMin, p%C_vShr_DMax, p%C_vShr_FMin, p%C_vShr_Exp) * p%k_vShr * max( lstar**2*abs(dvdr(j)) , lstar*(xd%Vx_wind_disk_filt(i-1) + Vx_wake_min ) )
                                                                     
      end do
      
         ! All of the a,b,c,d vectors use states at time increment [n]
         ! These need to be inside another radial loop because dvtdr depends on the j+1 and j-1 indices of vt()
      
      dvtdr(0) = 0.0_ReKi
      a(0)     = 0.0_ReKi
      b(0)     = p%dr * ( xd%Vx_wind_disk_filt(i-1) + xd%Vx_wake(i-1,0)  ) / dx + vt(0)/p%dr
      c(0)     = -vt(0)/p%dr
      c(p%NumRadii-1) = 0.0_ReKi
      d(0)     = (p%dr * (xd%Vx_wind_disk_filt(i-1) + xd%Vx_wake(i-1,0)) / dx - vt(0)/p%dr  ) * xd%Vx_wake(i,0) + ( vt(0)/p%dr ) * xd%Vx_wake(i,1) 
      
      do j = p%NumRadii-1, 1, -1 
         
         if (j <= p%NumRadii-2) then
            dvtdr(j) = ( vt(j+1) - vt(j-1) ) / (2*p%dr)
            c(j) = real(j,ReKi)*xd%Vr_wake(i-1,j)/4.0 - (1+2*real(j,ReKi))*vt(j)/(4.0*p%dr) - real(j,ReKi)*dvtdr(j)/4.0
            d(j) =    ( real(j,ReKi)*xd%Vr_wake(i-1,j)/4.0 - (1-2*real(j,ReKi))*vt(j)/(4.0*p%dr) - real(j,ReKi)*dvtdr(j)/4.0) * xd%Vx_wake(i,j-1) &
                    + ( p%r(j)*( xd%Vx_wind_disk_filt(i-1) + xd%Vx_wake(i-1,j)  )/dx -  real(j,ReKi)*vt(j)/p%dr  ) * xd%Vx_wake(i,j) &
                    + (-real(j,ReKi)*xd%Vr_wake(i-1,j)/4.0 + (1+2*real(j,ReKi))*vt(j)/(4.0*p%dr) + real(j,ReKi)*dvtdr(j)/4.0 ) * xd%Vx_wake(i,j+1)
             
         else
            dvtdr(j) = 0.0_ReKi
            d(j) = ( real(j,ReKi)*xd%Vr_wake(i-1,j)/4.0 - (1-2*real(j,ReKi))*vt(j)/(4.0*p%dr) - real(j,ReKi)*dvtdr(j)/4.0) * xd%Vx_wake(i,j-1) &
                    + ( p%r(j)*( xd%Vx_wind_disk_filt(i-1) + xd%Vx_wake(i-1,j)  )/dx -  real(j,ReKi)*vt(j)/p%dr  ) * xd%Vx_wake(i,j) 
                    
         end if  
         
         a(j) = -j*xd%Vr_wake(i-1,j)/4.0_ReKi + (1.0_ReKi-2.0*real(j,ReKi))*vt(j)/(4.0*p%dr) + real(j,ReKi)*dvtdr(j)/4.0 
         b(j) = p%r(j) * ( xd%Vx_wind_disk_filt(i-1) + xd%Vx_wake(i-1,j)  ) / dx + real(j,ReKi)*vt(j)/p%dr
         
        
      end do ! j = 1,p%NumRadii-1
   
         ! Update these states to [n+1]
      
      
      xd%x_plane        (  i) = xd%x_plane   (  i-1) + dx !dot_product(xd%xhat_plane(:,i-1),V_planeDT)
                              
      xd%p_plane        (:,i) = xd%p_plane   (:,i-1) + V_planeDT 
      xd%xhat_plane     (:,i) = xd%xhat_plane(:,i-1)
      xd%Vx_wind_disk_filt(i) = xd%Vx_wind_disk_filt(i-1)
      xd%TI_amb_filt      (i) = xd%TI_amb_filt(i-1)
      xd%D_rotor_filt     (i) = xd%D_rotor_filt(i-1)
      
         ! Update Vx_wake to [n+1]
      call ThomasAlgorithm(p%NumRadii, a, b, c, d, xd%Vx_wake(i,:))
      
   end do ! i = 1,p%NumPlanes-1 
   
      ! Downstream planes 
   do j = 1,p%NumRadii-1
      do i=p%NumPlanes-1,1,-1
            !  Vr_wake is for the                [n+1]        ,          [n]                  , and           [n]              increments
         xd%Vr_wake(i,j) = real(j-1,ReKi)*(  xd%Vr_wake(i,j-1)     + xd%Vr_wake(i-1,j-1) )/real(j,ReKi) - xd%Vr_wake(i-1,j)    &
            !  Vx_wake is for the                           [n+1]      ,      [n+1]        ,      [n]          , and  [n]             increments             
                           + real(2*j-1,ReKi)*p%dr * (  xd%Vx_wake(i,j) + xd%Vx_wake(i,j-1) - xd%Vx_wake(i-1,j) - xd%Vx_wake(i-1,j-1)    )  / dx
      end do  
   end do
  
      ! Update states at disk-plane to [n+1] 
      
   xd%xhat_plane     (:,0) =  xd%xhat_plane(:,0)*p%filtParam + u%xhat_disk(:)*(1-p%filtParam)  ! 2-step calculation for xhat_plane at disk
      xd%xhat_plane     (:,0) =  xd%xhat_plane(:,0) / norm2( xd%xhat_plane(:,0) ) 
   xd%p_plane        (:,0) =  xd%p_plane   (:,0)*p%filtParam + u%p_hub(:)*(1-p%filtParam)
   xd%Vx_wind_disk_filt(0) =  xd%Vx_wind_disk_filt(0)*p%filtParam + u%Vx_wind_disk*(1-p%filtParam)   
   xd%TI_amb_filt      (0) =  xd%TI_amb_filt(0)*p%filtParam + u%TI_amb*(1-p%filtParam)   
   xd%D_rotor_filt     (0) =  xd%D_rotor_filt(0)*p%filtParam + u%D_rotor*(1-p%filtParam)   
   xd%Vx_rel_disk_filt     =  xd%Vx_rel_disk_filt*p%filtParam + u%Vx_rel_disk*(1-p%filtParam)   
   
      !  filtered, azimuthally-averaged Ct values at each radial station
   
   do i=p%NumRadii-1,0,-1  
      
      xd%Ct_azavg_filt (i) =  xd%Ct_azavg_filt(i)*p%filtParam + u%Ct_azavg(i)*(1-p%filtParam) 
      
         ! compute a using the [n+1] values of Ct_azavg_filt
      if ( xd%Ct_azavg_filt(i) > 2.0 ) then
         ! THROW ERROR because we are in the prop-brake region
      else if ( xd%Ct_azavg_filt(i) >= 24./25. ) then
         a(i) = (2 + 3*sqrt(14*xd%Ct_azavg_filt(i)-12))/14.0
      else
         a(i) = 0.5 - 0.5*sqrt(1-xd%Ct_azavg_filt(i))
      end if
      
      if ( a(i) >= 1.0_ReKi / p%C_NearWake ) then
         ! THROW ERROR
      end if
      
      
   end do
   
   ! NOTE: We need another loop over NumRadii because we need the complete a() vector so that we can interpolate into a() using the value of r_wake
   
      ! Use the [n+1] version of xd%Vx_rel_disk_filt to determine the [n+1] version of Vx_wake(0,:)
   xd%Vx_wake(0,0) = - xd%Vx_rel_disk_filt*p%C_Nearwake*a(0)
   r_wake = 0.0_ReKi
   
   do i=1, p%NumRadii-1     
      r_wake = sqrt(r_wake**2 + p%dr*( ((1 - a(i))*p%r(i)) / (1-p%C_NearWake*a(i)) + ((1 - a(i-1))*p%r(i-1)) / (1-p%C_NearWake*a(i-1)) ) )
         ! given r_wake and a at p%dr increments, find value of a(r_wake) using interpolation 
      a_interp = InterpBin( r_wake, p%r, a, ILo, p%NumRadii ) !( XVal, XAry, YAry, ILo, AryLen )
      
         !                 [n+1] 
      xd%Vx_wake(0,i) = xd%Vx_rel_disk_filt*p%C_NearWake*a_interp
      xd%Vr_wake(0,i) = 0.0_ReKi
   end do
   
   
   
   call Cleanup()
   
contains
   subroutine Cleanup()
      call WD_DestroyInput( uInterp, errStat2, errMsg2)
   end subroutine Cleanup
end subroutine WD_UpdateStates
!----------------------------------------------------------------------------------------------------------------------------------
!> Routine for computing outputs, used in both loose and tight coupling.
!! This subroutine is used to compute the output channels (motions and loads) and place them in the WriteOutput() array.
!! The descriptions of the output channels are not given here. Please see the included OutListParameters.xlsx sheet for
!! for a complete description of each output parameter.
subroutine WD_CalcOutput( t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
! NOTE: no matter how many channels are selected for output, all of the outputs are calcalated
! All of the calculated output channels are placed into the m%AllOuts(:), while the channels selected for outputs are
! placed in the y%WriteOutput(:) array.
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )  :: t           !< Current simulation time in seconds
   TYPE(WD_InputType),           INTENT(IN   )  :: u           !< Inputs at Time t
   TYPE(WD_ParameterType),       INTENT(IN   )  :: p           !< Parameters
   TYPE(WD_ContinuousStateType), INTENT(IN   )  :: x           !< Continuous states at t
   TYPE(WD_DiscreteStateType),   INTENT(IN   )  :: xd          !< Discrete states at t
   TYPE(WD_ConstraintStateType), INTENT(IN   )  :: z           !< Constraint states at t
   TYPE(WD_OtherStateType),      INTENT(IN   )  :: OtherState  !< Other states at t
   TYPE(WD_OutputType),          INTENT(INOUT)  :: y           !< Outputs computed at t (Input only so that mesh con-
                                                               !!   nectivity information does not have to be recalculated)
   type(WD_MiscVarType),         intent(inout)  :: m           !< Misc/optimization variables
   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   integer, parameter                           :: indx = 1  
   integer(intKi)                               :: i, j
   integer(intKi)                               :: ErrStat2
   character(ErrMsgLen)                         :: ErrMsg2
   character(*), parameter                      :: RoutineName = 'WD_CalcOutput'
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

      ! Check if we are fully initialized
   if ( m%firstPass ) then
     ! call InitStatesWithInputs(p%NumPlanes, p%NumRadii, u, xd, errStat, errMsg)
      
      do i = 0, p%NumPlanes - 1
         y%p_plane     (:,i) = u%p_hub(:)
         y%xhat_plane  (:,i) = u%xhat_disk(:)
         do j = 0, p%NumRadii - 1
            y%Vx_wake(j,i) = 0.0_ReKi
            y%Vr_wake(j,i) = 0.0_ReKi
         end do          
      end do
              
      do i = 0, p%NumPlanes - 1
         y%D_wake(i)  =  WakeDiam( p%Mod_WakeDiam, p%NumRadii, p%dr, p%r, y%Vx_wake(:,i), u%Vx_wind_disk, u%D_rotor, p%C_WakeDiam)
      end do
      return
   else
      y%p_plane    = xd%p_plane
      y%xhat_plane = xd%xhat_plane
      y%Vx_wake    = xd%Vx_wake
      y%Vr_wake    = xd%Vr_wake
      do i = 0, p%NumPlanes - 1
         y%D_wake(i)  =  WakeDiam( p%Mod_WakeDiam, p%NumRadii, p%dr, p%r, xd%Vx_wake(:,i), xd%Vx_wind_disk_filt(i), xd%D_rotor_filt(i), p%C_WakeDiam)
      end do
   end if
   
  ! call SetInputs(p, u, m, indx, errStat2, errMsg2)      
  !    call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)
            
 
   !-------------------------------------------------------   
   !     get values to output to file:  
   !-------------------------------------------------------   
 !  if (p%NumOuts > 0) then
#ifdef DBG_OUTS
  !    call Calc_WriteDbgOutput( p, u, m, y, ErrStat2, ErrMsg2 ) 
#else
 !     call Calc_WriteOutput( p, u, m, y, indx, ErrStat2, ErrMsg2 )   
#endif   
!      call SetErrStat(ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName)      
   
      !...............................................................................................................................   
      ! Place the selected output channels into the WriteOutput(:) array with the proper sign:
      !...............................................................................................................................   

  !    do i = 1,p%NumOuts  ! Loop through all selected output channels
#ifdef DBG_OUTS
  !       y%WriteOutput(i) = m%AllOuts( i )
#else
  !       y%WriteOutput(i) = p%OutParam(i)%SignM * m%AllOuts( p%OutParam(i)%Indx )
#endif

  !    end do             ! i - All selected output channels
      
 !  end if
   
   
   
end subroutine WD_CalcOutput
!----------------------------------------------------------------------------------------------------------------------------------
!> Tight coupling routine for solving for the residual of the constraint state equations
subroutine WD_CalcConstrStateResidual( Time, u, p, x, xd, z, OtherState, m, z_residual, ErrStat, ErrMsg )
!..................................................................................................................................

   REAL(DbKi),                   INTENT(IN   )   :: Time        !< Current simulation time in seconds
   TYPE(WD_InputType),           INTENT(IN   )   :: u           !< Inputs at Time
   TYPE(WD_ParameterType),       INTENT(IN   )   :: p           !< Parameters
   TYPE(WD_ContinuousStateType), INTENT(IN   )   :: x           !< Continuous states at Time
   TYPE(WD_DiscreteStateType),   INTENT(IN   )   :: xd          !< Discrete states at Time
   TYPE(WD_ConstraintStateType), INTENT(IN   )   :: z           !< Constraint states at Time (possibly a guess)
   TYPE(WD_OtherStateType),      INTENT(IN   )   :: OtherState  !< Other states at Time
   TYPE(WD_MiscVarType),         INTENT(INOUT)   :: m           !< Misc/optimization variables
   TYPE(WD_ConstraintStateType), INTENT(INOUT)   :: Z_residual  !< Residual of the constraint state equations using
                                                                !!     the input values described above
   INTEGER(IntKi),               INTENT(  OUT)   :: ErrStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)   :: ErrMsg      !< Error message if ErrStat /= ErrID_None


   
      ! Local variables   
   integer, parameter                            :: indx = 1  
   integer(intKi)                                :: ErrStat2
   character(ErrMsgLen)                          :: ErrMsg2
   character(*), parameter                       :: RoutineName = 'WD_CalcConstrStateResidual'
   
   
   
   ErrStat = ErrID_None
   ErrMsg  = ""

         
   
   
end subroutine WD_CalcConstrStateResidual

subroutine InitStatesWithInputs(numPlanes, numRadii, u, xd, errStat, errMsg)

   integer(IntKi),               intent(in   )   :: numPlanes
   integer(IntKi),               intent(in   )   :: numRadii
   TYPE(WD_InputType),           INTENT(IN   )   :: u           !< Inputs at Time
   TYPE(WD_DiscreteStateType),   INTENT(INOUT)   :: xd          !< Discrete states at Time
   INTEGER(IntKi),               INTENT(  OUT)   :: errStat     !< Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)   :: errMsg      !< Error message if ErrStat /= ErrID_None

   integer(IntKi) :: i
   
   ! Note, all of these states will have been set to zero in the WD_Init routine
   
     
   do i = 0, numPlanes - 1
      xd%p_plane     (:,i) = u%p_hub(:)
      xd%xhat_plane  (:,i) = u%xhat_disk(:)
     
         ! Vx_wake and Vr_wake are already initialized to zero so we do not need to do that here.
      !xd%Vx_wake(j,i)  = 0.0_ReKi
      !xd%Vr_wake(j,i)  = 0.0_ReKi
   end do
   
      ! Only need to set the 0 index because everythin else is already initialized to zero.
   xd%Vx_wind_disk_filt(0) = u%Vx_wind_disk
   xd%TI_amb_filt      (0) = u%TI_amb
   xd%D_rotor_filt     (0) = u%D_rotor
   xd%Vx_rel_disk_filt     = u%Vx_rel_disk 
   
   do i = 0, numRadii - 1
      xd%Ct_azavg_filt (i) = u%Ct_azavg(i) 
   end do
  
end subroutine InitStatesWithInputs
   
!----------------------------------------------------------------------------------------------------------------------------------
!> This routine validates the inputs from the WakeDynamics input files.
SUBROUTINE ValidateInitInputData( DT, InitInp, InputFileData, ErrStat, ErrMsg )
!..................................................................................................................................
      
      ! Passed variables:
   real(DbKi),               intent(in   )  :: DT                                !< requested simulation time step size (s)
   type(WD_InitInputType),   intent(in   )  :: InitInp                           !< Input data for initialization routine
   type(WD_InputFileType),   intent(in)     :: InputFileData                     !< All the data in the WakeDynamics input file
   integer(IntKi),           intent(out)    :: ErrStat                           !< Error status
   character(*),             intent(out)    :: ErrMsg                            !< Error message

   
      ! local variables
   integer(IntKi)                           :: k                                 ! Blade number
   integer(IntKi)                           :: j                                 ! node number
   character(*), parameter                  :: RoutineName = 'ValidateInitInputData'
   
   ErrStat = ErrID_None
   ErrMsg  = ""
   
   if (len(trim(InitInp%RootName)) == 0) call SetErrStat ( ErrID_Fatal, 'Rootname must contain at least one character.', ErrStat, ErrMsg, RoutineName )  
   
   
   !if (NumBl > MaxBl .or. NumBl < 1) call SetErrStat( ErrID_Fatal, 'Number of blades must be between 1 and '//trim(num2lstr(MaxBl))//'.', ErrSTat, ErrMsg, RoutineName )
   if (  DT          <=  0.0)  call SetErrStat ( ErrID_Fatal, 'DT must be greater than zero.', ErrStat, ErrMsg, RoutineName )  
   if (  InputFileData%NumPlanes   <   2  )  call SetErrStat ( ErrID_Fatal, 'Number of wake planes must be greater than one.', ErrSTat, ErrMsg, RoutineName )
   if (  InputFileData%NumRadii    <   2  )  call SetErrStat ( ErrID_Fatal, 'Number of radii in the radial finite-difference grid must be greater than one.', ErrSTat, ErrMsg, RoutineName )
   if (  InputFileData%dr          <=  0.0)  call SetErrStat ( ErrID_Fatal, 'dr must be greater than zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%f_c         <=  0.0)  call SetErrStat ( ErrID_Fatal, 'f_c must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_NearWake  <= -1.0)  call SetErrStat ( ErrID_Fatal, 'C_NearWake must be greater than -1.0.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%k_vAmb      <   0.0)  call SetErrStat ( ErrID_Fatal, 'k_vAmb must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%k_vShr      <   0.0)  call SetErrStat ( ErrID_Fatal, 'k_vShr must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vAmb_DMin <   0.0)  call SetErrStat ( ErrID_Fatal, 'C_vAmb_DMin must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vAmb_DMax <=  InputFileData%C_vAmb_DMin)  call SetErrStat ( ErrID_Fatal, 'C_vAmb_DMax must be greater than C_vAmb_DMin.', ErrStat, ErrMsg, RoutineName ) 
   if ( (InputFileData%C_vAmb_FMin <   0.0)  .or. (InputFileData%C_vAmb_FMin > 1.0) ) call SetErrStat ( ErrID_Fatal, 'C_vAmb_FMin must be greater than or equal to zero and less than or equal to 1.0.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vAmb_Exp  <=  0.0)  call SetErrStat ( ErrID_Fatal, 'C_vAmb_Exp must be greater than zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vShr_DMin <   0.0)  call SetErrStat ( ErrID_Fatal, 'C_vShr_DMin must be greater than or equal to zero.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vShr_DMax <=  InputFileData%C_vShr_DMin)  call SetErrStat ( ErrID_Fatal, 'C_vShr_DMax must be greater than C_vShr_DMin.', ErrStat, ErrMsg, RoutineName ) 
   if ( (InputFileData%C_vShr_FMin <   0.0)  .or. (InputFileData%C_vShr_FMin > 1.0) ) call SetErrStat ( ErrID_Fatal, 'C_vShr_FMin must be greater than or equal to zero and less than or equal to 1.0.', ErrStat, ErrMsg, RoutineName ) 
   if (  InputFileData%C_vShr_Exp <=   0.0)  call SetErrStat ( ErrID_Fatal, 'C_vShr_Exp must be greater than zero.', ErrStat, ErrMsg, RoutineName ) 
   if (.not. ((InputFileData%Mod_WakeDiam == 1) .or. (InputFileData%Mod_WakeDiam == 2) .or. (InputFileData%Mod_WakeDiam == 3) ) ) call SetErrStat ( ErrID_Fatal, 'Mod_WakeDiam must be equal to 1, 2, or 3.', ErrStat, ErrMsg, RoutineName ) 
   if ( (InputFileData%C_WakeDiam <=   0.0)  .or. (InputFileData%C_WakeDiam >= 1.0) ) call SetErrStat ( ErrID_Fatal, 'C_WakeDiam must be greater than zero and less than 1.0.', ErrStat, ErrMsg, RoutineName ) 
   
END SUBROUTINE ValidateInitInputData



!=======================================================================
! Unit Tests
!=======================================================================

subroutine WD_TEST_Init_BadData(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(WD_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(WD_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(WD_ParameterType)       :: p             !< Parameters
   type(WD_ContinuousStateType) :: x             !< Initial continuous states
   type(WD_DiscreteStateType)   :: xd            !< Initial discrete states
   type(WD_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(WD_OtherStateType)      :: OtherState    !< Initial other states
   type(WD_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(WD_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(WD_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
                                                                        
   

   
      ! Set up the initialization inputs
   
   InitInp%RootName       = ''   
   interval               = 0.0_DbKi
   InitInp%InputFileData%NumPlanes      = 0
   InitInp%InputFileData%NumRadii       = 0
   InitInp%InputFileData%dr             = 0.0_ReKi
   InitInp%InputFileData%f_c            = 0.0
   InitInp%InputFileData%C_HWkDfl_O     = -2.9
   InitInp%InputFileData%C_HWkDfl_OY    = -.24*D2R
   InitInp%InputFileData%C_HWkDfl_x     = -0.0054
   InitInp%InputFileData%C_HWkDfl_xY    = 0.00039*D2R
   InitInp%InputFileData%C_NearWake     = -1.001
   InitInp%InputFileData%C_vAmb_DMin    = -0.01
   InitInp%InputFileData%C_vAmb_DMax    = -3
   InitInp%InputFileData%C_vAmb_FMin    = 2
   InitInp%InputFileData%C_vAmb_Exp     = 0
   InitInp%InputFileData%C_vShr_DMin    = -.1
   InitInp%InputFileData%C_vShr_DMax    = -3
   InitInp%InputFileData%C_vShr_FMin    = -1
   InitInp%InputFileData%C_vShr_Exp     = -1
   InitInp%InputFileData%k_vAmb         = -2
   InitInp%InputFileData%k_vShr         = -2
   InitInp%InputFileData%Mod_WakeDiam   = 0
   InitInp%InputFileData%C_WakeDiam     = 0
   
   
   call WD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, Interval, InitOut, ErrStat, ErrMsg )
   
   return
   
end subroutine WD_TEST_Init_BadData
subroutine WD_TEST_SetGoodInitInpData(interval, InitInp)
   real(DbKi)            , intent(out)       :: interval
   type(WD_InitInputType), intent(out)       :: InitInp       !< Input data for initialization routine

   InitInp%RootName       = 'GoodData'   
   interval               = 0.1_DbKi
   InitInp%InputFileData%NumPlanes      = 2
   InitInp%InputFileData%NumRadii       = 2
   InitInp%InputFileData%dr             = 0.1_ReKi
   InitInp%InputFileData%f_c            = 0.03333333333333
   InitInp%InputFileData%C_HWkDfl_O     = -2.9
   InitInp%InputFileData%C_HWkDfl_OY    = -.24*D2R
   InitInp%InputFileData%C_HWkDfl_x     = -0.0054
   InitInp%InputFileData%C_HWkDfl_xY    = 0.00039*D2R
   InitInp%InputFileData%C_NearWake     = 2
   InitInp%InputFileData%C_vAmb_DMin    = 0
   InitInp%InputFileData%C_vAmb_DMax    = 2
   InitInp%InputFileData%C_vAmb_FMin    = 0
   InitInp%InputFileData%C_vAmb_Exp     = 1
   InitInp%InputFileData%C_vShr_DMin    = 2
   InitInp%InputFileData%C_vShr_DMax    = 11
   InitInp%InputFileData%C_vShr_FMin    = .035
   InitInp%InputFileData%C_vShr_Exp     = .4
   InitInp%InputFileData%k_vAmb         = .07
   InitInp%InputFileData%k_vShr         = .0178
   InitInp%InputFileData%Mod_WakeDiam   = 1
   InitInp%InputFileData%C_WakeDiam     = .95

end subroutine WD_TEST_SetGoodInitInpData


subroutine WD_TEST_Init_GoodData(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(WD_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(WD_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(WD_ParameterType)       :: p             !< Parameters
   type(WD_ContinuousStateType) :: x             !< Initial continuous states
   type(WD_DiscreteStateType)   :: xd            !< Initial discrete states
   type(WD_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(WD_OtherStateType)      :: OtherState    !< Initial other states
   type(WD_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(WD_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(WD_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
                                                                        
   

   
      ! Set up the initialization inputs
   call WD_TEST_SetGoodInitInpData(interval, InitInp)

   call WD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, interval, InitOut, ErrStat, ErrMsg )
   
   return
   
end subroutine WD_TEST_Init_GoodData

subroutine WD_TEST_UpdateStates(errStat, errMsg)

   integer(IntKi),           intent(out)    :: errStat                           !< Error status
   character(*),             intent(out)    :: errMsg                            !< Error message


   type(WD_InitInputType)       :: InitInp       !< Input data for initialization routine
   type(WD_InputType)           :: u             !< An initial guess for the input; input mesh must be defined
   type(WD_ParameterType)       :: p             !< Parameters
   type(WD_ContinuousStateType) :: x             !< Initial continuous states
   type(WD_DiscreteStateType)   :: xd            !< Initial discrete states
   type(WD_ConstraintStateType) :: z             !< Initial guess of the constraint states
   type(WD_OtherStateType)      :: OtherState    !< Initial other states
   type(WD_OutputType)          :: y             !< Initial system outputs (outputs are not calculated;
                                                 !!   only the output mesh is initialized)
   type(WD_MiscVarType)         :: m             !< Initial misc/optimization variables
   real(DbKi)                   :: interval      !< Coupling interval in seconds: the rate that
   
   type(WD_InitOutputType)      :: initOut                         !< Input data for initialization routine
   
   integer(IntKi)  :: i
   real(DbKi) :: t
   
      ! Set up the initialization inputs
   call WD_TEST_SetGoodInitInpData(interval, InitInp)

   call WD_Init( InitInp, u, p, x, xd, z, OtherState, y, m, interval, InitOut, ErrStat, ErrMsg )
   
      ! Set up the inputs
   u%xhat_disk(1) = 1.0_ReKi
   u%p_hub(1)     = 10.0_ReKi
   u%p_hub(2)     = 20.0_ReKi
   u%p_hub(3)     = 50.0_ReKi
   do i=0,p%NumPlanes-1
      u%V_plane(1,i)  = 5 + .5*i
      u%V_plane(2,i)  = 1 + .1*i
      u%V_plane(3,i)  = 1 + .3*i
   end do
   u%Vx_wind_disk = 6.0
   u%TI_amb   = .15    
   u%D_rotor  = 30.0      
   u%Vx_rel_disk = 5.0
   do i=0,p%NumRadii-1
      u%Ct_azavg(i) = 0.2 + 0.05*i
   end do
   t = 0.0_DbKi
   call WD_UpdateStates(t, 1, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   
   
   call WD_CalcOutput(t, u, p, x, xd, z, OtherState, y, m, ErrStat, ErrMsg )
   
      ! Verify that xd and y are the same
   
   if (errStat == ErrID_None) then
      call WD_UpdateStates(0.0_DbKi, 1, u, p, x, xd, z, OtherState, m, errStat, errMsg )
   end if
   
   return


end subroutine WD_TEST_UpdateStates
END MODULE WakeDynamics
