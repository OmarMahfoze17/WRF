 MODULE module_cpl_mui
   
#ifdef key_cpp_mui
   !!======================================================================
   !!                    ***  MODULE cpl_oasis  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS3-MCT
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!   cpl_mui_init     : initialization of coupled mode communication
   !!   cpl_mui_define   : definition of grid and fields
   !!   cpl_mui_snd      : send out fields in coupled mode
   !!   cpl_mui_rcv      : receive fields in coupled mode
   !!   cpl_mui_finaliz : finalize the coupled mode communication
   !!----------------------------------------------------------------------
   USE module_domain     
   USE module_driver_constants, ONLY : max_domains, max_cplfld, max_extdomains
   ! USE mod_oasis              ! OASIS3-MCT module

   use iso_c_binding
   use mui_3d_f
   use mui_general_f

   ! USE module_driver_constants
   ! USE module_nesting
   ! USE module_configure
   ! USE module_timing
   ! USE module_utility

   IMPLICIT NONE
   PRIVATE

   TYPE ::   FLD_CPL                                  ! Coupling field information
      CHARACTER(len = 64) ::   clname                    ! Name of the coupling field, jpeighty defined in oasis
      INTEGER             ::   nid                       ! Id of the field
#if ( RWORDSIZE == 8 )
      REAL        , POINTER, DIMENSION(:,:) ::   dbl2d   ! 2d array to store received field 
#else
      REAL(kind=8), POINTER, DIMENSION(:,:) ::   dbl2d   ! 2d array to store received field 
#endif
   END TYPE FLD_CPL
   TYPE(FLD_CPL), DIMENSION(max_domains,max_extdomains,max_cplfld) :: srcv, ssnd   ! Coupling fields informations
   INTEGER :: ndm_comm                                                         ! MPI communicator between the computing nodes
   INTEGER :: ncomp_id                                                         ! id returned by oasis_init_comp
   INTEGER :: nlevdbg  = 1                                                     ! verbosity level
   INTEGER :: nlevdbg2 = 10                                                    ! verbosity level
   CHARACTER(len = 256) :: cltxt                                               ! messages or debug string
   integer :: mui_initiated=0;
   !! Routine accessibility
   PUBLIC cpl_mui_init
   PUBLIC cpl_mui_def_dmcomm
   PUBLIC cpl_mui_define
   PUBLIC cpl_mui_toreceive
   PUBLIC cpl_mui_tosend
   PUBLIC cpl_mui_snd
   PUBLIC cpl_mui_rcv
   PUBLIC cpl_mui_finalize
   PUBLIC cpl_mui_abort
   !Define pointer type to collcet uniface pointers for multi-domian function
   ! type ptr_typ_3d
   !    type(c_ptr) :: ptr
   ! end type ptr_typ_3d

   !Create an allocatable array to collect MUI uniface pointers with the type of
   ! ptr_typ_3d for multi-domian function
   ! type(ptr_typ_3d), target, save, allocatable :: uniface_pointers_3d(:)
   character(len=1024) :: domainName, interfaceName
   character(len=1024) :: sptlSmpType,tmpSmpType
   character(len=1024) :: arg_interface_count
   integer(c_int) :: interface_count=1,nForget
   character(:), allocatable, target :: interfaces3d(:)
   character(:), allocatable :: domain3d
   integer :: interfaceDirection(200), interfaceLocation(200)
   integer,allocatable,dimension(:) :: ifsDir,ifsLoc
   integer :: groupNumb  ! number of data groups to be pushed. The group is defined as a box 
   integer :: groupVort(200,6)
   !The domain is alwyes start at (0,0,0), but there are cases where 
   !data is needed to be pushed/fitched to/from domain that is shifted from 0  
   real(c_double) :: dataOrgShft(3) = 0.0_c_double

   type(c_ptr), target :: uniface_3d=c_null_ptr
   type(c_ptr), target :: spatial_sampler=c_null_ptr,spatial_sync=c_null_ptr
   type(c_ptr), target :: temporal_sampler=c_null_ptr,temporal_sync=c_null_ptr
   real(c_double) :: tolerance=1e-37_c_double
   ! integer(c_int) :: MUI_COMM_WORLD
   
   integer :: MUI_x1_ID,MUI_xn_ID,MUI_y1_ID,MUI_yn_ID,MUI_z1_ID,MUI_zn_ID   
   integer :: MUIBC_ID(6),sendReceiveMode,pointIndxMode,coupledDomain
   integer :: smartSendReceive=0,iSync=0,nSyncAhead=1

   ! RBF filter varaibles
   real(c_double) :: rSpatialSamp = 1.0_c_double,sigmaSpatialSamp = 1.0_c_double
   real(c_double) :: rTempSamp = 1.0_c_double,sigmaTempSamp = 1.0_c_double
   real(c_double) :: tempMeanSampLower,tempMeanSampUpper    !! used by the  temporal mean sampler 
   real           :: lat_ref,lon_ref,x1,x2,y1,y2,z1,z2,boxThickness
   INTEGER :: MUI_COMM_WORLD




#if ( defined( DM_PARALLEL ) && ( ! defined( STUBMPI ) ) )
   INCLUDE 'mpif.h'               ! only for MPI_COMM_NULL
#else
   INTEGER :: MPI_COMM_NULL = -1  ! define a fake (and not used) MPI_COMM_NULL, so it is compiling 
#endif

CONTAINS

   SUBROUTINE cpl_mui_init( kl_comm ) 
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_mui_init  ***
      !!
      !! ** Purpose :   Initialize coupled mode communication for WRF
      !!--------------------------------------------------------------------
      INTEGER, INTENT(OUT) :: kl_comm       ! local communicator of the model
      character(len=1024) :: numberSuffix
      integer :: ifsIndx , ierr
      character(len=50) :: myMessage,arg
      TYPE (domain)          :: grid
      NAMELIST/MUICoupling/domainName,interfaceName,interface_count, &
      interfaceDirection, interfacelocation,MUIBC_ID, groupNumb, &
      groupVort, dataOrgShft,tolerance,sendReceiveMode,smartSendReceive,iSync,nSyncAhead

      NAMELIST/MUISamplers/ sptlSmpType,tmpSmpType,rSpatialSamp,sigmaSpatialSamp,rTempSamp, &
      sigmaTempSamp,nForget,tempMeanSampLower,tempMeanSampUpper

      NAMELIST/MUIWRFDomainChar/coupledDomain, lat_ref,lon_ref,pointIndxMode
      NAMELIST/MUISendDomain/x1,x2,y1,y2,z1,z2,boxThickness
      !!!!! TEST TEST TEST !!!!

      ! Get the first command-line argument
      call get_command_argument(1, arg)

      ! Trim the argument and check if it is equal to "mui_coupled"
      if (trim(arg) == "mui_coupled") then
         CALL wrf_message(" MUI initialisation started ")
         ! IF (.FALSE.) kl_comm = -1 ! to avoid compilation warning
         ! CALL mpi_init ( ierr )
         call mui_mpi_split_by_app_f(kl_comm)
         MUI_COMM_WORLD = kl_comm
         open(10, file="MUI_CPL.input")
         read(10, nml=MUICoupling); rewind(10)  
         read(10, nml=MUISamplers); rewind(10)
         read(10, nml=MUIWRFDomainChar); rewind(10)

         read(10, nml=MUISendDomain); rewind(10)

         allocate(character(len_trim(interfaceName)+5) :: interfaces3d(interface_count))
         !For multi-domain function, "uniface_pointers_1d" should be used to collect the array of
         ! MUI uniface pointers. It is decleared in the MUI FORTRAN wrapper.
         allocate(uniface_pointers_3d(interface_count),&
                  ifsDir(interface_count),ifsLoc(interface_count))

         do ifsIndx = 1, interface_count
            !Generate character type of number suffix
            if (ifsIndx < 10) then
               write (numberSuffix, "(I1)") ifsIndx
            else if ((ifsIndx < 100) .and. (ifsIndx > 9)) then
               write (numberSuffix, "(I2)") ifsIndx
            else if ((ifsIndx < 1000) .and. (ifsIndx > 99)) then
               write (numberSuffix, "(I3)") ifsIndx
            else
               write (numberSuffix, "(I4)") ifsIndx

            endif

            !Create and collect interface names
            interfaces3d(ifsIndx) = trim(interfaceName) ! // "_" // trim(numberSuffix)
            ifsDir(ifsIndx)=interfaceDirection(ifsIndx)
            ifsLoc(ifsIndx)=interfaceLocation(ifsIndx)
         end do 

         call create_and_get_uniface_multi_3d_f(uniface_pointers_3d, trim(domainName), interfaces3d, interface_count)         

         call mui_create_sampler_exact_3d_f(spatial_sampler, tolerance)
         call mui_create_temporal_sampler_exact_3d_f(temporal_sampler, tolerance)
         
         mui_initiated=1   ! Set it to 1 if MUI initialisation is successful. 
         CALL wrf_message( " ************* MUI initialisation finished *****************")

      else
         CALL mpi_init ( ierr )
         kl_comm = MPI_COMM_WORLD
      end if

   END SUBROUTINE cpl_mui_init


   SUBROUTINE cpl_mui_def_dmcomm( kdm_comm ) 
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_mui_def_dmcomm  ***
      !!
      !! ** Purpose :   define ndm_comm: the MPI communicator between the computing nodes
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdm_comm       ! computing nodes communicator
      !!--------------------------------------------------------------------
      IF (.FALSE.) WRITE(*,*) kdm_comm ! to avoid compilation warning
      ! ndm_comm = kdm_comm   ! store it to used it in cpl_mui_define

      ! WRITE(cltxt,*) 'cpl_mui_def_dmcomm : ', kdm_comm
      ! CALL wrf_debug(nlevdbg, cltxt)
      ! CALL wrf_debug(nlevdbg, '~~~~~~~~~~~~~~~~~~~~~~~')

   END SUBROUTINE cpl_mui_def_dmcomm


   SUBROUTINE cpl_mui_define( cdsndname, cdrcvname, pgrid )
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_mui_define  ***
      !!
      !! ** Purpose :   Define grid and coupling field information for WRF
      !!--------------------------------------------------------------------
      CHARACTER(*), INTENT(IN), DIMENSION(:,:,:)  :: cdsndname, cdrcvname   ! coupling field names
      TYPE(domain), INTENT(IN), OPTIONAL, POINTER :: pgrid                  ! grid structure
      !
      IF (.FALSE.) WRITE(*,*) cdsndname, cdrcvname, pgrid%id  ! to avoid compilation warning

    
      
   END SUBROUTINE cpl_mui_define
   
   
   SUBROUTINE cpl_mui_enddef()
      !!-------------------------------------------------------------------
      !!             ***  ROUTINE cpl_mui_enddef  ***
      !!
      !! ** Purpose :   tells to OASIS that exchanged field definition is finished
      !!--------------------------------------------------------------------
      INTEGER :: ierror         ! return error code
      ! CALL oasis_enddef(ierror)
      ! IF( ierror /= OASIS_Ok )   CALL cpl_mui_abort( 'cpl_mui_define', 'Failure in oasis_enddef')
      
   END SUBROUTINE cpl_mui_enddef
   
   
   FUNCTION cpl_mui_toreceive( kdomwrf, kdomext, kfldid )
      !!-------------------------------------------------------------------
      !!             ***  FUNCTION cpl_mui_toreceive  ***
      !!
      !! ** Purpose :   send back a logical to tell if a variable is received or not
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kdomext   ! external model domain index
      INTEGER, INTENT(IN) :: kfldid    ! field index
      !
      LOGICAL :: cpl_mui_toreceive
      !!--------------------------------------------------------------------
      
      ! cpl_mui_toreceive = srcv(kdomwrf,kdomext,kfldid)%nid /= -1

   END FUNCTION cpl_mui_toreceive


   FUNCTION cpl_mui_tosend( kdomwrf, kdomext, kfldid )
      !!-------------------------------------------------------------------
      !!             ***  FUNCTION cpl_mui_tosend  ***
      !!
      !! ** Purpose :   send back a logical to tell if a variable is tosend or not
      !!--------------------------------------------------------------------
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kdomext   ! external model domain index
      INTEGER, INTENT(IN) :: kfldid    ! field index
      !
      LOGICAL :: cpl_mui_tosend
      !!--------------------------------------------------------------------
      
      ! cpl_mui_tosend = ssnd(kdomwrf,kdomext,kfldid)%nid /= -1

   END FUNCTION cpl_mui_tosend


   ! SUBROUTINE cpl_mui_snd( kdomwrf, kdomext, kfldid, ksec, pdata )
   SUBROUTINE cpl_mui_snd(  grid, max_edom                &
      &                     , ids,ide,jds,jde,kds,kde &
      &                     , ims,ime,jms,jme,kms,kme &
      &                     , ips,ipe,jps,jpe,kps,kpe, nsecrun)
      TYPE(domain), INTENT(IN), POINTER :: grid
      INTEGER,      INTENT(IN)          :: max_edom    ! max number of external model domains
      INTEGER,      INTENT(IN)          :: ids,ide,jds,jde,kds,kde
      INTEGER,      INTENT(IN)          :: ims,ime,jms,jme,kms,kme
      INTEGER,      INTENT(IN)          :: ips,ipe,jps,jpe,kps,kpe
      !
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: cplsnd
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: u_uo
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: v_vo
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: wspd
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: taut
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: x_u,y_u,x_v,y_v,x,y
      REAL(8)                                :: point_x,point_y,point_z,temp(3),TempValue
      REAL :: tStartMUISync,tEndMUISync
      INTEGER :: icnt, i,j,k,sendig_peers_size
      INTEGER :: ifldid,nrank,ierr,nproc
      INTEGER :: nsecrun
      integer :: unit,synchronised=1
      TYPE(WRFU_TimeInterval) :: simTimeFromStart,time_from_start,time_from_simulationStart
      TYPE(WRFU_Time) :: start_time,current_time,simulationStartTime
      REAL(8) :: pushTime,startTime,xMin  ,yMin , zMin, xMax ,yMax, zMax, veryLarg = 10e34
      LOGICAL,DIMENSION(max_edom) :: lltosend
      character(len=50) :: fileName
      call mui_create_temporal_sampler_exact_3d_f(temporal_sync, tolerance)
      call mui_create_sampler_exact_3d_f(spatial_sync, tolerance)
      !!--------------------------------------------------------------------

      !!--------------------------------------------------------------------
      ! note : data order is XZY , If it is different then this routine should be adjusted accordingly
      !!--------------------------------------------------------------------
   !    tmpTimeInterval  = domain_get_current_time ( grid ) - domain_get_sim_start_time ( grid )
   ! tmpTimeInterval2 = domain_get_current_time ( grid ) - domain_get_start_time ( grid )
   ! curr_secs  = real_time(tmpTimeInterval)
   ! curr_secs2 = real_time(tmpTimeInterval2)
      if (grid%id == coupledDomain ) then
         write (*,*) "WRF commetind the test "
         call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, real(1,c_double), real(2,c_double),real(0,c_double), 0.1234_c_double)
         call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, real(1.0,c_double))
         ! write(*,*) "check ref ll" , grid%latc_loc,grid%lat_ll_u,grid%scm_lat,grid%cen_lat, &
         ! grid%truelat1,grid%truelat2,grid%moad_cen_lat,grid%pole_lat
         if (pointIndxMode ==0) then  !! define data location at lat-long
            x_u(ips:ipe,jps:jpe) = grid%xlong_u(ips:ipe,jps:jpe)
            y_u(ips:ipe,jps:jpe) = grid%xlat_u(ips:ipe,jps:jpe)

            x_v(ips:ipe,jps:jpe) = grid%xlong_v(ips:ipe,jps:jpe)
            y_v(ips:ipe,jps:jpe) = grid%xlat_v(ips:ipe,jps:jpe)

            x(ips:ipe,jps:jpe) = grid%xlong(ips:ipe,jps:jpe)
            y(ips:ipe,jps:jpe) = grid%xlat(ips:ipe,jps:jpe)
         else if (pointIndxMode ==1) then !! define data location at X-Y
            call ll_to_xy(grid%xlat_u(ips:ipe,jps:jpe),grid%xlong_u(ips:ipe,jps:jpe), &
            lat_ref  ,lon_ref  ,ips,ipe,jps,jpe,x_u(ips:ipe,jps:jpe),y_u(ips:ipe,jps:jpe)) 
            call ll_to_xy(grid%xlat_v(ips:ipe,jps:jpe),grid%xlong_v(ips:ipe,jps:jpe), &
            lat_ref  ,lon_ref  ,ips,ipe,jps,jpe,x_v(ips:ipe,jps:jpe),y_v(ips:ipe,jps:jpe)) 
            call ll_to_xy(grid%xlat(ips:ipe,jps:jpe),grid%xlong(ips:ipe,jps:jpe), &
            lat_ref  ,lon_ref  ,ips,ipe,jps,jpe,x(ips:ipe,jps:jpe),y(ips:ipe,jps:jpe)) 
         else if (pointIndxMode ==3) then !! define data location at indix i-j
            ! x_u(ips:ipe,jps:jpe) = [ips:ipe,jps:jpe]

         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!! Time is given at date, hr, minutes, sec
         current_time = domain_get_current_time ( grid )  
         start_time   = domain_get_start_time ( grid )! start time of the current run
         simulationStartTime = domain_get_sim_start_time( grid ) ! start time of the simulation including restart time

         time_from_start = current_time - start_time ! get the value as interval

         

         ! simStartTime=domain_get_sim_start_time( grid )
         simTimeFromStart = domain_get_time_since_sim_start(grid)


         ! startTime = real(simStartTime%basetime%s,8) + &
         ! real(simStartTime%basetime%sn,8) / real(simStartTime%basetime%sd,8) - grid%dt

         !! get the value in seconds
         startTime = real(time_from_start%basetime%s,8) + &
         real(time_from_start%basetime%sn,8) / real(time_from_start%basetime%sd,8) - grid%dt

         pushTime = real(simTimeFromStart%basetime%s,8) + &
         real(simTimeFromStart%basetime%sn,8) / real(simTimeFromStart%basetime%sd,8) - grid%dt

         
         !    ll_to_xy(lat        ,lon         ,lat_ref,lon_ref,ips,ipe,jps,jpe,x  ,y  )
         ! write(*,*) "===============================0========================================"
         ! write(*,*) simTimeFromStart%basetime%s, simTimeFromStart%basetime%sn , simTimeFromStart%basetime%sd
         ! write(*,*) "indx d" , ids,ide,jds,jde,kds,kde 
         ! write(*,*) "indx m" , ims,ime,jms,jme,kms,kme ,lat_ref  ,lon_ref 
         ! write(*,*) "indx p" ,ips,ipe-1,jps,jpe-1,kps,kpe-1
         ! write(*,*) "Grid lat long", grid%xlat_u(ips,jps),grid%xlat_u(ipe-1,jpe-1),grid%xlong_u(ips,jps),grid%xlong_u(ipe-1,jpe-1)
         ! write(*,*) "Grid X y",pushTime, x_u(ips,jps),y_u(ips,jps),grid%z(ips,1,jps)
         ! write(*,*) "Grid X y",  pushTime, x_u(ipe,jpe),y_u(ipe,jpe),grid%z(ipe,1,jpe)

         

         if (abs(startTime - grid%dt) .lt. 0.00001 .and. smartSendReceive .ne. 0) then 
            write(*,*) "================================================"
            write(*,*) "MUI: initialising smart send and recieve ...... "
            xMin = veryLarg  ;yMin = veryLarg ; zMin = -100000
            xMax = -veryLarg ;yMax = -veryLarg; zMax =  100000
            do i=ips,ipe            
               do j=jps,jpe
                  if (x_u(i,j) .ge. x1 .and. x_u(i,j).le. x2 .and. y_v(i,j).ge.y1 .and. y_v(i,j).le.y2 ) then

                     if(x_u(i,j) .lt. xMin) xMin = x_u(i,j)
                     if(y_v(i,j) .lt. yMin) yMin = y_v(i,j)

                     if(x_u(i,j) .gt. xMax) xMax = x_u(i,j)
                     if(y_v(i,j) .gt. yMax) yMax = y_v(i,j)
                  endif
               enddo
            enddo


            if (veryLarg - xMin < 10) xMin = 0.0+veryLarg
            if (veryLarg - yMin < 10) yMin = 0.0+veryLarg
            if (veryLarg - abs(xMax) < 10) xMax = 0.0+veryLarg
            if (veryLarg - abs(yMax) < 10) yMax = 0.0+veryLarg


            write(*,*) "Smart Send Box is ", xMin,xMax, yMin,yMax 

            call mui_announce_recv_span_3d_box_f(uniface_pointers_3d(1)%ptr, &
            xMin,zMin,yMin,xMax, zMax,yMax,real(0.0,c_double),real(1000000.0, c_double),synchronised)
      
            call mui_announce_send_span_3d_box_f(uniface_pointers_3d(1)%ptr, &
            xMin,zMin, yMin,xMax,zMax,yMax &
            ,real(0.0,c_double),real(1000000.0, c_double),synchronised)
           
            call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, real(0.0,8))
            call mui_sendig_peers_size_3d_f(uniface_pointers_3d(1)%ptr, sendig_peers_size)

            write(*,*) "MUI: Finished initialising smart send and recieve. sendig_peers_size ", sendig_peers_size
            write(*,*) "================================================"
         endif


         ! note : data order is XZY , If it is different then this routine should be adjusted accordingly
         do i=ips,ipe            
            do j=jps,jpe

               ! temp = 0 
               ! point_x = real(x(i,j),c_double) ! real(i,c_double)! grid%xlong_u(i,j)
               ! point_y = real(y(i,j),c_double) ! real(j,c_double) !grid%xlat_u(i,j)
               ! call mui_push_3d_vector(uniface_pointers_3d(1)%ptr, "velocity"//c_null_char, point_x, &
               ! REAL(0.0,8),point_y, temp)

               ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, x_u(i,j), &
               ! grid%z(i,k,j),y_u(i,j), 0.0_c_double)
               ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, x_v(i,j), &
               ! grid%z(i,k,j),y_v(i,j), 0.0_c_double)
               ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, x_u(i,j), &
               ! grid%z_at_w(i,k,j),y_w(i,j), 0.0_c_double)
               if (x_u(i,j) .ge. x1 .and. x_u(i,j).le. x2 .and. y_v(i,j).ge.y1 .and. y_v(i,j).le.y2 ) then 
                  if (x_u(i,j) .le. x1+boxThickness .or. x_u(i,j).ge. x2-boxThickness &
                   .or. y_v(i,j).le.y1+boxThickness .or. y_v(i,j).ge.y2-boxThickness ) then 

                     do k=kps,kpe !kps,kpe      
                        
                        if (grid%z(i,k,j).ge.z1 .and. grid%z(i,k,j).le.z2 ) then
                        
                           ! NOTE : this part is writen to fit XCOMPACT3D coordinates 
                           ! where ux is streamwise, uy is vertical , and uz is spanwise velocities. 

                           ! point_x = real(x_u(i,j),c_double) 
                           ! point_y = real(y_u(i,j),c_double) 
                           ! point_z = real(grid%z(i,k,j),c_double)
                           ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "ux"//c_null_char, point_x, &
                           ! point_z,point_y, REAL(grid%u_2(i,k,j),c_double))

                           ! point_x = real(x_v(i,j),c_double) 
                           ! point_y = real(y_v(i,j),c_double) 
                           ! point_z = real(grid%z(i,k,j),c_double)
                           ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uz"//c_null_char, point_x, &
                           ! point_z,point_y, REAL(grid%v_2(i,k,j),c_double))

                           ! point_x = real(x(i,j),c_double) 
                           ! point_y = real(y(i,j),c_double)
                           ! point_z = real(grid%z_at_w(i,k,j),c_double)
                           ! call mui_push_3d_f(uniface_pointers_3d(1)%ptr, "uy"//c_null_char, point_x, &
                           ! point_z,point_y, REAL(grid%w_2(i,k,j),c_double))


                           ! As the grid is fully stagared, we interpolate at w
                           temp(1) =  REAL((grid%u_2(i-1,k,j-1)+grid%u_2(i,k,j-1) &
                                           +grid%u_2(i-1,k,j)+grid%u_2(i,k,j))/4.0,c_double)
                           temp(2) =  REAL(grid%w_2(i,k,j),c_double)
                           temp(3) =  REAL((grid%v_2(i-1,k,j-1)+grid%v_2(i,k,j-1) &
                           +grid%v_2(i-1,k,j)+grid%v_2(i,k,j))/4.0,c_double) 


                           ! temp(1) =  REAL(grid%u_2(i,k,j),c_double)
                           ! temp(2) =  REAL(grid%w_2(i,k,j),c_double)
                           ! temp(3) =  REAL(grid%v_2(i,k,j),c_double) 

                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           point_x = real(x(i,j),c_double) 
                           point_y = real(y(i,j),c_double) 
                           point_z = real(grid%z(i,k,j),c_double)

                           call mui_push_3d_vector(uniface_pointers_3d(1)%ptr, "velocity"//c_null_char, point_x,  &
                           point_z,point_y, temp)   
                        endif

                     enddo
                  endif
               endif
            enddo
         enddo

         

         call mui_commit_3d_f(uniface_pointers_3d(1)%ptr, pushTime)  

         if (iSync ==1 .and. startTime >= nSyncAhead*grid%dt ) then 
            call cpu_time(tStartMUISync)
            write(*,*) "---------- waiting for xcom to sync ----------------"
            call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(1)%ptr, &
             "tempvalue"//c_null_char, real(0.0,8),real(0.0,8), real(0.0,8) &
             , pushTime - nSyncAhead*grid%dt  , spatial_sync, temporal_sync, TempValue)
             call cpu_time(tEndMUISync)
             write (*,*) "WRF sync time is ", tEndMUISync-tStartMUISync

             write(*,*) "---------- Done waiting for xcom to sync ----------------"
            
            ! call mui_barrier_3d_f(uniface_pointers_3d(1)%ptr, real((t-dt*nForget)* timeScale,c_double))
   
            ! call mui_forget_upper_3d_f(uniface_pointers_3d(1)%ptr, &
            !       real((t-dt*nSyncAhead)* timeScale,c_double),reset_log)
   
   
         endif


         write (*,*) pushTime,grid%dt,grid%time_step,"WRF Commited Velocity field to MUI interface" 
         ! if (iSync ==1 .and. itrCount > nSyncAhead) then 
         !    call mui_fetch_exact_exact_3d_f(uniface_pointers_3d(1)%ptr, &
         !    "tempvalue"//c_null_char, REAL(0.0,8),REAL(0.0,8), REAL(0.0,8) &
         !    , real(T+subT-nSyncAhead,c_double),spatial_sampler, temporal_sampler, TempValue)
         ! endif
      endif

      
    END SUBROUTINE cpl_mui_snd

    SUBROUTINE ll_to_xy(lat,lon,lat_ref,lon_ref,ips,ipe,jps,jpe,x,y)
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE ll_to_xy  ***
      !!
      !! ** Purpose : Convert lat_long to X, Y coordinates
      !!---------------------------------------------------------------------
      INTEGER,              INTENT(IN   ) :: ips,ipe,jps,jpe
      REAL, DIMENSION( ips:ipe, jps:jpe ) :: lat, lon
      REAL, DIMENSION( ips:ipe, jps:jpe )::  x, y
      INTEGER                             :: i,j
      REAL                 :: lon_ref,lat_ref,lat1,lon1,dlon,dlat,lat_ref_rad,lon_ref_rad
      REAL                 :: lon2,lat2
      Real                 :: pi=3.141592653589793238462643383279502884197


      !! Convert lat and Long to radian 

      lat_ref_rad = lat_ref*pi/180
      lon_ref_rad = lon_ref*pi/180
      
      do i = ips,ipe
         do j = jps,jpe

            lat1 = lat(i,j)*pi/180.0;
            lat2 =  lat_ref_rad

            lon1 = lon_ref_rad;
            lon2 = lon_ref_rad;

            dlon = lon2 - lon1 ;
            dlat = lat2 - lat1;

            Y(i, j) = 2.0 * 6371000.0 * ASIN(SQRT(SIN(dlat / 2.0)**2 + COS(lat1) &
            * COS(lat2) * SIN(dlon / 2.0)**2)) * SIGN(1.0,-dlat)
           
            lat1 = lat_ref_rad;
            lon1 = lon(i,j)*pi/180.0;
            
            dlon = lon_ref_rad - lon1 ;
            dlat = lat_ref_rad - lat1;

            X(i, j) = 2.0 * 6371000.0 *ASIN(SQRT(SIN(dlat / 2.0)**2 + COS(lat1) * COS(lat_ref_rad) &
            * SIN(dlon / 2.0)**2)) *  SIGN(1.0,-dlon)
         enddo
      enddo


 
    
    END SUBROUTINE  ll_to_xy


   SUBROUTINE cpl_mui_rcv( kdomwrf, kdomext, kfldid, ksec, pcplrcv )

      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_mui_rcv  ***
      !!
      !! ** Purpose : - At each coupling time-step, this routine check if it is the good time
      !!                to receive field from the coupler
      !!----------------------------------------------------------------------
      INTEGER,              INTENT(IN   ) :: kdomwrf   ! wrf domain index
      INTEGER,              INTENT(IN   ) :: kdomext   ! external model domain index
      INTEGER,              INTENT(IN   ) :: kfldid    ! variable index
      INTEGER,              INTENT(IN   ) :: ksec      ! number of seconds since the last restart
      REAL, DIMENSION(:,:), INTENT(  OUT) :: pcplrcv   ! output data
      !!
!       INTEGER :: info              ! OASIS3 info argument
!       LOGICAL :: llaction          ! true if we received data from the coupler
!       !!--------------------------------------------------------------------
!       !
!       WRITE(cltxt,*) 'OASIS_GET in: kdomwrf, kdomext, kfldid, name, ksec', &
!          kdomwrf, kdomext, kfldid, ' ', TRIM(srcv(kdomwrf,kdomext,kfldid)%clname), ksec
!       CALL wrf_debug(nlevdbg, cltxt)

!       CALL oasis_get( srcv(kdomwrf,kdomext,kfldid)%nid, ksec, srcv(kdomwrf,kdomext,kfldid)%dbl2d, info )
! #if ( RWORDSIZE == 8 )
!       pcplrcv(:,:) =      srcv(kdomwrf,kdomext,kfldid)%dbl2d
! #else
!       pcplrcv(:,:) = REAL(srcv(kdomwrf,kdomext,kfldid)%dbl2d, kind=4)
! #endif

!       WRITE(cltxt,*) 'OASIS_GET out: info', info   ;   CALL wrf_debug(nlevdbg, cltxt)

!       ! llaction = info == OASIS_Recvd   .OR. info == OASIS_FromRest .OR.   &
!       !    &       info == OASIS_RecvOut .OR. info == OASIS_FromRestOut 

!       WRITE(cltxt,*) "llaction : ", llaction       ;   CALL wrf_debug(nlevdbg, cltxt)

!       IF( llaction ) THEN
!             WRITE(cltxt,*) '****************'                                                   ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) 'oasis_get: Incoming ', TRIM(srcv(kdomwrf,kdomext,kfldid)%clname)    ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) 'oasis_get:    varid ', srcv(kdomwrf,kdomext,kfldid)%nid             ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) 'oasis_get:     ksec ', ksec                                         ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) 'oasis_get:     info ', info                                         ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - shape         ', SHAPE(pcplrcv)                               ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local shape   ', SHAPE(srcv(kdomwrf, kdomext,kfldid)%dbl2d)   ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local minimum ', MINVAL(pcplrcv)                              ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local maximum ', MAXVAL(pcplrcv)                              ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local sum     ', SUM(pcplrcv)                                 ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '****************'                                                   ;   CALL wrf_debug(nlevdbg, cltxt)
!       ELSE
!             WRITE(cltxt,*) '****************'                                                   ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) 'oasis_get: field not received as info = ', info                     ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local minimum ', MINVAL(pcplrcv)                              ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local maximum ', MAXVAL(pcplrcv)                              ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '    - local sum     ', SUM(pcplrcv)                                 ;   CALL wrf_debug(nlevdbg, cltxt)
!             WRITE(cltxt,*) '****************'                                                   ;   CALL wrf_debug(nlevdbg, cltxt)
!       ENDIF

   END SUBROUTINE cpl_mui_rcv


   SUBROUTINE cpl_mui_finalize()
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_mui_finalize  ***
      !!
      !! ** Purpose : - Finalizes the coupling. If MPI_init has not been
      !!      called explicitly before cpl_mui_init it will also close
      !!      MPI communication.
      !!----------------------------------------------------------------------
      INTEGER :: ierror         ! return error code
      INTEGER :: jw,je,jf       ! local loop indicees
      !!--------------------------------------------------------------------
      ! DO jf = 1, max_cplfld
      !    DO je = 1, max_extdomains
      !       DO jw = 1, max_domains
      !          ierror = 0
      !          IF ( ASSOCIATED(srcv(jw,je,jf)%dbl2d) ) DEALLOCATE( srcv(jw,je,jf)%dbl2d, stat = ierror )
      !          IF (ierror > 0) THEN
      !             CALL cpl_mui_abort( 'cpl_mui_finalize', 'Failure in deallocating ')
      !             RETURN
      !          ENDIF
      !       END DO
      !    END DO
      ! END DO
      ! CALL oasis_terminate ( ierror )         

   END SUBROUTINE cpl_mui_finalize


   SUBROUTINE cpl_mui_abort( cdroutine, cdtxt )
      !!---------------------------------------------------------------------
      !!              ***  ROUTINE cpl_mui_abort  ***
      !!
      !! ** Purpose :   abort coupling simulation
      !!----------------------------------------------------------------------
      CHARACTER(*), INTENT(IN) :: cdroutine   ! name of the subroutine calling cpl_mui_abort
      CHARACTER(*), INTENT(IN) :: cdtxt       ! aborting text
      !!--------------------------------------------------------------------

      ! CALL wrf_message( '     ==== ABORTING ====' )
      ! CALL wrf_message( 'cpl_abort called by '//TRIM(cdroutine) )
      ! CALL wrf_message( '   ==> '//TRIM(cdtxt) )
      ! CALL oasis_abort( ncomp_id, cdroutine, cdtxt )
      
   END SUBROUTINE cpl_mui_abort
   
#else
   !!----------------------------------------------------------------------
   !!   Dummy modules just for compilation...
   !!----------------------------------------------------------------------
   USE module_domain, ONLY : domain
   use iso_c_binding
   use mui_3d_f
   use mui_general_f

   IMPLICIT NONE
   PRIVATE
   PUBLIC cpl_mui_init
   PUBLIC cpl_mui_def_dmcomm
   PUBLIC cpl_mui_define
   PUBLIC cpl_mui_toreceive
   PUBLIC cpl_mui_tosend
   PUBLIC cpl_mui_snd
   PUBLIC cpl_mui_rcv
   PUBLIC cpl_mui_finalize
   PUBLIC cpl_mui_abort



CONTAINS

   SUBROUTINE cpl_mui_init( kl_comm ) 
      INTEGER, INTENT(OUT) :: kl_comm       ! local communicator of the model
      IF (.FALSE.) kl_comm = -1 ! to avoid compilation warning
   END SUBROUTINE cpl_mui_init

   SUBROUTINE cpl_mui_def_dmcomm( kdm_comm ) 
      INTEGER, INTENT(IN) :: kdm_comm       ! computing nodes communicator
      IF (.FALSE.) WRITE(*,*) kdm_comm ! to avoid compilation warning
   END SUBROUTINE cpl_mui_def_dmcomm

   SUBROUTINE cpl_mui_define( cdsndname, cdrcvname, pgrid )
      CHARACTER(*), INTENT(IN), DIMENSION(:,:,:)  :: cdsndname, cdrcvname   ! coupling field names
      TYPE(domain), INTENT(IN), OPTIONAL, POINTER :: pgrid                  ! grid structure
      IF (.FALSE.) WRITE(*,*) cdsndname, cdrcvname, pgrid%id  ! to avoid compilation warning
   END SUBROUTINE cpl_mui_define

   FUNCTION cpl_mui_toreceive( kdomwrf, kdomext, kfldid )
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kdomext   ! external model domain index
      INTEGER, INTENT(IN) :: kfldid    ! field index
      LOGICAL :: cpl_mui_toreceive
      IF (.FALSE.) WRITE(*,*) kdomwrf, kdomext, kfldid  ! to avoid compilation warning
      IF (.FALSE.) cpl_mui_toreceive = .false.  ! to avoid compilation warning
   END FUNCTION cpl_mui_toreceive

   FUNCTION cpl_mui_tosend( kdomwrf, kdomext, kfldid )
      INTEGER, INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER, INTENT(IN) :: kdomext   ! external model domain index
      INTEGER, INTENT(IN) :: kfldid    ! field index
      LOGICAL :: cpl_mui_tosend
      IF (.FALSE.) WRITE(*,*) kdomwrf, kdomext, kfldid  ! to avoid compilation warning
      IF (.FALSE.) cpl_mui_tosend = .false.  ! to avoid compilation warning
   END FUNCTION cpl_mui_tosend

   SUBROUTINE cpl_mui_snd( kdomwrf, kdomext, kfldid, ksec, pdata )
      !!----------------------------------------------------------------------
      INTEGER,              INTENT(IN) :: kdomwrf   ! wrf domain index
      INTEGER,              INTENT(IN) :: kdomext   ! external model domain index
      INTEGER,              INTENT(IN) :: kfldid    ! field index
      INTEGER,              INTENT(IN) :: ksec      ! time-step in seconds
      REAL, DIMENSION(:,:), INTENT(IN) :: pdata     ! data to be sent
      IF (.FALSE.) WRITE(*,*) kdomwrf, kdomext, kfldid, ksec, pdata ! to avoid compilation warning
   END SUBROUTINE cpl_mui_snd

   SUBROUTINE cpl_mui_rcv( kdomwrf, kdomext, kfldid, ksec, pcplrcv )
      INTEGER,              INTENT(IN   ) :: kdomwrf   ! wrf domain index
      INTEGER,              INTENT(IN   ) :: kdomext   ! external model domain index
      INTEGER,              INTENT(IN   ) :: kfldid    ! variable index
      INTEGER,              INTENT(IN   ) :: ksec      ! number of seconds since the last restart
      REAL, DIMENSION(:,:), INTENT(  OUT) :: pcplrcv   ! output data
      IF (.FALSE.) WRITE(*,*) kdomwrf, kdomext, kfldid, ksec ! to avoid compilation warning
      IF (.FALSE.) pcplrcv(:,:) = -1. ! to avoid compilation warning
   END SUBROUTINE cpl_mui_rcv

   SUBROUTINE cpl_mui_finalize()
      IF (.FALSE.) WRITE(*,*) 'You should not be there...'
   END SUBROUTINE cpl_mui_finalize

   SUBROUTINE cpl_mui_abort( cdroutine, cdtxt )
      CHARACTER(*), INTENT(IN) :: cdroutine   ! name of the subroutine calling cpl_mui_abort
      CHARACTER(*), INTENT(IN) :: cdtxt       ! aborting text
      IF (.FALSE.) WRITE(*,*) cdroutine, cdtxt   ! to avoid compilation warning
   END SUBROUTINE cpl_mui_abort
#endif

END MODULE module_cpl_mui
