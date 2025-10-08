!the program of cal the grid system      
Program output
  ! Contact Haopeng Wang at haopengwang2022@outlook.com.
  ! For more details of the RBF interpolation used in this code, refer to
  ! Wang et al, 2022. Implicit solar coronal magnetohydrodynamic (MHD) modeling with a low-dissipation 
  ! hybridized AUSM-HLL Riemann solver. ApJ 935(1), 46. doi:10.3847/1538-4357/ac78e0.
  implicit double precision(a-h,o-z)

  !>>for 30Rs_lvl5
  !INTEGER(4),PARAMETER::N_th_p=50, N_ph_p=N_th_p*2, Rindex=72
  !For interface of COCONUT and euhforia
  INTEGER(4),PARAMETER::N_th_p=180, N_ph_p=N_th_p*2, Rindex=72
  INTEGER(4),PARAMETER::Nodstart=1,Nodend=192150,Cellstart=1,Cellend=378880,R_IndMax=74
  !>>for BCindependent6
  !INTEGER(4),PARAMETER::N_th_p=100, N_ph_p=N_th_p*2, Rindex=71
  !For interface of COCONUT and euhforia
  !INTEGER(4),PARAMETER::N_th_p=90, N_ph_p=N_th_p*2, Rindex=71
  !INTEGER(4),PARAMETER::Nodstart=1,Nodend=757908,Cellstart=1,Cellend=1495040,R_IndMax=73
  !<<
  
  integer,parameter :: nvar=6 
  real*8,dimension(:,:,:,:),allocatable :: utotal,utotal_ori
  character (len=60)  :: filename,filename_out
  real*8,dimension(:),  allocatable :: theta11,phi11
  real*8,dimension(:,:,:),  allocatable :: usph,Rusph_interp,Lousph_interp
  
  real*8,dimension(:,:),allocatable :: Nodal_Pos,Cell_Pos,Cell_state
  integer,dimension(:,:),allocatable :: TriPrismIndex,InletTriPrismIndex,OutletTriPrismIndex
  integer :: ComLinStart,ComLinEnd,CurLin,PerNodLayer,PerCellLayer,Rindex_
  !!>> CentralPos(1:3,1:2) corresponds to inner and outer triangle face
  real*8 :: VertexPos(1:3,1:6),CentralPos(1:3,1:2),r
  real*8 :: r_max,temp1,temp2,temp
  double precision :: kB,mu_cor,mH,Rhos
  real*8 :: muT2Gs,PlasmaBeta,BEXP2,dabsB
  real*8 :: radiu(0:98),radiusph_interp(1:17,1:93),HDLSandLDHS(1:4,1:93)
  real*8 :: HD,LS,LD,HS,eps,time,time0,time_begine,Time_cur,cadence
  real*8 :: x,y,z,rho,u,v,w,Bx,By,Bz,p,B_sph(3),phi_L1,phi_L1R,xxx
  real*8 :: Area_1,Area_2,Cross(1:3),AB(1:3),AC(1:3)
  real*8 :: FR_cent,Long0,Long,Longitude_shift,halfLength,rad_temp,minute_drop
  integer :: ii,i,j,k,Seq_,Seq_Start,Seq_End,int_time,index_j,index_k
  integer :: timescale,index_,time_step
  integer :: layer_index(1:5,0:Rindex,0:N_th_p+1,0:N_ph_p+1)
  integer :: mmddhh,year,hh_increased,seconds,second,minute,hours,days,months
  logical :: justnearest,from_uniform,writeBxyzonly,write_CFmesh,write_r
  logical :: unicostheta,from_ori,Time_sequence,writereformated_dat
  logical :: CH_trace,SortorNot,meiridians,pannels,rdial_profile
  logical :: interface_for_euhforia,update,Zerod1AU_only,Mesh_30Rs_lvl5,Mesh_BCindependent6
  character (len=2)  :: Char_mm,Char_dd,Char_hh,Char_minute,Char_ss
  character (len=4)  :: Char_year
  !character(len=20) :: fmt
  character(len=100) :: line
  character(len=:),allocatable :: ctime,ctimeo,ctime_2,central_meridian
  allocate(character(len=10)::ctime,ctimeo,ctime_2)
  allocate(character(len=10)::central_meridian)
    
  Time_sequence= .true. !.false. ! ! choose false when generating cfmesh file
  writereformated_dat = .false.
  !Please define "Time_sequence= .false." before seting "SortorNot = .true."
  SortorNot = .false. !.true.  !
  !-------------------------------------
  from_uniform= .false.  !goto 114
  from_ori = .false.
  unicostheta = .false.
  justnearest = .false.
  !-------------------------------------
  CH_trace= .false.
  meiridians= .false. !.true. ! 
  pannels= .true. !.false. !
  Zerod1AU_only= .true. 
  interface_for_euhforia = .false. !.true. !.true. ! 
  update= .true. !.false. 
  rdial_profile= .false.
  writeBxyzonly= .false. !.true. !
  write_CFmesh = .false. !.true. ! 
  write_r=.true.
  if(Time_sequence .eqv. .true.)then
    write_r=.false.
	write_CFmesh = .false.
  end if
  Mesh_30Rs_lvl5= .true.
  Mesh_BCindependent6= .false.
  !--------------------------------------
  year=2019
  mmddhh=62911
  minute_drop=14.d0
  Long0=300.d0 
  time_begine=0.d0  
  time_step=600
  Seq_Start=3762 !-------------------------------------------------------------?
  if(Time_sequence .eqv. .false.)then
      Seq_End=Seq_Start
  else
      Seq_End=4200 !
  end if

  time0=0.00 !0.32 !------------------------------------------------------------?
  time=time0
  
  pi=4.*atan(1.0)
  gamma=1.667
  Ts = 1.5d6 
  nrho= 1.d8  !/cm^3
  Rhos = nrho*1.d6*1.67d-27     ! Kg/m^3 rho at surface
  VS = 480363.085276 !m/s
  dmu0 = 4.d-7*pi 
  kB = 1.3806503e-23 ! Boltzmann constant in SI
  mu_cor = 1.27 ! Molecular weight in the corona 
  mH = 1.6733e-27 ! kg
  Ps=0.03851 !rhos*VS*VS   Pa
  Bs = 2.2d-4 !dsqrt(dmu0*rhos*VS*VS) Tesla
  muT2Gs=10.d0
  nrt=Rindex
  write(*,*)Bs  
  allocate ( &        
  theta11(0:N_th_p+1),&      
  phi11(0:N_ph_p+1),&
  Cell_Pos(1:3,Cellstart:Cellend)) 

  dth=pi/real(N_th_P)
  dph = dth
  DO i = 0,N_th_p+1
      THETA11(i) = i*dth-dth/2.
  ENDDO
  DO j = 0,N_ph_p+1
      PHI11(j) = j*dph-dph/2. 
  ENDDO  
  int_time=int(time)
  call int2str(int_time,ctime)
  write(*,*)'ctime=',ctime,' len(ctime)=',len(ctime)
  ctime=TRIM(ctime)
  write(*,*)'ctime=',ctime,' len(ctime)=',len(ctime)
  index_=0  
  timescale=100 !1 ! ---------------------------------------------------------------?
  !!>> read un-uniform structured grid sysem
  do Seq_=Seq_Start*timescale,Seq_End*timescale ! debuging
      index_=index_+1
      !if((index_>4) .and. (mod(index_,2) .gt. 0))cycle
	  if(mod(index_-1,time_step) .gt. 0)cycle
         deallocate(central_meridian)
         allocate(character(len=10)::central_meridian)
         allocate ( & 
             utotal_ori(1:11,0:Rindex,0:N_th_p+1,0:N_ph_p+1),& 
             utotal(1:11,0:nrt,0:N_th_p+1,0:N_ph_p+1),&
             Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1),&
             Lousph_interp(1:17,1:nrt,0:N_th_P+1),&
             usph(1:17,0:N_th_p+1,0:N_ph_p+1))   
   
  time=time0*timescale+real(Seq_)
  eps=1.d-3
  int_time=int(time+eps)
  write(*,*)"int_time=",int_time
  deallocate(ctime,ctimeo,ctime_2)
  allocate(character(len=10)::ctime,ctimeo,ctime_2)
  !call int2str(int_time,ctime)
  !call int2strwithonesuffix(int_time,ctime)
  call int2strwithsuffix(int_time,timescale,ctime,ctimeo,ctime_2)
  ctime=TRIM(ctime)
  ctimeo=TRIM(ctimeo)
  ctime_2=TRIM(ctime_2)
  write(*,*)'int_time=',int_time,'ctime=',ctime,' len(ctime)=',len(ctime),'ctimeo=',ctimeo,' len(ctimeo)=',len(ctimeo)   
  time=time/real(timescale)
  Time_cur=time-time_begine
  
  if(from_uniform .eqv. .false.)then
!!>>Step1: Reformat COCONUT data to structured data
  !!>> Start read raw data
        if(from_ori .eqv. .true.)then  
      filename='Radian.dat'
      open(21,file=filename)
      read(21,*)
      do i=0,nrt    
          read(21,*) radiu(i) 
      end do  
      close(21) 
      !!>> Read the initial reformated COCONUT data 
      filename='corona-time_'//ctime//'.dat'    
      open(27,file=filename)    
      read(27,*)  
      do ii=0,Rindex       
          do i = 0,N_th_p+1           
              do j = 0,N_ph_p+1                
                  read(27,*) utotal_ori(1,ii,i,j),utotal_ori(2,ii,i,j),utotal_ori(3,ii,i,j),&				
                  utotal_ori(4,ii,i,j),utotal_ori(5,ii,i,j),utotal_ori(6,ii,i,j),utotal_ori(7,ii,i,j),&				
                  utotal_ori(8,ii,i,j),utotal_ori(9,ii,i,j),utotal_ori(10,ii,i,j),utotal_ori(11,ii,i,j)            
              end do
          end do  
      end do
      close(27) 
        else
  !!>>1:3 corresponds to x,y,z in global coordinate system
  allocate(Nodal_Pos(1:3,Nodstart:Nodend))
  !!>>Index of 6 nodal in the Triangle prism cell.
  allocate(TriPrismIndex(1:6,Cellstart:Cellend))
  !!>>1:9 corresponds to rho,u,v,w,Bx,By,Bz,p,Phi
  allocate(Cell_state(1:9,Cellstart:Cellend))
  PerNodLayer=(Nodend-Nodstart+1)/(R_IndMax+1)
  PerCellLayer=(Cellend-Cellstart+1)/R_IndMax  
  filename='./corona-time_'//ctime//'.plt'
  open(27,file=filename)
  ComLinStart=1
  ComLinEnd=3
  do i=ComLinStart,ComLinEnd
      read(27,*)
  end do
  CurLin=ComLinEnd-ComLinStart+1
  do j=1,3
      do i=Nodstart,Nodend
          read(27,*)Nodal_Pos(j,i)
      end do
      i=Nodend
      CurLin=CurLin+(Nodend-Nodstart+1)
  end do
  do j=1,9
      do i=Cellstart,Cellend
          read(27,*)Cell_state(j,i)
      end do
      i=Cellend
      CurLin=CurLin+(Cellend-Cellstart+1)
      write(*,*)'Cell_state(j,i)=',Cell_state(j,i),'j=',j,'i=',i,'CurLin=',CurLin
  end do
    
  do i=Cellstart,Cellend
      read(27,*)TriPrismIndex(1,i),TriPrismIndex(2,i),TriPrismIndex(3,i),TriPrismIndex(3,i),&
      TriPrismIndex(4,i),TriPrismIndex(5,i),TriPrismIndex(6,i),TriPrismIndex(6,i)
  end do      
  close(27)
  write(*,*)'Seq_=',Seq_  
  !!<< End read raw   
  
  if (write_CFmesh .eqv. .true.) then
  deallocate (TriPrismIndex,Nodal_Pos)
  allocate (TriPrismIndex(1:7,Cellstart:Cellend),&
          InletTriPrismIndex(1:6,1:PerCellLayer),&
		 OutletTriPrismIndex(1:6,1:PerCellLayer),&
		          Nodal_Pos(1:3,Nodstart:Nodend)&
  )
  filename='./corona-time_0.CFmesh'
  open(27,file=filename)
  ComLinStart=1
  ComLinEnd=15
  do i=ComLinStart,ComLinEnd
      read(27,*)
  end do
  CurLin=ComLinEnd-ComLinStart+1
  do i=Cellstart,Cellend
     read(27,*)TriPrismIndex(1,i),TriPrismIndex(2,i),TriPrismIndex(3,i),TriPrismIndex(4,i),&
	           TriPrismIndex(5,i),TriPrismIndex(6,i),TriPrismIndex(7,i)
	CurLin=CurLin+1		   
  end do
  write(*,*)'CellIndexEnd:CurLin=',CurLin
  ComLinStart=1
  ComLinEnd=6
  do i=ComLinStart,ComLinEnd
      read(27,*)
  end do
  CurLin=CurLin+ComLinEnd-ComLinStart+1
  do i=1,PerCellLayer
      read(27,*)InletTriPrismIndex(1,i),InletTriPrismIndex(2,i),InletTriPrismIndex(3,i),InletTriPrismIndex(4,i),&
	           InletTriPrismIndex(5,i),InletTriPrismIndex(6,i)
	CurLin=CurLin+1		   
  end do
  write(*,*)'InletCellIndexEnd:CurLin=',CurLin
  ComLinStart=1
  ComLinEnd=5
  do i=ComLinStart,ComLinEnd
      read(27,*)
  end do
  CurLin=CurLin+ComLinEnd-ComLinStart+1
  do i=1,PerCellLayer
      read(27,*)OutletTriPrismIndex(1,i),OutletTriPrismIndex(2,i),OutletTriPrismIndex(3,i),OutletTriPrismIndex(4,i),&
	           OutletTriPrismIndex(5,i),OutletTriPrismIndex(6,i)
	CurLin=CurLin+1		   
  end do
  write(*,*)'OutletCellIndexEnd:CurLin=',CurLin
  ComLinStart=1
  ComLinEnd=2
  do i=ComLinStart,ComLinEnd
      read(27,*)
  end do
  CurLin=CurLin+ComLinEnd-ComLinStart+1
  do i=Nodstart,Nodend
      read(27,*)Nodal_Pos(1,i),Nodal_Pos(2,i),Nodal_Pos(3,i)
	CurLin=CurLin+1		   
  end do	
write(*,*)'NodalPosIndexEnd:CurLin=',CurLin
close(27)
write(*,*)'Start Writing CFcase of t = ',ctime,' hour'
  filename='./LayoutandCFmesh/corona-time_'//ctime//'.CFmesh'
  open(27,file=filename)
  write(27,*)'!COOLFLUID_VERSION 2013.9'
  write(27,*)'!CFMESH_FORMAT_VERSION 1.3'
  write(27,*)'!NB_DIM 3'
  write(27,*)'!NB_EQ 9'
 if(Mesh_BCindependent6 .eqv. .true.)then      
  write(27,*)'!NB_NODES 757908 0'
  write(27,*)'!NB_STATES 1495040 0'
  write(27,*)'!NB_ELEM 1495040'
  write(27,*)'!NB_ELEM_TYPES 1'
  write(27,*)'!GEOM_POLYORDER 1'
  write(27,*)'!SOL_POLYORDER 0'
  write(27,*)'!ELEM_TYPES Prism'
  write(27,*)'!NB_ELEM_PER_TYPE 1495040'
  write(27,*)'!NB_NODES_PER_TYPE 6'
 else if(Mesh_30Rs_lvl5 .eqv. .true.)THEN
  write(27,*)'!NB_NODES 192150 0'
  write(27,*)'!NB_STATES 378880 0'
  write(27,*)'!NB_ELEM 378880'
  write(27,*)'!NB_ELEM_TYPES 1'
  write(27,*)'!GEOM_POLYORDER 1'
  write(27,*)'!SOL_POLYORDER 0'
  write(27,*)'!ELEM_TYPES Prism'
  write(27,*)'!NB_ELEM_PER_TYPE 378880'
  write(27,*)'!NB_NODES_PER_TYPE 6'
 end if
 
  write(27,*)'!NB_STATES_PER_TYPE 1'
  write(27,*)'!LIST_ELEM'  
  do i=Cellstart,Cellend
     write(27,*)TriPrismIndex(1,i),TriPrismIndex(2,i),TriPrismIndex(3,i),TriPrismIndex(4,i),&
	           TriPrismIndex(5,i),TriPrismIndex(6,i),TriPrismIndex(7,i)	   
  end do  
  write(27,*)'!NB_TRSs 2'
  write(27,*)'!TRS_NAME Inlet'
  write(27,*)'!NB_TRs 1'
  if(Mesh_BCindependent6 .eqv. .true.)then 
    write(27,*)'!NB_GEOM_ENTS 20480'
  else if(Mesh_30Rs_lvl5 .eqv. .true.)THEN
    write(27,*)'!NB_GEOM_ENTS 5120'
  end if
  write(27,*)'!GEOM_TYPE Face'
  write(27,*)'!LIST_GEOM_ENT'
  do i=1,PerCellLayer
      write(27,*)InletTriPrismIndex(1,i),InletTriPrismIndex(2,i),InletTriPrismIndex(3,i),InletTriPrismIndex(4,i),&
	           InletTriPrismIndex(5,i),InletTriPrismIndex(6,i)	   
  end do
  write(27,*)'!TRS_NAME Outlet'
  write(27,*)'!NB_TRs 1'
  if(Mesh_BCindependent6 .eqv. .true.)then
    write(27,*)'!NB_GEOM_ENTS 20480'
  else if(Mesh_30Rs_lvl5 .eqv. .true.)THEN
    write(27,*)'!NB_GEOM_ENTS 5120'
  end if
  
  write(27,*)'!GEOM_TYPE Face'
  write(27,*)'!LIST_GEOM_ENT'
  do i=1,PerCellLayer
      write(27,*)OutletTriPrismIndex(1,i),OutletTriPrismIndex(2,i),OutletTriPrismIndex(3,i),OutletTriPrismIndex(4,i),&
	           OutletTriPrismIndex(5,i),OutletTriPrismIndex(6,i)	   
  end do
  write(27,*)'!EXTRA_VARS' 
  write(27,*)'!LIST_NODE' 
  do i=Nodstart,Nodend
      write(27,*)Nodal_Pos(1,i),Nodal_Pos(2,i),Nodal_Pos(3,i)	   
  end do 
  write(27,*)'!LIST_STATE 1' 
  do i=Cellstart,Cellend
      !write(27,trim(fmt)//"e24.14e2")Cell_state(1,i),Cell_state(2,i),Cell_state(3,i),Cell_state(4,i),Cell_state(5,i),&
	  !    Cell_state(6,i),Cell_state(7,i),Cell_state(8,i),Cell_state(9,i)
	  write(27,'(9(E24.14))')Cell_state(1,i),Cell_state(2,i),Cell_state(3,i),Cell_state(4,i),Cell_state(5,i),&
	      Cell_state(6,i),Cell_state(7,i),Cell_state(8,i),Cell_state(9,i)
  end do
  write(27,*)'!END'
  close(27)  
  goto 112
  end if

   !!>> Calculate Cell-centroid 
   if(Seq_ == Seq_Start*timescale)then
  do i=Cellstart,Cellend
      do j=1,6
          VertexPos(1:3,j)=Nodal_Pos(1:3,TriPrismIndex(j,i))
      end do
	  AB(1:3)=VertexPos(1:3,2)-VertexPos(1:3,1)
	  AC(1:3)=VertexPos(1:3,3)-VertexPos(1:3,1)
	  call Cross_product(AB(1:3),AC(1:3),Cross(1:3))
	  Area_1=0.5*dsqrt(DOT_PRODUCT(Cross(1:3),Cross(1:3)))
	  AB(1:3)=VertexPos(1:3,5)-VertexPos(1:3,4)
	  AC(1:3)=VertexPos(1:3,6)-VertexPos(1:3,4)
	  call Cross_product(AB(1:3),AC(1:3),Cross(1:3))
	  Area_2=0.5*dsqrt(DOT_PRODUCT(Cross(1:3),Cross(1:3)))
      !temp1=DOT_PRODUCT(VertexPos(1:3,1)-VertexPos(1:3,2),VertexPos(1:3,1)-VertexPos(1:3,2))+&
      !        DOT_PRODUCT(VertexPos(1:3,1)-VertexPos(1:3,3),VertexPos(1:3,1)-VertexPos(1:3,3))+&
      !        DOT_PRODUCT(VertexPos(1:3,2)-VertexPos(1:3,3),VertexPos(1:3,2)-VertexPos(1:3,3))
      !temp2=DOT_PRODUCT(VertexPos(1:3,4)-VertexPos(1:3,5),VertexPos(1:3,4)-VertexPos(1:3,5))+&
      !        DOT_PRODUCT(VertexPos(1:3,4)-VertexPos(1:3,6),VertexPos(1:3,4)-VertexPos(1:3,6))+&
      !        DOT_PRODUCT(VertexPos(1:3,5)-VertexPos(1:3,6),VertexPos(1:3,5)-VertexPos(1:3,6))
      do ii=1,3
          CentralPos(ii,1:2)=0.d0
          do j=1,3
              CentralPos(ii,1)=CentralPos(ii,1)+VertexPos(ii,j)
              CentralPos(ii,2)=CentralPos(ii,2)+VertexPos(ii,j+3)
          end do
          CentralPos(ii,1)=CentralPos(ii,1)/3.0
          CentralPos(ii,2)=CentralPos(ii,2)/3.0
		  Cell_Pos(ii,i)=(CentralPos(ii,1)*Area_1+CentralPos(ii,2)*Area_2)/(Area_1+Area_2)	
          !if(temp1<temp2)then
          !    Cell_Pos(ii,i)=CentralPos(ii,1)/3.0+CentralPos(ii,2)*2.0/3.0
          !else
          !    Cell_Pos(ii,i)=CentralPos(ii,1)*2.0/3.0+CentralPos(ii,2)/3.0
          !end if
      end do
  end do
  !!<< finish calculating cell centroid 
    !!>> Calculate radial distance of each layer
    do i=Cellstart,Cellend-PerCellLayer,PerCellLayer
        temp1=0.0
        Rindex_=int(real((i-1)/PerCellLayer))
        radiu(Rindex_)=0.0
        do j=0,PerCellLayer-1
            temp2=dsqrt(DOT_PRODUCT(Cell_Pos(1:3,i+j)-Cell_Pos(1:3,i+j+PerCellLayer),Cell_Pos(1:3,i+j)-&
			Cell_Pos(1:3,i+j+PerCellLayer)))
            r=dsqrt(DOT_PRODUCT(Cell_Pos(1:3,i+j),Cell_Pos(1:3,i+j)))
            radiu(Rindex_)=radiu(Rindex_)+r
            if(temp1<temp2)then
                temp1=temp2
                r_max=dsqrt(DOT_PRODUCT(Cell_Pos(1:3,i+j),Cell_Pos(1:3,i+j)))
            end if
        end do 
        write(*,*)'disRadial_max=',temp1,'Radian_max=',r_max,'i=',int(real((i-1)/PerCellLayer))
        radiu(Rindex_)=radiu(Rindex_)/real(PerCellLayer)
		write(*,*)"radiu(",Rindex_,")=",radiu(Rindex_)
        temp1=0.0
    end do
    if(write_r .eqv. .true.)then	
	  filename_out='./LayoutandCFmesh/Radian.dat'
      open(27,file=filename_out)
       write(27,*)'Rindex_start=',0,'Rindex_end=',Rindex
      do i=0,Rindex
          write(27,*)radiu(i)
      end do
      close(27)
	end if
   end if
    !!<< End calculating radial distance of each layer 
  
  !!>> Reformat unstructeured COCONUT data to structured one
  layer_index(1:5,0:Rindex,0:N_th_p+1,0:N_ph_p+1)=0
  if(SortorNot .eqv. .false.)then
      filename='./LayoutandCFmesh/layerindex.dat'
	  open(27,file=filename)
      do ii=0,Rindex
        do i = 0,N_th_p+1
            do j = 0,N_ph_p+1
                read(27,*)layer_index(1,ii,i,j),layer_index(2,ii,i,j),layer_index(3,ii,i,j),&
				layer_index(4,ii,i,j),layer_index(5,ii,i,j)
            end do
        end do  
    end do
    close(27)
  end if
  if(Zerod1AU_only .eqv. .true.)then
    do ii=66,71
        !call StructuredFormat(ii,radiu(ii),utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1))
		call UniformStructuredFormat(ii,SortorNot,layer_index(1:5,ii,0:N_th_p+1,0:N_ph_p+1),radiu(ii),&
		      utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1))
	    utotal(1:11,ii,0:N_th_p+1,0:N_ph_p+1)=utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1)
    end do
  else 
    do ii=0,Rindex
        !call StructuredFormat(ii,radiu(ii),utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1))
		call UniformStructuredFormat(ii,SortorNot,layer_index(1:5,ii,0:N_th_p+1,0:N_ph_p+1),radiu(ii),&
		      utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1))
	    utotal(1:11,ii,0:N_th_p+1,0:N_ph_p+1)=utotal_ori(1:11,ii,0:N_th_p+1,0:N_ph_p+1)
    end do
  end if
  if(SortorNot .eqv. .true.)then
      filename='./LayoutandCFmesh/layerindex.dat'
	  open(27,file=filename)
      do ii=0,Rindex
        do i = 0,N_th_p+1
            do j = 0,N_ph_p+1
                write(27,*)layer_index(1,ii,i,j),layer_index(2,ii,i,j),layer_index(3,ii,i,j),&
				layer_index(4,ii,i,j),layer_index(5,ii,i,j)
            end do
        end do  
    end do
    close(27)
  end if
  !>> Finish reformating unstructeured COCONUT data to structured one
        end if
!  deallocate(Nodal_Pos,TriPrismIndex, Cell_Pos,Cell_state)   
    
!!>>Step 2: Rearrange to uniform grid
   !!>>transffered to uniform grid system
!!   call transtouniformgrid(radiu(0:nrt),THETA11(0:N_th_P+1),PHI11(0:N_ph_p+1),&
!!   utotal_ori(1:11,0:Rindex,0:N_th_p+1,0:N_ph_p+1),utotal(1:11,0:nrt,0:N_th_p+1,0:N_ph_p+1))   
   if((Time_sequence .eqv. .false.) .and. (writereformated_dat .eqv. .true.))then
          filename_out='corona-time_'//ctime//'.dat'   
          open(27,file=filename_out)  
          write(27,*)'x',', y',', z',', rho',', u',', v',', w',', Bx',', By',', Bz',', p'   
          do ii=0,Rindex      
              do i = 0,N_th_p+1           
                  do j = 0,N_ph_p+1                
                      write(27,*) utotal_ori(1,ii,i,j),utotal_ori(2,ii,i,j),utotal_ori(3,ii,i,j),utotal_ori(4,ii,i,j),&				
                      utotal_ori(5,ii,i,j),utotal_ori(6,ii,i,j),utotal_ori(7,ii,i,j),utotal_ori(8,ii,i,j),&				
                      utotal_ori(9,ii,i,j),utotal_ori(10,ii,i,j),utotal_ori(11,ii,i,j)            
                  end do       
              end do    
          end do   
          close(27)
   
          filename_out='corona-uniform_t'//ctime//'.dat'	   
          open(27,file=filename_out)        
          write(27,*)'x',', y',', z',', rho',', u',', v',', w',', Bx',', By',', Bz',', p'       
          do ii=0,nrt         
              do i = 0,N_th_p+1           
                  do j = 0,N_ph_p+1               
                      write(27,*)utotal(1,ii,i,j),utotal(2,ii,i,j),utotal(3,ii,i,j),&				 
                      utotal(4,ii,i,j),utotal(5,ii,i,j),utotal(6,ii,i,j),utotal(7,ii,i,j),&				  
                      utotal(8,ii,i,j),utotal(9,ii,i,j),utotal(10,ii,i,j),utotal(11,ii,i,j)             
                  end do         
              end do         
          end do       
          close(27)
		  
          filename_out='corona-uniform_t'//ctime//'.plt'    
          open(28,file=filename_out)  
          WRITE(28,*) 'TITLE="SOLWIND MODELING DATA"'    
          WRITE(28,*) 'VARIABLES = "x0","x1","x2","rho","u","v","w","Bx","By","Bz","p"'    
          WRITE(28,*) 'ZONE T="ONLY1",I=',nrt+1,'J=',int(N_th_p+2),'K=',int(N_ph_p+2)    
          do k=0,N_ph_p+1          
              do j=0,N_th_p+1             
                  do i=0,Rindex                   
                      x=utotal(1,i,j,k)                 
                      y=utotal(2,i,j,k)                 
                      z=utotal(3,i,j,k)                 
                      rho=utotal(4,i,j,k)                
                      u=utotal(5,i,j,k)                 
                      v=utotal(6,i,j,k)                 
                      w=utotal(7,i,j,k)
                      Bx=utotal(8,i,j,k)
                      By=utotal(9,i,j,k)
                      Bz=utotal(10,i,j,k)
                      p=utotal(11,i,j,k)
                      WRITE(28,*)x,y,z,rho,u,v,w,Bx,By,Bz,p
                  end do
              end do
          end do
          close(28)		  
   end if	  
   !!<< end transffering to uniform grid system 
  else if(from_uniform .eqv. .true.)then
      !>>read uniform structured grid
    filename='corona-uniform_t'//ctime//'.dat'
	  open(27,file=filename)
      read(27,*)
      do ii=0,nrt
        do i = 0,N_th_p+1
            do j = 0,N_ph_p+1
                read(27,*)utotal(1,ii,i,j),utotal(2,ii,i,j),utotal(3,ii,i,j),&
				utotal(4,ii,i,j),utotal(5,ii,i,j),utotal(6,ii,i,j),utotal(7,ii,i,j),&
				utotal(8,ii,i,j),utotal(9,ii,i,j),utotal(10,ii,i,j),utotal(11,ii,i,j)
            end do
        end do  
    end do
    close(27)
  end if

    if(writeBxyzonly .eqv. .true.)then
          filename_out='BxyzOnly/corona_Bxyz'//ctimeo//'.plt'   
          open(28,file=filename_out)  
          WRITE(28,*) 'TITLE="SOLWIND MODELING DATA"'    
          WRITE(28,*) 'VARIABLES = "x0","x1","x2","Bx","By","Bz"'    
          WRITE(28,*) 'ZONE T="ONLY1",I=',69,'J=',int(N_th_p/2)+2,'K=',int(N_ph_p/2)+2, ',SOLUTIONTIME=',Time_cur    
          do k=0,int(N_ph_p/2)+1
		      index_k=min(k*2,N_ph_p+1)    
              do j=0,int(N_th_p/2)+1
			    index_j=min(j*2,N_th_p+1)
                  do i=0,68                   
                      x=utotal(1,i,index_j,index_k)                 
                      y=utotal(2,i,index_j,index_k)                 
                      z=utotal(3,i,index_j,index_k)                 
                      Bx=utotal(8,i,index_j,index_k)
                      By=utotal(9,i,index_j,index_k)
                      Bz=utotal(10,i,index_j,index_k)
                      WRITE(28,*)x,y,z,Bx,By,Bz
                  end do
              end do
          end do
          close(28)
    end if
	 
    if((CH_trace .eqv. .true.) .and. (Time_sequence .eqv. .false.))then
    filename_out='B_sph'//ctimeo//'.dat'
    open(1,file=filename_out)
	write(1,*)'Br',', Bth',', Bph'
      do ii=0,nrt
        do i = 0,N_th_p+1
            do j = 0,N_ph_p+1
              call xyz_rtp(THETA11(i),PHI11(j),utotal(8,ii,i,j),utotal(9,ii,i,j),&
			      utotal(10,ii,i,j),B_sph(1),B_sph(2),B_sph(3))
              write(1,*) B_sph(1),B_sph(2),B_sph(3) 
              end do
          end do
      end do
      close(1)	
    end if
	
      !!<<	  
      if(rdial_profile .eqv. .true.)then
      N_radius=68
      !The variation of the number density N and radial speed vr from 1Rs to 20Rs  interpolatedbyRBF    ! mark 2021 01 18   
      THETA=95.d0*pi/180.d0
      PHI=250.d0*pi/180.d0
      call sph_interpradial(PHI,N_radius,radiu(0:nrt),THETA,radiusph_interp(1:17,1:N_radius))
      filename_out='line_vr_N_HDLS'//ctimeo//'.plt'
      open(1,file=filename_out)
      !datafile='line_vr_N_HDLSExp'//ctime//'.plt'  
      !open(1,file=datafile)      
      write(1,*) 'VARIABLES = "X","N","Vr "'      
      write(1,*) 'ZONE T="ONLY1",i=',N_radius
      do i=1,N_radius
          radii=radiu(i)
          rho=radiusph_interp(4,i)
          vr=radiusph_interp(12,i)  
          HD=log10(rho*nrho)
          LS=vr*VS/1000.
          HDLSandLDHS(1,i)=HD
          HDLSandLDHS(2,i)=LS 
          write(1,*) radii,HD,LS                       
      end do
      close(1)
!!      write(*,*) "phi=",PHI/pi*180.,"theta=",THETA/pi*180. 
      
      THETA=160.d0*pi/180.d0
      PHI=250.d0*pi/180.d0
      call sph_interpradial(PHI,N_radius,radiu(0:nrt),THETA,radiusph_interp(1:17,1:N_radius)) 
      filename_out='line_vr_N_LDHS'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'VARIABLES = "X","N","Vr "'      
      write(1,*) 'ZONE T="ONLY1",i=',N_radius
      do i=1,N_radius
          radii=radiu(i)
          rho=radiusph_interp(4,i)
          vr=radiusph_interp(12,i)                      
          LD=log10(rho*nrho)
          HS=vr*VS/1000.
          HDLSandLDHS(3,i)=LD
          HDLSandLDHS(4,i)=HS 
          write(1,*) radii,LD,HS              
      end do
      close(1)
!!      write(*,*) "phi=",PHI/pi*180.,"theta=",THETA/pi*180.  
      filename_out='N-Vr_HDLSandLDHS_sph'//ctimeo//'.plt'
      open(22,file=filename_out)
      !open(22,file='N-Vr_HDLSandLDHS_sphExp'//ctime//'.dat')
      do i=1,N_radius
          write(22,*) radiu(i),HDLSandLDHS(1,i),HDLSandLDHS(2,i),HDLSandLDHS(3,i),HDLSandLDHS(4,i)
      end do
      close(22)
    !goto 112  
    end if
      
      
      !do i=r_indexLbd,nrt
      !    call sph_interp(i,radiu(i),usph_interp(i,:,:,:))
      !end do     
      !Synoptic maps at 3.0Rs--------------
	  if(pannels .eqv. .true.)then
	  if(Zerod1AU_only .eqv. .true.)goto 5121
      usph=0.d0
	  if(justnearest .eqv. .false.)then
          call sph_interp(46,radiu(46),Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1))
	  else
	      call sph(46,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)) 
	  end if
      usph(1:17,0:N_th_p+1,0:N_ph_p+1)=Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)
      

      filename_out='./Pannel/Synoptic_maps3.0Rs'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'TITLE="2.5sph DATA"'
      write(1,*) 'VARIABLES = "Heliolongitude (deg)","Heliolatitude (deg)","Vr (Km/s)","Vth (Km/s)","Vph (Km/s)",&
	        "N (10<sup>5</sup> cm<sup>-3</sup>)","Br (G)","Bth (G)","Bph (G)","T (10<sup>5</sup> K)","p (Pa)",&
            "Vx (Km/s)","Vy (Km/s)","Vz (Km/s)","Bx (G)","By (G)","Bz (G)","PlasmaBeta","|<b>B</b>| (G)"'
      write(1,*) 'ZONE T="ONLY1",i=',N_ph_p  ,',k=',N_th_P, ',SOLUTIONTIME=',Time_cur  
  
      do i1=1,N_th_P 
          th1=90.-THETA11(i1)/pi*180.
          do j1=1,N_ph_p 
        
              ph1=phi11(j1)/pi*180.
              rho=usph(4,i1,j1)*(Rhos/1.672E-27)/1.E6  !number density: X/cm^3
              vr =usph(12,i1,j1)*Vs/1000.              !velocity: km/s   
              vt =usph(13,i1,j1)*Vs/1000.  
              vp =usph(14,i1,j1)*Vs/1000.
              Br =usph(15,i1,j1)*Bs*1.E4               !Gauss 1T=10000G=1.E9nT
              Bt =usph(16,i1,j1)*Bs*1.E4  
              Bp =usph(17,i1,j1)*Bs*1.E4
              p_=usph(11,i1,j1)*Ps
              vx=usph(5,i1,j1)*Vs/1000.
              vy=usph(6,i1,j1)*Vs/1000.
              vz=usph(7,i1,j1)*Vs/1000.
              Bx=usph(8,i1,j1)*Bs*1.E4
              By=usph(9,i1,j1)*Bs*1.E4
              Bz=usph(10,i1,j1)*Bs*1.E4  !Gauss 1T=10000G=1.E9nT
              !T  =gamma*p_/(usph(4,i1,j1)*Rhos)*1.E-5              !Temperature:10^5 K
              T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))*1.E-5              !Temperature:10^5 K
              
              BEXP2=Br**2.d0+Bt**2.d0+Bp**2.d0
              dabsB=dsqrt(BEXP2)
              PlasmaBeta=p_/(BEXP2/1.E4/1.E4/2.0/dmu0)

              write(1,*)ph1,th1,vr,vt,vp,rho/(1.d5),Br,Bt,Bp,T,p_,&
              vx,vy,vz,Bx,By,Bz,PlasmaBeta,dabsB
          end do
      end do
      close(1)  
      
      !Synoptic maps at 21.5Rs--------------
5121      usph=0.d0
	  if(justnearest .eqv. .false.)then
	      rad_temp=21.5
          call sph_interp(69,rad_temp,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1))
      else
          call sph(68,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)) 
	  end if
	  usph(1:17,0:N_th_p+1,0:N_ph_p+1)=Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)
      filename_out='./Pannel/Synoptic_maps21d5Rs'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'TITLE="2.5sph DATA"'
      write(1,*) 'VARIABLES = "Heliolongitude (deg)","Heliolatitude (deg)","Vr (Km/s)","Vth (Km/s)","Vph (Km/s)",&
	        "N (10<sup>3</sup> cm<sup>-3</sup>)","Br (G)","Bth (G)","Bph (G)","T (10<sup>5</sup> K)","p (Pa)",&
            "Vx (Km/s)","Vy (Km/s)","Vz (Km/s)","Bx (G)","By (G)","Bz (G)","PlasmaBeta","|<b>B</b>| (G)"'
      write(1,*) 'ZONE T="ONLY1",i=',N_ph_p  ,',k=',N_th_P, ',SOLUTIONTIME=',Time_cur  
  
      do i1=1,N_th_P 
          th1=90.-THETA11(i1)/pi*180.
          do j1=1,N_ph_p 
        
              ph1=phi11(j1)/pi*180.
              rho=usph(4,i1,j1)*(Rhos/1.672E-27)/1.E6  !number density: X/cm^3
              vr =usph(12,i1,j1)*Vs/1000.              !velocity: km/s   
              vt =usph(13,i1,j1)*Vs/1000.  
              vp =usph(14,i1,j1)*Vs/1000.
              Br =usph(15,i1,j1)*Bs*1.E4               !Gauss 1T=10000G=1.E9nT
              Bt =usph(16,i1,j1)*Bs*1.E4  
              Bp =usph(17,i1,j1)*Bs*1.E4
              p_=usph(11,i1,j1)*Ps
              vx=usph(5,i1,j1)*Vs/1000.
              vy=usph(6,i1,j1)*Vs/1000.
              vz=usph(7,i1,j1)*Vs/1000.
              Bx=usph(8,i1,j1)*Bs*1.E4
              By=usph(9,i1,j1)*Bs*1.E4
              Bz=usph(10,i1,j1)*Bs*1.E4  !Gauss 1T=10000G=1.E9nT
              !T  =gamma*p_/(usph(4,i1,j1)*Rhos)*1.E-5              !Temperature:10^5 K
              T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))*1.E-5              !Temperature:10^5 K
              
              BEXP2=Br**2.d0+Bt**2.d0+Bp**2.d0
              dabsB=dsqrt(BEXP2)
              PlasmaBeta=p_/(BEXP2/1.E4/1.E4/2.0/dmu0)

              write(1,*)ph1,th1,vr,vt,vp,rho/(1.d3),Br,Bt,Bp,T,p_,&
              vx,vy,vz,Bx,By,Bz,PlasmaBeta,dabsB
          end do
      end do
      close(1)
      end if
	  
      goto 113	  

      
      !Synoptic maps at 1.0Rs--------------
      usph=0.d0
	  if(justnearest .eqv. .false.)then
	      call sph_interp(1,radiu(1),Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1))
	  else
	      call sph(1,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)) 
	  end if
      usph(1:17,0:N_th_p+1,0:N_ph_p+1)=Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)
      
      filename_out='Pannel/Synoptic_maps1.0Rs'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'TITLE="1.015sph DATA"'
      write(1,*) 'VARIABLES = "Heliolongitude (deg)","Heliolatitude (deg)","Vr (Km/s)","Vth (Km/s)","Vph (Km/s)",&
	        "N (10<sup>5</sup> cm<sup>-3</sup>)","Br (G)","Bth (G)","Bph (G)","T (10<sup>5</sup> K)","p (Pa)",&
            "Vx (Km/s)","Vy (Km/s)","Vz (Km/s)","Bx (G)","By (G)","Bz (G)","PlasmaBeta","|<b>B</b>| (G)"'
      write(1,*) 'ZONE T="ONLY1",i=',N_ph_p  ,',k=',N_th_P, ',SOLUTIONTIME=',Time_cur  
  
      do i1=1,N_th_P 
          th1=90.-THETA11(i1)/pi*180.
          do j1=1,N_ph_p 
        
              ph1=phi11(j1)/pi*180.
              rho=usph(4,i1,j1)*(Rhos/1.672E-27)/1.E6  !number density: X/cm^3
              vr =usph(12,i1,j1)*Vs/1000.              !velocity: km/s   
              vt =usph(13,i1,j1)*Vs/1000.  
              vp =usph(14,i1,j1)*Vs/1000.
              Br =usph(15,i1,j1)*Bs*1.E4               !Gauss 1T=10000G=1.E9nT
              Bt =usph(16,i1,j1)*Bs*1.E4  
              Bp =usph(17,i1,j1)*Bs*1.E4
              p_=usph(11,i1,j1)*Ps
              vx=usph(5,i1,j1)*Vs/1000.
              vy=usph(6,i1,j1)*Vs/1000.
              vz=usph(7,i1,j1)*Vs/1000.
              Bx=usph(8,i1,j1)*Bs*1.E4
              By=usph(9,i1,j1)*Bs*1.E4
              Bz=usph(10,i1,j1)*Bs*1.E4  !Gauss 1T=10000G=1.E9nT
              !T  =gamma*p_/(usph(4,i1,j1)*Rhos)*1.E-5              !Temperature:10^5 K
              T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))*1.E-5              !Temperature:10^5 K
              
              BEXP2=Br**2.d0+Bt**2.d0+Bp**2.d0
              dabsB=dsqrt(BEXP2)
              PlasmaBeta=p_/(BEXP2/1.E4/1.E4/2.0/dmu0)

              write(1,*)ph1,th1,vr,vt,vp,rho/(1.d5),Br,Bt,Bp,T,p_,&
              vx,vy,vz,Bx,By,Bz,PlasmaBeta,dabsB
          end do
      end do
      close(1)    
	  
!!      end if
!!	    goto 113

      usph=0.d0
	  if(justnearest .eqv. .false.)then
          call sph_interp(14,radiu(14),Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1))
	  else
          call sph(14,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)) 
	  end if
      usph(1:17,0:N_th_p+1,0:N_ph_p+1)=Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)
	  
      filename_out='Synoptic_maps1.035Rs'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'TITLE="1.015sph DATA"'
      write(1,*) 'VARIABLES = "Heliolongitude (deg)","Heliolatitude (deg)","Vr (Km/s)","Vth (Km/s)","Vph (Km/s)",&
	        "N (10<sup>5</sup> cm<sup>-3</sup>)","Br (G)","Bth (G)","Bph (G)","T (10<sup>5</sup> K)","p (Pa)",&
            "Vx (Km/s)","Vy (Km/s)","Vz (Km/s)","Bx (G)","By (G)","Bz (G)","PlasmaBeta","|<b>B</b>| (G)"'
      write(1,*) 'ZONE T="ONLY1",i=',N_ph_p  ,',k=',N_th_P, ',SOLUTIONTIME=',Time_cur  
  
      do i1=1,N_th_P 
          th1=90.-THETA11(i1)/pi*180.
          do j1=1,N_ph_p 
        
              ph1=phi11(j1)/pi*180.
              rho=usph(4,i1,j1)*(Rhos/1.672E-27)/1.E6  !number density: X/cm^3
              vr =usph(12,i1,j1)*Vs/1000.              !velocity: km/s   
              vt =usph(13,i1,j1)*Vs/1000.  
              vp =usph(14,i1,j1)*Vs/1000.
              Br =usph(15,i1,j1)*Bs*1.E4               !Gauss 1T=10000G=1.E9nT
              Bt =usph(16,i1,j1)*Bs*1.E4  
              Bp =usph(17,i1,j1)*Bs*1.E4
              p_=usph(11,i1,j1)*Ps
              vx=usph(5,i1,j1)*Vs/1000.
              vy=usph(6,i1,j1)*Vs/1000.
              vz=usph(7,i1,j1)*Vs/1000.
              Bx=usph(8,i1,j1)*Bs*1.E4
              By=usph(9,i1,j1)*Bs*1.E4
              Bz=usph(10,i1,j1)*Bs*1.E4  !Gauss 1T=10000G=1.E9nT
              !T  =gamma*p_/(usph(4,i1,j1)*Rhos)*1.E-5              !Temperature:10^5 K
              T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))*1.E-5   !Temperature:10^5 K
              
              BEXP2=Br**2.d0+Bt**2.d0+Bp**2.d0
              dabsB=dsqrt(BEXP2)
              PlasmaBeta=p_/(BEXP2/1.E4/1.E4/2.0/dmu0)

              write(1,*)ph1,th1,vr,vt,vp,rho/(1.d5),Br,Bt,Bp,T,p_,&
              vx,vy,vz,Bx,By,Bz,PlasmaBeta,dabsB
          end do
      end do
      close(1)  
      !goto 112
	  
113	  if(interface_for_euhforia .eqv. .true.)then 
      !Synoptic maps at 21Rs--------------
      usph=0.d0
	  if(justnearest .eqv. .false.)then
	      rad_temp=21.5
          call sph_interp(69,rad_temp,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1))
      else
          call sph(69,Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)) 
	  end if
	  usph(1:17,0:N_th_p+1,0:N_ph_p+1)=Rusph_interp(1:17,0:N_th_p+1,0:N_ph_p+1)
      filename_out='foreuhforia/plt/21d5Rs'//ctimeo//'.plt'
      open(1,file=filename_out)
      write(1,*) 'TITLE="2.5sph DATA"'
      write(1,*) 'VARIABLES = "Heliolongitude (deg)","Heliolatitude (deg)","Vr (Km/s)","Vth (Km/s)","Vph (Km/s)",&
	        "N (10<sup>3</sup> cm<sup>-3</sup>)","Br (G)","Bth (G)","Bph (G)","T (10<sup>5</sup> K)","p (Pa)",&
            "Vx (Km/s)","Vy (Km/s)","Vz (Km/s)","Bx (G)","By (G)","Bz (G)","PlasmaBeta","|<b>B</b>| (G)"'
      write(1,*) 'ZONE T="ONLY1",i=',N_ph_p  ,',k=',N_th_P, ',SOLUTIONTIME=',Time_cur  
  
      do i1=1,N_th_P 
          th1=90.-THETA11(i1)/pi*180.
          do j1=1,N_ph_p 
        
              ph1=phi11(j1)/pi*180.
              rho=usph(4,i1,j1)*(Rhos/1.672E-27)/1.E6  !number density: X/cm^3
              vr =usph(12,i1,j1)*Vs/1000.              !velocity: km/s   
              vt =usph(13,i1,j1)*Vs/1000.  
              vp =usph(14,i1,j1)*Vs/1000.
              Br =usph(15,i1,j1)*Bs*1.E4               !Gauss 1T=10000G=1.E9nT
              Bt =usph(16,i1,j1)*Bs*1.E4  
              Bp =usph(17,i1,j1)*Bs*1.E4
              p_=usph(11,i1,j1)*Ps
              vx=usph(5,i1,j1)*Vs/1000.
              vy=usph(6,i1,j1)*Vs/1000.
              vz=usph(7,i1,j1)*Vs/1000.
              Bx=usph(8,i1,j1)*Bs*1.E4
              By=usph(9,i1,j1)*Bs*1.E4
              Bz=usph(10,i1,j1)*Bs*1.E4  !Gauss 1T=10000G=1.E9nT
              !T  =gamma*p_/(usph(4,i1,j1)*Rhos)*1.E-5              !Temperature:10^5 K
              T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))*1.E-5              !Temperature:10^5 K
              
              BEXP2=Br**2.d0+Bt**2.d0+Bp**2.d0
              dabsB=dsqrt(BEXP2)
              PlasmaBeta=p_/(BEXP2/1.E4/1.E4/2.0/dmu0)

              write(1,*)ph1,th1,vr,vt,vp,rho/(1.d3),Br,Bt,Bp,T,p_,&
              vx,vy,vz,Bx,By,Bz,PlasmaBeta,dabsB
          end do
      end do
      close(1)
	  hh_increased=int(Time_cur)
	  seconds=int((Time_cur-real(hh_increased))*3600.d0+minute_drop*60.d0)
	  second=mod(seconds,60)
	  minute=int(seconds/60)
	  call hours2yymmddhh(mmddhh,hh_increased,hours,days,months,year)
      write(Char_year ,"(i4.4)") year
	  write(Char_mm ,"(i2.2)") months
	  write(Char_dd ,"(i2.2)") days
	  write(Char_hh ,"(i2.2)") hours
	  write(Char_minute,"(i2.2)") minute
	  write(Char_ss,"(i2.2)") second
	  
	  
	  if(update .eqv. .true.)then
	  filename_out='./foreuhforia/dat/solar_wind_boundary_'//ctime_2//'.dat'
	  open(1,file=filename_out)
	  write(1,'(A)')"Time:"
	  !write(formatted_value, '(E22.19E3)') Time_cur
	  !write(*,*)"formatted_value=",formatted_value,'Time_cur=',Time_cur
	  !write(1, '(A)') trim(adjustl(formatted_value))
	  write(1,'(A)')Char_year//"-"//Char_mm//"-"//Char_dd//"T"//Char_hh//":"//Char_minute//":"//Char_ss
	  write(1,'(A)')"Radius of sphere:"
	  write(1,'(A)')"14959787070.0"
	  write(1,'(A)')"Number of colatitude grid points:"
	  write(1,'(I0)')N_th_P
	  write(1,'(A)')"Colatitude grid points:"
	  !fmt='(1L(E15.6))'
	  do i1=1,N_th_P
	     line = ''
         write(line, '(es30.19)') THETA11(i1)
         !write(*,*)ADJUSTL(line)
         write(1, '(A)')ADJUSTL(line)
	    !write(formatted_value, '(E22.19E3)') THETA11(i1)  ! Format the value
        !write(1, '(A)') trim(adjustl(formatted_value))  ! Write to file
	    !write(1,'(1(E28.18))')THETA11(i1)
	  end do
	  write(1,'(A)')"Number of longitude grid points:"
	  write(1,'(I0)')N_ph_P
	  write(1,'(A)')"Longitude grid points:"
	  do j1=1,N_ph_P
	  	 line = ''
         write(line, '(es30.19)') PHI11(j1)
         !write(*,*)ADJUSTL(line)
         write(1, '(A)')ADJUSTL(line)
	    !write(formatted_value, '(E22.19E3)') PHI11(j1)
	    !write(1, '(A)') trim(adjustl(formatted_value))
	    !write(1,'(1(E28.18))')PHI11(j1)
	  end do
	  write(1,'(A)')"vr"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  vr =usph(12,i1,j1)*Vs
		  line = ''
          write(line, '(es30.19)') vr
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') vr
	      !write(1, '(A)') trim(adjustl(formatted_value))
		  !write(1,'(1(E28.18))')vr
		end do
	  end do
	  write(1,'(A)')"vp"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  vp =usph(14,i1,j1)*Vs
		  line = ''
          write(line, '(es30.19)') vp
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') vp
	      !write(1, '(A)') trim(adjustl(formatted_value))
		  !write(1,'(1(E28.18))')vp
		end do
	  end do
	  write(1,'(A)')"vt"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  vt =usph(13,i1,j1)*Vs
		  line = ''
          write(line, '(es30.19)') vt
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') vt
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')vt
		end do
	  end do
	  write(1,'(A)')"number_density"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P	  
	      rho=usph(4,i1,j1)*(Rhos/1.672E-27)  !number density: X/m^3
		  line = ''
          write(line, '(es30.19)') rho
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') rho
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')rho
		end do
	  end do
	  write(1,'(A)')"temperature"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P	  
	      T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))
		  line = ''
          write(line, '(es30.19)') T
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') T
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')T
		end do
	  end do
      write(1,'(A)')"Br"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  Br=usph(15,i1,j1)*Bs
		  line = ''
          write(line, '(es30.19)') Br
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') Br
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')Br
		end do
	  end do
	  write(1,'(A)')"Bp"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  Bp=usph(17,i1,j1)*Bs
		  line = ''
          write(line, '(es30.19)') Bp
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)')Bp
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')Bp
		end do
	  end do
	  write(1,'(A)')"Bt"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  Bt=usph(16,i1,j1)*Bs
		  line = ''
          write(line, '(es30.19)') Bt
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)
		  !write(formatted_value, '(E22.19E3)') Bt
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')Bt
		end do
	  end do
      close(1)
      end if	
      if(update .eqv. .false.)then
	  filename_out='./foreuhforia/dat_ori/solar_wind_boundary_'//ctime_2//'.dat'
	  open(1,file=filename_out)
	  write(1,'(A)')"Time:"
	  !write(formatted_value, '(E22.19E3)') Time_cur
	  !write(*,*)"formatted_value=",formatted_value,'Time_cur=',Time_cur
	  !write(1, '(A)') trim(adjustl(formatted_value))
	  write(1,'(A)')Char_year//"-"//Char_mm//"-"//Char_dd//"T"//Char_hh//":"//Char_minute//":"//Char_ss
	  write(1,'(A)')"Radius of sphere:"
	  write(1,'(A)')"14959787070.0"
	  write(1,'(A)')"Number of colatitude grid points:"
	  write(1,'(I0)')N_th_P
	  write(1,'(A)')"Colatitude grid points:"
	  !fmt='(1L(E15.6))'
	  do i1=1,N_th_P
	    !write(formatted_value, '(E22.19E3)') THETA11(i1)  ! Format the value
        !write(1, '(A)') trim(adjustl(formatted_value))  ! Write to file
		 line = ''
         write(line, '(es30.19)') THETA11(i1)
         !write(*,*)ADJUSTL(line)
         write(1, '(A)')ADJUSTL(line)
	     !write(1,'(1(E28.18))')THETA11(i1)
	  end do
	  write(1,'(A)')"Number of longitude grid points:"
	  write(1,'(I0)')N_ph_P
	  write(1,'(A)')"Longitude grid points:"
	  do j1=1,N_ph_P
		 line = ''
         write(line, '(es30.19)') PHI11(j1)
         !write(*,*)ADJUSTL(line)
         write(1, '(A)')ADJUSTL(line)
	    !write(formatted_value, '(E22.19E3)') PHI11(j1)
	    !write(1, '(A)') trim(adjustl(formatted_value))
	    !write(1,'(1(E28.18))')PHI11(j1)
	  end do
	  write(1,'(A)')"vr"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  vr =usph(12,i1,j1)*Vs
		  line = ''
          write(line, '(es30.19)') vr
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)		  
		  !write(formatted_value, '(E22.19E3)') vr
	      !write(1, '(A)') trim(adjustl(formatted_value))
		  !write(1,'(1(E28.18))')vr
		end do
	  end do
      write(1,'(A)')"Br"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P
		  Br=usph(15,i1,j1)*Bs
		  line = ''
          write(line, '(es30.19)') Br
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)	
		  !write(formatted_value, '(E22.19E3)') Br
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')Br
		end do
	  end do
	  write(1,'(A)')"number_density"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P	  
	      rho=usph(4,i1,j1)*(Rhos/1.672E-27)  !number density: X/m^3
		  line = ''
          write(line, '(es30.19)') rho
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)	
		  !write(formatted_value, '(E22.19E3)') rho
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')rho
		end do
	  end do
	  write(1,'(A)')"temperature"
	  do i1= 1, N_th_P
	    do j1 = 1, N_ph_P	  
	      T = p_/(2.0*usph(4,i1,j1)*Rhos*kB/(mu_cor*mH))
		  line = ''
          write(line, '(es30.19)') T
          !write(*,*)ADJUSTL(line)
          write(1, '(A)')ADJUSTL(line)	
		  !write(formatted_value, '(E22.19E3)') T
	      !write(1, '(A)') trim(adjustl(formatted_value))		  
		  !write(1,'(1(E28.18))')T
		end do
	  end do
      close(1)	  
      end if	  
	  end if
	  
	  if(meiridians .eqv. .true.)then
      FR_cent=250.d0 !76.d0 !61.d0 ! 235.d0 !
      !Long=256.d0 !274.d0 !281.0 !
	  Long=Long0-360.d0/653.d0*Time_cur
	  halfLength=0.0
      Longitude_shift=FR_cent-Long !180.d0
	  if(Longitude_shift<0)Longitude_shift=Longitude_shift+2.d0*180.0
      if(Longitude_shift>2.d0*180.0)Longitude_shift=Longitude_shift-2.d0*180.0
	  
	  N_radius=68 !57 ! 20Rs
      eps=1.d-20
	  phi_L1=Longitude_shift+halfLength
      !Earth View
        !phi_L1=(201.252+time*58.466/61.0/24.0)  !*pi/180.d0  !143.d0*pi/180.d0
        if(phi_L1>2.d0*180.0)phi_L1=phi_L1-2.d0*180.0
		if(phi_L1<0.0)phi_L1=phi_L1+2.d0*180.0
!        write(*,*)'phi_L1=',phi_L1,'int(phi_L1)=',int(phi_L1)
        deallocate(central_meridian)
        allocate(character(len=10)::central_meridian)
        call int2str(int(phi_L1),central_meridian)
!        write(*,*)'len(central_meridian)=',len(central_meridian)
        central_meridian=TRIM(central_meridian)
        write(*,*)' len(central_meridian)=',len(central_meridian)
        write(*,*)'central_meridian=',central_meridian
        phi_L1=phi_L1*pi/180.d0
        filename_out='./Meridian/Earth/SCMEC_20_'//ctimeo//'_'//central_meridian//'.plt'
        !filename_out='SCMEC_20(111-291)'//ctime//'.plt'
        open(1,file=filename_out)
        write(1,*) 'TITLE="SOLWIND MODELING DATA"'
        write(1,*) 'VARIABLES = "X","Z","N*(r/Rs)<sup>2</sup>","Vx (Km/s)","Vz (Km/s)","Bx (G)","Bz (G)",&
		"Vr (Km/s)","Br (G)","log<sub>10</sub>N","p (Pa)","Bth (G)","Bph (G)","|<b>V</b>|/V<sub>A</sub>"'
 	   
        write(1,*) 'ZONE T="ONLY1",i=',N_radius,',j=',N_th_P+2, ',SOLUTIONTIME=',Time_cur                
        PHI0=phi_L1+0.5d0*pi
        if(PHI0>2.d0*pi)PHI0=PHI0-2.d0*pi
        PHI=PHI0  !PHI11(K) !k=52
        call sph_interplong(PHI,N_radius,radiu(0:nrt),THETA11(0:N_th_P+1),&
		 Lousph_interp(1:17,1:N_radius,0:N_th_P+1))
        do j=0,N_th_P+1
           do i=1,N_radius
               !call coordinate_rtp_xyz(x0,y0,z0,radiu(i),THETA11(j),PHI)  
        x0 =radiu(i)*dsin(THETA11(j))*dcos(PHI)           
        y0 =radiu(i)*dsin(THETA11(j))*dsin(PHI)      
        z0 =radiu(i)*dcos(THETA11(j))               
              x=Lousph_interp(1,i,j)
              y=Lousph_interp(2,i,j)
		      z=Lousph_interp(3,i,j) 
!!              if((j==N_th_P/2 .or. j==N_th_P/2-1 .or. j==N_th_P/2+1) .and. &
!!			  abs(real(i-N_radius))<2.0)print *,"i,j=",i,j,"dr=",&
!!			  dsqrt((x-x0)**2.d0+(y-y0)**2.d0+(z-z0)**2.d0),"THETA11(j)=",THETA11(j) !,"radiu(i)=",radiu(i)
               call XYZ_ROTATE(X,Y,Z,X0,Y0,Z0,PHI0)
               
              Bx=Lousph_interp(8,i,j)
              By=Lousph_interp(9,i,j)
		      Bz=Lousph_interp(10,i,j)     
               call XYZ_ROTATE(Bx,By,Bz,Bx0,By0,Bz0,PHI0)
               
              vx=Lousph_interp(5,i,j)
              vy=Lousph_interp(6,i,j)
		      vz=Lousph_interp(7,i,j)    
               call XYZ_ROTATE(vx,vy,vz,vx0,vy0,vz0,PHI0) 

              rho=Lousph_interp(4,i,j)
              vx=vx0
              vy=vy0
              vz=vz0
              vv=dsqrt(vx**2+vy**2+vz**2)
              vr=Lousph_interp(12,i,j)
              P=Lousph_interp(11,i,j)
              Bx=Bx0
              By=By0
		      Bz=Bz0
              Br=Lousph_interp(15,i,j)
              Bth=Lousph_interp(16,i,j)
              Bph=Lousph_interp(17,i,j)
              bb=dsqrt(br**2+bth**2+bph**2)
                            
              write(1,*)x0,z0,radiu(i)**2.0*rho*(Rhos/1.672E-27)/1.E6,vx*vs/1000.,vz*vs/1000.,bx*Bs*1.d4,bz*Bs*1.d4,vr*vs/1000.,&
			  br*Bs*1.d4,log10(rho*(Rhos/1.672E-27)/1.E6),P*Ps,Bth*Bs*1.d4,Bph*Bs*1.d4,&
			  dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)/dsqrt((Bx**2.d0+By**2.d0+Bz**2.d0)/rho) !,dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)*vs/dsqrt((Bx**2.d0+By**2.d0+Bz**2.d0)*(Bs*1.E4)**2.d0/rho/Rhos/8.d0/pi)   !
           end do
        end do
       
	    write(1,*) 'ZONE T="ONLY2",i=',N_radius,',j=',N_th_P+2
        PHI=phi_L1-0.5d0*pi
        if(PHI<0.d0*pi)PHI=PHI+2.d0*pi        !k=12    
        
        call sph_interplong(PHI,N_radius,radiu(0:nrt),THETA11(0:N_th_P+1),&
		  Lousph_interp(1:17,1:N_radius,0:N_th_P+1))
        do j=0,N_th_P+1
           do i=1,N_radius
               !call coordinate_rtp_xyz(x,y,z,radiu(i),THETA11(j),PHI11(k))                          
        x0 =radiu(i)*dsin(THETA11(j))*dcos(PHI)           
        y0 =radiu(i)*dsin(THETA11(j))*dsin(PHI)      
        z0 =radiu(i)*dcos(THETA11(j))               
              x=Lousph_interp(1,i,j)
              y=Lousph_interp(2,i,j)
		      z=Lousph_interp(3,i,j) 
!!              if((j==N_th_P/2 .or. j==N_th_P/2-1 .or. j==N_th_P/2+1) .and. &
!!			    abs(real(i-N_radius))<2.0)print *,"i,j=",i,j,"dr=",&
!!			  dsqrt((x-x0)**2.d0+(y-y0)**2.d0+(z-z0)**2.d0),"THETA11(j)=",THETA11(j)
               call XYZ_ROTATE(X,Y,Z,X0,Y0,Z0,PHI0)
               
              Bx=Lousph_interp(8,i,j)
              By=Lousph_interp(9,i,j)
		      Bz=Lousph_interp(10,i,j)     
               call XYZ_ROTATE(Bx,By,Bz,Bx0,By0,Bz0,PHI0)
               
              vx=Lousph_interp(5,i,j)
              vy=Lousph_interp(6,i,j)
		      vz=Lousph_interp(7,i,j)    
               call XYZ_ROTATE(vx,vy,vz,vx0,vy0,vz0,PHI0) 
              
              rho=Lousph_interp(4,i,j)
              vx=vx0
              vy=vy0
              vz=vz0
              vv=dsqrt(vx**2+vy**2+vz**2)
              vr=Lousph_interp(12,i,j)
              P=Lousph_interp(11,i,j)
              Bx=Bx0
              By=By0
		      Bz=Bz0
              Br=Lousph_interp(15,i,j)
              Bth=Lousph_interp(16,i,j)
              Bph=Lousph_interp(17,i,j)
              bb=dsqrt(br**2+bth**2+bph**2)
                            
              write(1,*)x0,z0,radiu(i)**2.0*rho*(Rhos/1.672E-27)/1.E6,vx*vs/1000.,vz*vs/1000.,bx*Bs*1.d4,bz*Bs*1.d4,vr*vs/1000.,&
			  br*Bs*1.d4,log10(rho*(Rhos/1.672E-27)/1.E6),P*Ps,Bth*Bs*1.d4,Bph*Bs*1.d4,&
			  dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)/dsqrt((Bx**2.d0+By**2.d0+Bz**2.d0)/rho)
           end do
        end do
        close(1)
        !Stereo-A View
        phi_L1R=phi_L1-0.5d0*pi
        !phi_L1R=(110.503+time*63.64/61.0/24.0) !*pi/180.d0
        if(phi_L1R<0.d0)phi_L1R=phi_L1R+2.d0*pi
        deallocate(central_meridian)
        allocate(character(len=10)::central_meridian)
        call int2str(int(phi_L1R*180.0/pi),central_meridian)
        central_meridian=TRIM(central_meridian)
        phi_L1R=phi_L1R
        !write(*,*)'phi_L1R=',phi_L1R,'int(phi_L1R)=',int(phi_L1R),'meridian_SteA=',meridian_SteA
        filename_out='./Meridian/SterA/SCMEC_20_'//ctimeo//'_'//central_meridian//'.plt'
      open(1,file=filename_out)
        write(1,*) 'TITLE="SOLWIND MODELING DATA"'     
        write(1,*) 'VARIABLES = "X","Z","N*(r/Rs)<sup>2</sup>","Vx (Km/s)","Vz (Km/s)","Bx (G)","Bz (G)",&
		"Vr (Km/s)","Br (G)","log<sub>10</sub>N","p (Pa)","Bth (G)","Bph (G)","|<b>V</b>|/V<sub>A</sub>"'
        write(1,*) 'ZONE T="ONLY1",i=',N_radius,',j=',N_th_P+2, ',SOLUTIONTIME=',Time_cur
        PHI0=phi_L1R+0.5d0*pi
        if(PHI0>2.d0*pi)PHI0=PHI0-2.d0*pi
        PHI=PHI0  !PHI11(K) k=72
        call sph_interplong(PHI,N_radius,radiu(0:nrt),THETA11(0:N_th_P+1),&
		Lousph_interp(1:17,1:N_radius,0:N_th_P+1))  
        do j=0,N_th_P+1
           do i=1,N_radius
               !call coordinate_rtp_xyz(x,y,z,radiu(i),THETA11(j),PHI11(k)) 
        x0 =radiu(i)*dsin(THETA11(j))*dcos(PHI)           
        y0 =radiu(i)*dsin(THETA11(j))*dsin(PHI)      
        z0 =radiu(i)*dcos(THETA11(j))               
              x=Lousph_interp(1,i,j)
              y=Lousph_interp(2,i,j)
		      z=Lousph_interp(3,i,j) 
!!              if((j==N_th_P/2 .or. j==N_th_P/2-1 .or. j==N_th_P/2+1) .and. &
!!			  abs(real(i-N_radius))<2.0)print *,"i,j=",i,j,"dr=",&
!!			  dsqrt((x-x0)**2.d0+(y-y0)**2.d0+(z-z0)**2.d0),"THETA11(j)=",THETA11(j)
               call XYZ_ROTATE(X,Y,Z,X0,Y0,Z0,PHI0)
               
              Bx=Lousph_interp(8,i,j)
              By=Lousph_interp(9,i,j)
		      Bz=Lousph_interp(10,i,j)     
               call XYZ_ROTATE(Bx,By,Bz,Bx0,By0,Bz0,PHI0)
               
              vx=Lousph_interp(5,i,j)
              vy=Lousph_interp(6,i,j)
		      vz=Lousph_interp(7,i,j)    
               call XYZ_ROTATE(vx,vy,vz,vx0,vy0,vz0,PHI0) 
              
              rho=Lousph_interp(4,i,j)
              vx=vx0
              vy=vy0
              vz=vz0
              vv=dsqrt(vx**2+vy**2+vz**2)
              vr=Lousph_interp(12,i,j)
              P=Lousph_interp(11,i,j)
              Bx=Bx0
              By=By0
		      Bz=Bz0
              Br=Lousph_interp(15,i,j)
              Bth=Lousph_interp(16,i,j)
              Bph=Lousph_interp(17,i,j)
              bb=dsqrt(br**2+bth**2+bph**2)
                            
              write(1,*) x0,z0,radiu(i)**2.0*rho*(Rhos/1.672E-27)/1.E6,vx*vs/1000.,vz*vs/1000.,bx*Bs*1.d4,bz*Bs*1.d4,vr*vs/1000.,&
			  br*Bs*1.d4,log10(rho*(Rhos/1.672E-27)/1.E6),P*Ps,Bth*Bs*1.d4,Bph*Bs*1.d4,&
			  dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)/dsqrt((Bx**2.d0+By**2.d0+Bz**2.d0)/rho)
           end do
        end do
       
	    write(1,*) 'ZONE T="ONLY2",i=',N_radius,',j=',N_th_P+2
        PHI=phi_L1R-0.5d0*pi
        if(PHI<0.d0*pi)PHI=PHI+2.d0*pi        !k=32        
        
        call sph_interplong(PHI,N_radius,radiu(0:nrt),THETA11(0:N_th_P+1),&
		Lousph_interp(1:17,1:N_radius,0:N_th_P+1))
        do j=0,N_th_P+1
           do i=1,N_radius
               !call coordinate_rtp_xyz(x,y,z,radiu(i),THETA11(j),PHI11(k)) 
        x0 =radiu(i)*dsin(THETA11(j))*dcos(PHI)           
        y0 =radiu(i)*dsin(THETA11(j))*dsin(PHI)      
        z0 =radiu(i)*dcos(THETA11(j))               
              x=Lousph_interp(1,i,j)
              y=Lousph_interp(2,i,j)
		      z=Lousph_interp(3,i,j) 
!!              if((j==N_th_P/2 .or. j==N_th_P/2-1 .or. j==N_th_P/2+1) .and. &
!!			  abs(real(i-N_radius))<2.0)print *,"i,j=",i,j,&
!!			  "dr=",dsqrt((x-x0)**2.d0+(y-y0)**2.d0+(z-z0)**2.d0),"THETA11(j)=",THETA11(j)
               call XYZ_ROTATE(X,Y,Z,X0,Y0,Z0,PHI0)
               
              Bx=Lousph_interp(8,i,j)
              By=Lousph_interp(9,i,j)
		      Bz=Lousph_interp(10,i,j)     
               call XYZ_ROTATE(Bx,By,Bz,Bx0,By0,Bz0,PHI0)
               
              vx=Lousph_interp(5,i,j)
              vy=Lousph_interp(6,i,j)
		      vz=Lousph_interp(7,i,j)    
               call XYZ_ROTATE(vx,vy,vz,vx0,vy0,vz0,PHI0) 
              
              rho=Lousph_interp(4,i,j)
              vx=vx0
              vy=vy0
              vz=vz0
              vv=dsqrt(vx**2+vy**2+vz**2)
              vr=Lousph_interp(12,i,j)
              P=Lousph_interp(11,i,j)
              Bx=Bx0
              By=By0
		      Bz=Bz0
              Br=Lousph_interp(15,i,j)
              Bth=Lousph_interp(16,i,j)
              Bph=Lousph_interp(17,i,j)
              bb=dsqrt(br**2+bth**2+bph**2)
                            
              write(1,*) x0,z0,radiu(i)**2.0*rho*(Rhos/1.672E-27)/1.E6,vx*vs/1000.,vz*vs/1000.,bx*Bs*1.d4,bz*Bs*1.d4,vr*vs/1000.,&
			  br*Bs*1.d4,log10(rho*(Rhos/1.672E-27)/1.E6),P*Ps,Bth*Bs*1.d4,Bph*Bs*1.d4,&
			  dsqrt(vx**2.d0+vy**2.d0+vz**2.d0)/dsqrt((Bx**2.d0+By**2.d0+Bz**2.d0)/rho)
           end do
        end do
        close(1)  
      end if		
        deallocate(Nodal_Pos,TriPrismIndex,  Cell_state)
        deallocate(utotal_ori,utotal, Rusph_interp,Lousph_interp,usph) 
     end do    
112    print *, "Hello, over"
    contains

    subroutine int2str(int_,str)
    integer,intent(in) :: int_
    character(len=10),intent(inout) :: str
    if(int_ .lt. 10)then
        write(str,'(i1.1)') int_
    else if(int_ .ge. 10 .and. int_ .lt. 100)then
        write(str,'(i2.2)') int_
    else if(int_ .ge. 100 .and. int_ .lt. 1000)then
        write(str,'(i3.3)') int_
    else if(int_ .ge. 1000 .and. int_ .lt. 10000)then
        write(str,'(i4.4)') int_
    else if(int_ .ge. 10000 .and. int_ .lt. 100000)then
        write(str,'(i5.5)') int_
    else if(int_ .ge. 100000 .and. int_ .lt. 1000000)then
        write(str,'(i6.6)') int_
    end if
    end subroutine
	
	subroutine int2strwithonesuffix(int_,str)
    integer,intent(in) :: int_
    character(len=10),intent(inout) :: str
	character(len=10) :: str_int
	character(len=1) :: str_mod
	!character(len=:),allocatable :: str_temp
	integer :: int_temp,mod_temp
	int_temp= floor(real(int_/10))
	mod_temp= mod(int_,10)
    if(int_temp .lt. 10)then
        write(str_int,'(i1.1)') int_temp
    else if(int_temp .ge. 10 .and. int_temp .lt. 100)then
        write(str_int,'(i2.2)') int_temp
    else if(int_temp .ge. 100 .and. int_temp .lt. 1000)then
        write(str_int,'(i3.3)') int_temp
    else if(int_temp .ge. 1000 .and. int_temp .lt. 10000)then
        write(str_int,'(i4.4)') int_temp
    else if(int_temp .ge. 10000 .and. int_temp .lt. 100000)then
        write(str_int,'(i5.5)') int_temp
    else if(int_temp .ge. 100000 .and. int_temp .lt. 1000000)then
        write(str_int,'(i6.6)') int_temp
    end if
	if(mod_temp>0)then
	    write(str_mod,'(i1.1)') mod_temp
	end if
	!length=len_trim(str_int)+length_mod
	!allocate(character(len=length):: str_temp)
	if(mod_temp>0)then
	   str=trim(str_int)//'_'//str_mod
	else
	   str=trim(str_int)
	end if   
	end subroutine int2strwithonesuffix
	
	subroutine int2strwithsuffix(int_,decimal_space,str,str2,str3)
    integer,intent(in) :: int_,decimal_space
    character(len=10),intent(inout) :: str,str2,str3
	character(len=10) :: str_int,str_mod 
	character(len=2) :: str_mod2
	integer :: int_temp,mod_temp
	int_temp= floor(real(int_/decimal_space))
	mod_temp= mod(int_,decimal_space)
    if(int_temp .lt. 10)then
        write(str_int,'(i1.1)') int_temp
    else if(int_temp .ge. 10 .and. int_temp .lt. 100)then
        write(str_int,'(i2.2)') int_temp
    else if(int_temp .ge. 100 .and. int_temp .lt. 1000)then
        write(str_int,'(i3.3)') int_temp
    else if(int_temp .ge. 1000 .and. int_temp .lt. 10000)then
        write(str_int,'(i4.4)') int_temp
    else if(int_temp .ge. 10000 .and. int_temp .lt. 100000)then
        write(str_int,'(i5.5)') int_temp
    else if(int_temp .ge. 100000 .and. int_temp .lt. 1000000)then
        write(str_int,'(i6.6)') int_temp
    end if
	
	if(mod_temp>0)then
	  if(mod(mod_temp,10) .eq. 0)then
	    write(str_mod,'(i1.1)') mod_temp/10
		write(str_mod2,'(i2.2)') mod_temp
	  else
	    write(str_mod,'(i2.2)') mod_temp
		write(str_mod2,'(i2.2)') mod_temp
	  end if
	   str=trim(str_int)//'_'//trim(str_mod)
	   str2=trim(str_int)//'_'//trim(str_mod2)
	   str3=trim(str_int)//trim(str_mod2)
	else
	   str=trim(str_int)
	   str2=trim(str_int)//'_'//'00'
	   str3=trim(str_int)//'00'
	end if   
	end subroutine int2strwithsuffix
	
	subroutine Cross_product(AA,BB,CC)
	implicit none
    real*8,intent(in) :: AA(3),BB(3)
    real*8,intent(out) :: CC(3)
	CC(1)=AA(2)*BB(3)-AA(3)*BB(2)
	CC(2)=AA(3)*BB(1)-AA(1)*BB(3)
	CC(3)=AA(1)*BB(2)-AA(2)*BB(1)
	end subroutine Cross_product

    subroutine coordinate_xyz_rtp(r,theta,phi,x,y,z)!        ֱ       ת  
        implicit double precision(a-h,o-z)
        pi = 4.0*datan(1.d0)
        r = dsqrt(x**2.0+y**2.0+z**2.0)
        rxy = dsqrt(x**2.0+y**2.0)
        if(y.gt.0.d0) then
           phi = dacos(x/rxy)
        else
           phi = -dacos(x/rxy)
        end if
        theta=dacos(z/r) 
     end subroutine coordinate_xyz_rtp

     subroutine coordinate_rtp_xyz(x,y,z,r,theta,phi)  
       implicit double precision(a-h,o-z)
       x = r*dsin(theta)*dcos(phi)
       y = r*dsin(theta)*dsin(phi)
       z = r*dcos(theta)
     end subroutine coordinate_rtp_xyz

     subroutine xyz_rtp(th,ph,vx,vy,vz,vr,vt,vp)
       implicit double precision(a-h,o-z)
       vr=vx*dsin(th)*dcos(ph)+vy*dsin(th)*dsin(ph)+vz*dcos(th)
       vt=vx*dcos(th)*dcos(ph)+vy*dcos(th)*dsin(ph)-vz*dsin(th)
       vp=-vx*dsin(ph)+vy*dcos(ph)
     end subroutine xyz_rtp
    
    subroutine rtp_xyz(th,ph,vr,vt,vp,vx,vy,vz)
      implicit double precision(a-h,o-z)
      vx=dsin(th)*dcos(Ph)*vr+dcos(th)*dcos(PH)*vt-dsin(PH)*vp
      vy=dsin(th)*dsin(Ph)*vr+dcos(th)*dsin(PH)*vt+dcos(PH)*vp
      vz=dcos(th)*vr-dsin(th)*vt
    end subroutine rtp_xyz 
 
    SUBROUTINE XYZ_ROTATE(X,Y,Z,X0,Y0,Z0,PHI)
    IMPLICIT NONE
    REAL*8,INTENT(IN)::X,Y,Z,PHI
    REAL*8,INTENT(OUT)::X0,Y0,Z0

    X0=X*DCOS(PHI)+Y*DSIN(PHI)
    Y0=Y*DCOS(PHI)-X*DSIN(PHI)
    Z0=Z      
    END SUBROUTINE XYZ_ROTATE      
    
    subroutine Lagrange_3D(rad,theta,phi,r,th,ph,varis,vari)
  implicit none
  real*8,intent(in) :: rad(2),theta(2),phi(2),r,th,ph,varis(2,2,2)
  real*8,intent(out) :: vari
  integer :: i,j,k
  real*8  :: COF1,COF2,COF3
  real*8  :: CCOF(2,2,2)
  
				COF1=(r-rad(1))/(rad(2)-rad(1))
				COF2=(th-theta(1))/(theta(2)-theta(1))
				COF3=(ph-phi(1))/(phi(2)-phi(1))

				CCOF(1,1,1)=(1.-COF1)*(1.-COF2)*(1.-COF3)
				CCOF(1,1,2)=(1.-COF1)*(1.-COF2)*COF3
				CCOF(1,2,1)=(1.-COF1)*COF2*(1.-COF3)
				CCOF(1,2,2)=(1.-COF1)*COF2*COF3
				CCOF(2,1,1)=COF1*(1.-COF2)*(1.-COF3)
				CCOF(2,1,2)=COF1*(1.-COF2)*COF3
				CCOF(2,2,1)=COF1*COF2*(1.-COF3)
				CCOF(2,2,2)=COF1*COF2*COF3

				vari=varis(1,1,1)*CCOF(1,1,1)+varis(2,1,1)*CCOF(2,1,1)+varis(1,2,1)*CCOF(1,2,1)&
     				+varis(1,1,2)*CCOF(1,1,2)+varis(1,2,2)*CCOF(1,2,2)+varis(2,2,1)*CCOF(2,2,1)&
     			    +varis(2,1,2)*CCOF(2,1,2)+varis(2,2,2)*CCOF(2,2,2)	 
    end subroutine
    
    subroutine Lagrange_2D(theta,phi,th,ph,varis,vari)
  implicit none
  real*8,intent(in) :: theta(2),phi(2),th,ph,varis(2,2)
  real*8,intent(out) :: vari
  real*8  :: COF2,COF3
  real*8  :: CCOF(2,2)
  
				COF2=(th-theta(1))/(theta(2)-theta(1))
				COF3=(ph-phi(1))/(phi(2)-phi(1))

				CCOF(1,1)=(1.-COF2)*(1.-COF3)
				CCOF(1,2)=(1.-COF2)*COF3
				CCOF(2,1)=COF2*(1.-COF3)
				CCOF(2,2)=COF2*COF3

				vari=varis(1,1)*CCOF(1,1)+varis(2,1)*CCOF(2,1)&
     				+varis(1,2)*CCOF(1,2)+varis(2,2)*CCOF(2,2)	 
    end subroutine

    SUBROUTINE dsvbksb(u,w,v,m,n,mp,np,b,x)
    implicit none
    INTEGER m,mp,n,np,NMAX
    REAL*8 :: b(mp),u(mp,np),v(np,np),w(np),x(np)
    PARAMETER (NMAX=500) !Maximum anticipated value of n.
    !Solves A*X = B for a vector X, where A is specified by the arrays u, w, v as returned by
    !svdcmp. m and n are the logical dimensions of a, and will be equal for square matrices. mp
    !and np are the physical dimensions of a. b(1:m) is the input right-hand side. x(1:n) is
    !the output solution vector. No input quantities are destroyed, so the routine may be called
    !sequentially with different b's.
    INTEGER i,j,jj
    REAL*8 :: s,tmp(NMAX),eps
    eps=1.d-10
    do j=1,n     !Calculate U^T*B.
        s=0.d0
        if(w(j).ne.0.d0)then    !Nonzero result only if wj is nonzero.
        ! if(dabs(w(j)).gt.eps)then
            do i=1,m
                s=s+u(i,j)*b(i)
            enddo 
            s=s/w(j)     !This is the divide by wj .
        endif
        tmp(j)=s
    enddo 
    do j=1,n    !Matrix multiply by V to get answer.
        s=0.d0
        do jj=1,n
            s=s+v(j,jj)*tmp(jj)
        enddo 
        x(j)=s
    enddo 
    return
    END SUBROUTINE dsvbksb
    
    subroutine SVD_solver(m,n,A,U,S,V)
        integer,intent(in) :: m,n
        real*8,intent(in) :: A(m,n)
        real*8,intent(out) :: U(m,n),S(n),V(n,n)
        !-------
        integer,parameter :: LWMAX = 1000
        integer :: LDA, LDU, LDVT
        integer :: INFO, LWORK
        real*8 :: VT(n,n),WORK(LWMAX),eps
        external DGESVD
        eps=1.d-8
        LDA=m
        LDU=m
        LDVT=n
    
        LWORK = -1
        CALL DGESVD( 'S', 'S', m, n, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
        LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
        CALL DGESVD( 'S', 'S', m, n, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )
        IF( INFO.GT.0 ) THEN
        !IF( INFO.GT.eps ) THEN
            WRITE(*,*)'The algorithm computing SVD failed to converge.'
            STOP
        END IF
    
        V=TRANSPOSE(VT)
    
    end subroutine SVD_solver    
    
    subroutine initial_reconstruction_B_cell(U,aB,wB,vB,agspt_weightB_point,Uxyz)          !rbfƬƬ  ų ʧ    ֵ      
    implicit none
    real*8,intent(in) :: U(3,7),aB(21,21),wB(21),vB(21,21),agspt_weightB_point(3,21)
    real*8,intent(out) :: Uxyz(3)
    !---------------
    integer :: im,jm
    real*8 :: b(21)
    integer :: m,n,mp,np
    real*8 :: gx(21),vari_max,vari_min
    
    m=21
    n=21
    mp=21
    np=21
    
    do jm=1,7             
    b((jm-1)*3+1)=U(1,jm)-U(1,1) !Phi_1
    b((jm-1)*3+2)=U(2,jm)-U(2,1) !Phi_2
    b((jm-1)*3+3)=U(3,jm)-U(3,1) !Phi_3
    end do
        gx=0.d0
        call dsvbksb(aB,wB,vB,m,n,mp,np,b,gx) !Now we can backsubstitute.
        !=============================================
        do im=1,3
		Uxyz(im)=U(im,1)		
		do jm=1,21		
            Uxyz(im)=Uxyz(im)+agspt_weightB_point(im,jm)*gx(jm)
        end do
        end do
!        do im=1,3
!            vari_max=max(u(im,1),u(im,2),u(im,3),u(im,4),u(im,5),u(im,6),u(im,7))
!            vari_min=min(u(im,1),u(im,2),u(im,3),u(im,4),u(im,5),u(im,6),u(im,7))
!            if(Uxyz(im).lt.vari_min)then
!                Uxyz(im)=vari_min
!            else
!                Uxyz(im)=min(vari_max,Uxyz(im))
!            end if
!        end do        
    end subroutine initial_reconstruction_B_cell     
    
    subroutine recon_cellrbf(W_,a,w,v,agspt_weight7_point,Wxyz)                           !rbfƬƬ       ֵ 
    implicit none
    real*8,intent(in) :: W_(7),a(7,7),w(7),v(7,7),agspt_weight7_point(7)
    real*8,intent(out) :: Wxyz
    !---------------
    integer :: irbf,jrbf
    real*8 :: b(7)
    integer :: m,n,mp,np,im
    real*8 :: gx(7)
    real*8 :: W_max,W_min,vari 
    real*8 :: delta_p1,delta_p2,delta_m
    real*8 :: Limit(6),Limit_Venkata,Tay 
    real*8 :: radi,theta,phi,vr,vt,vp
    
    m=7
    n=7
    mp=7
    np=7
    
        do jrbf=1,7
            b(jrbf)=W_(jrbf)-W_(1)
        end do
        gx=0.d0
        call dsvbksb(a,w,v,m,n,mp,np,b,gx) !Now we can backsubstitute.
        !=============================================
		Wxyz=W_(1)		
		do jrbf=1,7		
            Wxyz=Wxyz+agspt_weight7_point(jrbf)*gx(jrbf)
        end do                           
           
!        W_max=max(W_(1),W_(2),W_(3),W_(4),W_(5),W_(6),W_(7))          
!        W_min=min(W_(1),W_(2),W_(3),W_(4),W_(5),W_(6),W_(7))           
!        if(Wxyz .lt. W_min)then           
!            Wxyz=W_min          
!        else           
!            Wxyz=min(W_max,Wxyz)          
!        end if
        
    end subroutine recon_cellrbf   
    subroutine rbf_limiter(Wi,vary)
	implicit none
    real*8,intent(in) :: Wi(7)
    real*8,intent(inout) :: vary
	real*8 :: W_max,W_min
	W_max=max(Wi(1),Wi(2),Wi(3),Wi(4),Wi(5),Wi(6),Wi(7))                   
	W_min=min(Wi(1),Wi(2),Wi(3),Wi(4),Wi(5),Wi(6),Wi(7))  
	if(vary.lt.W_min)then            
        vary=W_min          
    else            
        vary=min(W_max,vary)          
    end if  
    end subroutine rbf_limiter
    
    subroutine sph(i,usph)
      implicit double precision(a-h,o-z)  
      real*8,dimension(17,0:N_th_p+1,0:N_ph_p+1) :: usph
	  integer,intent(in)::i
      j=0
      k=0
      do i1=0,N_th_P+1
         th1=THETA11(i1)
         do j1=0, N_ph_p+1
            ph1=phi11(j1)
            x1 =dsin(th1)*dcos(ph1)
            y1 =dsin(th1)*dsin(ph1)
            z1 =dcos(th1)
            detmin=100.   
            DO j=2,N_th_p-1
               DO k=2,N_ph_p-1
                  ead1=dsqrt(utotal(1,i,j,k)**2+utotal(2,i,j,k)**2+utotal(3,i,j,k)**2)
                  x_yin=utotal(1,i,j,k)/ead1
                  y_yin=utotal(2,i,j,k)/ead1
                  z_yin=utotal(3,i,j,k)/ead1
                  radd=dsqrt((x1-x_yin)**2+(y1-y_yin)**2+(z1-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     jBOU=j
                     kBOU=k
                  ENDIF
               ENDDO
            ENDDO
           j=jBOU
           k=kBOU
              x_yin =  utotal(1,i,j,k) 
              y_yin =  utotal(2,i,j,k) 
              z_yin =  utotal(3,i,j,k) 
              
              radii = dsqrt(x_yin**2.0+y_yin**2.0+z_yin**2.0)
              rxy = dsqrt(x_yin**2+y_yin**2)
              if(y_yin.gt.0.d0) then
                ph0 = dacos(x_yin/rxy)
              else
                ph0 = 2.d0*pi-dacos(x_yin/rxy)
              end if
              th0 = dacos(z_yin/radii)
              
              rho=utotal(4,i,j,k)
              vx=utotal(5,i,j,k)
              vy=utotal(6,i,j,k)
              vz=utotal(7,i,j,k)
			  pressure=utotal(8,i,j,k)
              Bx=utotal(9,i,j,k) 
              By=utotal(10,i,j,k)
              Bz=utotal(11,i,j,k)
              call xyz_rtp(th0,ph0,vx,vy,vz,vr,vth,vph)
              call xyz_rtp(th0,ph0,Bx,By,Bz,Br,Bth,Bph)
              usph(1,i1,j1)=x_yin
              usph(2,i1,j1)=y_yin
              usph(3,i1,j1)=z_yin
              usph(4,i1,j1)=rho
              usph(5,i1,j1)=vx
              usph(6,i1,j1)=vy
              usph(7,i1,j1)=vz
              usph(8,i1,j1)=Bx
              usph(9,i1,j1)=By
              usph(10,i1,j1)=Bz
              usph(11,i1,j1)=pressure
              usph(12,i1,j1)=vr
              usph(13,i1,j1)=vth
              usph(14,i1,j1)=vph
              usph(15,i1,j1)=Br
              usph(16,i1,j1)=Bth
              usph(17,i1,j1)=Bph
        enddo
     enddo

    end subroutine sph
 
    subroutine sph_interp(i,radius,usphinterp)
      implicit double precision(a-h,o-z)  
      real*8,intent(out) :: usphinterp(17,0:N_th_p+1,0:N_ph_p+1)
	  integer,intent(in)::i
      real*8,intent(in)  :: radius 
      integer :: ii,iBOU,jBOU,kBOU,ivari,jj,m,n,i1,j1
      real*8 :: pos_local(3),cen_nei(3,7),Vari_output(8,7),a_divB_rbf81(21,21),a_divB_rbf7(7,7)
	  real*8 :: aB(21,21),wB(21),vB(21,21),a(7,7),w(7),v(7,7)
      real*8 :: agspt_weightB_point(3,21),agspt_weight7_point(7),Vari_interp(11),bbrbf7_(7,7)
	  real*8 :: wrbf7(7),vrbf7(7,7),bbrbf_(21,21),wrbf(21),vrbf(21,21)  
      real*8 :: wg(7),wg_tot
      logical :: average,R_inverse
      R_inverse=.true.
      c_rbf=10.d0
      j=0
      k=0
      do i1=0,N_th_P+1
         th1=THETA11(i1)
         do j1=0, N_ph_p+1
            ph1=phi11(j1)

            x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)
            !call coordinate_rtp_xyz(x1,y1,z1,radius,th1,ph1)
            pos_local(1)=x1   
            pos_local(2)=y1 
            pos_local(3)=z1
            detmin=100.
            DO j=2,N_th_p-1
               DO k=2,N_ph_p-1
                   do ii=i-1,i+1
                  x_yin=utotal(1,ii,j,k)
                  y_yin=utotal(2,ii,j,k)
                  z_yin=utotal(3,ii,j,k)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     iBOU=ii
                     jBOU=j
                     kBOU=k
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                     !rad=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
                     !print *,"rad=",rad,"nyy=",nyy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU
                  ENDIF
                   end do                  
               ENDDO
            ENDDO             
            x_inter=pos_local(1)  
            y_inter=pos_local(2)
            z_inter=pos_local(3)            
           !radd=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
           !raxy=dsqrt(x_local**2.d0+y_local**2.d0)
           rad=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0+pos_local(3)**2.d0)
           raxy=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0)
           eps=1.d-8
           if(rad<eps .or. raxy<eps)then
               print *,"rad=",rad,"raxy=",raxy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU !,"rad=",rad
           end if
           !print *,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU,"detmin=",detmin
           ii=iBOU
           j=jBOU
           k=kBOU           
           th0=th1
           ph0=ph1
           !radii = dsqrt(x1**2.0+y1**2.0+z1**2.0)             
           !rxy = dsqrt(x1**2+y1**2)            
           !if(y1.gt.0.d0) then           
           !    ph0 = dacos(x1/rxy)            
           !else            
           !    ph0 = 2.d0*pi-dacos(x1/rxy)            
           !end if            
           !th0 = dacos(z1/radii)  
               cen_nei(1:3,1)=utotal(1:3,ii,j,k)
               Vari_output(1:8,1)=utotal(4:11,ii,j,k)
               average= .false.
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(ii+2,Rindex)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(ii+3,Rindex)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.
                       end if
                   end if
               else
                   ivari=min(ii+1,Rindex)
               end if
               cen_nei(1:3,2)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,2)=utotal(4:11,ivari,j,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(ii-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(ii-3,0)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.
                       end if
                   end if
               else
                   ivari=max(ii-1,0)
               end if
               cen_nei(1:3,3)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,3)=utotal(4:11,ivari,j,k)
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(j+2,N_th_p+1)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(j+3,N_th_p+1)                                         
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then                       
                           average = .true.  
                       end if
                   end if
               else
                   ivari=min(j+1,N_th_p+1)
               end if
               cen_nei(1:3,4)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,4)=utotal(4:11,ii,ivari,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(j-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(j-3,0)  
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then                      
                           average = .true.                        
                       end if                            
                   end if                 
               else
                   ivari=max(j-1,0)
               end if
               cen_nei(1:3,5)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,5)=utotal(4:11,ii,ivari,k)              
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(k+2,N_ph_p+1)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(k+3,N_ph_p+1)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.   
                       end if
                   end if 
               else
                   ivari=min(k+1,N_ph_p+1)
               end if
               cen_nei(1:3,6)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,6)=utotal(4:11,ii,j,ivari)      
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(k-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(k-3,0)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.   
                       end if
                   end if                    
               else
                   ivari=max(k-1,0)
               end if
               cen_nei(1:3,7)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,7)=utotal(4:11,ii,j,ivari)  
               !average=.true.  
               if(average)then
                   if(R_inverse)then
                       wg_tot=0.d0
                       do ivari=1,7
                           temp=(x1-cen_nei(1,ivari))**2.d0+(y1-cen_nei(2,ivari))**2.d0+&
                                (z1-cen_nei(3,ivari))**2.d0
                           wg(ivari)=1.d0/dsqrt(temp+eps)
                           wg_tot=wg_tot+wg(ivari)
                       end do
                       Vari_interp(4:11)=0.d0
                       do ivari=1,7
                           Vari_interp(4:11)=Vari_interp(4:11)+Vari_output(1:8,ivari)*wg(ivari)/wg_tot
                       end do
                   else
                   Vari_interp(4:11)=0.d0
                   do ivari=1,7
                       Vari_interp(4:11)=Vari_interp(4:11)+Vari_output(1:8,ivari)
                   end do
                   Vari_interp(4:11)=Vari_interp(4:11)/7.d0 
                   end if                   
               else               
			   davg=0.d0
			   do ivari=2,7
			   davg=davg+dsqrt((cen_nei(1,1)-cen_nei(1,ivari))**2.d0+(cen_nei(2,1)-cen_nei(2,ivari))**2.d0+&
			   (cen_nei(3,1)-cen_nei(3,ivari))**2.d0)
			   end do
			   davg=davg/6.d0
               a_divB_rbf81=0.d0
               do ivari=1,7  
                   do jj=1,7                	                   
                    basefun=c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-&
					cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0
                                              
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+1)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					(((cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					davg**4.d0)*basefun**(-1.5d0))                                                    
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+2)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					  (cen_nei(1,ivari)-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+2)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					(cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
                                     
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+2)=-((cen_nei(2,ivari)-cen_nei(2,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+3)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))	                                
                   end do
               end do   
               m=21 
               n=21 
               call SVD_solver(m,n,a_divB_rbf81,bbrbf_,wrbf,vrbf)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf(jj).gt.wmax)wmax=wrbf(jj)
               end do 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf(jj).lt.wmin)wrbf(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               end do                              
               aB=bbrbf_
               wB=wrbf
               vB=vrbf  
               
               a_divB_rbf7=0.d0
			   do ivari=1,7                 			  
                   do jj=1,7              
                       a_divB_rbf7(ivari,jj)=dsqrt(c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					   (cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0)             			   
                   end do
               end do
               m=7 
               n=7 
               call SVD_solver(m,n,a_divB_rbf7,bbrbf7_,wrbf7,vrbf7)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf7(jj).gt.wmax)wmax=wrbf7(jj)
               enddo 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf7(jj).lt.wmin)wrbf7(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               enddo                  
               a=bbrbf7_
               w=wrbf7
               v=vrbf7        
                              
               do jj=1,7 
               agspt_weight7_point(jj)=dsqrt(c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			     (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0) 
               
               basefun=c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0+&
			     (z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0
               
               agspt_weightB_point(1,(jj-1)*3+1)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			    (((y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(1,(jj-1)*3+2)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)         
               agspt_weightB_point(1,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(x_inter-cen_nei(1,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                            
               agspt_weightB_point(2,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(2,(jj-1)*3+2)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)        
               agspt_weightB_point(2,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                              
               agspt_weightB_point(3,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(3,(jj-1)*3+2)=-((y_inter-cen_nei(2,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)      
               agspt_weightB_point(3,(jj-1)*3+3)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)                              
              end do  
          
      call initial_reconstruction_B_cell(Vari_output(5:7,1:7),aB,wB,vB,agspt_weightB_point(1:3,1:21),Vari_interp(8:10))
	  do ivari=8,10
	      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))
	  end do
	  
      do ivari=4,7  !pri_variables   
          call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari))   
          call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))			  
      end do
      ivari=11
      call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari))   
      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))	
               end if
      !Vari_interp(4:11)=Vari_output(1:8,1)
      rho=Vari_interp(4)
      vx=Vari_interp(5)             
      vy=Vari_interp(6)             
      vz=Vari_interp(7)           
      Bx=Vari_interp(8)           
      By=Vari_interp(9)          
      Bz=Vari_interp(10)                   
      pressure=Vari_interp(11)
      call xyz_rtp(th0,ph0,vx,vy,vz,vr,vth,vph)
      call xyz_rtp(th0,ph0,Bx,By,Bz,Br,Bth,Bph)            
            
              usphinterp(1,i1,j1)=x1
              usphinterp(2,i1,j1)=y1
              usphinterp(3,i1,j1)=z1
              usphinterp(4,i1,j1)=rho
              usphinterp(5,i1,j1)=vx
              usphinterp(6,i1,j1)=vy
              usphinterp(7,i1,j1)=vz
              usphinterp(8,i1,j1)=Bx
              usphinterp(9,i1,j1)=By
              usphinterp(10,i1,j1)=Bz
              usphinterp(11,i1,j1)=pressure
              usphinterp(12,i1,j1)=vr
              usphinterp(13,i1,j1)=vth
              usphinterp(14,i1,j1)=vph
              usphinterp(15,i1,j1)=Br
              usphinterp(16,i1,j1)=Bth
              usphinterp(17,i1,j1)=Bph
        enddo
     enddo

    end subroutine sph_interp
  
    subroutine sph_interplong(ph1,Nradius,radi,theta,usphinterplong)
      implicit double precision(a-h,o-z)  
      real*8,intent(in)::ph1
	  integer,intent(in):: Nradius
      real*8,intent(out) :: usphinterplong(1:17,1:Nradius,0:N_th_p+1)
      real*8,intent(in)  :: radi(0:nrt),theta(0:N_th_p+1)
      integer :: ii,iBOU,jBOU,kBOU,ivari,jj,m,n,i1,j1
      real*8 :: radius,pos_local(3),cen_nei(3,7),Vari_output(8,7),a_divB_rbf81(21,21)
	  real*8 :: a_divB_rbf7(7,7),aB(21,21),wB(21),vB(21,21),a(7,7),w(7),v(7,7)
      real*8 :: agspt_weightB_point(3,21),agspt_weight7_point(7),Vari_interp(11)
	  real*8 :: bbrbf7_(7,7),wrbf7(7),vrbf7(7,7),bbrbf_(21,21),wrbf(21),vrbf(21,21)
	  real*8 :: eps
      real*8 :: wg(7),wg_tot
      logical :: writemultiPoint,average,R_inverse
      R_inverse=.true.
	  writemultiPoint=.true.
      c_rbf=10.d0
      !ph1=phi11(kk)
      do j1=0,N_th_P+1
         th1=theta(j1)
         do i1=1, Nradius
            radius=radi(i1)  
			x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)
            !call coordinate_rtp_xyz(x1,y1,z1,radius,th1,ph1)
            pos_local(1)=x1   
            pos_local(2)=y1 
            pos_local(3)=z1
            detmin=100.
            DO j=2,N_th_p-1
               DO k=2,N_ph_p-1
                   do ii=i1-1,i1+1
                  x_yin=utotal(1,ii,j,k)
                  y_yin=utotal(2,ii,j,k)
                  z_yin=utotal(3,ii,j,k)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     iBOU=ii
                     jBOU=j
                     kBOU=k
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                     !rad=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
                     !print *,"rad=",rad,"nyy=",nyy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU
                  ENDIF
                   end do                  
               ENDDO
            ENDDO             
            x_inter=pos_local(1)  
            y_inter=pos_local(2)
            z_inter=pos_local(3)            
           !radd=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
           !raxy=dsqrt(x_local**2.d0+y_local**2.d0)
           rad=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0+pos_local(3)**2.d0)
           raxy=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0)
           eps=1.d-8
           if(rad<eps .or. raxy<eps)then
               print *,"rad=",rad,"raxy=",raxy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU !,"rad=",rad
           end if	   
!!		   if(abs(real(iBOU-Nradius))<2.0 .and. detmin>1.d-1)then
!!               print *,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU,"detmin=",detmin
!!		   end if
           ii=iBOU
           j=jBOU
           k=kBOU           
           th0=th1
           ph0=ph1
           !radii = dsqrt(x1**2.0+y1**2.0+z1**2.0)             
           !rxy = dsqrt(x1**2+y1**2)            
           !if(y1.gt.0.d0) then           
           !    ph0 = dacos(x1/rxy)            
           !else            
           !    ph0 = 2.d0*pi-dacos(x1/rxy)            
           !end if            
           !th0 = dacos(z1/radii)  
               cen_nei(1:3,1)=utotal(1:3,ii,j,k)
               Vari_output(1:8,1)=utotal(4:11,ii,j,k)
    if(radius<1.0d0)then
		       Vari_interp(4:11)=Vari_output(1:8,1)
    else
               average=.false.
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(ii+2,Rindex)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(ii+3,Rindex)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.
                           if(writemultiPoint)write(*,*)'ii+3=',ivari,'ii=',ii,'j=',j,'k=',k
                       end if
                   end if
               else
                   ivari=min(ii+1,Rindex)
               end if
               cen_nei(1:3,2)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,2)=utotal(4:11,ivari,j,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(ii-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(ii-3,0)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ivari,j,k)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.
                           if(writemultiPoint)write(*,*)'ii-3=',ivari,'ii=',ii,'j=',j,'k=',k
                       end if
                   end if
               else
                   ivari=max(ii-1,0)
               end if
               cen_nei(1:3,3)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,3)=utotal(4:11,ivari,j,k)
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(j+2,N_th_p+1)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(j+3,N_th_p+1)                                         
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then                       
                           average = .true.                       
                           if(writemultiPoint)write(*,*)'j+3=',ivari,'ii=',ii,'j=',j,'k=',k
                       end if
                   end if
               else
                   ivari=min(j+1,N_th_p+1)
               end if
               cen_nei(1:3,4)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,4)=utotal(4:11,ii,ivari,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(j-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(j-3,0)  
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,ivari,k)-utotal(1:3,ii,j,k)))<eps)then                      
                           average = .true.                     
                           if(writemultiPoint)write(*,*)'j-3=',ivari,'ii=',ii,'j=',j,'k=',k                  
                       end if                            
                   end if                 
               else
                   ivari=max(j-1,0)
               end if
               cen_nei(1:3,5)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,5)=utotal(4:11,ii,ivari,k)              
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(k+2,N_ph_p+1)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=min(k+3,N_ph_p+1)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.   
                           if(writemultiPoint)write(*,*)'k+3=',ivari,'ii=',ii,'j=',j,'k=',k
                       end if
                   end if 
               else
                   ivari=min(k+1,N_ph_p+1)
               end if
               cen_nei(1:3,6)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,6)=utotal(4:11,ii,j,ivari)      
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(k-2,0)
                   if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                   utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                       ivari=max(k-3,0)
                       if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k),&
                       utotal(1:3,ii,j,ivari)-utotal(1:3,ii,j,k)))<eps)then
                           average = .true.   
                           if(writemultiPoint)write(*,*)'k-3=',ivari,'ii=',ii,'j=',j,'k=',k
                       end if
                   end if                    
               else
                   ivari=max(k-1,0)
               end if
               cen_nei(1:3,7)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,7)=utotal(4:11,ii,j,ivari) 
               !average=.true.
               if(average)then
                   if(R_inverse)then
                       wg_tot=0.d0
                       do ivari=1,7
                           temp=(x1-cen_nei(1,ivari))**2.d0+(y1-cen_nei(2,ivari))**2.d0+&
                                (z1-cen_nei(3,ivari))**2.d0
                           wg(ivari)=1.d0/dsqrt(temp+eps)
                           wg_tot=wg_tot+wg(ivari)
                       end do
                       Vari_interp(4:11)=0.d0
                       do ivari=1,7
                           Vari_interp(4:11)=Vari_interp(4:11)+Vari_output(1:8,ivari)*wg(ivari)/wg_tot
                       end do
                   else
                   Vari_interp(4:11)=0.d0
                   do ivari=1,7
                       Vari_interp(4:11)=Vari_interp(4:11)+Vari_output(1:8,ivari)
                   end do
                   Vari_interp(4:11)=Vari_interp(4:11)/7.d0 
                   end if
               else
			   davg=0.d0
			   do ivari=2,7
			   davg=davg+dsqrt((cen_nei(1,1)-cen_nei(1,ivari))**2.d0+&
			   (cen_nei(2,1)-cen_nei(2,ivari))**2.d0+(cen_nei(3,1)-cen_nei(3,ivari))**2.d0)
			   end do
			   davg=davg/6.d0			 
               a_divB_rbf81=0.d0
               do ivari=1,7  
                   do jj=1,7                	                   
                    basefun=c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-&
					cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0
                                              
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+1)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))                                                    
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+2)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					  (cen_nei(1,ivari)-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+2)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
                                     
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+2)=-((cen_nei(2,ivari)-cen_nei(2,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+3)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))	                                
                   end do
               end do   
               m=21 
               n=21 
               call SVD_solver(m,n,a_divB_rbf81,bbrbf_,wrbf,vrbf)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf(jj).gt.wmax)wmax=wrbf(jj)
               end do 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf(jj).lt.wmin)wrbf(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               end do                              
               aB=bbrbf_
               wB=wrbf
               vB=vrbf  
               
               a_divB_rbf7=0.d0
			   do ivari=1,7                 			  
                   do jj=1,7              
                       a_divB_rbf7(ivari,jj)=dsqrt(c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					     (cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0)             			   
                   end do
               end do
               m=7 
               n=7 
               call SVD_solver(m,n,a_divB_rbf7,bbrbf7_,wrbf7,vrbf7)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf7(jj).gt.wmax)wmax=wrbf7(jj)
               enddo 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf7(jj).lt.wmin)wrbf7(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               enddo                  
               a=bbrbf7_
               w=wrbf7
               v=vrbf7        
                              
               do jj=1,7 
               agspt_weight7_point(jj)=dsqrt(c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			     (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0) 
               
               basefun=c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			     (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0
               
               agspt_weightB_point(1,(jj-1)*3+1)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(1,(jj-1)*3+2)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)         
               agspt_weightB_point(1,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(x_inter-cen_nei(1,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                            
               agspt_weightB_point(2,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(2,(jj-1)*3+2)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)        
               agspt_weightB_point(2,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                              
               agspt_weightB_point(3,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(3,(jj-1)*3+2)=-((y_inter-cen_nei(2,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)      
               agspt_weightB_point(3,(jj-1)*3+3)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)                              
              end do  
          
      call initial_reconstruction_B_cell(Vari_output(5:7,1:7),aB,wB,vB,agspt_weightB_point(1:3,1:21),Vari_interp(8:10))
	  do ivari=8,10
	      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))
	  end do
      do ivari=4,7  !pri_variables   
          call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari)) 
          call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))		  
      end do
      ivari=11
      call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari))  
      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))	  
               end if 
    end if	  
      !Vari_interp(4:11)=Vari_output(1:8,1)
      rho=Vari_interp(4)
      vx=Vari_interp(5)             
      vy=Vari_interp(6)             
      vz=Vari_interp(7)           
      Bx=Vari_interp(8)           
      By=Vari_interp(9)          
      Bz=Vari_interp(10)                   
      pressure=Vari_interp(11)
      call xyz_rtp(th0,ph0,vx,vy,vz,vr,vth,vph)
      call xyz_rtp(th0,ph0,Bx,By,Bz,Br,Bth,Bph)
              

              usphinterplong(1,i1,j1)=x1
              usphinterplong(2,i1,j1)=y1
              usphinterplong(3,i1,j1)=z1
              usphinterplong(4,i1,j1)=rho
              usphinterplong(5,i1,j1)=vx
              usphinterplong(6,i1,j1)=vy
              usphinterplong(7,i1,j1)=vz
              usphinterplong(8,i1,j1)=Bx
              usphinterplong(9,i1,j1)=By
              usphinterplong(10,i1,j1)=Bz
              usphinterplong(11,i1,j1)=pressure
              usphinterplong(12,i1,j1)=vr
              usphinterplong(13,i1,j1)=vth
              usphinterplong(14,i1,j1)=vph
              usphinterplong(15,i1,j1)=Br
              usphinterplong(16,i1,j1)=Bth
              usphinterplong(17,i1,j1)=Bph
        enddo
     enddo

    end subroutine sph_interplong
  
subroutine sph_interpradial(ph1,Nradius,radi,th1,radiusphinterp)
      implicit double precision(a-h,o-z)  
      real*8,intent(in)::ph1,th1
	  integer,intent(in) :: Nradius
      real*8,intent(out) :: radiusphinterp(1:17,1:Nradius)
      real*8,intent(in)  :: radi(0:nrt)
      integer :: ii,iBOU,jBOU,kBOU,ivari,jj,m,n,i1
      real*8 :: pos_local(3),cen_nei(3,7),Vari_output(8,7),a_divB_rbf81(21,21)
	  real*8 :: a_divB_rbf7(7,7),aB(21,21),wB(21),vB(21,21),a(7,7),w(7),v(7,7)
      real*8 :: agspt_weightB_point(3,21),agspt_weight7_point(7),Vari_interp(11)
	  real*8 :: bbrbf7_(7,7),wrbf7(7),vrbf7(7,7),bbrbf_(21,21),wrbf(21),vrbf(21,21)
      c_rbf=10.d0
         do i1=1, Nradius
            radius=radi(i1)
            x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)
            !call coordinate_rtp_xyz(x1,y1,z1,radius,th1,ph1)             
            pos_local(1)=x1   
            pos_local(2)=y1 
            pos_local(3)=z1
            detmin=100.
            DO j=2,N_th_p-1
               DO k=2,N_ph_p-1
                   do ii=i1-1,i1+1
                  x_yin=utotal(1,ii,j,k)
                  y_yin=utotal(2,ii,j,k)
                  z_yin=utotal(3,ii,j,k)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     iBOU=ii
                     jBOU=j
                     kBOU=k
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                     !rad=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
                     !print *,"rad=",rad,"nyy=",nyy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU
                  ENDIF
                   end do                  
               ENDDO
            ENDDO             
            x_inter=pos_local(1)  
            y_inter=pos_local(2)
            z_inter=pos_local(3)
           !radd=dsqrt(x_local**2.d0+y_local**2.d0+z_local**2.d0)
           !raxy=dsqrt(x_local**2.d0+y_local**2.d0)
           rad=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0+pos_local(3)**2.d0)
           raxy=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0)
           eps=1.d-8
           if(rad<eps .or. raxy<eps)then
               print *,"rad=",rad,"raxy=",raxy,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU !,"rad=",rad
           end if
           !print *,"iBOU=",iBOU,"jBOU=",jBOU,"kBOU=",kBOU,"detmin=",detmin
           
           ii=iBOU
           j=jBOU
           k=kBOU           
           th0=th1
           ph0=ph1
           !radii = dsqrt(x1**2.0+y1**2.0+z1**2.0)             
           !rxy = dsqrt(x1**2+y1**2)            
           !if(y1.gt.0.d0) then           
           !    ph0 = dacos(x1/rxy)            
           !else            
           !    ph0 = 2.d0*pi-dacos(x1/rxy)            
           !end if            
           !th0 = dacos(z1/radii)                   
           
               cen_nei(1:3,1)=utotal(1:3,ii,j,k)
               Vari_output(1:8,1)=utotal(4:11,ii,j,k)
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii+1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(ii+2,Rindex)
               else
                   ivari=min(ii+1,Rindex)
               end if
               cen_nei(1:3,2)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,2)=utotal(4:11,ivari,j,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii-1,j,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(ii-2,0)
               else
                   ivari=max(ii-1,0)
               end if
               cen_nei(1:3,3)=utotal(1:3,ivari,j,k)
               Vari_output(1:8,3)=utotal(4:11,ivari,j,k)
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(j+2,N_th_p+1)
               else
                   ivari=min(j+1,N_th_p+1)
               end if
               cen_nei(1:3,4)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,4)=utotal(4:11,ii,ivari,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(j-2,0)
               else
                   ivari=max(j-1,0)
               end if
               cen_nei(1:3,5)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,5)=utotal(4:11,ii,ivari,k)              
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(k+2,N_ph_p+1)
               else
                   ivari=min(k+1,N_ph_p+1)
               end if
               cen_nei(1:3,6)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,6)=utotal(4:11,ii,j,ivari)      
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k),&
			     utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(k-2,0)
               else
                   ivari=max(k-1,0)
               end if
               cen_nei(1:3,7)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,7)=utotal(4:11,ii,j,ivari)     

			   davg=0.d0
			   do ivari=2,7
			   davg=davg+dsqrt((cen_nei(1,1)-cen_nei(1,ivari))**2.d0+(cen_nei(2,1)-&
			     cen_nei(2,ivari))**2.d0+(cen_nei(3,1)-cen_nei(3,ivari))**2.d0)
			   end do
			   davg=davg/6.d0
               a_divB_rbf81=0.d0
               do ivari=1,7  
                   do jj=1,7                	                   
                    basefun=c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0
                                              
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+1)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					 davg**4.d0)*basefun**(-1.5d0))                                                    
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+2)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					 (cen_nei(1,ivari)-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+2)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					 davg**4.d0)*basefun**(-1.5d0))                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
                                     
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+2)=-((cen_nei(2,ivari)-cen_nei(2,jj))*&
					 (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+3)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))	                                
                   end do
               end do   
               m=21 
               n=21 
               call SVD_solver(m,n,a_divB_rbf81,bbrbf_,wrbf,vrbf)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf(jj).gt.wmax)wmax=wrbf(jj)
               end do 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf(jj).lt.wmin)wrbf(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               end do                              
               aB=bbrbf_
               wB=wrbf
               vB=vrbf  
               
               a_divB_rbf7=0.d0
			   do ivari=1,7                 			  
                   do jj=1,7              
                       a_divB_rbf7(ivari,jj)=dsqrt(c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					   (cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0)             			   
                   end do
               end do
               m=7 
               n=7 
               call SVD_solver(m,n,a_divB_rbf7,bbrbf7_,wrbf7,vrbf7)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf7(jj).gt.wmax)wmax=wrbf7(jj)
               enddo 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf7(jj).lt.wmin)wrbf7(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               enddo                  
               a=bbrbf7_
               w=wrbf7
               v=vrbf7        
                              
               do jj=1,7 
               agspt_weight7_point(jj)=dsqrt(c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			   (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0) 
               
               basefun=c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0+&
			   (z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0
               
               agspt_weightB_point(1,(jj-1)*3+1)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+(((y_inter-cen_nei(2,jj))**2.d0&
				                                       +(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(1,(jj-1)*3+2)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)         
               agspt_weightB_point(1,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(x_inter-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0) 
                            
               agspt_weightB_point(2,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(2,(jj-1)*3+2)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+(((x_inter-cen_nei(1,jj))**2.d0&
				                                       +(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)        
               agspt_weightB_point(2,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(y_inter-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0) 
                              
               agspt_weightB_point(3,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(z_inter-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(3,(jj-1)*3+2)=-((y_inter-cen_nei(2,jj))*(z_inter-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)      
               agspt_weightB_point(3,(jj-1)*3+3)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+(((x_inter-cen_nei(1,jj))**2.d0&
				                                       +(y_inter-cen_nei(2,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)                              
              end do  
          
      call initial_reconstruction_B_cell(Vari_output(5:7,1:7),aB,wB,vB,agspt_weightB_point(1:3,1:21),Vari_interp(8:10))
	  do ivari=8,10
	      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))
	  end do
      do ivari=4,7  !pri_variables   
          call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari)) 
          call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))			  
      end do
      ivari=11
      call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari)) 
      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))		  
      !Vari_interp(4:11)=Vari_output(1:8,1)
      rho=Vari_interp(4)
      vx=Vari_interp(5)             
      vy=Vari_interp(6)             
      vz=Vari_interp(7)           
      Bx=Vari_interp(8)           
      By=Vari_interp(9)          
      Bz=Vari_interp(10)                  
      pressure=Vari_interp(11)
      call xyz_rtp(th0,ph0,vx,vy,vz,vr,vth,vph)
      call xyz_rtp(th0,ph0,Bx,By,Bz,Br,Bth,Bph)
              radiusphinterp(1,i1)=x1
              radiusphinterp(2,i1)=y1
              radiusphinterp(3,i1)=z1      
              radiusphinterp(4,i1)=rho
              radiusphinterp(5,i1)=vx
              radiusphinterp(6,i1)=vy
              radiusphinterp(7,i1)=vz
              radiusphinterp(8,i1)=Bx
              radiusphinterp(9,i1)=By
              radiusphinterp(10,i1)=Bz 
              radiusphinterp(11,i1)=pressure
              radiusphinterp(12,i1)=vr
              radiusphinterp(13,i1)=vth
              radiusphinterp(14,i1)=vph
              radiusphinterp(15,i1)=Br
              radiusphinterp(16,i1)=Bth
              radiusphinterp(17,i1)=Bph
         enddo
    end subroutine sph_interpradial
  
subroutine transtouniformgrid(radi,theta,phi,usph_ori,usph_out)
      implicit double precision(a-h,o-z)  
      real*8,intent(in)::usph_ori(1:11,0:Rindex,0:N_th_p+1,0:N_ph_p+1)
      real*8,intent(in)  :: radi(0:nrt),theta(0:N_th_p+1),phi(0,N_ph_p+1)
	  real*8,intent(out) :: usph_out(1:11,0:nrt,0:N_th_p+1,0:N_ph_p+1)
      integer :: ii,iBOU,jBOU,kBOU,ivari,jj,m,n,i1,j1,k1
      real*8 :: radius,pos_local(3),cen_nei(3,7),Vari_output(8,7),a_divB_rbf81(21,21)
	  real*8 :: a_divB_rbf7(7,7),aB(21,21),wB(21),vB(21,21),a(7,7),w(7),v(7,7)
      real*8 :: agspt_weightB_point(3,21),agspt_weight7_point(7),Vari_interp(11)
	  real*8 :: bbrbf7_(7,7),wrbf7(7),vrbf7(7,7),bbrbf_(21,21),wrbf(21),vrbf(21,21)
	  real*8 :: eps
      c_rbf=10.d0
      do j1=0,N_th_P+1
	    th1=theta(j1)
	    do k1=0, N_ph_p+1
         ph1=phi11(k1)
         do i1=0,nrt
            radius=radi(i1)  
			x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)
            !call coordinate_rtp_xyz(x1,y1,z1,radius,th1,ph1)
            pos_local(1)=x1   
            pos_local(2)=y1 
            pos_local(3)=z1
            detmin=100.
			j_begin=max(j1-2,1)
			j_end=min(j1+2,N_th_p)
			k_begin=max(k1-2,1)
			k_end=min(k1+2,N_ph_p)	
			if(j1<2 .or. j1>N_th_P-1)then
			  k_begin=1
			  k_end=N_ph_p
			end if
            DO j=j_begin,j_end
               DO k=k_begin,k_end		   
                   do ii=max(i1-1,0),min(i1+1,Rindex)
                  x_yin=usph_ori(1,ii,j,k)
                  y_yin=usph_ori(2,ii,j,k)
                  z_yin=usph_ori(3,ii,j,k)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     iBOU=ii
                     jBOU=j
                     kBOU=k
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                  ENDIF
				  !if(detmin<0.75*dabs(radi(min(i1+1,nrt))-radi(i1)))exit
                   end do                  
               ENDDO
            ENDDO   
!!          if(j1==int(N_th_P/3) .and. k1==int(N_ph_p/3))write(*,*)'detmin=',detmin,'i1=',i1,'j1=',j1,'k1=',k1
          !write(*,*)'detmin=',detmin,'i1=',i1,'j1=',j1,'k1=',k1		  
            x_inter=pos_local(1)  
            y_inter=pos_local(2)
            z_inter=pos_local(3)            
           rad=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0+pos_local(3)**2.d0)
           raxy=dsqrt(pos_local(1)**2.d0+pos_local(2)**2.d0)
           eps=1.d-8
           ii=iBOU
           j=jBOU
           k=kBOU           
           th0=th1
           ph0=ph1 
               cen_nei(1:3,1)=usph_ori(1:3,ii,j,k)
               Vari_output(1:8,1)=usph_ori(4:11,ii,j,k)
    if(radius<1.0d0 .or. ii>=Rindex-1)then
		       Vari_interp(4:11)=Vari_output(1:8,1)
	else
               if(dsqrt(DOT_PRODUCT(usph_ori(1:3,ii+1,j,k)-usph_ori(1:3,ii,j,k),&
			   usph_ori(1:3,ii+1,j,k)-usph_ori(1:3,ii,j,k)))<eps)then
                   ivari=min(ii+2,Rindex)
               else
                   ivari=min(ii+1,Rindex)
               end if
               cen_nei(1:3,2)=usph_ori(1:3,ivari,j,k)
               Vari_output(1:8,2)=usph_ori(4:11,ivari,j,k)
               if(dsqrt(DOT_PRODUCT(usph_ori(1:3,ii-1,j,k)-usph_ori(1:3,ii,j,k),&
			   usph_ori(1:3,ii-1,j,k)-usph_ori(1:3,ii,j,k)))<eps)then
                   ivari=max(ii-2,0)
               else
                   ivari=max(ii-1,0)
               end if
               cen_nei(1:3,3)=usph_ori(1:3,ivari,j,k)
               Vari_output(1:8,3)=usph_ori(4:11,ivari,j,k)
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii,j+1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(j+2,N_th_p+1)
               else
                   ivari=min(j+1,N_th_p+1)
               end if
               cen_nei(1:3,4)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,4)=utotal(4:11,ii,ivari,k)
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii,j-1,k)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(j-2,0)
               else
                   ivari=max(j-1,0)
               end if
               cen_nei(1:3,5)=utotal(1:3,ii,ivari,k)
               Vari_output(1:8,5)=utotal(4:11,ii,ivari,k)              
               
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii,j,k+1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=min(k+2,N_ph_p+1)
               else
                   ivari=min(k+1,N_ph_p+1)
               end if
               cen_nei(1:3,6)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,6)=utotal(4:11,ii,j,ivari)      
               if(dsqrt(DOT_PRODUCT(utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k),&
			   utotal(1:3,ii,j,k-1)-utotal(1:3,ii,j,k)))<eps)then
                   ivari=max(k-2,0)
               else
                   ivari=max(k-1,0)
               end if
               cen_nei(1:3,7)=utotal(1:3,ii,j,ivari)
               Vari_output(1:8,7)=utotal(4:11,ii,j,ivari)        

			   davg=0.d0
			   do ivari=2,7
			   davg=davg+dsqrt((cen_nei(1,1)-cen_nei(1,ivari))**2.d0+(cen_nei(2,1)-&
			   cen_nei(2,ivari))**2.d0+(cen_nei(3,1)-cen_nei(3,ivari))**2.d0)
			   end do
			   davg=davg/6.d0			 
               a_divB_rbf81=0.d0
               do ivari=1,7  
                   do jj=1,7                	                   
                    basefun=c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-&
					cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0
                                              
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+1)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					 davg**4.d0)*basefun**(-1.5d0))                                                    
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+2)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					 (cen_nei(1,ivari)-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+2)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					 davg**4.d0)*basefun**(-1.5d0))                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					 (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
                                     
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					 (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+2)=-((cen_nei(2,ivari)-cen_nei(2,jj))*&
					 (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+3)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					 (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0)/&
					 davg**4.d0)*basefun**(-1.5d0))	                                
                   end do
               end do   
               m=21 
               n=21 
               call SVD_solver(m,n,a_divB_rbf81,bbrbf_,wrbf,vrbf)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf(jj).gt.wmax)wmax=wrbf(jj)
               end do 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf(jj).lt.wmin)wrbf(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               end do                              
               aB=bbrbf_
               wB=wrbf
               vB=vrbf  
               
               a_divB_rbf7=0.d0
			   do ivari=1,7                 			  
                   do jj=1,7              
                       a_divB_rbf7(ivari,jj)=dsqrt(c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					   (cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0)             			   
                   end do
               end do
               m=7 
               n=7 
               call SVD_solver(m,n,a_divB_rbf7,bbrbf7_,wrbf7,vrbf7)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf7(jj).gt.wmax)wmax=wrbf7(jj)
               enddo 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf7(jj).lt.wmin)wrbf7(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               enddo                  
               a=bbrbf7_
               w=wrbf7
               v=vrbf7        
                              
               do jj=1,7 
               agspt_weight7_point(jj)=dsqrt(c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			    (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0) 
               
               basefun=c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0+&
			    (z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0
               
               agspt_weightB_point(1,(jj-1)*3+1)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			    (((y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(1,(jj-1)*3+2)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			    davg**4.d0)*basefun**(-1.5d0)         
               agspt_weightB_point(1,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(x_inter-cen_nei(1,jj))/&
			    davg**4.d0)*basefun**(-1.5d0) 
                            
               agspt_weightB_point(2,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			    davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(2,(jj-1)*3+2)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			    (((x_inter-cen_nei(1,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)        
               agspt_weightB_point(2,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(y_inter-cen_nei(2,jj))/&
			    davg**4.d0)*basefun**(-1.5d0) 
                              
               agspt_weightB_point(3,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(z_inter-cen_nei(3,jj))/&
			    davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(3,(jj-1)*3+2)=-((y_inter-cen_nei(2,jj))*(z_inter-cen_nei(3,jj))/&
			    davg**4.d0)*basefun**(-1.5d0)      
               agspt_weightB_point(3,(jj-1)*3+3)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			    (((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)                              
              end do  
          
      call initial_reconstruction_B_cell(Vari_output(5:7,1:7),aB,wB,vB,agspt_weightB_point(1:3,1:21),Vari_interp(8:10))
	  do ivari=8,10
	      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))
	  end do
      do ivari=4,7  !pri_variables   
          call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari)) 
          call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))		  
      end do
      ivari=11
      call recon_cellrbf(Vari_output(ivari-3,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari))  
      call rbf_limiter(Vari_output(ivari-3,1:7),Vari_interp(ivari))	  
    end if	  
	  usph_out(1,i1,j1,k1)=x1
	  usph_out(2,i1,j1,k1)=y1
	  usph_out(3,i1,j1,k1)=z1
	  usph_out(4:11,i1,j1,k1)=Vari_interp(4:11)
        enddo
     enddo
	end do
end subroutine transtouniformgrid
  
    subroutine StructuredFormat(i,radius,usphinterp)
      implicit double precision(a-h,o-z)  
      real*8,intent(out) :: usphinterp(11,0:N_th_p+1,0:N_ph_p+1)
	  integer,intent(in)::i
      real*8,intent(in)  :: radius
      integer :: iBOU,i1,j1,ii1,jj1,index
      real*8 :: pos_local(3),th1,ph1,x1,y1,z1,x_local,y_local,z_local,detmin,x_yin,y_yin,z_yin,radd
      do i1=0,N_th_P+1
         th1=THETA11(i1)
         do j1=0, N_ph_p+1
            ph1=phi11(j1)
            x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)              
            pos_local(1)=x1    
            pos_local(2)=y1
            pos_local(3)=z1
            x_local=0.0
            y_local=0.0
            z_local=0.0
            iBOU=0
            detmin=100.
            DO jj1=1,PerCellLayer
                do ii1=max(i-1,0),min(i+1,R_IndMax)
                  index=ii1*PerCellLayer+jj1
                  x_yin=Cell_Pos(1,index)
                  y_yin=Cell_Pos(2,index)
                  z_yin=Cell_Pos(3,index)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
                     iBOU=index
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                  ENDIF
                end do                  
            ENDDO
            temp=dsqrt(x_local**2.0+y_local**2.0+z_local**2.0)
            if(temp<0.9)write(*,*)'iBOU=',iBOU,'index=',index
            usphinterp(1,i1,j1)=x_local
            usphinterp(2,i1,j1)=y_local
            usphinterp(3,i1,j1)=z_local
            usphinterp(4:11,i1,j1)=Cell_state(1:8,iBOU)
         end do
      end do
    end subroutine StructuredFormat
	
    subroutine R_inverse_Avg(x1,y1,z1,cennei,Varinei,Variinterp)	
	implicit none
    real*8,intent(in) :: x1,y1,z1
	real*8,intent(in) :: cennei(:,:),Varinei(:,:)
    real*8,intent(out) :: Variinterp(:)
	real*8 :: wg_tot,temp,eps
	real*8,allocatable :: wg(:)
	integer :: n_vari,n_nei,n_dim,ivari,jvari
	eps=1.d-12
	n_dim=size(cennei,1)
	n_nei=size(cennei,2)
	n_vari=size(Varinei,1)
	allocate(wg(1:n_nei))
	if(n_dim .ne. 3)write(*,*)"Warning! size function error"
	wg_tot=0.d0
    do ivari=1,n_nei
        temp=(x1-cennei(1,ivari))**2.d0+(y1-cennei(2,ivari))**2.d0+&
              (z1-cennei(3,ivari))**2.d0
        wg(ivari)=1.d0/dsqrt(temp+eps)
        wg_tot=wg_tot+wg(ivari)
    end do
    Variinterp(1:n_vari)=0.d0
    do ivari=1,n_vari
	    do jvari=1,n_nei
            Variinterp(ivari)=Variinterp(ivari)+Varinei(ivari,jvari)*wg(jvari)/wg_tot
		end do
    end do	
	deallocate(wg)
	end subroutine R_inverse_Avg
	
    subroutine UniformStructuredFormat(i,SortorNot,layerindex,radius,usphinterp)
      implicit double precision(a-h,o-z)  
	  integer,intent(in):: i
	  logical,intent(in) :: SortorNot
	  integer,intent(inout) :: layerindex(1:5,0:N_th_p+1,0:N_ph_p+1)
      real*8,intent(in)  :: radius
	  real*8,intent(out) :: usphinterp(11,0:N_th_p+1,0:N_ph_p+1)
      integer :: iBOU(7),rBOU,ii1,jj1,index
      real*8 :: th1,ph1,x1,y1,z1,x_local,y_local,z_local,detmin,x_yin,y_yin,z_yin,radd
	  integer :: ii,jBOU,kBOU,ivari,jj,m,n,i1,j1
      real*8 :: pos_local(3),cen_nei(3,7),Vari_output(8,7),a_divB_rbf81(21,21)
	  real*8 :: a_divB_rbf7(7,7),aB(21,21),wB(21),vB(21,21),a(7,7),w(7),v(7,7)
      real*8 :: agspt_weightB_point(3,21),agspt_weight7_point(7),Vari_interp(11)
	  real*8 :: bbrbf7_(7,7),wrbf7(7),vrbf7(7,7),bbrbf_(21,21),wrbf(21),vrbf(21,21)
	  c_rbf=10.d0
      do i1=0,N_th_P+1
         th1=THETA11(i1)
         do j1=0, N_ph_p+1
            ph1=phi11(j1)
            x1 =radius*dsin(th1)*dcos(ph1)
            y1 =radius*dsin(th1)*dsin(ph1)
            z1 =radius*dcos(th1)              
            pos_local(1)=x1    
            pos_local(2)=y1
            pos_local(3)=z1
		  if(SortorNot .eqv. .false.)then
		    iBOU(1:5)=layerindex(1:5,i1,j1)
		  else			
            x_local=0.0
            y_local=0.0
            z_local=0.0
            iBOU(1:7)=0
			rBOU=0
            detmin=100.
            DO jj1=1,PerCellLayer
                do ii1=max(i-1,0),min(i+1,R_IndMax)
                  index=ii1*PerCellLayer+jj1
                  x_yin=Cell_Pos(1,index)
                  y_yin=Cell_Pos(2,index)
                  z_yin=Cell_Pos(3,index)
                  radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
                  if(radd<detmin) then
                     detmin=radd
					 !central cell of stencil
                     iBOU(1)=index
					 rBOU=ii1
                     x_local=x_yin
                     y_local=y_yin
                     z_local=z_yin
                  ENDIF
                end do                  
            ENDDO
			ii1=rBOU
			 !1th neighbor
			detmin=100.
			do jj1=1,PerCellLayer
			    index=ii1*PerCellLayer+jj1
			    x_yin=Cell_Pos(1,index)
                y_yin=Cell_Pos(2,index)
                z_yin=Cell_Pos(3,index)
                radd=dsqrt((x_local-x_yin)**2+(y_local-y_yin)**2+(z_local-z_yin)**2)
				if((radd<detmin) .and. (index .ne. iBOU(1))) then
				  detmin=radd
				  iBOU(2)=index
				end if
			end do
			!2nd neighbor
			detmin=100.
			do jj1=1,PerCellLayer
			    index=ii1*PerCellLayer+jj1
			    x_yin=Cell_Pos(1,index)
                y_yin=Cell_Pos(2,index)
                z_yin=Cell_Pos(3,index)
                radd=dsqrt((x_local-x_yin)**2+(y_local-y_yin)**2+(z_local-z_yin)**2)
				if((radd<detmin) .and. (index .ne. iBOU(1)) .and. (index .ne. iBOU(2))) then
				  detmin=radd
				  iBOU(3)=index
				end if
			end do
			!3rd neighbor
			detmin=100.
			do jj1=1,PerCellLayer
			    index=ii1*PerCellLayer+jj1
			    x_yin=Cell_Pos(1,index)
                y_yin=Cell_Pos(2,index)
                z_yin=Cell_Pos(3,index)
                radd=dsqrt((x_local-x_yin)**2+(y_local-y_yin)**2+(z_local-z_yin)**2)
				if((radd<detmin) .and. (index .ne. iBOU(1)) .and. (index .ne. iBOU(2)).and. &
				  (index .ne. iBOU(3))) then
				  detmin=radd
				  iBOU(4)=index
				end if
			end do
			!4th neighbor
			detmin=100.
			do jj1=1,PerCellLayer
			    index=ii1*PerCellLayer+jj1
			    x_yin=Cell_Pos(1,index)
                y_yin=Cell_Pos(2,index)
                z_yin=Cell_Pos(3,index)
                radd=dsqrt((pos_local(1)-x_yin)**2+(pos_local(2)-y_yin)**2+(pos_local(3)-z_yin)**2)
				if((radd<detmin) .and. (index .ne. iBOU(1)) .and. (index .ne. iBOU(2)).and. &
				  (index .ne. iBOU(3)) .and. (index .ne. iBOU(4))) then
				  detmin=radd
				  iBOU(5)=index
				end if
			end do	
			layerindex(1:5,i1,j1)=iBOU(1:5)
		  end if
		    cen_nei(1:3,1:5)=Cell_Pos(1:3,iBOU(1:5))
			Vari_output(1:8,1:5)=Cell_state(1:8,iBOU(1:5))
		  
		    x_local=cen_nei(1,1)
			y_local=cen_nei(2,1)
			z_local=cen_nei(3,1)
			x_inter=pos_local(1)
			y_inter=pos_local(2)
			z_inter=pos_local(3)
            temp=dsqrt(x_local**2.0+y_local**2.0+z_local**2.0)
            if((temp<1.d0) .or. (i .le. 1) .or. (i .ge. R_IndMax-1))then
			  !write(*,*)'ii1=',ii1,'index=',iBOU(1)
			  call R_inverse_Avg(x_inter,y_inter,z_inter,cen_nei(1:3,1:5),Vari_output(1:8,1:5),Vari_interp(1:8))
			else
			  !5th neighbor
			  iBOU(6)=iBOU(1)+PerCellLayer
			  cen_nei(1:3,6)=Cell_Pos(1:3,iBOU(6))
			  Vari_output(1:8,6)=Cell_state(1:8,iBOU(6))
			  !6th neighbor
			  iBOU(7)=iBOU(1)-PerCellLayer
			  cen_nei(1:3,7)=Cell_Pos(1:3,iBOU(7))
			  Vari_output(1:8,7)=Cell_state(1:8,iBOU(7))
			  
               davg=0.d0
			   do ivari=2,7
			   davg=davg+dsqrt((cen_nei(1,1)-cen_nei(1,ivari))**2.d0+&
			   (cen_nei(2,1)-cen_nei(2,ivari))**2.d0+(cen_nei(3,1)-cen_nei(3,ivari))**2.d0)
			   end do
			   davg=davg/6.d0	
               if(davg<eps)	write(*,*)'davg=',davg,'ii,jj=',ii,jj		   
               a_divB_rbf81=0.d0
               do ivari=1,7  
                   do jj=1,7                	                   
                    basefun=c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-&
					cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0
                                              
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+1)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))                                                    
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+2)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                    a_divB_rbf81((ivari-1)*3+1,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					  (cen_nei(1,ivari)-cen_nei(1,jj))/davg**4.d0)*basefun**(-1.5d0)                 
                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)                
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+2)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))                   
                    a_divB_rbf81((ivari-1)*3+2,(jj-1)*3+3)=-((cen_nei(3,ivari)-cen_nei(3,jj))*&
					  (cen_nei(2,ivari)-cen_nei(2,jj))/davg**4.d0)*basefun**(-1.5d0)
                                     
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+1)=-((cen_nei(1,ivari)-cen_nei(1,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+2)=-((cen_nei(2,ivari)-cen_nei(2,jj))*&
					  (cen_nei(3,ivari)-cen_nei(3,jj))/davg**4.d0)*basefun**(-1.5d0)                                   
                    a_divB_rbf81((ivari-1)*3+3,(jj-1)*3+3)=(-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
					  (((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+(cen_nei(2,ivari)-cen_nei(2,jj))**2.d0)/&
					  davg**4.d0)*basefun**(-1.5d0))	                                
                   end do
               end do   
               m=21 
               n=21 
               call SVD_solver(m,n,a_divB_rbf81,bbrbf_,wrbf,vrbf)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf(jj).gt.wmax)wmax=wrbf(jj)
               end do 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf(jj).lt.wmin)wrbf(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               end do                              
               aB=bbrbf_
               wB=wrbf
               vB=vrbf  
               
               a_divB_rbf7=0.d0
			   do ivari=1,7                 			  
                   do jj=1,7              
                       a_divB_rbf7(ivari,jj)=dsqrt(c_rbf/davg+((cen_nei(1,ivari)-cen_nei(1,jj))**2.d0+&
					     (cen_nei(2,ivari)-cen_nei(2,jj))**2.d0+(cen_nei(3,ivari)-cen_nei(3,jj))**2.d0)/davg**2.d0)             			   
                   end do
               end do
               m=7 
               n=7 
               call SVD_solver(m,n,a_divB_rbf7,bbrbf7_,wrbf7,vrbf7)
               !---------------------------------------    
               wmax=0.d0     !Will be the maximum singular value obtained.
               do jj=1,n
                   if(wrbf7(jj).gt.wmax)wmax=wrbf7(jj)
               enddo 
               wmin=wmax*1.d-10                !This is where we set the threshold for singular values
               do jj=1,n                         !allowed to be nonzero. The constant is typical,
                   if(wrbf7(jj).lt.wmin)wrbf7(jj)=0.d0     !but not universal. You have to experiment with your own application.     
               enddo                  
               a=bbrbf7_
               w=wrbf7
               v=vrbf7        
                              
               do jj=1,7 
               agspt_weight7_point(jj)=dsqrt(c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			     (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0) 
               
               basefun=c_rbf/davg+((x_inter-cen_nei(1,jj))**2.d0+&
			     (y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**2.d0
               
               agspt_weightB_point(1,(jj-1)*3+1)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((y_inter-cen_nei(2,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(1,(jj-1)*3+2)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)         
               agspt_weightB_point(1,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(x_inter-cen_nei(1,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                            
               agspt_weightB_point(2,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(2,(jj-1)*3+2)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(z_inter-cen_nei(3,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)        
               agspt_weightB_point(2,(jj-1)*3+3)=-((z_inter-cen_nei(3,jj))*(y_inter-cen_nei(2,jj))/&
			     davg**4.d0)*basefun**(-1.5d0) 
                              
               agspt_weightB_point(3,(jj-1)*3+1)=-((x_inter-cen_nei(1,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)
               agspt_weightB_point(3,(jj-1)*3+2)=-((y_inter-cen_nei(2,jj))*(z_inter-cen_nei(3,jj))/&
			     davg**4.d0)*basefun**(-1.5d0)      
               agspt_weightB_point(3,(jj-1)*3+3)=-(2.d0/davg**2.d0)*basefun**(-0.5d0)+&
			     (((x_inter-cen_nei(1,jj))**2.d0+(y_inter-cen_nei(2,jj))**2.d0)/davg**4.d0)*basefun**(-1.5d0)                              
              end do  
          
              call initial_reconstruction_B_cell(Vari_output(5:7,1:7),aB,wB,vB,agspt_weightB_point(1:3,1:21),Vari_interp(5:7))
	          do ivari=5,7
	              call rbf_limiter(Vari_output(ivari,1:7),Vari_interp(ivari))
	          end do
              do ivari=1,4  !pri_variables   
                  call recon_cellrbf(Vari_output(ivari,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari)) 
                  call rbf_limiter(Vari_output(ivari,1:7),Vari_interp(ivari))		  
              end do
              ivari=8
              call recon_cellrbf(Vari_output(ivari,1:7),a,w,v,agspt_weight7_point(1:7),Vari_interp(ivari))  
              call rbf_limiter(Vari_output(ivari,1:7),Vari_interp(ivari))				   
			end if
			  
            usphinterp(1,i1,j1)=x_inter
            usphinterp(2,i1,j1)=y_inter
            usphinterp(3,i1,j1)=z_inter
            usphinterp(4:11,i1,j1)=Vari_interp(1:8)
         end do
      end do
    end subroutine UniformStructuredFormat	

    subroutine hours2yymmddhh(N_CARRINGTON_NO,i_cur,hours,days,months,year)	
	integer,intent(in):: N_CARRINGTON_NO,i_cur
	integer,intent(inout):: hours,days,months,year
	integer :: hour,day,month,day_M,days_Feb
	integer :: current_data,N_CARRINGTON_NOs,i
	current_data=N_CARRINGTON_NO
  do i=1,i_cur
    N_CARRINGTON_NOs=current_data
	hour=mod(N_CARRINGTON_NOs,100)+1
	day=mod(int(N_CARRINGTON_NOs/100),100)
	month=int(N_CARRINGTON_NOs/10000)
	
	hours=mod(hour,24)
	day=day+int(hour/24)
	if(month==1 .or. month==3 .or. month==5 .or. month==7 .or. month==8 .or. month==10 .or. month==12)then
	    day_M=31
	else if(month==2)then
	    day_M=days_Feb
	else
	    day_M=30
	end if
	days=mod(day,day_M)
	month=month+int(day/day_M)
	if(days==0)then
	  days=day_M
	  month=month-1
	end if
	if(month>12)then
	  year=year+int(month/12)
      days_Feb=28
	   if(mod(year,4)==0)then
	     if(mod(year,100)==0)then
		   if(mod(year,400)==0)then
		    days_Feb=29
		   end if
		 else
		   days_Feb=29
		 end if
	   end if   
	end if
	months=mod(month,12)
	if(months==0)months=12
	N_CARRINGTON_NOs=months*10000+days*100+hours
	current_data=N_CARRINGTON_NOs
  end do
	end subroutine hours2yymmddhh
    
end Program output


 





