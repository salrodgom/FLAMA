module mod_random
! module for pseudo random numbers
 implicit none
 private
 public init_random_seed, randint, r4_uniform
 contains
!
 subroutine init_random_seed( )
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer              :: i, n, un, istat, dt(8), pid
  integer(int64)       :: t
  call random_seed(size = n)
  allocate(seed(n))
  write(6,'(a)')"Random seed:"
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
   read(un) seed
   close(un)
   write(6,'(a)')"OS provides a random number generator"
  else
   ! Fallback to XOR:ing the current time and pid. The PID is
   ! useful in case one launches multiple instances of the same
   ! program in parallel.
   call system_clock(t)
   if (t == 0) then
    call date_and_time(values=dt)
    t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
      + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
      + dt(3) * 24_int64 * 60 * 60 * 1000 &
      + dt(5) * 60 * 60 * 1000 &
      + dt(6) * 60 * 1000 + dt(7) * 1000 &
      + dt(8)
   end if
   pid = getpid()
   t = ieor(t, int(pid, kind(t)))
   do i = 1, n
      seed(i) = lcg(t)
   end do
   write(6,'(a)')"Fallback to the current time and pid."
  end if
  call random_seed(put=seed)
  write(6,*) seed
 contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
   integer :: lcg
   integer(int64) :: s
   if (s == 0) then
    s = 104729
   else
    s = mod(s, 4294967296_int64)
   end if
   s = mod(s * 279470273_int64, 4294967291_int64)
   lcg = int(mod(s, int(huge(0), int64)), kind(0))
  end function lcg
 end subroutine init_random_seed
!
 integer function randint(i,j)
  integer,intent(in)    :: i,j
  real                  :: r
  call random_number(r)
  randint = i + floor((j+1-i)*r)
 end function randint
!
 real function r4_uniform(x,y)
  implicit none
  real,intent(in)       :: x,y
  real                  :: r
  call random_number(r)
  r4_uniform=(r*(y-x))+x
  return
 end function r4_uniform
end module mod_random
!
module types
 implicit none
 integer,parameter   :: max_atom_number = 100
 integer,parameter   :: maxnp = 80, maxcompounds = 1, maxdata = 500 !100
 integer,parameter   :: maxlinelength=maxnp*32
 type   CIFFile
  character(len=100) :: filename
  integer            :: n_atoms
  real               :: rv(1:3,1:3)
  real               :: vr(1:3,1:3)
  real               :: cell_0(1:6)
  real               :: atom_xcrystal(1:3,1:max_atom_number)
  real               :: atom_charge(1:max_atom_number)
  character(len=4)   :: atom_label(1:max_atom_number)
  character(len=2)   :: type_symbol(1:max_atom_number)
  real               :: obs_energy
  real               :: cal_energy
  real               :: obs_energy_weight
  !integer            :: composition        ! 0,1,2,3, etc. (not implemented yet)
 end type
 type                          :: typ_ga
  character(len=maxlinelength) :: genotype
  real                         :: phenotype(1:maxnp)
  real                         :: fitness
 end type
end module types
!
module GeometricProperties
 implicit none
 private
 public cell,uncell,output_gulp,output_lmp,cellnormal2lammps,Clen,Clen_trim ! output_gulp_fit
 contains
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3))
 call inverse(rv,vr,3)
 RETURN
 END SUBROUTINE cell
!
 !SEIJAS   Typical wierd parameters in LAMMPS.
 subroutine cellnormal2lammps(cell_0,xlo_bound,ylo_bound,zlo_bound,&
                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
  implicit none
  real,intent(in)  :: cell_0(6)
  real,intent(out) :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real,intent(out) :: xy,xz,yz
  real,parameter :: pi=acos(-1.0)
  real,parameter :: radtodeg = 180.0/pi
  real,parameter :: degtorad = pi/180.0
  real :: xlo,xhi,ylo,yhi,zlo,zhi,lx,ly,lz,cosa,cosb,cosg
  real :: alp,bet,gam,sing
  IF(cell_0(4) == 90.0) THEN
    cosa = 0.0
  ELSE
    ALP=cell_0(4)*degtorad
    COSA=cos(ALP)
  ENDIF
  IF(cell_0(5) == 90.0) THEN
    cosb = 0.0
  ELSE
    bet = cell_0(5)*degtorad
    cosb = cos(bet)
  ENDIF
  IF(cell_0(6) == 90.0) then
    sing = 1.0
    cosg = 0.0
  ELSE
    gam = cell_0(6)*degtorad
    sing = sin(gam)
    cosg = cos(gam)
  ENDIF
  lx=cell_0(1)
  xy=cell_0(2)*cosg
  xz=cell_0(3)*cosb
  ly=sqrt(cell_0(2)*cell_0(2)-xy*xy)
  yz=(cell_0(2)*cell_0(3)*cosa-xy*xz)/ly
  lz=sqrt(cell_0(3)*cell_0(3)-xz*xz-yz*yz)
  xlo_bound=0.0
  ylo_bound=0.0
  zlo_bound=0.0
  xhi_bound=lx
  yhi_bound=ly
  zhi_bound=lz
  return
 end subroutine cellnormal2lammps
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0)
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    end do
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
 subroutine output_gulp(CIFFiles,GULPFilename)
  use types
  implicit none
  type(CIFfile),intent(inout)     :: CIFFiles
  character(len=100),intent(out)  :: GULPFilename
  integer                         :: u
  integer                         :: i,k
  real                            :: mmm,rrr
  integer                         :: zzz
  character(len=2)                :: zlz
  character(len=4)                :: extension=".gin"
  GULPFilename=CIFFiles%filename(1:Clen_trim(CIFFiles%filename)-4)//extension
  GULPFilename=adjustl(GULPfilename)
  open(newunit=u,file=GULPFilename)
  write(u,'(a)')'single conv molecule'
  write(u,'(A)')'cell'
  write(u,'(6(f9.5,1x))') (CIFFiles%cell_0(i) , i=1,6)
  write(u,'(A)')'fractional'
  do i=1,CIFFiles%n_atoms
   write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')CIFFiles%atom_label(i),&
    (CIFFiles%atom_xcrystal(k,i),k=1,3),CIFFiles%atom_charge(i)
  end do
  write(u,'(a)')'supercell 1 1 1'
  write(u,'(A)')'library peros'
  close(u)
 end subroutine output_gulp
!
 !SEIJAS 
 subroutine output_lmp(CIFFiles,LAMMPSFilename_inp,LAMMPSFilename_dat,&
            xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
  use types
  implicit none
  type(CIFfile),intent(inout)     :: CIFFiles
  character(len=100),intent(out)  :: LAMMPSFilename_inp, LAMMPSFilename_dat 
  character(len=500)  :: line = " "
  character(len=100)              :: LAMMPSFilename_sav
  real,intent(in) :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real,intent(in) :: xy,xz,yz
  integer                         :: i,k
  character(len=4)                :: extension1=".inp"
  character(len=4)                :: extension2=".dat"
  character(len=4)                :: extension3=".sav"
  real                            :: real_coord(1:3)
  LAMMPSFilename_inp=CIFFiles%filename(1:Clen_trim(CIFFiles%filename)-4)//extension1
  LAMMPSFilename_inp=adjustl(LAMMPSfilename_inp)
  LAMMPSFilename_dat=CIFFiles%filename(1:Clen_trim(CIFFiles%filename)-4)//extension2
  LAMMPSFilename_dat=adjustl(LAMMPSfilename_dat)
  LAMMPSFilename_sav=CIFFiles%filename(1:Clen_trim(CIFFiles%filename)-4)//extension3
  LAMMPSFilename_sav=adjustl(LAMMPSfilename_sav)

  !Copy the files.
  write(line,*)"cp template.inp ",LAMMPSFilename_inp
  call system(line)
  write(line,*)"cp template.dat ",LAMMPSFilename_dat
  call system(line)

  !Replace the read_data in the input file.
  !(The sed command gets confused with the / simbol, that is why I had to do it that way)
  write(line,'(a,a,a,a,a,a,a,a)')"sed -i 's/.*read_data.*/read_data    ",LAMMPSFilename_dat(1:5),"\",&
            LAMMPSFilename_dat(6:len(LAMMPSFilename_dat)),"/'i ",LAMMPSFilename_inp(1:5),"\",&
            LAMMPSFilename_inp(6:len(LAMMPSFilename_inp))
  call system(line)

  !Replace cell shape in the data file.
  write(line,*)"sed -i 's/.*xlo xhi.*/",xlo_bound, xhi_bound,&
               " xlo xhi/' ", LAMMPSFilename_dat
  call system(line)
  write(line,*)"sed -i 's/.*ylo yhi.*/",ylo_bound, yhi_bound,&
               " ylo yhi/' ", LAMMPSFilename_dat
  call system(line)
  write(line,*)"sed -i 's/.*zlo zhi.*/",zlo_bound, zhi_bound,&
               " zlo zhi/' ", LAMMPSFilename_dat
  call system(line)
  write(line,*)"sed -i 's/.*xy xz yz.*/",xy, xz, yz,&
               " xy xz yz/' ", LAMMPSFilename_dat
  call system(line)

  !Replace coordinates. It works with 99 atoms or less.
  !Atoms have to be in the same order in both data files (CIF and LAMMPS).
  do i=1,CIFFiles%n_atoms
   !Fractional to cartesian
   do k=1,3
    real_coord(k) = CIFFiles%atom_xcrystal(1,i)*CIFFiles%rv(k,1)&
                   +CIFFiles%atom_xcrystal(2,i)*CIFFiles%rv(k,2)&
                   +CIFFiles%atom_xcrystal(3,i)*CIFFiles%rv(k,3)
   end do
   if (i.le.9) then
    write(line,'(a,i1,a,F12.8,F12.8,F12.8,a,a)')"sed -i 's/Atom_",&
           i,"_xyz/",(real_coord(k),k=1,3),"/g' ", LAMMPSFilename_dat
   else
    write(line,'(a,i2,a,F12.8,F12.8,F12.8,a,a)')"sed -i 's/Atom_",&
           i,"_xyz/",(real_coord(k),k=1,3),"/g' ", LAMMPSFilename_dat
   end if
   call system(line)
  end do

  !Create save files.
  write(line,*)"cp ",LAMMPSFilename_dat," ",LAMMPSFilename_sav
  call system(line)
 end subroutine output_lmp
!
! subroutine output_gulp_fit(n_files,CIFFiles)
!  use types
!  implicit none
!  integer,intent(in)            :: n_files
!  integer                       :: u = 444
!  type(CIFfile),intent(inout)   :: CIFFiles(n_files)
!  character(len=100)            :: GULPFilename="fit.gin"
!  integer                       :: i, j, k
!  real                          :: mmm,rrr,obs_energy_min
!  integer                       :: zzz
!  character(len=2)              :: zlz
!  obs_energy_min=minval(CIFFiles%obs_energy)
!  open(newunit = u, file = GULPFilename)
!  write(u,'(a)')'fit conv molecule'
!  do i=1,n_files
!   write(u,'(A)')'cell'
!   write(u,'(6(f9.5,1x))') (CIFFiles(i)%cell_0(k) , k=1,6)
!   write(u,'(A)')'fractional'
!   do j=1,CIFFiles(i)%n_atoms
!    write(u,'(a4,1x,3(f14.7,1x),1x,f14.7)')CIFFiles(i)%atom_label(j),&
!    (CIFFiles(i)%atom_xcrystal(k,j),k=1,3),CIFFiles(i)%atom_charge(j)
!   end do
!   write(u,'(a)')'supercell 1 1 1'
!   write(u,'(a)')'observable'
!   write(u,*)'energy eV'
!   write(u,*)CIFFiles(i)%obs_energy - obs_energy_min, 100.0
!   write(u,*)'end'
!  end do
!  !# Tell the program to fit the overall shift
!  write(u,'(A)')'vary'
!  write(u,'(A)')' shift'
!  write(u,'(A)')'end'
!  write(u,'(A)')'library peros'
!  close(u)
! end subroutine output_gulp_fit
end module GeometricProperties
!
module get_structures
 implicit none
 private
 public GenerateCIFFileList, ReadListOfCIFFiles, ReadCIFFiles, WriteEnergies
 contains
!
  subroutine GenerateCIFFileList()
   implicit none
   character(len=1000)  :: string=" "
   write(string,'(a,1x,a)')"if [ -f list_debug ] ; then rm list_debug ; touch list_debug ; fi",&
    "; ls struc/*.cif > list_debug"
   call system(string)
  end subroutine GenerateCIFFileList
!
  subroutine ReadListOfCIFFiles(n_files)
   implicit none
   character(len=100)  :: line = " "
   integer             :: ierr = 0
   integer,intent(out) :: n_files
   n_files = 0
   open(111,file="list",iostat=ierr)
   read_list: do
    read(111,'(a)',iostat=ierr) line
    if(ierr/=0) exit read_list
    n_files=n_files+1
   end do read_list
   rewind(111)
   close(111)
   return
  end subroutine
!
  subroutine Recalculate_Energies_and_Weights(n_files, CIFFiles, total_Shift)
   use types
   implicit none
   integer,intent(in)            :: n_files
   type(CIFfile),intent(inout)   :: CIFFiles(n_files)
   integer                       :: i
   real                          :: Z_partition
   real, intent(out)             :: total_Shift
   real,parameter                :: kBT = 0.05170368566 ! (eV)
   total_Shift = minval(CIFFiles(1:n_files)%obs_energy)
   do i=1,n_files
    CIFFiles(i)%obs_energy = CIFFiles(i)%obs_energy - total_Shift
   end do
   Z_partition = sum( exp( - CIFFiles(1:n_files)%obs_energy / kBT ) )
   do i=1,n_files
    CIFFiles(i)%obs_energy_weight = exp( - CIFFiles(i)%obs_energy / kBT ) / Z_partition
   end do
   write(6,'(a,1x,f14.7)')'Weight sum = ',sum( CIFFiles(1:n_files)%obs_energy_weight )
   return
  end subroutine
!
  subroutine ReadCIFFiles(n_files,CIFFiles)
   use types
   use GeometricProperties
   implicit none
   real                              :: energy = 0.0
   integer                           :: i,j,k,l
   integer                           :: ierr = 0
   integer,intent(in)                :: n_files
   type(CIFfile),intent(inout)       :: CIFFiles(n_files)
   character(len=200)               :: line = " "
   character(len=100)                :: filename = " "
   character(len=100)                :: filename_inp = " "
   character(len=100)                :: filename_dat = " "
   character(len=20)                 :: spam
   !string_stop_head depends on how are written the CIF files. It's the last line before the atoms list.
   !character(len=80)                 :: string_stop_head= "_atom_site_fract_z"
   character(len=80)                 :: string_stop_head= "   _atom_site_type_symbol"
   real                              :: infinite = 3.4028e38, total_Shift
   !SEIJAS
   real    :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
   real    :: xy,xz,yz
   integer :: number_of_lines, ios
   real    :: potential_coeff
   open(111,file="list",iostat=ierr)
   if(ierr/=0)stop
   do i=1,n_files
    read(111,'(a)')line
    read(line(1:49),'(a)') CIFFiles(i)%filename
    read(line(50:),*)      CIFFiles(i)%obs_energy
    write(6,'(a)')trim(CIFFiles(i)%filename)
    open(100,file=trim(CIFFiles(i)%filename),status='old',iostat=ierr)
    if(ierr/=0) stop 'CIFFile does not found'
    read_cif: do
     read(100,'(a)',iostat=ierr) line
     if(ierr/=0) exit read_cif
     if(line(1:14)=="_cell_length_a")then
      read(line,*)spam,CIFFiles(i)%cell_0(1)
      cycle read_cif
     end if
     if(line(1:14)=="_cell_length_b")then
      read(line,*)spam,CIFFiles(i)%cell_0(2)
      cycle read_cif
     end if
     if(line(1:14)=="_cell_length_c")then
      read(line,*)spam,CIFFiles(i)%cell_0(3)
      cycle read_cif
     end if
     if(line(1:17)=="_cell_angle_alpha")then
      read(line,*)spam,CIFFiles(i)%cell_0(4)
      cycle read_cif
     end if
     if(line(1:16)=="_cell_angle_beta")then
      read(line,*)spam,CIFFiles(i)%cell_0(5)
      cycle read_cif
     end if
     if(line(1:17)=="_cell_angle_gamma")then
      read(line,*)spam,CIFFiles(i)%cell_0(6)
      cycle read_cif
     end if
     if(line(1:)==string_stop_head) exit read_cif
    end do read_cif
    call cell(CIFFiles(i)%rv,CIFFiles(i)%vr,CIFFiles(i)%cell_0)
    CIFFiles(i)%n_atoms=0
    read_natoms: do
     read(100,'(a)',iostat=ierr) line
     if(ierr/=0) exit read_natoms
     CIFFiles(i)%n_atoms=CIFFiles(i)%n_atoms+1
    end do read_natoms
    rewind(100)
    write(6,'(80a)')('=',k=1,80)
    write(6,*)'Observable, energy:',CIFFiles(i)%obs_energy
    write(6,*)'Atoms:', CIFFiles(i)%n_atoms, CIFFiles(i)%filename
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(1,j), j=1,3 )
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(2,j), j=1,3 )
    write(6,'(3(f14.7,1x))') ( CIFFiles(i)%rv(3,j), j=1,3 )
    do
     read(100,'(a)') line
     if (line(1:)==string_stop_head) exit
    end do
    do j=1,CIFFiles(i)%n_atoms
     read(100,'(a)') line
     !!Salva CIF files. SEIJAS
     !read(line,*) CIFFiles(i)%atom_label(j),CIFFiles(i)%type_symbol(j),&
     !            (CIFFiles(i)%atom_xcrystal(k,j),k=1,3)
     !CIFFiles(i)%atom_charge(j)=0.0
     !Bipasa CIF files.
     read(line,*) CIFFiles(i)%atom_label(j),CIFFiles(i)%atom_charge(j),&
                 (CIFFiles(i)%atom_xcrystal(k,j),k=1,3)
     CIFFiles(i)%atom_charge(j)=0.0
     CIFFiles(i)%type_symbol(j)=''
     !------
     write(6,'(a2,1x,3(f14.7,1x))')CIFFiles(i)%type_symbol(j),(CIFFiles(i)%atom_xcrystal(k,j),k=1,3)
    end do
    close(100)
    !SEIJAS  Call LAMMPS.
!R!    LAMMPS_or_GULP: if ( flag_lmp ) then
     call cellnormal2lammps(CIFFiles(i)%cell_0,xlo_bound,ylo_bound,zlo_bound,&
                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)    
     call output_lmp(CIFFiles(i),filename_inp(1:Clen_trim(filename_inp)),&
          filename_dat(1:Clen_trim(filename_dat)),&
          xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
     !Write in the data file the values of the pair coefs. It can not be done in the Output_lmp routine
     !because of incompatibilities with the use of this routine in the openmp calculations.
     !Replace actual values of PairIJ Coeffs. 
     !Number of lines in the phenotype.dat file.
     OPEN (UNIT=19, FILE="phenotype.dat", STATUS='OLD', ACTION='READ')
     number_of_lines = 0
     DO
      READ (19,*, END=10)
      number_of_lines = number_of_lines + 1
     END DO
     10 CLOSE (19)

     OPEN (UNIT=19, FILE="phenotype.dat", STATUS='OLD', ACTION='READ')
     do l = 1,number_of_lines
      read(19,*) potential_coeff 
      if (l.le.9) then
       write(line,'(a,i1,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",potential_coeff,"/g' ",filename_dat
      else
       write(line,'(a,i2,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",potential_coeff,"/g' ",filename_dat
      end if
      call system(line)
     end do
     CLOSE (19)
     !Write pair coefs until here.
     write(6,*)'LAMMPS file:',filename_dat(1:Clen_trim(filename_dat))
     line="lmp -in "//filename_inp(1:Clen_trim(filename_inp))//" > tmp "
     call system(line)
     line="grep -A1 'TotEng' tmp | awk '{print $1}' > c"
     call system(line)
     open(456,file="c")
     read(456,*)
     read(456,'(a)')line
     if(line(1:4)=="-nan")then
      CIFFiles(i)%cal_energy=infinite
     else
      read(line,*) CIFFiles(i)%cal_energy
     end if
     close(456)
     line="rm c tmp"
     call system(line)
!     line="rm ",filename_inp(1:Clen_trim(filename_inp))
!     call system(line)
!     line="rm ",filename_dat(1:Clen_trim(filename_dat))
!     call system(line)
!R!    else
!!!--------------------------------
!R!     call output_gulp(CIFFiles(i),filename(1:Clen_trim(filename)))
!R!     write(6,*)'GULP file:',filename(1:Clen_trim(filename))
!R!     line="~/bin/gulp < "//filename(1:Clen_trim(filename))//" > tmp "
!R!     call system(line)
!R!     line="grep 'Total lattice energy       =' tmp | grep 'eV' | awk '{print $5}' > c"
!R!     call system(line)
!R!     open(456,file="c")
!R!     read(456,'(a)')line
!R!     if(line(1:20)=="********************")then
!R!      CIFFiles(i)%cal_energy=infinite
!R!     else
!R!      read(line,*) CIFFiles(i)%cal_energy
!R!     end if
!R!     close(456)
!R!     line="rm c tmp "//filename(1:Clen_trim(filename))//" "
!R!     call system(line)
!R!    end if LAMMPS_or_GULP
!!!---------------------------------
    write(6,*)'Calculated, energy', CIFFiles(i)%cal_energy
   end do
   close(111)
   call Recalculate_Energies_and_Weights(n_files, CIFFiles, total_Shift)
   write(6,'(a)') 'Situation:'
   write(6,'(a)') '================================================================='
   write(6,'(a,1x,f14.7)') 'Shift:', total_Shift
   write(6,'(4(a,4x))') 'Obs Energy/ eV','Calc. Energy / eV',' weight / -',' filename'
   do i=1,n_files
    write(6,*)CIFFiles(i)%obs_energy, CIFFiles(i)%cal_energy, &
              CIFFiles(i)%obs_energy_weight, trim(CIFFiles(i)%filename)
   end do
   return
  end subroutine ReadCIFFiles
!
  subroutine WriteEnergies(n_files,CIFFiles,add)
   use types
   implicit none
   integer, intent(in)         :: n_files
   character(len=3),intent(in) :: add
   type(CIFfile),intent(inout) :: CIFFiles(n_files)
   integer                     :: i,u=123,ierr=0
   character(len=20)           :: filename = " "
   filename="fitness_"//add(1:3)//".txt"
   open(u,file=filename,iostat=ierr)
   if(ierr/=0) stop "fitness.txt can not be open"
   write(u,'(a)')"# struc.;  CellSize / A ;  Energy Obs. / eV; Energy Cal. / eV"
   do i=1,n_files
    write(u,*)i,CIFFiles(i)%cell_0(1),CIFFiles(i)%obs_energy,CIFFiles(i)%cal_energy
   end do
   close(u)
  end subroutine WriteEnergies
!
end module get_structures
!
module qsort_c_module
! Recursive Fortran 95 quicksort routine sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms, 1997 printing
! Made F conformant by Walt Brainerd
 implicit none
 public  :: QsortC
 private :: Partition
 contains
!
 recursive subroutine QsortC(A)
  real(16), intent(in out), dimension(:) :: A
  integer                                :: iq
  if(size(A) > 1) then
   call Partition(A, iq)
   call QsortC(A(:iq-1))
   call QsortC(A(iq:))
  endif
 end subroutine QsortC
!
 subroutine Partition(A, marker)
  real(16), intent(in out), dimension(:) :: A
  integer, intent(out)               :: marker
  integer                            :: i, j
  real(16)                           :: temp
  real(16)                           :: x
  x = A(1)
  i = 0
  j = size(A) + 1
  do
   j = j-1
   do
    if (A(j) <= x) exit
    j = j-1
   end do
   i=i+1
   do
    if (A(i) >= x) exit
    i=i+1
   end do
   if (i < j) then
    temp = A(i)
    A(i) = A(j)
    A(j) = temp
   elseif (i == j) then
    marker = i+1
    return
   else
    marker = i
    return
   endif
  end do
 end subroutine Partition
!
end module qsort_c_module
!
module flama_globals
 use types
 use mod_random
 use get_structures
 implicit none
 integer                    :: npar,i
 integer,allocatable        :: np(:)
 integer                    :: err_apertura,ii,intervalos,j,n_refits
 real,allocatable           :: param(:,:)
 real,target                :: datas(2,maxcompounds,maxdata),f(maxdata)
 real,pointer               :: x(:),y(:),alldat(:,:)
 character(15),allocatable  :: ajuste(:,:)
 character(32*maxnp),allocatable :: string_IEEE(:)
 character(100)             :: line,string
 character(5)               :: inpt
 logical                    :: flag = .true., FlagFire = .false.,seed_flag=.true.
 logical                    :: physical_constrains = .false., range_flag =.true.
 logical                    :: refit_flag = .false., flag_shift=.false.
!R! logical                    :: flag_lmp = .false. !SEIJAS   To use LAMMPS instead of GULP.
 real,parameter             :: R = 0.008314472 ! kJ / mol / K
 real                       :: T = 298.0
 contains
  subroutine read_input()
  implicit none
  character(len=32)  :: chain
  allocate(np(1))
  npar = 0
  read_input_do: do
   read(5,'(A)',iostat=err_apertura)line
   if ( err_apertura /= 0 ) exit read_input_do
   if(line(1:1)=='#') cycle read_input_do
   if(line(1:7)=="shifted") then
    write(6,*)'[WARN] Shifting ( n_par -> n_par + 1 )'
    read(line,*)inpt,flag_shift
    if( flag_shift ) then
     npar=npar+1
    end if
    cycle read_input_do
   end if
   if(line(1:5)=='n_par')then
    read(line,*)inpt,np(1)
    npar=npar+np(1)
    allocate(ajuste(1:2,1:npar))
    allocate(param(1,0:npar-1))
    np(1) = npar
    param(1,:) = 0.0
    write(6,'(a,1x,i2)')'Parameters to fit:', npar
    do i=1,npar
     read(5,'(a)') line
     read(line( 1:10),'(a)') ajuste(1,i)
     read(line(12:),'(a)')   ajuste(2,i)
     write(6,'(a10,1x,a10)') ajuste(1,i),ajuste(2,i)
    end do
    cycle read_input_do
   end if
   !SEIJAS   To call LAMMPS instead of GULP.
!R!   if(line(1:6)=='LAMMPS')then
!R!    read(line,*)inpt,flag_lmp
!R!    cycle read_input_do
!R!   end if
   if(line(1:5)=='ffit?')then
    read(line,*)inpt,flag
    cycle read_input_do
   end if
   if(line(1:5)=='refit') then
    write(6,*)'[WARN] Refitting parameters'
    read(line,*)inpt,refit_flag
    if(refit_flag)then
     read(5,*) n_refits
     allocate(string_IEEE(n_refits))
     do i=1,n_refits
      do j=1,npar
       read(5,'(a32)') chain(1:32)
       write(6,'(a32)')chain(1:32)
       string_IEEE(i)(32*(j-1)+1:32*j)=chain(1:32)
      end do
      write(6,'(a)')string_IEEE(i)(1:32*npar)
     end do
    else
     allocate(string_IEEE(1))
     string_IEEE(i) = ' '
    end if
   end if
   if(line(1:22)=='physically_constrained') then
    physical_constrains=.true.
    write(6,'(a)') '[WARN] The fits are physically constrained'
   end if
   if(err_apertura/=0) exit read_input_do
  end do read_input_do
 end subroutine read_input
!
 subroutine ReadObservables(n_files,CIFFiles)
  implicit none
  integer                        :: i
  integer,intent(in)             :: n_files
  type(CIFfile),intent(inout)    :: CIFFiles(n_files)
  do i=1,n_files
   datas(1,1,i)=real(i)
   datas(2,1,i)=CIFFiles(i)%obs_energy
   write(6,'(2f20.5)') datas(1,1,i),datas(2,1,i)
  end do
 end subroutine ReadObservables
end module flama_globals
!
module mod_genetic
 use types
 use mod_random
 use flama_globals
 use qsort_c_module
 use get_structures
 use GeometricProperties
 implicit none
 public fit
 integer,parameter             :: ga_size         = 800 ! number of GA agents
 real,parameter                :: ga_mutationrate = 0.6 !0.3333 !2000/real(ga_size)
 real,parameter                :: ga_eliterate= 0.01, GA_DisasterRate = 0.0000001
 integer,parameter             :: ga_elitists = int( ga_size * ga_eliterate)
 type(typ_ga), pointer         :: parents(:)
 type(typ_ga), pointer         :: children(:)
 type(typ_ga), target          :: pop_alpha( ga_size )
 type(typ_ga), target          :: pop_beta( ga_size )
 contains
!
 !type(type_ga) function NewCitizenNoConstrain(compound,n_files_CIFFiles)
 ! implicit None
 ! integer                     :: i,j,k
 ! integer,intent(in)          :: n_files,compound
 ! type(CIFFile),intent(inout) :: CIFFiles(n_files)
 ! real                        :: infinite = 3.4028e38
 ! character(len=15)           :: funk = " "
 ! ! Initialise variables
 ! NewCitizenNoConstrain%fitness = infinite
 ! NewCitizenNoConstrain%genotype = ' '
 ! ! Initialise in binary and transform to real:
 ! do i = 1,32*np(compound)
 !  ! random array of 0 and 1
 !  NewCitizenNoConstrain%genotype(i:i) = achar(randint(48,49))
 ! end do
 ! do i = 1,np(compound)
 !  read(NewCitizenNoConstrain%genotype(32*(i-1)+1:32*i),'(b32.32)') &
 !   NewCitizenNoConstrain%phenotype(i)
 ! end do
 ! NewCitizenNoConstrain%fitness = &
 !   fitness( NewCitizenNoConstrain%phenotype,compound,n_files,CIFFiles)
 ! return
 !end function NewCitizenNoConstrain
!
 type(typ_ga) function NewCitizen(compound,n_files,CIFFiles)
  implicit none
  integer                     :: i,j,k
  integer,intent(in)          :: n_files,compound
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  real                        :: infinite = 3.4028e38
  character(len=15)           :: funk = " "
  ! Initialise variables
  NewCitizen%fitness = infinite
  NewCitizen%genotype = ' '
  ! Initialise in real, in a range, and transform to binary:
  make_new_values_under_constrains: do i = 1,np(compound)
   funk = " "
   phys_constrains_make_values: if ( physical_constrains ) then
    funk=ajuste(2,i)(1:Clen_trim(ajuste(2,i)))
    funk=adjustl(funk)
    select case(funk)
     case("A_buck")
      NewCitizen%phenotype(i) = r4_uniform( 0.1 , 5e+6 ) 
     case("A_lj")
      NewCitizen%phenotype(i) = r4_uniform( 1e-5, 1e10 )
     case("epsilon")
      NewCitizen%phenotype(i) = r4_uniform( 0.0002, 0.01 )
     case("sigma")
      NewCitizen%phenotype(i) = r4_uniform( 1.0, 4.0 )
     case("rho_buck")
      NewCitizen%phenotype(i) = r4_uniform( 1e-1, 0.6 )
     case("C_buck")
      NewCitizen%phenotype(i) = r4_uniform( 0.0, 5e4 )
     !SEIJAS   To fit also the Ks of the organic molecule.
     case("K_bond")
      NewCitizen%phenotype(i) = r4_uniform( 1.0, 300.0 )
     case("K_angle")
      NewCitizen%phenotype(i) = r4_uniform( 1.0, 4.0 )
     case("K_dihed")
      NewCitizen%phenotype(i) = r4_uniform( 0.01, 5.0 )
    end select
   end if phys_constrains_make_values
   NewCitizen%genotype(32*(i-1)+1:32*i)=real2bin( NewCitizen%phenotype(i))
  end do make_new_values_under_constrains
  NewCitizen%fitness = fitness( NewCitizen%phenotype,compound,n_files,CIFFiles)
  return
 end function NewCitizen
!
 subroutine UpdateCitizen( axolotl ,compound , n_files, CIFFiles)
  implicit none
  integer,intent(in)          :: compound,n_files
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  integer                     :: i,GA_ELITISTS
  real                        :: infinite = 0.0
  type(typ_ga), intent(inout) :: axolotl
  do i = 1,np(compound)
   read(axolotl%genotype(32*(i-1)+1:32*i),'(b32.32)') axolotl%phenotype(i)
  end do
  axolotl%fitness = fitness( axolotl%phenotype,compound,n_files,CIFFiles)
  return
 end subroutine UpdateCitizen
!
 integer*4 function get_file_unit (lu_max)
! get_file_unit returns a unit number that is not in use
  integer*4  :: lu_max, lu, m, iostat
  logical    :: opened
  m = lu_max  ;  if (m < 1) m = 97
  do lu = m,1,-1
   inquire (unit=lu, opened=opened, iostat=iostat)
   if (iostat.ne.0) cycle
   if (.not.opened) exit
  end do
  get_file_unit = lu
  return
end function get_file_unit
!
 subroutine WriteLib(compound,phenotype)
  implicit none
  real, intent(in)    :: phenotype(maxnp)
  integer,intent(in)  :: compound
  character(len=200)  :: line
  integer             :: i
  !SEIJAS changed to include phenotype.dat
! YA NO VALE AL CAMBIAR LA ESCRITURA DE LOS PAIR COEFFS
!  call system("rm phenotype.dat")
!  OPEN (19,file="phenotype.dat")
  call system("cp peros_input.lib peros.lib")
  referw: do i = 1,np(compound)
   write(line,'(a,a,a,e20.15,a)')"sed -i 's/",&
    trim(ajuste(1,i)(1:Clen_trim(ajuste(1,i)))),"/",phenotype(i),"/g' peros.lib"
   call system(line)
!   write(19,*) phenotype(i)
  end do referw
  CLOSE(19, STATUS='KEEP')
 end subroutine WriteLib
!
 real function Fitness(phenotype,compound,n_files,CIFFiles)
!bueno  use OMP_LIB
  implicit none
  real, intent(in)    :: phenotype(maxnp)
  integer,intent(in)  :: n_files
  type(CIFfile),intent(inout) :: CIFFIles(n_files)
  integer             :: i, j, jjj,compound,k = 0, u,ii,np_real,l
  real                :: a(0:np(compound)-1),xx,yy,penalty
  real                :: obs_energy_min, cal_energy_min, obs_energy_max
  real                :: partition = 0.0
  character(len=100)  :: funk = " ",filename(n_files) !,filename_inp(n_files),filename_dat(n_files)
  integer,dimension(8) :: values           !! SEIJAS BASURA holaa
  integer              :: ierr, number_of_lines  !se puede quitar cuando se arregle tambien
  character(len=200)  :: line
  character(len=1000) :: script
  logical             :: flagzero = .false.
  real                :: infinite = 3.4028e38
  real :: xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound
  real :: xy,xz,yz
  call system("cp peros_input.lib peros.lib")
  penalty=0.0
  fitness = 0.0
  ! For each parameter, we add a penalty if the parameter is out of range
!XR  call system("rm phenotype.dat")
!XR  OPEN (19,file="phenotype.dat") !SEIJAS
  refer: do i = 1,np(compound)
   funk = " "
   phys_constrains: if ( physical_constrains ) then
    funk=ajuste(2,i)(1:Clen_trim(ajuste(2,i)))
    funk=adjustl(funk)
    select case(funk)
     case("A_buck")
      if (phenotype(i)<=0.1.or.isnan(phenotype(i)).or.phenotype(i)>=5e6)then
       penalty = infinite
       exit refer
      end if
     case("A_lj")
      if (phenotype(i)<=1e-7.or.isnan(phenotype(i)).or.phenotype(i)>=1e10)then
       penalty = infinite
       exit refer
      end if
     case("epsilon")
      if (phenotype(i)<0.0002.or.isnan(phenotype(i)).or.phenotype(i)>0.01)then
       penalty = infinite
       exit refer
      end if
     case("sigma")
      if (phenotype(i)<1.0.or.isnan(phenotype(i)).or.phenotype(i)>4.0)then
       penalty = infinite
       exit refer
      end if
     case("rho_buck")
      if (phenotype(i)<1.0e-1.or.phenotype(i)>0.6.or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
     case("C_buck")
      if(phenotype(i)<0.0.or.phenotype(i)>1e7.or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
     case("E_shift")
      if(phenotype(i)>=maxval(CIFFiles%obs_energy).or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
     !SEIJAS   To fit also the Ks of the organic molecule.
     case("K_bond")
      if(phenotype(i)<1.0.or.phenotype(i)>300.0.or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
     case("K_angle")
      if(phenotype(i)<1.0.or.phenotype(i)>4.0.or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
     case("K_dihed")
      if(phenotype(i)<0.01.or.phenotype(i)>5.0.or.isnan(phenotype(i)))then
       penalty = infinite
       exit refer
      end if
    end select
   end if phys_constrains
   !SEIJAS    not needed for LAMMPS
   !write(line,'(a,a,a,e20.15,a)')&
   ! "sed -i 's/",trim(ajuste(1,i)(1:Clen_trim(ajuste(1,i)))),"/",phenotype(i),"/g' peros.lib"
   !call system(line)
  end do refer
  CLOSE(19)
  ! We proced with the GULP interface, in order to calculate the fitness:
  ! Interface with GULP code
  calgulp: if ( penalty < 1.0 ) then
!..!$omp parallel default(private) shared(n_files, CIFFiles, filename, line, filename_inp, filename_dat, phenotype)
!..!$omp do
!..   scan_: do i=1,n_files
!NN    do i=1,n_files
    ! get a free unit for read/write for each CPU
    !u = get_file_unit(444)
    ! the functionality newunit is for Fortran 2008:
    !!!!$omp critical
    !SEIJAS  Call LAMMPS.
!!!---------------------------------
!R!    LMP_or_GULP: if ( flag_lmp ) then
!NN     call cellnormal2lammps(CIFFiles(i)%cell_0,xlo_bound,ylo_bound,zlo_bound,&
!NN                              xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
!NN     call output_lmp(CIFFiles(i),filename_inp(i)(1:Clen_trim(filename_inp(i))),&
!NN          filename_dat(i)(1:Clen_trim(filename_dat(i))),&
!NN          xlo_bound,ylo_bound,zlo_bound,xhi_bound,yhi_bound,zhi_bound,xy,xz,yz)
!NN     !Use actual phenotype.
!NN     do l = 1,np(compound)
!NN      if (l.le.9) then
!NN       write(line,'(a,i1,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",phenotype(l),"/g' ",filename_dat(i)
!NN      else
!NN       write(line,'(a,i2,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",phenotype(l),"/g' ",filename_dat(i)
!NN      end if
!NN      call system(line)
!NN     end do
!NN    end do
!call date_and_time(VALUES=values)
!write(*,*) values
    do i=1,n_files
     line="cp "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".sav "//&
          CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".dat"
     call system(line)
     !
     do l = 1,np(compound)
      if (l.le.9) then
       write(line,'(a,i1,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",phenotype(l),"/g' ",&
                  CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".dat"
      else
       write(line,'(a,i2,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",phenotype(l),"/g' ",&
                  CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".dat"
      end if
      call system(line)
     end do
    end do
!call date_and_time(VALUES=values)
!write(*,*) values
     !
!$omp parallel default(private) shared(n_files, CIFFiles, filename) !!, filename_inp, filename_dat)
!$omp do
    do i=1,n_files
     line="lmp -in "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".inp"//" > "&
           //CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".out"
     call system(line)
     line="grep -A1 'TotEng' "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".out | awk '{print $1}' > "&
          //CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".tmp "
     call system(line)
!
!
! Si esto no funciona borrar esto de leer los archivos dentro del paralelo y descomentar los temp.
     ierr=0
     number_of_lines = 0
     open(newunit = u,file=CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".tmp",iostat=ierr)
     DO
      READ (u,*,iostat=ierr)
      if(ierr/=0) exit 
      number_of_lines = number_of_lines + 1
     END DO
     if (number_of_lines == 2) then
      rewind(u)
      read(u,*)
      read(u,*)line
      if(line(1:4)=="-nan")then
       CIFFiles(i)%cal_energy=infinite
       write(*,*) ':::::::::   WARNING: -nan obtained in LAMMPS   '
      else
       read(line,*) CIFFiles(i)%cal_energy
      end if
     else
      !Creo que es una chapuza pero no se me ocurre como hacer que esto no de mas fallos.
      write(*,*) ':::::::::   WARNING: Something wrong with the files ',&
                 CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)       
      close(u)
      write(*,*) ':::::::::   Trying again... '
      line="lmp -in "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".inp"//" > "&
            //CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".out"
      call system(line)
      line="grep -A1 'TotEng' "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".out | awk '{print $1}' > "&
           //CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".tmp "
      call system(line)
      ierr=0
      number_of_lines = 0
      open(newunit = u,file=CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".tmp",iostat=ierr)
      DO
       READ (u,*,iostat=ierr)
       if(ierr/=0) exit
       number_of_lines = number_of_lines + 1
      END DO
      if (number_of_lines == 2) then
       rewind(u)
       read(u,*)
       read(u,*)line
       if(line(1:4)=="-nan")then
        CIFFiles(i)%cal_energy=infinite
        write(*,*) ':::::::::   WARNING: -nan obtained in LAMMPS   '
       else
        write(*,*) ':::::::::   FIXED!!'
        read(line,*) CIFFiles(i)%cal_energy
       end if
      else
       write(*,*) ':::::::::   WARNING: Something wrong AGAIN '  
       line="cp "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".dat fallo.dat"
       call system(line)
       line="cp "//CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".inp fallo.inp"
       call system(line)
       CIFFiles(i)%cal_energy=infinite
      end if
     end if
     close(u)
!
!
!
    end do
!$omp end do
!$omp end parallel
!:......creo que este comando jode el paralelo. el otro ya tiene un synchronization implicito. !$omp barrier
!Nr    do i=1,n_files
!Nr     open(22,file=CIFFiles(i)%filename(1:Clen_trim(CIFFiles(i)%filename)-4)//".tmp")
!     open(newunit = u,file=filename_inp(i)(1:Clen_trim(filename_inp(i))-4)//".tmp")
!Nr     read(22,*)
!     read(u,*)
!Nr     read(22,'(a)')line
!     read(u,'(a)')line
!Nr     close(22)
  !Nr   if(line(1:4)=="-nan")then
!Nr      CIFFiles(i)%cal_energy=infinite
!Nr     else
!Nr      read(line,*) CIFFiles(i)%cal_energy
!Nr     end if
!R!    else
!!!---------------------------------
!R!     call output_gulp(CIFFiles(i),filename(i))
!R!     ! here we interact with GULP using the call system, very expensive in call mem
!R!     write(script,'(a,a,a,a,a,a,a,a,a,a,a,a,a,a,a)')"~/bin/gulp < ",&
!R!      filename(i)(1:Clen_trim(filename(i)))," > ",filename(i)(1:Clen_trim(filename(i))),".gout ; ",&
!R!      "grep 'Total lattice energy       =' ",filename(i)(1:Clen_trim(filename(i))),&
!R!      ".gout | grep 'eV' | awk '{print $5}' > ",filename(i)(1:Clen_trim(filename(i))),".tmp ;",&
!R!      "grep 'ERROR' ",filename(i)(1:Clen_trim(filename(i))),&
!R!      ".gout | wc -l | awk '{print $1}' >> ",filename(i)(1:Clen_trim(filename(i))),".tmp "
!R!     call system(script)
!R!     open(newunit = u,file=filename(i)(1:Clen_trim(filename(i)))//".tmp")
!R!     read(u,'(a)')line
!R!    if(line(1:20)=="********************")then
!R!      CIFFiles(i)%cal_energy=infinite
!R!     else
!R!      read(line,*) CIFFiles(i)%cal_energy
!R!     end if
!R!    end if LMP_or_GULP 
    ! CHECK: PROBLEMS IN PARALELL
    !read(u,*) jjj
    !if (jjj>0) then
    ! CIFFiles(i)%cal_energy = infinite
    ! write(6,*)'There are errors in output', jjj
    ! stop '#'
    !end if
!temp    close(u)
    ! debug:
    !write(6,*) "file:", i, CIFFiles(i)%filename, OMP_GET_THREAD_NUM(), u, CIFFiles(i)%cal_energy
!Nr   end do !..scan_
!..!$omp end do
!..!$omp end parallel
!..!$omp barrier
   cal_energy_min=minval(CIFFiles%cal_energy)
   do i=1,n_files
    CIFFiles(i)%cal_energy = CIFFiles(i)%cal_energy - cal_energy_min
    fitness = fitness + &
    0.5*CIFFiles(i)%obs_energy_weight*abs(CIFFiles(i)%obs_energy-(CIFFiles(i)%cal_energy))**2/real(n_files)
   end do
  else
   fitness = fitness + penalty
  end if calgulp
  !write(6,*)'Rosenbluth:',(exp(-CIFFiles(j)%obs_energy/(10*obs_energy_min))/partition, j=1,n_files)
  !write(6,*)'Suma:',sum( exp(-CIFFiles(1:n_files)%obs_energy/10.0/obs_energy_min)/partition )
  !SEIJAS podria llamar rm struc/*inp struc/* out,....  call system("rm -rf *.tmp *.gout")
  return
 end function Fitness
!
 subroutine WriteCitizen(k,kk,kkk,compound,lod) !,vgh)
  implicit none
  integer,intent(in)             :: k,kk,compound,lod
  real,intent(in)                :: kkk ! biodiversity
  integer                        :: i
  character(len=100)             :: fmt_
  character(len=32*np(compound)) :: wnowaste
  real                           :: wnowasteparam(1:32*np(compound)),wfitness
  do i=1,32*np(compound)
   wnowaste(i:i)=' '
  end do
  wnowaste = parents(k)%genotype
  do i=1,np(compound)
   wnowasteparam(i) = parents(k)%phenotype(i)
  end do
  wfitness = parents(k)%fitness
  write(6,'(i2,1x,i5,1x,a32,1x,e25.12,1x,e25.12,1x,a,1x,a,e25.12,a,1x,a)') &
       k,kk,wnowaste(1:32),&
       wnowasteparam(1),wfitness,'[Fitness]'  !,'(',kkk,')','[Similarity]' !,k
  do i=2,np(compound)
   if(lod>0.and.i==3)then
    write(6,'(9x,a32,1x,e25.12,10x,a,1x,i2,a)')&
        wnowaste(1+32*(i-1):32*i),wnowasteparam(i),'Finishing:',lod,'/10'
   else
    write(6,'(9x,a32,1x,e25.12)')wnowaste(1+32*(i-1):32*i),wnowasteparam(i)
   end if
  end do
 end subroutine WriteCitizen
!
 subroutine SortByFitness()
  type(typ_ga)             ::  sorted(1:ga_size)
  integer                  ::  k,i
  real(16)                 ::  ftnss(1:ga_size)
  do k=1,ga_size
   ftnss(k)=dble(parents(k)%fitness)
   if(isnan(parents(k)%fitness)) ftnss(k) = 9999999999.d99
  end do
  call QsortC( ftnss )
  exter:do k = 1, ga_size ! <- ordered
   inter:do i = 1, ga_size
   if( dble(parents(i)%fitness) == ftnss(k))then
     sorted(k) = parents(i)
     cycle inter
   end if
   end do inter
  end do exter
  parents=sorted
  return
 end subroutine SortByFitness
!
 pure character(len=32) function real2bin(x)
  implicit none
  real,intent(in) :: x
  write(real2bin(1:32) ,'(b32.32)') x
  return
 end function real2bin
!
 real function Biodiversity(compound,agent)
! calculates the standard deviation of the fitness centered in the first agent
! lowest value in the list
  implicit None
  integer,intent(in)          :: Compound
  type(typ_ga), intent(inout) :: Agent(1:ga_size)
  integer                     :: k
  Biodiversity = sqrt(sum( ( agent(2:ga_size)%fitness - agent(1)%fitness)**2))
  Biodiversity = Biodiversity/real(ga_size -1)
  return
 end function Biodiversity
!
 subroutine Mutate( agent , compound )
  implicit none
  type(typ_ga), intent(inout) :: agent
  integer                     :: spos, compound
  real                        :: r
  r = r4_uniform(0.0,1.0)
  if ( r < real(np(compound)-1)/real(np(compound))) then ! (np-1)/np probability
   spos = randint(1, 32*np(compound) )
   if( agent%genotype(spos:spos) == achar(48) ) then ! 0 -> 1
    agent%genotype(spos:spos) = achar(49)
   else                                              ! 1 -> 0
    agent%genotype(spos:spos) = achar(48)
   end if
  else                                                   ! 1/np      probability
   do i = 1,np(compound)
    spos = randint(32*(i-1)+1,32*i)
    if( agent%genotype(spos:spos) == achar(48) ) then
     agent%genotype(spos:spos) = achar(49)
    else
     agent%genotype(spos:spos) = achar(48)
    end if
   end do
  end if
  return
 end subroutine Mutate
!
 subroutine NuclearDisaster(Compound,n_files,CIFFiles)
! all the poor people mutate, the elistists survive
  implicit none
  integer,intent(in)      ::  Compound,n_files
  type(CIFFile),intent(inout):: CIFFiles(n_files)
  integer                 :: k = 0, i, j
  real                    :: rrr
  do i = GA_ELITISTS + 1, GA_Size
   do j=1,32*np(compound)
    Children%genotype(j:j) = achar(randint(48,49))
   end do
   call UpdateCitizen(Children(i),Compound,n_files,CIFFiles)
  end do
  return
 end subroutine NuclearDisaster
!
 subroutine Swap()
  if (associated(parents, target=pop_alpha)) then
   parents => pop_beta
   children => pop_alpha
  else
   parents => pop_alpha
   children => pop_beta
  end if
  return
 end subroutine Swap
!
 subroutine Elitism()
  children(:GA_ELITISTS) = parents(:GA_ELITISTS)
  return
 end subroutine
!
 subroutine Mate(compound,n_files,CIFFiles)
  integer                     :: i, i1, i2, spos
  integer, intent(in)         :: compound,n_files
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  real                :: rrr
  call Elitism()
  do i = GA_ELITISTS + 1, ga_size
   ! Crossover:
   ! {{ eleccion random del primer 50% de la tabla
   call choose_randomly(i1,i2)
   ! }}
   ! {{ eleccion proporcionalmente a su fitness
   !call choose_propto_fitness(i1,i2)
   !write(6,*)i1,i2
   ! }}
   spos = randint(0, 32*np(compound) )
   children(i)%genotype = parents(i1)%genotype(:spos) // parents(i2)%genotype(spos+1:)
   ! Mutate and NuclearDisaster:
   rrr = r4_uniform(0.0,1.0)
   if ( rrr < GA_MUTATIONRATE) then
    call Mutate(children(i),compound)
   else if ( rrr >= GA_MutationRate .and. rrr <= GA_MutationRate + GA_DisasterRate ) then
    call NuclearDisaster(Compound,n_files,CIFFiles)
    return
   end if
  end do
  do i = 1, ga_size
   call UpdateCitizen(children(i),compound,n_files,CIFFiles)
  end do
  return
 end subroutine Mate
!
 subroutine choose_randomly(j1,j2)
  implicit none
  integer,intent(out) :: j1,j2
  j1  = randint(1, int(ga_size/2))
  j2  = randint(1, int(ga_size/2))
  do while ( j1 == j2 )
   j2 = randint(1, int(ga_size/2))
  end do
  return
 end subroutine choose_randomly
!
 subroutine choose_propto_fitness(j1,j2)
  implicit none
  integer,intent(out) :: j1,j2
  integer             :: i
  real                :: ftnss(ga_size),prop(0:ga_size)=0.0,rrr1,rrr2
  real                :: infinity = HUGE(2147483647)
  rrr1 = 0.0
  do i = 1, ga_size
   ftnss(i) = 1.0/parents(i)%fitness
   if ( isnan( parents(i)%fitness ) ) ftnss(i) = 0.0
   if ( parents(i)%fitness > infinity ) ftnss(i) = 0.0
   prop(i) = ftnss(i)
   if( ftnss(i) >= infinity ) then
     rrr1 = rrr1 + infinity
   else
     rrr1 = rrr1 + ftnss(i)
   end if
  end do
  prop = prop / rrr1
  ! select 1:
   rrr1 = r4_uniform(0.0,1.0)
   slct1: do i=1,ga_size
    if(rrr1<=prop(i-1).and.rrr1>prop(i))then
     j1 = i
    end if
   end do slct1
   ! select 2:
   rrr2 = r4_uniform(0.0,1.0)
   do while ( rrr1 == rrr2 )
    rrr2 = r4_uniform(0.0,1.0)
   end do
   slct2: do i=1,ga_size
    if(rrr2<=prop(i-1).and.rrr2>prop(i))then
     j2 = i
    end if
   end do slct2
  return
 end subroutine choose_propto_fitness
!
 subroutine Fit(Compound,n_files,CIFFiles)
  implicit none
  integer,intent(in)          :: Compound, n_files
  type(CIFFile),intent(inout) :: CIFFiles(n_files)
  integer,parameter           :: maxstep = 100, minstep = 10  !maxstep era 100
  integer                     :: kk, ii, i, j, k,vgh, l
  real                        :: fit0 = 0.0
  real                        :: eps
  character(len=100)          :: string, line
  write(6,'(a)')'Initialising GA World:'
  write(6,'(i5,1x,a,1x,i5,1x,a,f14.7,1x,a)')ga_size,'(elites:',int(ga_eliterate*ga_size),'mutate:',ga_mutationrate,')'
  write(6,'(a)')'[...]'
  if ( refit_flag ) then
   do i=1,n_refits
    write(6,'(a,1x,i4)')'Initial Agent from input file',i
    pop_alpha(i)%genotype=string_IEEE(i)
    call UpdateCitizen(pop_alpha(i),compound,n_files,CIFFiles)
    write(string,'("(a13, ", i4, "a, a1)" )') 32*np(Compound)
    write(6,string) '> genotype: [',pop_alpha(i)%genotype(1:32*np(Compound)),']'
    write(6,'(a,10(f14.7,1x))') '> phenotype: ',( pop_alpha(i)%phenotype(j), j=1,np(compound) )
    write(6,'(a,1x,e25.12)')'Fitness:',pop_alpha(i)%fitness
    write(6,'(a)')'--------------------------'
   end do
   finish_make: do i=n_refits+1,ga_size
    if(i>ga_size) exit finish_make
    write(6,'(a,1x,i4)') 'Initial agent (random generation)',i
    pop_alpha(i) = NewCitizen(compound,n_files,CIFFiles)
    !SEIJAS
    !The next UpdateCitizen was here, but I think is redundant.
    !call UpdateCitizen(pop_alpha(i),compound,n_files,CIFFiles)
    write(string,'("(a13, ", i4, "a, a1)" )') 32*np(Compound)
    write(6,string) '> genotype: [',pop_alpha(i)%genotype(1:32*np(Compound)),']'
    write(6,'(a,10(f14.7,1x))') '> phenotype:',( pop_alpha(i)%phenotype(j), j=1,np(compound) )
    write(6,'(a,1x,e25.12)')'Fitness:',pop_alpha(i)%fitness
    write(6,'(a)')'--------------------------'
   end do finish_make
   write(6,'(a)')'[...]'
   parents =>  pop_alpha
   children => pop_beta
   call SortByFitness()
  else
   pop_alpha = [(NewCitizen(compound,n_files,CIFFiles), i = 1,ga_size)]
   parents =>  pop_alpha
   children => pop_beta
  end if
  ii = 0 ; kk = 0
  write(6,'(a)') 'Go Down the Rabbit Hole > '
  write(6,'(a)') '[Only showing the elite ...]'
  converge: do while ( .true. )
   ii=ii+1
   call SortByFitness()
   fit0 = fitness( parents(1)%phenotype,Compound,n_files,CIFFiles)
   call WriteEnergies(n_files,CIFFiles,"res")
   eps=0.0 !antes ponia esto, pero no hace falta la desv_stnd !eps  = Biodiversity( compound, children)
   do i = 1, GA_ELITISTS
    call WriteCitizen(i,ii,eps,compound,kk) !int(0.5*ga_size*ga_size-ga_size) )
   end do
   !TOCADO EL 26 DIC
  !SEIJAS To write file.dat ready to run MD in LAMMPS.
  line="cp "//CIFFiles(1)%filename(1:Clen_trim(CIFFiles(1)%filename)-4)//".sav for_MD_LAMMPS.dat"
  call system(line)
  do l = 1,np(compound)
   if (l.le.9) then
    write(line,'(a,i1,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",parents(1)%phenotype(l),"/g' for_MD_LAMMPS.dat"
   else
    write(line,'(a,i2,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",parents(1)%phenotype(l),"/g' for_MD_LAMMPS.dat"
   end if
   call system(line)
  end do
  !HASTA AQUI 26 DIC
   if( ii>=minstep .and. parents(1)%fitness <= 0.1 .and. abs(parents(1)%fitness-fit0) <= 1e-4)then
    kk = kk + 1
   else
    kk = 0
   end if
   if ( ii >= maxstep .or. kk >= 10 ) exit converge
   call Mate(compound,n_files,CIFFiles)
   call Swap()
  end do converge
  call SortByFitness()
  !SEIJAS lo he comentado pporque no me vale y asi no dara problemas. call WriteLib(compound,children(1)%phenotype)
  do i = 0, np( compound )-1
   param( compound,i ) = children(1)%phenotype(i+1)
  end do
  write(6,*)'#',(param(compound,i),i=0,np(compound )-1)
  write(6,*)'#','Fitness:',fit0  !,'Similarity:',eps
  !SEIJAS To write file.dat ready to run MD in LAMMPS.
  line="cp "//CIFFiles(1)%filename(1:Clen_trim(CIFFiles(1)%filename)-4)//".sav for_MD_LAMMPS.dat"
  call system(line)
  do l = 1,np(compound)
   if (l.le.9) then
    write(line,'(a,i1,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",param(compound,l-1),"/g' for_MD_LAMMPS.dat"
   else
    write(line,'(a,i2,a,e20.15,a,a)')"sed -i 's/pairIJ_",l," /",param(compound,l-1),"/g' for_MD_LAMMPS.dat"
   end if
   call system(line)
  end do
  return
 end subroutine fit
end module mod_genetic
!
!module mod_simplex
!use mod_random
!use flama_globals
!use qsort_c_module
!use mod_genetic
!implicit none
!private
!public  :: fit_simplex
!contains
! ======================================================================
! Interface between program and module:
!subroutine fit_simplex()
!implicit none
!real               :: e = 1.0e-4, scale = 1.0
!integer            :: iprint = 0
!integer            :: i,j
!real               :: sp(0:np(compound)-1)
!type(typ_ga)       :: axolotl
!write(6,'(a)')' '
!write(6,'(a)')'Minimising the cost function using the Nelder-Mead SIMPLEX method:'
!write(6,*)'# Compound:',compound, np(compound)
!do i = 0,np(compound)-1
! axolotl%phenotype(i+1) = param(compound,i)
! sp(i) = axolotl%phenotype(i+1)
!end do
!call simplex(sp,compound,np(compound),e,scale,iprint)
!do i= 0, np(compound)-1
! axolotl%phenotype(i+1)=sp(i)
! param(compound,i)=sp(i)
!end do
!write(111,*)'#',(param(compound,i),i=0,np(compound )-1)
!write(111,*)'#','Fitness:',func(np(compound),sp,compound)
!end subroutine fit_simplex

! ======================================================================
! This is the function to be minimized
!real function func(n,x,compound) result(rosen)
! implicit none
! integer,intent(in)   :: n
! real,   intent (in)  :: x(0:n-1)
! integer,intent(in)   :: compound
! real                 :: sp(1:np(compound))
! integer              :: i
! sp = 0.0
! do i=0,np(compound)-1
!  sp(i+1) = x(i)
! end do
! rosen = fitness(sp,compound)
! return
!end function func
!! This is the simplex routine
!! Michael F. Hutt
!subroutine simplex(start, compound, n, EPSILON, scale, iprint)
! implicit none
! integer, intent (in)                   :: n, iprint, compound
! real, intent (inout), dimension(0:n-1) :: start
! real, intent (in)                      :: EPSILON, scale
! integer, parameter :: MAX_IT = 1000000
! real, parameter    :: ALPHA=1.0
! real, parameter    :: BETA=0.5
! real, parameter    :: GAMMA=2.0
! ======================================================================
! Variable Definitions
! vs = vertex with the smallest value
! vh = vertex with next smallest value
! vg = vertex with largest value
! i,j,m,row
! k = track the number of function evaluations
! itr = track the number of iterations
! v = holds vertices of simplex
! pn,qn = values used to create initial simplex
! f = value of function at each vertex
! fr = value of function at reflection point
! fe = value of function at expansion point
! fc = value of function at contraction point
! vr = reflection - coordinates
! ve = expansion - coordinates
! vc = contraction - coordinates
! vm = centroid - coordinates
! min
! fsum,favg,s,cent
! vtmp = temporary array passed to FUNC
! ======================================================================
! Integer :: vs,vh,vg
! Integer :: i,j,k,itr,m,row
! real, dimension(:,:), allocatable :: v
! real, dimension(:), allocatable  :: f
! real, dimension(:), allocatable :: vr
! real, dimension(:), allocatable :: ve
! real, dimension(:), allocatable :: vc
! real, dimension(:), allocatable :: vm
! real, dimension(:), allocatable :: vtmp
! real :: pn,qn
! real :: fr,fe,fc
! real :: min,fsum,favg,cent,s

! allocate (v(0:n,0:n-1))
! allocate (f(0:n))
! allocate (vr(0:n-1))
! allocate (ve(0:n-1))
! allocate (vc(0:n-1))
! allocate (vm(0:n-1))
! allocate (vtmp(0:n-1))

! create the initial simplex
! assume one of the vertices is 0.0

! pn = scale*(sqrt(n+1.)-1.+n)/(n*sqrt(2.))
! qn = scale*(sqrt(n+1.)-1.)/(n*sqrt(2.))

! DO i=0,n-1
!   v(0,i) = start(i)
! END DO

! DO i=1,n
!   DO j=0,n-1
!     IF (i-1 == j) THEN
!       v(i,j) = pn + start(j)
!     ELSE
!       v(i,j) = qn + start(j)
!     END IF
!   END DO
! END DO


! find the initial function values

! DO j=0,n
! put coordinates into single dimension array
! to pass it to FUNC
!   DO m=0,n-1
!     vtmp(m) = v(j,m)
!   END DO
!   f(j)=FUNC(n,vtmp,compound)
! END DO

! Print out the initial simplex
! Print out the initial function values
! find the index of the smallest value for printing
! IF (iprint == 0) THEN
!  vs=0
!  DO j=0,n
!   If (f(j) .LT. f(vs)) Then
!     vs = j
!   END IF
!  END DO
! print out the value at each iteration
!  Write(6,'(a)') "Initial Values from genetic algorithm:"
!  Write(6,*) (v(vs,j),j=0,n-1),'Fit:',f(vs)
! END IF
! k = n+1
! begin main loop of the minimization

!DO itr=1,MAX_IT
! find the index of the largest value
! vg = 0
! DO j=0,n
!   IF (f(j) .GT. f(vg)) THEN
!     vg = j
!   END IF
! END DO

! find the index of the smallest value
! vs = 0
! DO j=0,n
!   If (f(j) .LT. f(vs)) Then
!     vs = j
!   END IF
! END DO

! find the index of the second largest value
! vh = vs
! Do j=0,n
!   If ((f(j) .GT. f(vh)) .AND. (f(j) .LT. f(vg))) Then
!     vh = j
!   END IF
! END DO

! calculate the centroid
! DO j=0,n-1
! cent = 0.0
!   DO m=0,n
!     If (m .NE. vg) Then
!       cent = cent + v(m,j)
!     END IF
!   END DO
!   vm(j) = cent/n
! END DO

! reflect vg to new vertex vr
! DO j=0,n-1
!   vr(j) = (1+ALPHA)*vm(j) - ALPHA*v(vg,j)
! END DO
! fr = FUNC(n,vr,compound)
! k = k+1

! If ((fr .LE. f(vh)) .AND. (fr .GT. f(vs))) Then
!   DO j=0,n-1
!     v(vg,j) = vr(j)
!   END DO
!   f(vg) = fr
! END IF

! investigate a step further in this direction
! If (fr .LE. f(vs)) Then
!   DO j=0,n-1
!     ve(j) = GAMMA*vr(j) + (1-GAMMA)*vm(j)
!   END DO
!   fe = FUNC(n,ve,compound)
!   k = k+1

! by making fe < fr as opposed to fe < f(vs), Rosenbrocks function
! takes 62 iterations as opposed to 64.

!   If (fe .LT. fr) Then
!     DO j=0,n-1
!       v(vg,j) = ve(j)
!     END DO
!     f(vg) = fe
!   Else
!     DO j=0,n-1
!       v(vg,j) = vr(j)
!     END DO
!     f(vg) = fr
!   END IF
! END IF

! check to see if a contraction is necessary
! If (fr .GT. f(vh)) Then
!   DO j=0,n-1
!     vc(j) = BETA*v(vg,j) + (1-BETA)*vm(j)
!   END DO
!   fc = FUNC(n,vc,compound)
!   k = k+1
!   If (fc .LT. f(vg)) Then
!     DO j=0,n-1
!       v(vg,j) = vc(j)
!     END DO
!   f(vg) = fc

! at this point the contraction is not successful,
! we must halve the distance from vs to all the
! vertices of the simplex and then continue.
! 10/31/97 - modified C program to account for
! all vertices.

! Else
!   DO row=0,n
!     If (row .NE. vs) Then
!       DO j=0,n-1
!         v(row,j) = v(vs,j)+(v(row,j)-v(vs,j))/2.0
!       END DO
!     END IF
!   END DO
!   DO m=0,n-1
!     vtmp(m) = v(vg,m)
!   END DO
!   f(vg) = FUNC(n,vtmp,compound)
!   k = k+1

!   DO m=0,n-1
!     vtmp(m) = v(vh,m)
!   END DO
!   f(vh) = FUNC(n,vtmp,compound)
!   k = k+1
!   END IF
! END IF
! find the index of the smallest value for printing
! !vs=0
! !DO j=0,n
! !  If (f(j) .LT. f(vs)) Then
! !    vs = j
! !  END IF
! !END DO
! print out the value at each iteration
! !IF (iprint == 0) THEN
! !  Write(6,*) "Iteration:",itr,(v(vs,j),j=0,n-1),'Value:',f(vs)
! !END IF
! test for convergence
! fsum = 0.0
! DO j=0,n
!   fsum = fsum + f(j)
! END DO
! favg = fsum/(n+1.)
! !s = 0.0
! !DO j=0,n
! !  s = s + ((f(j)-favg)**2.)/n
! !END DO
! !s = sqrt(s)
! If (favg .LT. EPSILON.or.itr==MAX_IT) Then

! print out the value at each iteration
!  Write(6,'(a,1x,i6)') "Final Values:", itr
!  Write(6,*) (v(vs,j),j=0,n-1),'Fit:',f(vs)
!  IF(itr/=MAX_IT)then
!   write(6,'(a,1x,f14.7,1x,a,1x,f14.7)')'Nelder-Mead has converged:',favg,'<',epsilon
!  else
!   write(6,'(a,1x,i6,1x,a,1x,i6)')'Maximun number of steps:',itr,'=',MAX_IT
!  end if
!  EXIT ! Nelder Mead has converged - exit main loop
! END IF
!END DO
! end main loop of the minimization
! find the index of the smallest value
! vs = 0
! DO j=0,n
!   If (f(j) .LT. f(vs)) Then
!     vs = j
!   END IF
! END DO
!  print out the minimum
! DO m=0,n-1
!   vtmp(m) = v(vs,m)
! END DO
! min = FUNC(n,vtmp,compound)
! !write(6,*)'The minimum was found at ',(v(vs,k),k=0,n-1)
! !write(6,*)'The value at the minimum is ',min
! DO i=0,n-1
!  start(i)=v(vs,i)
! END DO
!250  FORMAT(A29,F7.4)
!300  FORMAT(F11.6,F11.6,F11.6)
! return
! end subroutine simplex
! ======================================================================
!end module mod_simplex
!
program flama
 use,intrinsic :: iso_fortran_env
 use types
 use mod_random
 use get_structures
 use flama_globals
 use mod_genetic
 !use mod_simplex
 implicit none
 type(CIFfile),allocatable       :: CIFFiles(:)
 integer                         :: n_files = 0
 write(6,'(4a)')'This file was compiled by ',compiler_version(),&
                ' using the options ',compiler_options()
 call read_input()
 call init_random_seed()
 !call GenerateCIFFileList()
 !First run for debuging
 call ReadListOfCIFFiles(n_files)
 allocate( CIFFiles(1:n_files) )
 call ReadCIFFiles(n_files,CIFFiles)
 call WriteEnergies(n_files,CIFFiles,"ini")
 call ReadObservables(n_files,CIFFiles)
 if (flag) then
  call Fit(1,n_files,CIFFiles)
  !call fit_simplex()
 end if
 call WriteEnergies(n_files,CIFFiles,"end")
 deallocate(CIFFiles)
end program flama
