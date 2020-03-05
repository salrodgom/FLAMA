program bin2real
 implicit none
 integer                     :: i
 character(len=32)           :: genotype
 real                        :: phenotype
 do
  read(5,'(a32)') genotype
  read(genotype(1:32),'(b32.32)')phenotype
  write(6,'(f20.10)') phenotype
 end do
end program bin2real
