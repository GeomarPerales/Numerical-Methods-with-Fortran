program integral_rectangulos

implicit none

real::a,b,h,VI
real,dimension(:),allocatable::F
integer::i,n

write(*,*)'limite inferior: '
read(*,*) a
write(*,*)'limite superior: '
read(*,*) b
write(*,*)'numero de rectangulos: '
read(*,*) n

h=(b-a)/n

VI=0 !valor inicial de la integral

allocate(F(n))

do i=1,n
	F(i)=(a+h*i)**2 !funcion a integrar
	VI=VI+f(i)*h !sumatoria de area	
end do

write(*,*) VI

end program
