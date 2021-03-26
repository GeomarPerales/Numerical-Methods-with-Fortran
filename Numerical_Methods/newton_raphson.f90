program newton_raphson

implicit none

real,dimension(:),allocatable::X,F,DF,Err
real,dimension(:,:),allocatable::TAB
integer::i,j,imax,n
real::vi,Error

Error=0.01


write(*,*)'numero de iteraciones: '
read(*,*) n		!numero de iteraciones 
write(*,*)'valor inicial: '
read(*,*) vi	!valor inicial

write(*,*) 'error de aproximación'
write(*,*) Error

allocate(X(n))
allocate(F(n-1))
allocate(DF(n-1))
allocate(Err(n-1))

X(1)=vi

do i=1,n-1

F(i)=X(i)**3-X(i)-1	 !funcion
DF(i)=3*X(i)**2-1 !derivada de la funcion
X(i+1)=X(i)-F(i)/DF(i) !metodo numerico Newton raphson
Err(i)=abs(X(i+1)-X(i))

if (Err(i)<Error) then
	write(*,*)'El valor de la raiz es:'
	write(*,*) X(i)
	imax=i
	write(*,*)'El numero maximo de iteraciones fue: '
	write(*,*) imax
	exit
end if

end do

allocate(TAB(n-1,2))

do i=1,imax+1
	TAB(i,1)=X(i)
	TAB(i,2)=Err(i)
end do
											
write(*,*) '----------------------------------- '
write(*,*) 'tabla de valores '
write(*,*) '----------------------------------- '
write(*,*) ' ------- raiz -------- error ------ '

do i=1,imax+1
write(*,*) (TAB(i,j),j=1,2)
end do

end program

