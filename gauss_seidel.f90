program gauss_seidel
implicit none

real,dimension(:,:),allocatable::MAT,A
real,dimension(:),allocatable::B
real,dimension(:,:),allocatable::X,e
integer::n,i,j,k
real::S1,S2,tol


write(*,*)'ingresar el numero de incognitas:'
read(*,*)n

write(*,*)'-------------------------------------'
write(*,*)' - tolerancia -'
write(*,*)' ejemplo: 0.1, 0.01, etc'
write(*,*)'-------------------------------------'
write(*,*)'ingresar la tolerancia del error:'
read(*,*)tol


!leyendo matriz almacenada en archivo matrix txt
open(unit=25,file='matrix.txt',status='old',action='read')
allocate(MAT(n,n+2))

do i=1,n
read(25,*) (MAT(i,j),j=1,n+2)
end do
close(25)

!Metodo doolitle
allocate(A(n,n))
allocate(B(n))
allocate(X(n,n))
allocate(e(n,n))


do i=1,n 
	do j=1,n
	A(i,j)=MAT(i,j)
	end do
end do

do i=1,n 
B(i)=MAT(i,n+1)
end do

k=1
do i=1,n 
X(k,i)=MAT(i,n+2)
end do

!
10 do i=1,n
	s1=0
	do j=1,n
		if (j>i) then
		s1=s1+A(i,j)*X(k,j)
		end if
	end do
s2=0
do j=1,n
if (j<i) then
s2=s2+A(i,j)*X(k+1,j)
end if
end do
X(k+1,i)=(B(i)-s1-s2)/A(i,i)
end do

!
do i=1,n
e(k,i)=abs(X(k+1,i)-X(k,i))
	if (e(k,i)<=tol) then
	else
	k=k+1
	goto 10	
	end if
end do



open(21,file="res.txt",action='write',status='replace')
do i=1,k+1
write(21,*) (X(i,j),j=1,n)
end do
close(21)

open(22,file="error.txt",action='write',status='replace')
do i=1,k
write(22,*) (e(i,j),j=1,n)
end do
close(22)


end program
