program doolittle
implicit none

real,dimension(:,:),allocatable::MAT,A,L,U
real,dimension(:),allocatable::B
real,dimension(:),allocatable::Y,X
integer::n,i,j,k
real::S1,S2


write(*,*)'ingresar el numero de incognitas'
read(*,*)n

!leyendo matriz almacenada en archivo matrix txt
open(unit=11,file='matrix.txt',status='old',action='read')
allocate(MAT(n,n+1))

do i=1,n
read(11,*) (MAT(i,j),j=1,n+1)
end do
close(11)

!Metodo doolitle
allocate(A(n,n))
allocate(B(n))
allocate(L(n,n))
allocate(U(n,n))
allocate(X(n))
allocate(Y(n))


do i=1,n 
	do j=1,n
	A(i,j)=MAT(i,j)
	end do
end do

do i=1,n 
B(i)=MAT(i,n+1)
end do

!generando matrices vacías
do i=1,n
do j=1,n
	L(i,j)=0
	U(i,j)=0
	Y(i)=0
	X(i)=0
end do
end do

!metodo doolitle - LU  - nakamura
do j=1,n
do i=1,n
if (i<=j)then
	U(i,j)=A(i,j)
		do k=1,i-1
	U(i,j)=U(i,j)-L(i,k)*U(k,j)
		end do
end if

if (j<=i) then
	L(i,j)=A(i,j)
	do k=1,j-1
	L(i,j)=L(i,j)-L(i,k)*U(k,j)
	end do
L(i,j)=L(i,j)/U(j,j)
end if

end do
end do

!MATRIZ L - B
Y(1)=B(1)

do i=2,n
S1=0
	do j=1,i-1
	S1=S1+L(i,j)*Y(j)
	end do
	Y(i)=B(i)-S1
end do

do i=n-1,1,-1
S2=0
	do j=i+1,n
	S2=S2+U(i,j)*X(j)
	end do
	X(i)=(Y(i)-S2)/U(i,i)
end do



!imprimiendo en txt la matriz final aumentada
open(12,file="matL.txt",action='write',status='replace')
do i=1,n
write(12,*) (L(i,j),j=1,n)
end	do
close(12)

open(13,file="matU.txt",action='write',status='replace')
do i=1,n
write(13,*) (U(i,j),j=1,n)
end	do
close(13)

open(14,file="matX.txt",action='write',status='replace')
do i=1,n
write(14,*) X(i)
end	do
close(14)

end program
