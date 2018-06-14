program gauss_jordan_matrix
implicit none

real,dimension(:,:),allocatable::MAT
real,dimension(:,:),allocatable::A
real,dimension(:),allocatable::B,s
integer::n,i,j

write(*,*)'ingresar el numero de incognitas'
read(*,*)n

!leyendo matriz almacenada en archivo matrix txt
open(unit=11,file='matrix.txt',status='old',action='read')
allocate(MAT(n,n+1))

do i=1,n
read(11,*) (MAT(i,j),j=1,n+1)
end do
close(11)

allocate(A(1:n,1:n))
allocate(B(1:n))
allocate(S(1:n))

write(*,*)'----- Matriz A -----'
do i=1,n 
	do j=1,n
	A(i,j)=MAT(i,j)
	end do
end do

do i=1,n
write(*,*) (A(i,j),j=1,n)
end do

write(*,*)'----- Matriz B -----'
do i=1,n 
B(i)=MAT(i,n+1)
write(*,*) B(i)
end do

write(*,*)'----- Matriz Aumentada -----'
do i=1,n 
write(*,*) (MAT(i,j),j=1,n+1)
end do

!eliminacion de Gauss / Gauss - Jordan 
open(12,file="matpap.txt",action='write',status='replace')
do i=1,n
	MAT(i,:)=MAT(i,:)/MAT(i,i)
	do j=1,n
		if (i/=j) then
	write(12,*) MAT(j,:)
	MAT(j,:)=MAT(j,:)-MAT(i,:)*MAT(j,i)
		end if
		
	end do
end do

!vector solution
S=MAT(1:n,n+1)

!imprimiendo en txt la matriz final aumentada
open(13,file="matfinal.txt",action='write',status='replace')
do i=1,n
write(13,*) (MAT(i,j),j=1,n+1)
end	do
close(13)

!imprimiendo en txt el vector solucion
open(14,file="solution.txt",action='write',status='replace')
do i=1,n
write(14,*) S(i)
end	do
close(14)

end program
