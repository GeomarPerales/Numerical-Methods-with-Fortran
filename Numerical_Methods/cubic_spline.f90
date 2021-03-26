program interpolation_cubica

real,dimension(:,:),allocatable::MAT
real,dimension(:,:),allocatable::V,AUM,P,A,B,M,TAB
real::MULT,R,ve
real,dimension(:),allocatable::X,Y
integer::n,i,j,k,L

n=4;
write(*,*)'interpolacion de grado:'
write(*,*) n-1

!leyendo matriz almacenada en archivo matrix txt
open(unit=11,file='data.txt',status='old',action='read')
allocate(MAT(n,2))

do i=1,n
read(11,*) (MAT(i,j),j=1,n-2)
end do
close(11)

allocate(V(n,n))
allocate(AUM(n,2*n))
allocate(A(n,n))
allocate(B(n,1))
allocate(P(n,n))
allocate(M(n,1))

!matriz de vandermonde

 do i=1,n 
	do j=1,n
	V(i,j)=MAT(i,1)**(j-1)
	end do
end do

!aumentada
do i=1,n 
	do j=1,n
	AUM(i,j)=V(i,j)
	do k=n+1,2*n
			if (k==i+n) then 
			AUM(i,k)=1
			else
			AUM(i,k)=0
			end	if
		end do
	end do
end do

!eliminacion de Gauss / Gauss - Jordan 
do k=1,n-1		
	do i=k+1,n
	MULT=AUM(i,k)/AUM(k,k)
		do j=k,2*n
	AUM(i,j)=AUM(i,j)-MULT*AUM(k,j)	
		end do
	end do
end do

do i=1,n
	MULT=AUM(i,i)
	do j=i,2*n
		AUM(i,j)=AUM(i,j)/MULT
	end do
end do

do k=n-1,1,-1
	do i=1,k
	MULT=AUM(i,K+1)
		do j=k,2*n
		AUM(i,j)=AUM(i,j)-AUM(k+1,j)*MULT
		end do
	end do
end do

!asignado valores a la matriz A y B
do i=1,n
	do j=n+1,2*n
	A(i,j-n)=AUM(i,j)
	end do
end	do

do i=1,n
	do j=2,n-2
	B(i,1)=MAT(i,j)
	end do
end do

write(*,*)'--------------------------------------------------'

!producto de matrices
M=MATMUL(A,B)

write(*,*)'coeficientes del polinomio'
write(*,*) M

write(*,*)'--------------------------------------------------'

write(*,*)'ingresar valor evaluado'
read(*,*)ve
R=M(4,1)*ve**(n-1)+M(3,1)*ve**(n-2)+M(2,1)*ve**(n-3)+M(1,1)


write(*,*)'valor evaluado en la interpolación'

write(*,*) R
write(*,*)'--------------------------------------------------'
!serie de datos interpolados

write(*,*)'insertar el limite de datos a interpolar:'
read(*,*) L

allocate(X(L))
allocate(Y(L))
allocate(TAB(L,2))

do i=1,L
	X(i)=i
	Y(i)=M(4,1)*X(i)**(n-1)+M(3,1)*X(i)**(n-2)+M(2,1)*X(i)**(n-3)+M(1,1)
end do

do i=1,L
	TAB(i,1)=X(i)
	TAB(i,2)=Y(i)
end do
write(*,*)'--------------------------------------------------'
write(*,*)'-------- X ---------Yajustado---'
do i=1,L
write(*,*) (TAB(i,j),j=1,2)
end do

open(12,file="res_interpolacion.txt",action='write',status='replace')
write(12,*)'----- X ------ Yajustado ----------'
do i=1,L
write(12,*) (TAB(i,j),j=1,2)
end do
close(12)
end program
