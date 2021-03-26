program minimos_cuadrados
implicit none
real,dimension(:,:),allocatable::MAT
real,dimension(:,:),allocatable::TAB

real,dimension(:),allocatable::X
real,dimension(:),allocatable::Y,Ya,err
real::r,sx,sy,sxy,sx2,sy2,se,erm
real::a,b
integer::n,i,j

write(*,*)'ingresar el tamaño de los datos'
read(*,*)n

!leyendo matriz almacenada en archivo matrix txt
open(unit=12,file='data.txt',status='old',action='read')
allocate(MAT(n,2))
do i=1,n
read(12,*) (MAT(i,j),j=1,2)
end do
close(12)

write(*,*)'-------------------------------------'
write(*,*)'DATA IMPORTADA'
do i=1,n
write(*,*)(MAT(i,j),j=1,2)
end do
write(*,*)'-------------------------------------'

allocate(X(n))
allocate(Y(n))
allocate(Ya(n))
allocate(err(n))

!asignando valores de mat
do i=1,n
X(i)=MAT(i,1)
Y(i)=MAT(i,2)
end do

!inicializando sumatorias
sx=0
sy=0
sxy=0
sx2=0
sy2=0

do i=1,n
	sx=sx+X(i)
	sy=sy+Y(i)
	sx2=sx2+X(i)**2
	sy2=sy2+Y(i)**2
	sxy=sxy+X(i)*Y(i)
end do

a=(n*sxy-sx*sy)/(n*sx2-sx**2)

b=(sy-a*sx)/n

r=(n*sxy-sx*sy)/sqrt((n*sx2-sx**2)*(n*sy2-sy**2))

!ajuste de datos
do i=1,n
Ya(i)=a*X(i)+b
err(i)=abs(Ya(i)-Y(i))
end do

se=0
do i=1,n
se=se+err(i)
end do
erm=se/n

Write(*,*) 'El valor de a (Pendiente) es:'
write(*,*) a
Write(*,*) 'El valor de b es:'
write(*,*) b
Write(*,*) 'La correlacion es:'
write(*,*) r
write(*,*)'-----------------------------------------------------------------'

allocate(TAB(n,4))

do i=1,n
	TAB(i,1)=X(i)
	TAB(i,2)=Y(i)
	TAB(i,3)=Ya(i)
	TAB(i,4)=err(i)
end do
write(*,*)'Ajuste de Datos'
write(*,*)'---------- X ------------- Y -------- Yajustado ------ Error ----'
do i=1,n
write(*,*) (TAB(i,j),j=1,4)
end do
write(*,*)'-----------------------------------------------------------------'
write(*,*) 'El error maximo es:'
write(*,*) maxval(err)
write(*,*) 'El error medio es:'
write(*,*) erm

!almacenamiento de datos en un archivo txt
open(13,file="res_interpolation.txt",action='write',status='replace')
write(13,*)'---------- X ------------- Y -------- Yajustado ------ Error ----'
do i=1,n
write(13,*) (TAB(i,j),j=1,4)
end do
close(13)

end program
