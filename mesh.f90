! MODULE: mesh
! AUTHOR: Cesar Garcia / Victor Sanchez
! DESCRIPTION:
! Loading and manipulating meshes consisting of triangles and tetrahedra.
! Supports the msh-format exported by Gmsh and the neutral mesh format of Netgen.
! Contains also functions for splitting a mesh into submeshes and various mesh manipulation routines
!
MODULE mesh
	IMPLICIT NONE

	
	TYPE nodo
		real(KIND=8), DIMENSION(3)::rp !posilbe vector de ubicacion de nodo
		 
	END TYPE nodo

	TYPE parchest
		integer:: id
		integer, DIMENSION(3):: indice_nodos
	END TYPE parchest

	TYPE lineas 
		integer:: id
		integer, DIMENSION(2):: indice_nodos
	END TYPE lineas

	TYPE solidos
		integer:: id
		integer, DIMENSION(4):: indice_nodos
	END TYPE solidos

	TYPE contenido_mesh
		type(nodo), DIMENSION(:),allocatable :: nodos 
		type(parchest), DIMENSION(:), allocatable ::parchest
		type(lineas), DIMENSION(:), allocatable ::lineas
		type(solidos), DIMENSION(:), allocatable ::solidos
		integer ::nnodos
		integer ::nparchest
		integer ::nlineas
		integer ::nsolidos 
	END TYPE contenido_mesh

Contains
	FUNCTION cargar_mesh(archivoMesh) RESULT(mesh)
		character(LEN=*),INTENT(IN):: archivoMesh
		character(LEN=256):: linea
		type(contenido_mesh)::mesh
		
		integer :: fid = 10, iovar, nnodos, numero_nodo, n,&
         nelementos, numero_elemento, tipo_elemento, nparchest, nlineas, nsolidos, etiqueta
        integer, DIMENSION(10) :: datos_elementos !CAMBIA NOMBRE

        integer,DIMENSION(:,:),allocatable::tipos_elementos

        real (KIND=8), DIMENSION(3) :: np
    	character (LEN=3) :: mshVersion

    	OPEN(fid,file=TRIM(archivoMesh), action='READ', IOSTAT=iovar)

    	if (iovar>0) then
    		write(*,*)'No se pudo abrir el archivo', TRIM(archivoMesh),'!'	
    		stop
    	end if

    	iovar=0
   
    	do while (iovar==0)
    		READ(fid,*,IOSTAT=iovar) linea
    		if(linea=='$Elements') then
    			READ(fid,*,IOSTAT=iovar) nelementos
    			allocate(tipos_elementos(nelementos,1:2))
    			n=0
    			do n=1,nelementos
    				READ(fid,*,IOSTAT=iovar) numero_elemento, tipos_elementos(n,1), tipos_elementos(n,2)
    			end do

    			READ(fid,*,IOSTAT=iovar) linea

    			if(linea /= '$EndElements') then
    				write(*,*) 'No se encuentra especificado $EndElements'
    			end if
    		end if
    	end do
    	CLOSE(fid)

    	!Comenzar
    	OPEN(fid, file=TRIM(archivoMesh), action='READ', IOSTAT=iovar)

    	if (iovar>0) then
    		write(*,*)'No se pudo abrir el archivo', TRIM(archivoMesh), '!'
    	end if

    	iovar=0

    	do while (iovar==0)
    		READ(fid,*,IOSTAT=iovar) linea
    		if(linea=='$MeshFormat')then
    			READ(fid,'(A3)',IOSTAT=iovar) mshVersion
    			if(mshVersion /= '2.1' .and. mshVersion /= '2.2') then
    				write(*,*) 'Version de mesh no soportada'
    				CLOSE(fid)
    				stop
    			end if
    			write(*,'(A,A3)') 'Mesh version: ', mshVersion
    			READ(fid,*,IOSTAT=iovar) linea
    			if(linea /= '$EndMeshFormat') then
    				write(*,*) 'No se encuentra especificado $EndMeshFormat'
    				stop
    			end if
    		elseif(linea =='$Nodes') then
    			READ(fid,*,IOSTAT=iovar) nnodos
    			allocate(mesh%nodos(nnodos)) !crear nodos
    			n=0
    			do n=1, nnodos
    				READ(fid,*,IOSTAT=iovar) numero_nodo, np(1:3)
    				mesh%nodos(n)%rp=np
    			end do
    			mesh%nnodos=nnodos
    			READ(fid, *, IOSTAT=iovar) linea
          		if(linea/='$EndNodes') then
             		write(*,*) 'No se encuentra especificado $EndNodes'
             		stop
             	end if	
          		
          	elseif (linea =='$Elements') then
          		mesh%nparchest= COUNT(tipos_elementos(:,1)==2)
          		mesh%nlineas= COUNT(tipos_elementos(:,1)==1)
          		mesh%nsolidos= COUNT(tipos_elementos(:,1)==4)

          		allocate(mesh%parchest(1:mesh%nparchest))

          		if (mesh%nlineas/=0) then 
          			allocate(mesh%lineas(1:mesh%nlineas))	
          		end if
          		if (mesh%nsolidos/=0) then
          			allocate(mesh%solidos(1:mesh%nsolidos))
          		end if
          		READ(fid,*,IOSTAT=iovar) nelementos

          		nparchest=0
          		nlineas=0
          		nsolidos=0
          		
          		n=0
          		do n=1,nelementos
          			etiqueta=tipos_elementos(n,2)
          			if(tipos_elementos(n,1)==2) then
          				nparchest=nparchest +1
          				READ(fid,*,IOSTAT=iovar) datos_elementos(1:(etiqueta+ 6))
          				mesh%parchest(nparchest)%id=datos_elementos(4)
          				mesh%parchest(nparchest)%indice_nodos(1:3)=datos_elementos((etiqueta+ 4):(etiqueta+ 6))
          			elseif (tipos_elementos(n,1)==1) then
          				nlineas=nlineas +1
          				READ(fid,*,IOSTAT=iovar) datos_elementos(1:(etiqueta+ 5))
          				mesh%lineas(nlineas)%id=datos_elementos(4)
          				mesh%lineas(nlineas)%indice_nodos(1:2)=datos_elementos((etiqueta+ 4):(etiqueta+ 5))
          			elseif (tipos_elementos(n,2)==4) then
          				nsolidos = nsolidos +1
          				READ(fid,*,IOSTAT=iovar) datos_elementos(1:(etiqueta+ 7))
          				mesh%solidos(nsolidos)%id=datos_elementos(4)
          				mesh%solidos(nsolidos)%indice_nodos(1:4)=datos_elementos((etiqueta +4):(etiqueta+ 7))
          			else
          				READ(fid,*,IOSTAT=iovar) datos_elementos(1)
          			end if
          		end do

          		READ(fid, *, IOSTAT=iovar) linea
          		if(linea/='$EndElements') then
            		write(*,*) 'No se encuentra especificado $EndElements'
            		stop
          		end if

    		end if
    	
    	end do
    	CLOSE(fid)
    	!CERRAR ARCHIVO
END FUNCTION cargar_mesh
END MODULE mesh
