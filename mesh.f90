! MODULE: mesh
! AUTHOR: Cesar Garcia / Victor Sanchez
! DESCRIPTION:
! Loading and manipulating meshes consisting of triangles and tetrahedra.
! Supports the msh-format exported by Gmsh and the neutral mesh format of Netgen.
! Contains also functions for splitting a mesh into submeshes and various mesh manipulation routines

MODULE mesh
	IMPLICIT NONE

	integer k=2
	TYPE nodo
		real(KIND=dp), DIMENSION(3)::rp !posilbe vector de ubicacion de nodo
		integer:: parent_index 
	END TYPE nodo


	FUNCTION cargar_mesh(archivoMesh) RESULT(mesh)
		character(LEN=256),INTENT(IN):: archivoMesh
		character(LEN=256):: linea
		type(contenedor_mesh)::mesh
		
		integer :: fid = 10, iovar, nnodos, numero_nodo, n,&
         nelementos, numero_elemento, tipo_elemento, ccaras, clineas, csolidos, ntags
        integer, DIMENSION(10) :: element_data !CAMBIA NOMBRE

        integer,DIMENSION(:,:)::tipos_elementos

        real (KIND=dp), DIMENSION(3) :: np
    	character (LEN=3) :: mshVersion

    	OPEN(fid,file=TRIM(archivoMesh), action='READ', IOSTAT=iovar)

    	if (iovar>0) then
    		write(*,*)'No se pudo abrir el archivo', TRIM(archivoMesh),'!'
    		stop
    	end if

    	iovar=0
    	!Este comentario es para probar git hub
    	!Este es un segundo comentario
    	do while (iovar==0)
    		read(fid,*,IOSTAT=iovar) linea
    		if(linea=='$Elements') then
    			read(fid,*,IOSTAT=iovar) nelementos
    			allocate(tipos_elementos(nelementos,1:2))
    			do n=1:nelementos
    				read(fid,*,IOSTAT=iovar) numero_elemento, tipos_elementos(n,1), tipos_elementos(n,2)
    			end do
    			read(fid,*,IOSTAT=iovar) linea
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
    		read(fid,*,IOSTAT=iovar) linea
    		if(linea=='$MeshFormat')then
    			read(fid,'(A3)',IOSTAT=iovar) mshVersion
    			if(mshVersion /= '2.1' .and. mshVersion /= '2.2') then
    				write(*,*) 'Version de mesh no soportada'
    				CLOSE(fid)
    				stop
    			end if
    			write(*,'(A,A3)') 'Mesh version: ', mshVersion
    			read(fid,*,IOSTAT=iovar) linea
    			if(linea /= '$EndMeshFormat') then
    				write(*,*) 'No se encuentra especificado $EndMeshFormat'
    				stop
    			end if
    		else if(linea =='$Nodes')
    			read(fid,*,IOSTAT=iovar) nnodos
    			allocate(mesh%nodos(nnodos)) !crear nodos
    			do n=1: nnodos
    				read(fid,*,IOSTAT=iovar) numero_nodo, np(1:3)
    				mesh%nodos(n)%rp=np
    			end do
    			mesh%nnodos=nnodos
    			read(fid, *, IOSTAT=iovar) linea
          		if(linea/='$EndNodes') then
             		write(*,*) 'No se encuentra especificado $EndNodes'
             		stop
          		end if
          	else if (linea =='$Elements')

    		end if 		
    	end do

	END FUNCTION cargar_mesh
END MODULE mesh