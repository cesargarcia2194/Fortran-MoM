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
		integer, numero_nodo
		real(KIND=8), DIMENSION(3)::rp !posilbe vector de ubicacion de nodo
		 
	END TYPE nodo
	TYPE cara
    	integer:: id, indice_cara
    	integer, DIMENSION(3):: indice_nodos, indice_ladoscomunes
    	! los que agregue
    	real,  DIMENSION(3) :: normal, baricentro ! normal: vector normal
    	! ladoscomunes: vectores directores de los lados comunes
    	! normal_ladoscomunes: vecor normal a los lados comunes
    	real, DIMENSION(3,3) :: ladoscomunes, normal_ladoscomunes
    	real :: area, pd !! no se que era pd, pero es dot(normal,p1)
 	END TYPE cara

	TYPE lineas 
		integer:: id
		integer, DIMENSION(2):: indice_nodos
	END TYPE lineas

  	TYPE solidos
  	  integer:: id
  	  integer, DIMENSION(4):: indice_nodos
  	  INTEGER, DIMENSION(4) :: indices_cara_solido
  	END TYPE solidos


	TYPE cara_solido
     	INTEGER, DIMENSION(3) :: indices_nodo
     	INTEGER, DIMENSION(2) :: indices_solido, indices_bnodo
     	INTEGER :: indece_cara ! -1 if not a boundary
     	REAL (KIND=8) :: area
  	END TYPE cara_solido

	TYPE contenido_mesh
		type(nodo), DIMENSION(:),allocatable :: nodos 
		type(cara), DIMENSION(:), allocatable ::caras
		type(lineas), DIMENSION(:), allocatable ::lineas
    	type(solidos), DIMENSION(:), allocatable ::solidos
		TYPE(cara_solido), DIMENSION(:), ALLOCATABLE :: caras_solidos
		integer ::nnodos
		integer ::ncaras
		integer ::nlineas
		integer ::nsolidos
		integer:: ncaras_solido 
	END TYPE contenido_mesh

Contains
	!Si n>4 regresa 1, sino n
  	FUNCTION indexrot4(n) RESULT(res)
    	INTEGER, INTENT(IN) :: n
    	INTEGER :: res
	
    	res = n
	
    	IF(res>4) THEN
    	   res = res - ((res-1)/4)*4
    	END IF
END FUNCTION indexrot4
	!Si los valores triplet 1 estan en triplet 2 devuel TRUE (no importa el orden)
	FUNCTION cmp_triplets(triplet1, triplet2) RESULT(res)
	    INTEGER, DIMENSION(3), INTENT(IN) :: triplet1, triplet2
	    LOGICAL :: res
	    ! Revisar esta declaracion y como es la matriz de verdad
	    INTEGER, PARAMETER, DIMENSION(3,6) :: indices = reshape((/1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1/),(/3,6/))
	    !1 1 2 2 3 3
	    !2 3 1 3 1 2
	    !3 2 3 1 2 1
	    INTEGER :: i
	
	    res = .FALSE.
	
	    DO i=1,6
	       IF(triplet1(1)==triplet2(indices(1,i)) .AND.&
	            triplet1(2)==triplet2(indices(2,i)) .AND.&
	            triplet1(3)==triplet2(indices(3,i))) THEN
	          res = .TRUE.
	          RETURN
	       END IF
	    END DO
END FUNCTION cmp_triplets

	FUNCTION cargar_mesh(archivoMesh) RESULT(mesh)
		character(LEN=*),INTENT(IN):: archivoMesh
		character(LEN=256):: linea
		type(contenido_mesh)::mesh
		
		integer :: fid = 10, iovar, nnodos, numero_nodo, n,&
         nelementos, numero_elemento, tipo_elemento, ncaras, nlineas, nsolidos, etiqueta
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
    				mesh%nodos(n)%numero_nodo=numero_nodo
    			end do
    			mesh%nnodos=nnodos
    			READ(fid, *, IOSTAT=iovar) linea
          		if(linea/='$EndNodes') then
             		write(*,*) 'No se encuentra especificado $EndNodes'
             		stop
             	end if	
          		
          	elseif (linea =='$Elements') then
          		mesh%ncaras= COUNT(tipos_elementos(:,1)==2)
          		mesh%nlineas= COUNT(tipos_elementos(:,1)==1)
          		mesh%nsolidos= COUNT(tipos_elementos(:,1)==4)

          		allocate(mesh%caras(1:mesh%ncaras))

          		if (mesh%nlineas/=0) then 
          			allocate(mesh%lineas(1:mesh%nlineas))	
          		end if
          		if (mesh%nsolidos/=0) then
          			allocate(mesh%solidos(1:mesh%nsolidos))
          		end if
          		READ(fid,*,IOSTAT=iovar) nelementos

          		ncaras=0
          		nlineas=0
          		nsolidos=0
          		
          		n=0
          		do n=1,nelementos
          			etiqueta=tipos_elementos(n,2)
          			if(tipos_elementos(n,1)==2) then
          				ncaras=ncaras +1
          				READ(fid,*,IOSTAT=iovar) datos_elementos(1:(etiqueta+ 6))
          				mesh%caras(ncaras)%id=datos_elementos(4)
          				mesh%caras(ncaras)%indice_nodos(1:3)=datos_elementos((etiqueta+ 4):(etiqueta+ 6))
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



	function prod_cruz3(v1,v2) RESULT(v1xv2)
		real(KIND=8),DIMENSION(3),INTENT(IN)::v1
		real(KIND=8),DIMENSION(3),INTENT(IN)::v2
		real(KIND=8),DIMENSION(3)::v1xv2

		v1xv2(1)=v1(2)*v2(3)-v2(2)*v1(3)
		v1xv2(2)=v2(1)*v1(3)-v1(1)*v2(3)
		v1xv2(3)=v1(1)*v2(2)-v2(1)-v1(2)
	end function prod_cruz3

	function norma_vector3(v)RESULT(res)
		real(KIND=8),DIMENSION(3),INTENT(IN)::v
		real(KIND=8)::res
		res=SQRT(v(1)**2 +v(2)**2 +v(3)**2)
	end function norma_vector3

	function hace_vertor3(v1,v2)RESULT(v)
		real(KIND=8),DIMENSION(3),INTENT(IN)::v1
		real(KIND=8),DIMENSION(3),INTENT(IN)::v2
		real(KIND=8),DIMENSION(3)::v
		v(1)=v2(1)-v1(1)
		v(2)=v2(2)-v1(2)
		v(3)=v2(3)-v1(3)
	end function hace_vertor3

	SUBROUTINE hacer_caras(mesh)


		type(contenido_mesh),INTENT(INOUT):: mesh
		integer:: n

		real(KIND=8),DIMENSION(3)::e1, e2, e3, normal, normalu


		DO n=1,mesh%ncaras

			e1=hace_vertor3(mesh%nodos%(mesh%caras(n)%indice_nodos(1))%rp, mesh%nodos%(mesh%caras(n)%indice_nodos(2))%rp)
			e2=hace_vertor3(mesh%nodos%(mesh%caras(n)%indice_nodos(2))%rp, mesh%nodos%(mesh%caras(n)%indice_nodos(3))%rp)
			e3=hace_vertor3(mesh%nodos%(mesh%caras(n)%indice_nodos(3))%rp, mesh%nodos%(mesh%caras(n)%indice_nodos(1))%rp)

			normal = prod_cruz3(e1,e2)
			area = norma_vector3(normal)

			e1 = e1/norma_vector3(e1)
			e2 = e2/norma_vector3(e2)
			e3 = e3/norma_vector3(e3) 

			normalu = normal/norma_vector3(normal)	

			mesh%caras(n)%normal(:)=normalu(:)

			mesh%caras(n)%ladoscomunes(1,:)=e1(:)
			mesh%caras(n)%ladoscomunes(2,:)=e2(:)
			mesh%caras(n)%ladoscomunes(3,:)=e3(:)

			mesh%caras(n)%area=area

		END DO
END SUBROUTINE hacer_caras


!--------------------------------------------------------------------
SUBROUTINE hacer_caras_solidos(mesh)
    TYPE(contenido_mesh), INTENT(INOUT) :: mesh
    INTEGER :: n, m, l, ccara, ncaras, nbnds
    TYPE(cara_solido), DIMENSION(:), ALLOCATABLE :: tmpcaras
    LOGICAL:: si_cara
    INTEGER, DIMENSION(3):: triplet

    ! Allocate face reservoir with upper bound size.
    ncaras = SIZE(mesh%solidos)*4
    ALLOCATE(tmpcaras(1:ncaras))

    ! Declare node indices with undefined values.
    DO n=1,ncaras
       tmpcaras(n)%indices_nodo(:) = -1
       tmpcaras(n)%indices_bnodo(:) = -1
       tmpcaras(n)%indices_solido(:) = -1
    END DO

    ccara = 0

    ! Create unique faces.
    DO n=1,SIZE(mesh%solidos)
       DO m=1,4

          ! This determined the face node indexing.
          triplet = (/mesh%solidos(n)%indice_nodos(m),&
               mesh%solidos(n)%indice_nodos(indexrot4(m+1)),&
               mesh%solidos(n)%indice_nodos(indexrot4(m+2))/)

          ! Check if this face already exists in the list.
          si_cara = .FALSE.
          DO l=1,ccara
             IF(cmp_triplets(tmpcaras(l)%indices_nodo, triplet)) THEN
                si_cara = .TRUE.

                ! Add this face index to edge's face list and the second one.
                ! If an edge is shared by more than two faces, only two connections
                ! are recorded.
                tmpcaras(l)%indices_solido(2) = n
                tmpcaras(l)%indices_bnodo(2) = mesh%solidos(n)%indice_nodos(indexrot4(m+3))

                ! Add this face index to tetrahedra's face list.
                mesh%solidos(n)%indices_cara_solido(m) = l

                EXIT
             END IF
          END DO

          IF(si_cara.eqv..FALSE.)THEN
             ! Add new face.
             ccara = ccara + 1
             tmpcaras(ccara)%indices_nodo = triplet

             ! Add this face index to edge's face list as the first one.
             tmpcaras(ccara)%indices_solido(1) = n
             tmpcaras(ccara)%indices_bnodo(1) = mesh%solidos(n)%indice_nodos(indexrot4(m+3))

             ! Add this face index to tetrahedra's face list.
             mesh%solidos(n)%indices_cara_solido(m) = ccara
          END IF
       END DO
    END DO

    ! Trim face arrays.
    mesh%ncaras_solido = ccara
    ALLOCATE(mesh%caras_solidos(1:ccara))
    DO n=1,ccara
       mesh%caras_solidos(n)%indices_nodo(:) = tmpcaras(n)%indices_nodo(:)
       mesh%caras_solidos(n)%indices_bnodo(:) = tmpcaras(n)%indices_bnodo(:)
       mesh%caras_solidos(n)%indices_solido(:) = tmpcaras(n)%indices_solido(:)
       mesh%caras_solidos(n)%indece_cara = -1
    END DO

    ! Deallocate temporary arrays.
    DEALLOCATE(tmpcaras)

    ! Connect boundary solid faces to surface faces.
    DO n=1,mesh%ncaras_solido
       IF(mesh%caras_solidos(n)%indices_solido(1)==-1 .OR.&
            mesh%caras_solidos(n)%indices_solido(2)==-1) THEN
          DO m=1,mesh%ncaras
             IF(cmp_triplets(mesh%caras_solidos(n)%indices_nodo, mesh%caras(m)%indice_nodos)) THEN
                mesh%caras_solidos(n)%indece_cara = m
                EXIT
             END IF
          END DO

          IF(mesh%caras_solidos(n)%indece_cara==-1) THEN
             WRITE(*,*) 'Could not connect solid boundary face to surface face!'
             STOP
          END IF
       END IF
    END DO

    WRITE(*,'(A,I0,:,A)') ' - Created ', ccara, ' unique solid faces'

END SUBROUTINE hacer_caras_solidos
END MODULE mesh
