! MODULE: mesh
! AUTHOR: Cesar Garcia / Victor Sanchez
! DESCRIPTION:
! Loading and manipulating meshes consisting of triangles and tetrahedra.
! Supports the msh-format exported by Gmsh and the neutral mesh format of Netgen.
! Contains also functions for splitting a mesh into submeshes and various mesh manipulation routines
!
MODULE mesh
  use constantes
  use alineal
  IMPLICIT NONE

  
  TYPE nodo
    integer::numero_nodo
    real(KIND=dp), DIMENSION(3)::rp !vector pocision de nodo
     
  END TYPE nodo
  TYPE cara
      integer:: id, indice_cara
      integer, DIMENSION(3):: indice_nodos, indice_ladoscomunes
      ! los que agregue
      real(KIND=dp),  DIMENSION(3) :: normal, baricentro ! normal: vector normal
      ! ladoscomunes: vectores directores de los lados comunes
      ! normal_ladoscomunes: vecor normal a los lados comunes
      real(KIND=dp), DIMENSION(3,3) :: vdirec, normal_lados
      real :: area, pd !! no se que era pd, pero es dot(normal,p1)
   END TYPE cara

  TYPE lineas 
    integer:: id
    integer, DIMENSION(2):: indice_nodos
  END TYPE lineas

  !! Lados comunes entre las caras
  TYPE lado
     integer,DIMENSION(2)::indices_nodo, indices_cara
    real::longitud

   !INTEGER :: bnd
    !INTEGER :: parent_index
    !INTEGER :: couple_index
    !INTEGER, DIMENSION(:,:), ALLOCATABLE :: child_indices !
  END TYPE lado

  TYPE solido
    integer:: id
    integer, DIMENSION(4):: indices_nodo
    INTEGER, DIMENSION(4) :: indices_cara_solido
    real(KIND=dp)::volumen
  END TYPE solido


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
    type(lado), DIMENSION(:), ALLOCATABLE :: lados
    type(solido), DIMENSION(:), allocatable ::solidos
    TYPE(cara_solido), DIMENSION(:), ALLOCATABLE :: caras_solidos
    real(kind=dp),DIMENSION(:,:),ALLOCATABLE::centros
    integer ::nnodos
    integer ::ncaras
    integer ::nlados
    integer ::nlineas
    integer ::nsolidos
    integer:: ncaras_solido
      
  END TYPE contenido_mesh

Contains

  FUNCTION cargar_mesh(archivoMesh) RESULT(mesh)
    character(LEN=*),INTENT(IN):: archivoMesh
    character(LEN=256):: linea
    type(contenido_mesh)::mesh
    
    integer :: fid = 10, iovar, nnodos, numero_nodo, n,&
         nelementos, numero_elemento, tipo_elemento, ncaras, nlineas, nsolidos, etiqueta
        integer, DIMENSION(10) :: datos_elementos !CAMBIA NOMBRE

        integer,DIMENSION(:,:),allocatable::tipos_elementos

        real (KIND=dp), DIMENSION(3) :: np
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
                  mesh%solidos(nsolidos)%indices_nodo(1:4)=datos_elementos((etiqueta +4):(etiqueta+ 7))
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

SUBROUTINE hacer_lados(mesh) 
  type(contenido_mesh), INTENT(INOUT) :: mesh
  integer::n, m, l, nlados, clados 
  type(lado), DIMENSION(:), allocatable:: temp_lados 
  logical:: si_lado
  integer, DIMENSION(2):: par

  nlados=SIZE(mesh%caras)*3
  allocate(temp_lados(1:nlados))

  do n=1,nlados
    temp_lados(n)%indices_cara(:)=-1
    temp_lados(n)%indices_nodo(:)=-1
  end do
  clados=0

  do n=1,SIZE(mesh%caras)
    do m=1,3
      par = (/mesh%caras(n)%indice_nodos(m), mesh%caras(n)%indice_nodos(indexrot3(m+1))/)
      si_lado=.FALSE.
      do l=1,clados
        if(comp_par(temp_lados(l)%indices_nodo, par)) then
          si_lado=.TRUE.
          temp_lados(l)%indices_cara(2)=n
          mesh%caras(n)%indice_ladoscomunes(m)=clados

          EXIT
        end if
      end do

      if(si_lado .eqv. .FALSE.) then
        clados = clados + 1

        temp_lados(clados)%indices_nodo = par
        temp_lados(clados)%indices_cara(1)=n
        
        mesh%caras(n)%indice_ladoscomunes(m)=clados
        
      end if
    end do
  end do
  mesh%nlados=clados
  allocate(mesh%lados(1:clados))
  do n=1,clados
    mesh%lados(n)%indices_nodo(:) = temp_lados(n)%indices_nodo(:)
    mesh%lados(n)%indices_cara(:) = temp_lados(n)%indices_cara(:)
  end do
END SUBROUTINE hacer_lados


!--------------------------------------------------------------------
SUBROUTINE hacer_caras(mesh)
    TYPE(contenido_mesh), INTENT(INOUT) :: mesh
    INTEGER :: n, m, l, ccara, ncaras, nbnds
    TYPE(cara_solido), DIMENSION(:), ALLOCATABLE :: temp_caras
    LOGICAL:: si_cara
    INTEGER, DIMENSION(3):: triplet

    ! Allocate face reservoir with upper bound size.
    ncaras = SIZE(mesh%solidos)*4
    ALLOCATE(temp_caras(1:ncaras))

    ! Declare node indices with undefined values.
    DO n=1,ncaras
       temp_caras(n)%indices_nodo(:) = -1
       temp_caras(n)%indices_bnodo(:) = -1
       temp_caras(n)%indices_solido(:) = -1
    END DO

    ccara = 0

    ! Create unique faces.
    DO n=1,SIZE(mesh%solidos)
       DO m=1,4

          ! This determined the face node indexing.
          triplet = (/mesh%solidos(n)%indices_nodo(m),&
               mesh%solidos(n)%indices_nodo(indexrot4(m+1)),&
               mesh%solidos(n)%indices_nodo(indexrot4(m+2))/)

          ! Check if this face already exists in the list.
          si_cara = .FALSE.
          DO l=1,ccara
             IF(comp_trio(temp_caras(l)%indices_nodo, triplet)) THEN
                si_cara = .TRUE.

                ! Add this face index to edge's face list and the second one.
                ! If an edge is shared by more than two faces, only two connections
                ! are recorded.
                temp_caras(l)%indices_solido(2) = n
                temp_caras(l)%indices_bnodo(2) = mesh%solidos(n)%indices_nodo(indexrot4(m+3))

                ! Add this face index to tetrahedra's face list.
                mesh%solidos(n)%indices_cara_solido(m) = l

                EXIT
             END IF
          END DO

          IF(si_cara.eqv..FALSE.)THEN
             ! Add new face.
             ccara = ccara + 1
             temp_caras(ccara)%indices_nodo = triplet

             ! Add this face index to edge's face list as the first one.
             temp_caras(ccara)%indices_solido(1) = n
             temp_caras(ccara)%indices_bnodo(1) = mesh%solidos(n)%indices_nodo(indexrot4(m+3))

             ! Add this face index to tetrahedra's face list.
             mesh%solidos(n)%indices_cara_solido(m) = ccara
          END IF
       END DO
    END DO

    ! Trim face arrays.
    mesh%ncaras_solido = ccara
    ALLOCATE(mesh%caras_solidos(1:ccara))
    DO n=1,ccara
       mesh%caras_solidos(n)%indices_nodo(:) = temp_caras(n)%indices_nodo(:)
       mesh%caras_solidos(n)%indices_bnodo(:) = temp_caras(n)%indices_bnodo(:)
       mesh%caras_solidos(n)%indices_solido(:) = temp_caras(n)%indices_solido(:)
       mesh%caras_solidos(n)%indece_cara = -1
    END DO

    ! Deallocate temporary arrays.
    DEALLOCATE(temp_caras)

    ! Connect boundary solid faces to surface faces.
    DO n=1,mesh%ncaras_solido
       IF(mesh%caras_solidos(n)%indices_solido(1)==-1 .OR.&
            mesh%caras_solidos(n)%indices_solido(2)==-1) THEN
          DO m=1,mesh%ncaras
             IF(comp_trio(mesh%caras_solidos(n)%indices_nodo, mesh%caras(m)%indice_nodos)) THEN
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

    !! la rutina no pasa por aqui
    WRITE(*,'(A,I0,:,A)') ' - Created ', ccara, ' unique solid faces'


    WRITE(*,*) 'Mesh data built successfully'


END SUBROUTINE hacer_caras

!>  Calcula los datos para las funciones base
SUBROUTINE calcular_datos_base(mesh)
  type(contenido_mesh), INTENT(INOUT):: mesh
  integer :: n
    REAL (KIND=dp), DIMENSION(3) :: b, p1, p2, p3, p4, p21, p32, p13

    ALLOCATE(mesh%centros(3,1:mesh%nlados))

    do n=1,mesh%ncaras
      p1=mesh%nodos(mesh%caras(n)%indice_nodos(1))%rp
      p2=mesh%nodos(mesh%caras(n)%indice_nodos(2))%rp
      p3=mesh%nodos(mesh%caras(n)%indice_nodos(3))%rp



      b = prod_cruz((p2-p1),(p3-p1))
      mesh%caras(n)%normal=b/normav(b)

      mesh%caras(n)%pd = prod_escalar(mesh%caras(n)%normal, p1)

      p21 = p2 - p1
        p32 = p3 - p2
        p13 = p1 - p3
        mesh%caras(n)%vdirec(:,1) = p21/normav(p21)
        mesh%caras(n)%vdirec(:,2) = p32/normav(p32)
        mesh%caras(n)%vdirec(:,3) = p13/normav(p13)

        mesh%caras(n)%normal_lados(:,1) = prod_cruz(mesh%caras(n)%vdirec(:,1), mesh%caras(n)%normal)
        mesh%caras(n)%normal_lados(:,2) = prod_cruz(mesh%caras(n)%vdirec(:,2), mesh%caras(n)%normal)
        mesh%caras(n)%normal_lados(:,3) = prod_cruz(mesh%caras(n)%vdirec(:,3), mesh%caras(n)%normal)

        mesh%caras(n)%area = (0.5_dp)*normav(prod_cruz(-p13,p32))  !0,5_dp convertir
        mesh%caras(n)%baricentro = (p1 + p2 + p3)/3.0_dp !convertir
    end do

    do n=1,mesh%nlados
      mesh%lados(n)%longitud = normav(mesh%nodos(mesh%lados(n)%indices_nodo(2))%rp -&
        mesh%nodos(mesh%lados(n)%indices_nodo(2))%rp)
      
      mesh%centros(:,n)=(mesh%nodos(mesh%lados(n)%indices_nodo(1))%rp+mesh%nodos(mesh%lados(n)%indices_nodo(1))%rp)/2
    end do
    do n=1,mesh%nsolidos
       p1 = mesh%nodos(mesh%solidos(n)%indices_nodo(1))%rp
       p2 = mesh%nodos(mesh%solidos(n)%indices_nodo(2))%rp
       p3 = mesh%nodos(mesh%solidos(n)%indices_nodo(3))%rp
       p4 = mesh%nodos(mesh%solidos(n)%indices_nodo(4))%rp
       mesh%solidos(n)%volumen = ABS(prod_escalar((p1-p4),prod_cruz((p2-p4),(p3-p4))))/6.0_dp
    end do

    do n=1,mesh%ncaras_solido
       p1 = mesh%nodos(mesh%caras_solidos(n)%indices_nodo(1))%rp
       p2 = mesh%nodos(mesh%caras_solidos(n)%indices_nodo(2))%rp
       p3 = mesh%nodos(mesh%caras_solidos(n)%indices_nodo(3))%rp

       mesh%caras_solidos(n)%area = 0.5_dp*normav(prod_cruz(p3-p1,p3-p2))
    end do
END SUBROUTINE calcular_datos_base


END MODULE mesh
