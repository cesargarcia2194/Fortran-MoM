program prueba3
  use mesh
  use mlfma
  use modRWG

  implicit none
  ! prueba leyendo un archivo mesh .msh
  TYPE(mesh_container) :: meshs
  Character(len=40) :: file
  REAL (KIND=dp) :: scale
  integer :: i
  ! variables para los modulos
  integer ( kind = il ), allocatable, dimension(:,:) :: t_p
  integer ( kind = il ), allocatable, dimension(:,:) :: e_p
  integer ( kind = il ), allocatable, dimension(:,:) :: e_t
  integer ( kind = il ), allocatable, dimension(:,:) :: e_po
  real ( kind = dp ), allocatable, dimension(:,:) :: p_coord
  real (kind = dp ), allocatable, dimension(:) :: e_long
  real ( kind = dp ), allocatable, dimension(:) :: t_area
  real ( kind= dp ), allocatable, dimension(:,:) :: t_normal
  real ( kind = dp ), allocatable, dimension(:,:) :: e_centro
  real ( kind = dp ), allocatable, dimension(:,:) :: t_baric
  real ( kind = dp ), allocatable, dimension(:,:,:) :: t_baric_sub
  integer ( kind = il ) :: num_p, num_t


  print *, 'Nombre del archivo'
  read *, file
  meshs = load_mesh(file)

  print *, 'Ahora vamos a esportarlo, usando build mesh'

  !call export_mesh('nuevo.msh',meshs)
  ! meshs es la variables con los mesh
  scale = 1.0
  call build_mesh(meshs,scale)

  allocate(t_area(meshs%nfaces))
  allocate(t_baric(3,meshs%nfaces))
  allocate(p_coord(3,meshs%nnodes))
  do i=1,meshs%nnodes
    p_coord(:,i) = meshs%nodes(i)%p
    !print*, meshs%edges(i)%bnode_indices
  enddo
  allocate(t_p(3,size(meshs%faces)))
  do i=1,meshs%nfaces
    t_area = meshs%faces(i)%area
    t_baric(:,i) = meshs%faces(i)%cp
    t_p(:,i) = meshs%faces(i)%node_indices
  enddo
  allocate(e_p(2,size(meshs%edges)))
  allocate(e_long(meshs%nedges))
  allocate(e_po(2,size(meshs%edges)))
  allocate(e_t(2,size(meshs%edges)))
  do i=1,meshs%nedges
    e_long = meshs%edges(i)%length
    e_po(:,i) = meshs%edges(i)%bnode_indices
    e_p(:,i) = meshs%edges(i)%node_indices
    e_t(:,i) = meshs%edges(i)%face_indices
  enddo
  !print*, shape(t_p)
  num_p = meshs%nnodes
  num_t = meshs%nfaces
  
  call hallar_geometria_rwg(p_coord,t_p,meshs%nnodes,meshs%nfaces,t_area,t_baric,t_baric_sub)

  call inicializar_MLFMA(meshs%centers, meshs%nedges, e_t, &
    e_po, e_long, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p,&
          num_t, t_normal)


end program prueba3