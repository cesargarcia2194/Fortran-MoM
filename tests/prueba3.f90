program prueba3
  use mesh
  use mlfma
  use modRWG
  use utils

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
  integer ( kind = il ) :: num_p, num_t, num_e
  real (kind = dp) :: frequency


    complex (kind = dp), dimension(3) :: pol_onda
    complex (kind = dp) :: ctte_onda

    real (kind = dp) :: error_solver
    real (kind = dp), dimension(3) :: dir_onda
    
    character (len = 256) :: comm, val
    character (len = 3), parameter :: strdefault= '/%\'
    character (len = 256) :: autocommand

    integer (kind = il) :: ent, numtest

    systemDelimiter = '/'
    !crea la carpeta para asegurarnos de que exista
    !call system('mkdir -p --parents' // '..' // trim(systemDelimiter) // 'Results')

    !Lectura de argumentos de la linea de comandos
    conf%argument1(1:len(conf%argument1)) = ' '
    conf%argument2(1:len(conf%argument2)) = ' '
    call getarg(1, conf%argument1)
    call getarg(2, conf%argument2)
    if (streq(conf%argument2, 'gui')) then
        conf%optionGUI = .true.
    else
        conf%optionGUI = .false.
    end if
    !

    !Inicializacion

    ctte_onda = 1.
    lambda = 1.
    call setlambda(lambda)
    !skipkD = .false.

    conf%efieOldMode = 0
    conf%test_name = 'lasttest'!strdefault
    conf%msj = .false.
    conf%numThetaRCS = 1
    conf%numPhiRCS = 360
    conf%errorSolverCGS = 0.001
    conf%Lforz = -1
    conf%dirTheta = 90.
    conf%dirPhi = 180.
    conf%rcsMonoFMin = 10
    conf%rcsMonoFMax = 275
    conf%rcsMonoSamples = 75
    conf%usar_precond = 'ILUT'
    conf%prec_activo = .true.
    conf%fill_in_prc = 5.
    conf%drop_tolerance = 0.001
    conf%analisis_rcs = 'MONO'
    conf%n_pasos_rcs_mono = 180
    conf%tipo_polarizac = 'VV'
    conf%variable_espec = 'FREC'
    conf%file_name = strdefault
    conf%file_script = strdefault
    conf%criterio_ladomax = 0.5_dp
    conf%CFIEfactor = 1._dp
    conf%force_ladomax = .false.
    conf%mult_maxprecision = 6
    conf%mult_minprecision = 0
    conf%isurf_phase = .false.
    conf%save_incrementals = .true.
    conf%rcsmonostepi = 1
    conf%rcsmonostepf = conf%n_pasos_rcs_mono

    alphaCFIE = conf%CFIEfactor      

    !conf%polVH = 1
    error_solver = conf%errorSolverCGS
    dir_onda(1) = sin(conf%dirTheta*pi/180.)*cos(conf%dirPhi*pi/180.)
    dir_onda(2) = sin(conf%dirTheta*pi/180.)*sin(conf%dirPhi*pi/180.)
    dir_onda(3) = cos(conf%dirTheta*pi/180.)
    !call ajustar_polariz(dir_onda,pol_onda)

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
  allocate(t_normal(3,meshs%nfaces))
  do i=1,meshs%nnodes
    p_coord(:,i) = meshs%nodes(i)%p
    !print*, meshs%edges(i)%bnode_indices
  enddo
  allocate(t_p(3,size(meshs%faces)))
  do i=1,meshs%nfaces
    t_area = meshs%faces(i)%area
    t_baric(:,i) = meshs%faces(i)%cp
    t_p(:,i) = meshs%faces(i)%node_indices
    t_normal(:,i) = meshs%faces(i)%n
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