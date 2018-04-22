program prueba3
  use mesh
  ! prueba leyendo un archivo mesh .msh
  TYPE(mesh_container) :: meshs
  Character(len=40) :: file
  REAL (KIND=dp) :: scale

  print *, 'Nombre del archivo'
  read *, file
  meshs = load_mesh(file)

  print *, 'Ahora vamos a esportarlo, usando build mesh'

  !call export_mesh('nuevo.msh',meshs)
  ! meshs es la variables con los mesh
  scale = 1.0
  call build_mesh(meshs,scale)

  print *, 'Fino, ojala aiga servido'

  do n=1,SIZE(meshs%faces)
    print *, 'la n'
    write(*,'(F9.3)') meshs%faces(n)%n, meshs%faces(n)%cp
    print *, 'area', meshs%faces(n)%area
    print *, 'la s'

    do i=1,3
      print *, meshs%faces(n)%s(i,1:3)
    end do

  enddo

  do i=1,meshs%nedges
      print *, meshs%edges(i)%face_indices(:)
    end do


end program prueba3