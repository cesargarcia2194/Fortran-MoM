program prueba1
  use mesh
  implicit none
  type(contenido_mesh) :: cont_mesh
  integer :: n, i

  cont_mesh = cargar_mesh('sphere.msh')
  

  print *, 'Primero Constuimos el mesh'

  !call export_mesh('nuevo.msh',mes
  call hacer_lados(cont_mesh)

  CALL calcular_datos_base(cont_mesh)
  print *, 'Primera prueba, como el test'

  do n=1,SIZE(cont_mesh%caras)
    print *, 'la n'
    write(*,'(F9.3)') cont_mesh%caras(n)%normal, cont_mesh%caras(n)%baricentro

    print *, 'la s'
    print *, 'area ', cont_mesh%caras(n)%area
    do i=1,3
      print *, cont_mesh%caras(n)%vdirec(i,1:3)
    end do

  enddo

  do i=1,cont_mesh%nlados
      print *, cont_mesh%lados(i)%indices_cara(:)
  end do


  
end program prueba1