program main
  use mesh
  implicit none
  type(contenido_mesh) :: cont_mesh
  cont_mesh = cargar_mesh('cubito.msh')
  print*, 'Hola'
end program main