program main
  use mesh
  implicit none
  type(contenido_mesh) :: cont_mesh
  cont_mesh = cargar_mesh('sphere.msh')
  print*, 'Hola'
  print*, 'Numero de nodos', cont_mesh%nnodos
  print*, 'Numero de parches triangulares', cont_mesh%nparchest
  print*, 'numero de solidos', cont_mesh%nsolidos
end program main