program main
  use mesh
  implicit none
  type(contenido_mesh) :: cont_mesh
  cont_mesh = cargar_mesh('sphere.msh')

  call hacer_caras_solidos(cont_mesh)

  print*, 'Hola'
  print*, 'Numero de nodos', cont_mesh%nnodos
 
  print*, 'numero de solidos', cont_mesh%nsolidos
  
end program main