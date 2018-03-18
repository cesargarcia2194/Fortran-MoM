program main
  use mesh
  implicit none
  type(contenido_mesh) :: cont_mesh
  cont_mesh = cargar_mesh('sphere.msh')
  call hacer_lados(cont_mesh)


  
end program main