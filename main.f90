program main
  use mesh
  use method
  implicit none
  type(contenido_mesh) :: cont_mesh
  integer::n
  cont_mesh = cargar_mesh('sphere.msh')
  call hacer_lados(cont_mesh)
  call  calcular_datos_base(cont_mesh)
  do n=1, cont_mesh%nlados
    print*, cont_mesh%centros(:,n)
  end do
  
  call printconsole('Hola')
  print *, cont_mesh%lados

  call inicializar_MLFMA(cont_mesh%centro, cont_mesh%nlados, cont_mesh%lados(:)%indices_nodo, &
    ont_mesh%lados(:)%bnodes, cont_mesh%lados(:)%longitud, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p,&
          num_t, t_normal)


end program main