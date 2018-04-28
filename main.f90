program main
  use mesh
  !use method
  implicit none
  type(mesh_container) :: cont_mesh
  integer::n
  cont_mesh = load_mesh_gmsh('lamina.msh')
  call build_mesh(cont_mesh)
  call  compute_basis_data(cont_mesh)
  do n=1, cont_mesh%nedges
    print*, cont_mesh%centers(:,n)
  end do
  
  !call printconsole('Hola')


  ! call inicializar_MLFMA(cont_mesh%centro, cont_mesh%nlados, cont_mesh%lados(:)%indices_nodo, &
  !  ont_mesh%lados(:)%bnodes, cont_mesh%lados(:)%longitud, t_area, t_baric, t_baric_sub, p_coord, t_p, e_p, num_p,&
  !        num_t, t_normal)


end program main