subroutine test_legacy_cv_faces

  use legacy_cv_numbering
  use legacy_cv_shape_functions

  implicit none

  integer :: optelm, eletyp, seletyp, ngi, sngi, nloc, snloc, svngi
  integer, dimension(:,:), allocatable :: neiloc
  real, dimension(:,:), allocatable :: svn, svnlx, svnly, m
  real, dimension(:), allocatable :: svweigh
  logical :: d3, dcyl, redquad

  integer i

  redquad = .false.
  optelm = 1 ! linear triangles

  call setelm(eletyp, ngi, nloc, optelm, seletyp, sngi, snloc, svngi, d3, dcyl, redquad)

  allocate(neiloc(nloc, svngi), svn(nloc, svngi), svnlx(nloc, svngi), svnly(nloc, svngi),&
           svweigh(svngi), m(nloc, ngi))
  call shapesv(eletyp, neiloc, ngi, nloc, svngi, svn, svnlx, svnly, svweigh, m)

  write(0,*) 'p1p1'
  do i = 1, size(svn,2)
    write(0,*) 'gi = ', i
    write(0,*) 'cv_p1p1%n(:,gi) = ', svn(:,i)
    write(0,*) 'dim = 1'
    write(0,*) 'cv_p1p1%dn(:,gi,dim) = ', svnlx(:,i)
    write(0,*) 'cv_p1p1%quadrature%weight(gi) = ', svweigh(i)
  end do

  deallocate(neiloc, svn, svnlx, svnly, svweigh, m)

  optelm = 5 ! linear tets

  call setelm(eletyp, ngi, nloc, optelm, seletyp, sngi, snloc, svngi, d3, dcyl, redquad)

  allocate(neiloc(nloc, svngi), svn(nloc, svngi), svnlx(nloc, svngi), svnly(nloc, svngi),&
           svweigh(svngi), m(nloc, ngi))
  call shapesv(eletyp, neiloc, ngi, nloc, svngi, svn, svnlx, svnly, svweigh, m)

  write(0,*) 'p1p1_tet'
  do i = 1, size(svn,2)
    write(0,*) 'gi = ', i
    write(0,*) 'cv_p1p1_tet%n(:,gi) = ', svn(:,i)
    write(0,*) 'dim = 1'
    write(0,*) 'cv_p1p1_tet%dn(:,gi,dim) = ', svnlx(:,i)
    write(0,*) 'dim = 2'
    write(0,*) 'cv_p1p1_tet%dn(:,gi,dim) = ', svnly(:,i)
    write(0,*) 'cv_p1p1_tet%quadrature%weight(gi) = ', svweigh(i)
  end do

  optelm = 2 ! quadratic triangles

  call setelm(eletyp, ngi, nloc, optelm, seletyp, sngi, snloc, svngi, d3, dcyl, redquad)

  allocate(neiloc(nloc, svngi), svn(nloc, svngi), svnlx(nloc, svngi), svnly(nloc, svngi),&
           svweigh(svngi), m(nloc, ngi))
  call shapesv(eletyp, neiloc, ngi, nloc, svngi, svn, svnlx, svnly, svweigh, m)

  write(0,*) 'p2p2'
  do i = 1, size(svn,2)
    write(0,*) 'gi = ', i
    write(0,*) 'cv_p2p2%n(:,gi) = ', svn(:,i)
    write(0,*) 'dim = 1'
    write(0,*) 'cv_p2p2%dn(:,gi,dim) = ', svnlx(:,i)
    write(0,*) 'cv_p2p2%quadrature%weight(gi) = ', svweigh(i)
  end do

  deallocate(neiloc, svn, svnlx, svnly, svweigh, m)

  write(0,*) 'ending'

end subroutine test_legacy_cv_faces