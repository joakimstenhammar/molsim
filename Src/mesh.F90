! ... 'version 6.4.7, Sep 18, 2015',

!************************************************************************
!************************************************************************
!**                                                                    **
!**  Copyright 2014 - 2014                                             **
!**                                                                    **
!**  Per Linse                                                         **
!**  Physical Chemistry                                                **
!**  Department of Chemistry                                           **
!**  Lund University                                                   **
!**  Sweden                                                            **
!**                                                                    **
!**  All rights reserved. The code may not be modified or              **
!**  redistributed without the written conscent of the copyright       **
!**  owner. The copyright owner does not take any responsibility       **
!**  for any error in the code or its documentation.                   **
!**                                                                    **
!************************************************************************
!************************************************************************

!************************************************************************
!*                                                                      *
!*     MeshModule                                                       *
!*                                                                      *
!************************************************************************

! ... module for mesh

module MeshModule
    implicit none
    type Node                             ! node in the DOP-tree
        !private
        real(8) :: dop(6)                 ! bounding box in object coordinate system
        integer(4) :: c(2)                ! id:s of children
        logical :: leaf                   ! leaf node: children is triangles
    end type

    type TriMesh  ! triangle mesh, with DOP-tree
        !private
        real(8), allocatable :: c(:,:)    ! c(3,np) coordinates of triangle verticies
        integer(4), allocatable :: t(:,:) ! t(3,nt) triangles as index as verticies into c
        type(node), allocatable :: n(:)   ! nodes in tree
        integer(4) :: levels              ! levels of subdivisions of triangles
    end type

    type AffineTrans                      ! affine transformation
        real(8) :: trans(3)               ! location of object origin in lab system
        real(8) :: rot(3,3)               ! rotation matrix applied to object
        integer(4) :: sel(3,6)            ! for each directed axis in lab system which 3 axises in object system contribute positively
    end type

    integer :: dopc                       ! counter of DOP overlap checks
    integer :: dopac                      ! counter of DOP overlap checks, with positive results
    integer :: tric                       ! counter of triangle overlap checks

    real(8), parameter ::  eye3(3,3) = reshape([1,0,0,0,1,0,0,0,1],[3,3])
    integer, parameter :: indrot(5) = [ 1, 2, 3, 1, 2 ]

contains

!........................................................................

!      BuildSuperball
!           !
!      BuildDopTree
!           !
!      BuildNode ----
!                   !
!                 FindExtrema, idSort, Volume, Union, Swap, BuildNode
!
!
!      TestOverlap
!           !
!      OverlapMesh
!           !
!      OverlapTree
!           !
!      OverlapLeaf
!           !
!      OverlapDop

!************************************************************************
!*                                                                      *
!*     BuildSuperball                                                   *
!*                                                                      *
!************************************************************************

! ... mesh a superball

function BuildSuperBall(radius, m, nlevels, max_err, rms_err) result(mesh)
    real(8), intent(in)    :: radius          ! r in x^m + y^m + z^m = r^m
    real(8), intent(in)    :: m               ! m in the same
    integer(4), intent(in) :: nlevels         ! nbr of refinement levels of the mesh
    real(8), intent(out), optional :: max_err ! max distance in r coordinate between mesh and ball
    real(8), intent(out), optional :: rms_err ! rms of the same distance, over all triangles
    type(TriMesh) :: mesh                     ! resulting mesh
    real(8) :: r,max_err2, sum2_err, c(3,3), bary(3), proj(3), err2, errv(3)
    integer(4) :: it, np, nt, a, i1, ip, k, level, i, ntl
    integer(4) :: p(3), q(3)
    integer(4), allocatable :: neigh(:,:)
    integer(1), allocatable :: d(:,:)
    real(8), parameter :: O = 0.0d0

    if (m < 200d0) then  ! superball

       np = 6 + 4*(lshift(1,2*nlevels)-1)
       nt = 8*lshift(1,2*nlevels)
       allocate( mesh%c(3,np) )
       allocate( mesh%t(3,nt) )
       allocate( d(3,nt) )
       d = 0
       allocate( neigh(3,nt) )
       neigh = 0

       r = radius
       mesh%c(:,1:6) = reshape( [ &
           r, O, O,  -r, O, O,   O, r, O,   O, O, r,   O,-r, O,   O, O,-r  &
       ], [ 3, 6 ] )
       mesh%t(:,1:8) = reshape( [ &
           1, 3, 4,   1, 5, 4,   1, 5, 6,   1, 3, 6,   2, 3, 4,   2, 5, 4,   2, 5, 6,   2, 3, 6  &
       ], [ 3, 8 ] )
       d(:,1:8) = reshape( [ &
           0, 0, 0,   1, 1, 0,   1, 0, 1,   0, 1, 1,   1, 0, 1,   0, 1, 1,   0, 0, 0,   1, 1, 0  &
       ], [ 3, 8 ] )

       ip = 6
       it = 8
       do level=1,nlevels
           ntl = it
           neigh(:,1:ip) = 0
           do k=1,ntl
               p = mesh%t(:,k)
               do i = 1,3
                   i1 = indrot(i+1)
                   a = p( indrot(i+d(i,k)))
                   if (neigh(i,a) == 0) then
                       ip = ip + 1
                       mesh%c(:,ip) = MidPoint(mesh%c(:,p(i)), mesh%c(:,p(i1)), m, radius)
                       neigh(i,a) = ip
                   end if
                   q(i) = neigh(i,a)
               end do
               mesh%t(:,it+1) = [ p(1), q(1), q(3) ]
               mesh%t(:,it+2) = [ q(1), p(2), q(2) ]
               mesh%t(:,it+3) = [ q(3), q(2), p(3) ]
               mesh%t(:,   k) = [ q(2), q(3), q(1) ]
               d(:,it+1:it+3) = reshape([ (d(:,k),i=1,3) ], [ 3, 3 ])
               !d(:,it+1:it+3) = spread(d(:,k),2,3)
               d(:,k) =  1- d(:,k)
               it = it + 3
           end do
       end do

    else     ! cube

       np = 8
       nt = 12
       allocate( mesh%c(3,np) )
       allocate( mesh%t(3,nt) )

       r = radius
       mesh%c(:,1:np) = reshape( [ &
          -r,-r,-r,   r,-r,-r,  -r, r,-r,   r, r,-r,   -r,-r, r,   r,-r, r,  -r, r, r,   r, r, r    &
       ], [ 3, np ] )

       mesh%t(:,1:nt) = reshape( [ &
           1, 4, 2,   1, 3, 4,   1, 2, 6,   1, 6, 5,   1, 7, 3,   1, 5, 7,    &
           4, 2, 8,   8, 2, 6,   6, 5, 8,   8, 5, 7,   3, 7, 4,   4, 7, 8     &
       ], [ 3, nt ] )

    end if

    max_err2 = 0d0
    sum2_err = 0d0
    do k=1,nt
        p = mesh%t(:,k)
        c = mesh%c(:,p)
        bary = sum(c,2)/3
        proj = MidPoint(bary,bary,m,radius)
        errv = bary - proj
        err2 = dot_product(errv,errv)
        max_err2 = max(max_err2, err2)
        sum2_err = sum2_err + err2
    end do
    if (present(max_err)) max_err = sqrt(max_err2)
    if (present(rms_err)) rms_err = sqrt(sum2_err/nt)

    call BuildDopTree(mesh)

end function BuildSuperball

!************************************************************************
!*                                                                      *
!*     BuildDopTree                                                     *
!*                                                                      *
!************************************************************************

! ... build a DOP tree

subroutine BuildDopTree(m)
    type(TriMesh), intent(inout) :: m          ! mesh, for which DOP-tree will be built
    real(8) :: bary(3,size(m%t,2))
    real(8) :: dop(6,size(m%t,2)), mdop(6)
    integer(4) :: ids(size(m%t,2))
    integer(4) :: nt, p(3), nn, nsize, i
    real(8) :: c(3,3)
    nt = size(m%t,2)
    do i = 1,nt
        p = m%t(:,i)
        c = m%c(:,p)
        bary(:,i) = sum(c,2)/3
        dop(1:3,i) = minval(c,2)
        dop(4:6,i) = maxval(c,2)
    end do
    mdop(1:3) = minval(m%c,2)
    mdop(4:6) = maxval(m%c,2)
    nsize = int( nt * 2.3 )
    allocate( m%n(nsize) )
    nn = 1
    ids = [ (i,i=1,nt) ]
    i = BuildNode(m,mdop,bary,dop,ids,nn,1)
end subroutine BuildDopTree

!************************************************************************
!*                                                                      *
!*     BuildNode                                                        *
!*                                                                      *
!************************************************************************

! ... build a node

recursive function BuildNode(m, ndop, bary, dop, ids, nn, lvl) result(nid)
    type(TriMesh), intent(inout) :: m  ! mesh/dop-tree
    real(8), intent(in)  :: ndop(6)    ! bounding box
    real(8), intent(in)  :: bary(:,:)  ! barycentric point of each triangle
    real(8), intent(in)  :: dop(:,:)   ! DOP for each trianlge
    integer(4), intent(inout) :: ids(:)! ids of the triangles this dop should contain
    integer(4), intent(inout) :: nn    ! next available node id (updated)
    integer(4), intent(in) :: lvl      ! current level in tree
    integer(4) :: nid                  ! id of this node
    integer(4) :: i(2), id, ntri, d, choice
    real(8) :: cdop(6,2), idop(6,2), eval(3,2)
    real(8) :: vol(2), ivol(2), dvol(2)
    logical :: pick_back, add_back
    nid = nn
#define node m%n(nid)
    nn = nn + 1
    if (lvl > m%levels ) m%levels = lvl

    node%dop = ndop
    ntri = size(ids)
    if (ntri <= 1) then
        node%leaf = .true.
        node%c(:) = 0
        node%c(1:ntri) = ids(1:ntri)
        return
    end if
    node%leaf = .false.

    eval = findExtrema(bary,ids)
    d = maxloc(abs(eval(:,2) - eval(:,1)),1)
    call idSort(bary(d,:),ids)
    cdop = dop(:,[ ids(1), ids(ntri) ])
    vol = volume(cdop)
    i =  [ 2, ntri-1 ]
    do while (i(1) <= i(2))
        pick_back = i(1)-1 > ntri-i(2)
        id = ids(i(1+l2i(pick_back)))
        idop(:,1) = union(cdop(:,1),dop(:,id))
        idop(:,2) = union(cdop(:,2),dop(:,id))
        ivol = volume(idop)
        dvol = ivol - vol
        if (dvol(1) == dvol(2)) then
            add_back = pick_back
        else
            add_back = dvol(1) > dvol(2)
        end if

        if (add_back .neqv. pick_back) then
            call swap(ids(i(1)),ids(i(2)))
        end if

        if (add_back) then
            i(2) = i(2) - 1
        else
            i(1) = i(1) + 1
        end if
        choice = 1+l2i(add_back)
        cdop(:,choice) = idop(:,choice)
        vol(choice) = ivol(choice)
    end do

    node%c(1) = BuildNode(m,cdop(:,1),bary,dop,ids(   1:i(2)),nn,lvl+1)
    node%c(2) = BuildNode(m,cdop(:,2),bary,dop,ids(i(1):ntri),nn,lvl+1)
    if (ntri > 18) then
        !print*, i(2), ntri-i(2), i(2)/real(ntri)
        !print*,
    end if
#undef node
end function BuildNode

!************************************************************************
!*                                                                      *
!*     TriTri                                                           *
!*                                                                      *
!************************************************************************

! ... check if two triangles intersect
!     Tomas Moller, Journal of Graphics Tools, 25, 2(2), 1997.

function TriTri(v, u) result(overlap)
    real(8), intent(in) :: u(3,3)            ! vertex coordinates of first triangle
    real(8), intent(in) :: v(3,3)            ! vertex coordinates of secound triangle
    logical :: overlap                       ! if the triangles did overlap
    real(8) :: isu(2), isv(2)
    real(8) :: e1(3), e2(3), nv(3), nu(3)
    real(8) :: mv, mu, du(3), dv(3), d(3)
    integer(4) :: ind

    e1 = v(:,2) - v(:,1)
    e2 = v(:,3) - v(:,1)
    nv = Cross(e1,e2)
    mv = -dot_product(v(:,1),nv)

    du = matmul(transpose(u),nv) + mv
    if (du(1)*du(2) > 0 .and. du(1)*du(3) > 0) then
        overlap = .false.
        return
    end if

    e1 = u(:,2) - u(:,1)
    e2 = u(:,3) - u(:,1)
    nu = Cross(e1,e2)
    mu = -dot_product(u(:,1),nu)

    dv = matmul(transpose(v),nu) + mu
    if (dv(1)*dv(2) > 0 .and. dv(1)*dv(3) > 0) then
        overlap = .false.
        return
    end if

    d = Cross(nv,nu)
    ind = maxloc(abs(d),1)
    if (cint(v(ind,:),dv,isv) .and. cint(u(ind,:),du,isu)) then
        call Sort2(isv)
        call Sort2(isu)
        overlap = isu(1) <= isv(2) .and. isv(1) <= isu(2)
    else
        overlap = coplanar(v,u,nu)
    end if

end function TriTri

!************************************************************************
!*                                                                      *
!*     Cint                                                             *
!*                                                                      *
!************************************************************************

! ... tries to find interval on the line intersection between the triangle planes

function cint(v, d, is)
    real(8), intent(in) :: v(3)         ! coordinates of vertices of triangle a in some axis
    real(8), intent(in) :: d(3)         ! projection of the vertices on the normal axis of the other trangle
    real(8), intent(out) :: is(2)       ! the interval triangle occupies om the intersection line
    logical :: cint                     ! true, if interval is nonempty
    real(8) :: v1, d1, vv(2), dd(2)
    integer(4) :: i(3)

    if (d(1)*d(2) > 0) then
        i = [3,1,2]
    elseif (d(1)*d(3) > 0) then
        i = [2,1,3]
    elseif (d(2)*d(3) > 0 .or. d(1) /= 0) then
        i = [1,2,3]
    elseif (d(2) /= 0) then
        i = [2,1,3]
    elseif (d(3) /= 0) then
        i = [3,1,2]
    else
        cint = .false.
        return
    end if
    v1 = v(i(1)); vv = v(i(2:3))
    d1 = d(i(1)); dd = d(i(2:3))
    is = v1 + (vv - v1)*d1/(d1-dd)
    cint = .true.
end function cint

!************************************************************************
!*                                                                      *
!*     Coplanar                                                         *
!*                                                                      *
!************************************************************************

! ... check for overlap in coplanar triangles

function coplanar(ut, vt, n) result(overlap)
    real(8), intent(in) :: ut(3,3)     ! vertex coordinates of first triangle
    real(8), intent(in) :: vt(3,3)     ! vertex coordinates of secound triangle
    real(8), intent(in) :: n(3)        ! normal to first triangle
    real(8) :: u(2,3), v(2,3)
    real(8) :: a(2), b(2), c(2), d, e, f
    logical :: overlap
    integer(4) :: ix, i(2), k, k2, m, m2

    ix = maxloc(abs(n),1)
    i = [indrot(ix+1),indrot(ix+2)]
    u = ut(i,:)
    v = vt(i,:)
    do k = 1,3
        k2 = indrot(k+1)
        a = v(:,k2) - v(:,k)
        do m = 1,3
            m2 = indrot(m+1)
            b = u(:,m) - u(:,m2)
            c = v(:,k) - u(:,m)
            f = Cross2(a,b)
            d = Cross2(b,c)
            e = Cross2(c,a)
            if (f == 0) cycle
            if (f > 0) then
                overlap = (d>=0 .and. d<=f .and. e>=0 .and. e<= f)
            else
                overlap = (d<=0 .and. d>=f .and. e<=0 .and. e>= f)
            end if
            if(overlap) return
        end do
    end do
    overlap = PointInTriangle(u(:,1),v) .or. PointInTriangle(v(:,1),u)
end function coplanar

!************************************************************************
!*                                                                      *
!*     PointInTriangle                                                  *
!*                                                                      *
!************************************************************************

! ... check if point in a triangle

function PointInTriangle(p, u) result(overlap)
    real(8), intent(in) :: p(2)    ! 2d coordinate of point
    real(8), intent(in) :: u(2,3)  ! 2d coordinates of vertices
    logical :: overlap             ! if the point was within the triangle
    real(8) :: a(2), c
    integer(4) :: d(3), i, i2
    do i = 1,3
        i2 = indrot(i+1)
        a = u(:,i2) - u(:,i)
        c = -a(1)*u(1,i)+a(2)*u(2,i)
        d(i) = l2i(a(1)*p(1)-a(2)*p(2)>c)
    end do
    overlap = d(1)==d(2) .and. d(2)==d(3)
end function PointInTriangle

!************************************************************************
!*                                                                      *
!*     Cross2                                                           *
!*                                                                      *
!************************************************************************

! ... cross product of two 2d vectors

function Cross2(a, b) result(c)
    real(8), intent(in) :: a(2)
    real(8), intent(in) :: b(2)
    real(8):: c
    c = a(1)*b(2) - a(2)*b(1)
end function Cross2

!************************************************************************
!*                                                                      *
!*     Cross                                                            *
!*                                                                      *
!************************************************************************

! ... Cross product of two 3d vectors

function Cross(a, b)
    real(8), intent(in) :: a(3)
    real(8), intent(in) :: b(3)
    real(8) :: Cross(3)
    Cross(1) = a(2) * b(3) - a(3) * b(2)
    Cross(2) = a(3) * b(1) - a(1) * b(3)
    Cross(3) = a(1) * b(2) - a(2) * b(1)
end function Cross

!************************************************************************
!*                                                                      *
!*     Sort2                                                            *
!*                                                                      *
!************************************************************************

! ... sort two elements in ascending order

subroutine Sort2(v)
    real(8), intent(inout) :: v(2)
    real(8) :: temp
    if (v(2) < v(1)) then
        temp = v(2)
        v(2) = v(1)
        v(1) = temp
    end if
end subroutine Sort2

!************************************************************************
!*                                                                      *
!*     MidPoint                                                         *
!*                                                                      *
!************************************************************************

! ... find midpoint of two 3d arrays, midpoint may be modified

function MidPoint(a, b, m, rad)
    real(8), intent(in) :: a(3), b(3)
    real(8), intent(in) :: m
    real(8), intent(in) :: rad
    real(8) :: MidPoint(3)
    MidPoint = 0.5d0*(a+b)
    MidPoint = MidPoint * rad / sum(abs(MidPoint)**m)**(1/m)
end function MidPoint

!************************************************************************
!*                                                                      *
!*     l2i                                                              *
!*                                                                      *
!************************************************************************

! ... change a logic variable to an integer

function l2i(l) result(i)
    logical, intent(in) :: l   ! logical
    integer(4) :: i            ! 0 or 1 for false and true
    !i = l
    i = merge(1,0,l)
end function

!************************************************************************
!*                                                                      *
!*     findExtrema                                                      *
!*                                                                      *
!************************************************************************

! ...

function findExtrema(c, ids) result(eval)
    real(8),    intent(in) :: c(:,:)  ! all coordnates
    integer(4), intent(in) :: ids(:)  ! ids of coordinates to consider
    real(8) ::  eval(3,2)             ! extrema for each axis (forms bounding box)
    integer(4) :: i
    eval = c(:,[ids(1), ids(1)])
    do i=2,size(ids)
        where (eval(:,1) > c(:,ids(i)))
            eval(:,1) = c(:,ids(i))
        end where
        where (eval(:,2) < c(:,ids(i)))
            eval(:,2) = c(:,ids(i))
        end where
    end do
end function

!************************************************************************
!*                                                                      *
!*     union                                                            *
!*                                                                      *
!************************************************************************

! ...

function union(a, b) result(u)
    real(8), intent(in) :: a(6)    ! box a
    real(8), intent(in) :: b(6)    ! box b
    real(8) ::  u(6)               ! smallest box enclosing both
    u = merge(a,b, [a(1:3) < b(1:3),a(4:6) > b(4:6)] )
end function

!************************************************************************
!*                                                                      *
!*     volume                                                           *
!*                                                                      *
!************************************************************************

! ...

function volume(dop) result(vol)
    real(8), intent(in) :: dop(:,:)  ! two dops
    real(8) :: vol(2)                ! respective volume
    where (any( dop(1:3,:) >= dop(4:6,:), 1))
        vol = 0
    elsewhere
        vol = product(dop(4:6,:) - dop(1:3,:), 1)
    end where
end function

!************************************************************************
!*                                                                      *
!*     idSort                                                           *
!*                                                                      *
!************************************************************************

! ...

subroutine idSort(vec,ind)
   real(8),    intent(in)  :: vec(:)   ! unsorted array
   integer(4), intent(inout) :: ind(:) ! index array

   integer(4) :: n, l, i, j, ir,indt
   real(8) :: vect

   n = size(ind)

   if (n == 1) return
   l = n/2+1
   ir = n
10 continue
      if (l>1) then
         l = l-1
         indt = ind(l)
         vect = vec(indt)
      else
         indt = ind(ir)
         vect = vec(indt)
         ind(ir) = ind(1)
         ir = ir-1
         if (ir == 1) then
         ind(1) = indt
            return
         end if
      end if
      i = l
      j = l+l
   20 if (j <= ir) then
         if (j < ir) then
            if (vec(ind(j)) < vec(ind(j+1))) j = j+1
         end if
         if (vect < vec(ind(j))) then
            ind(i) = ind(j)
            i = j
            j = j+j
         else
            j = ir+1
         end if
      goto 20
      end if
      ind(i) = indt
   goto 10
end subroutine idSort

!************************************************************************
!*                                                                      *
!*     swap                                                             *
!*                                                                      *
!************************************************************************

! ...

subroutine swap(a, b)
   integer(4), intent (inout)   :: a, b    ! a and b are swapped
   integer(4) :: temp
   temp = a ; a = b ; b = temp
end subroutine swap

!************************************************************************
!*                                                                      *
!*     transformation                                                   *
!*                                                                      *
!************************************************************************

! ...

function transformation(rot, trans) result(tf)
    real(8), intent(in) :: rot(3,3) ! rotation matrix
    real(8), intent(in) :: trans(3) ! translation vector
    type(AffineTrans) :: tf         ! affine transformation with precalculated selections
    integer(4) :: a,b, i, j
    tf%rot = rot
    tf%trans = trans
    do i = 1,3
        do j = 1,3
            a = j
            b = j + 3
            if(rot(i,j) < 0) call swap(a,b)
            tf%sel(j,i) = a
            tf%sel(j,i+3) = b
        end do
    end do
end function transformation

!************************************************************************
!*                                                                      *
!*     rotation                                                         *
!*                                                                      *
!************************************************************************

! ...

function rotation(ang) result(rot)
    real(8), intent(in) :: ang(3)    ! rotational angles
    real(8) :: rot(3,3)              ! rotational matrix
    real(8) :: c ,s, ri(3,3)
    integer(4) :: i,k(2)
    rot = eye3
    do i=1,3
        ri = eye3
        k = [indrot(i+1),indrot(i+2)]
        c = cos(ang(i))
        s = sin(ang(i))
        ri(k,k) = reshape([ c, s, -s, c], [2,2])
        rot = matmul(rot,ri)
    end do
end function rotation

!************************************************************************
!*                                                                      *
!*     TestOverlap                                                      *
!*                                                                      *
!************************************************************************

! ...

function TestOverlap(m, ang, tr) result(overlap)
    type(TriMesh), intent(in) :: m     ! mesh
    real(8), intent(in) :: ang(3)      ! relative rotation
    real(8), intent(in) :: tr(3)       ! relative translaton
    logical :: overlap
    real(8) :: rot(3,3)
    rot = rotation(ang)
    overlap = OverlapMesh(m,tr,rot)
end function TestOverlap

!************************************************************************
!*                                                                      *
!*     OverlapMesh                                                      *
!*                                                                      *
!************************************************************************

! ... translate and rotate a mesh and check overlap of two trees

function OverlapMesh(m, trans, rot) result(overlap)
    type(TriMesh), intent(in) :: m    ! the mesh describing both particles
    real(8), intent(in) :: trans(3)   ! relative translation
    real(8), intent(in) :: rot(3,3)   ! relative rotation
    type(AffineTrans) :: tf
    logical :: overlap
    tf = transformation(rot,trans)
    overlap = OverlapTree(m,tf)
end function OverlapMesh

!************************************************************************
!*                                                                      *
!*     OverlapTree                                                      *
!*                                                                      *
!************************************************************************

! ... check overlap of two trees

function OverlapTree(m, tf) result(overlap)
    type(TriMesh) :: m         ! mesh
    type(AffineTrans) :: tf    ! transforation of b relative a
    logical :: overlap         ! if they overlap
    integer(4) :: ai, bi, i, j
    integer(4) :: st(2,4*m%levels+4), sp
    st(:,1) = [1,1] ! root nodes
    sp = 1
    do while (sp > 0)
       ai = st(1,sp); bi = st(2,sp)
       sp = sp - 1
#define a m%n(ai)
#define b m%n(bi)
       dopc = dopc + 1
       overlap = OverlapDop(tf,a%dop,b%dop)
       if (.not. overlap) cycle
       dopac = dopac + 1

       if (a%leaf) then
           if (b%leaf) then
               overlap = OverlapLeaf(m,tf,ai,bi)
               if (overlap) return
               cycle
           else
               call swap(ai,bi)
           end if
       end if

       do i=1,2
           if (b%leaf) then
               sp = sp + 1
               st(:,sp) = [a%c(i), bi]
           else
               do j=1,2
                   sp = sp + 1
                   st(:,sp) = [a%c(i), b%c(j)]
               end do
           end if
       end do
#undef a
#undef b
    end do
end function OverlapTree

!************************************************************************
!*                                                                      *
!*     OverlapDop                                                       *
!*                                                                      *
!************************************************************************

! ... check overlap between two bounding boxes

function OverlapDop(tf, d1, d2) result(overlap)

    type(AffineTrans), intent(in) :: tf      ! rotation of box 2 into box 1's system
    real(8), intent(in) :: d1(6)             ! bounding box 1 in its own system
    real(8), intent(in) :: d2(6)             ! bounding box 2 in its own system
    logical :: overlap                       ! whenever the enclosing box around box 2 overlaps with box 1
    real(8) :: h, l
    integer(4) :: i
    overlap = .false.
    do i = 1,3
        l = tf%trans(i) + dot_product(tf%rot(i,:),d2(tf%sel(:,i)))
        h = tf%trans(i) + dot_product(tf%rot(i,:),d2(tf%sel(:,i+3)))
        if (max(d1(i),l) > min(d1(i+3),h)) return
    end do
    overlap = .true.
end function OverlapDop

!************************************************************************
!*                                                                      *
!*     OverlapLeaf                                                      *
!*                                                                      *
!************************************************************************

! ... check overlap of two leafs

function OverlapLeaf(m, tf, ai, bi) result(overlap)
    type(TriMesh), intent(in) :: m           ! the mesh
    type(AffineTrans), intent(in) :: tf      ! rotation of leaf node b
    integer(4), intent(in) :: ai             ! idx of leaf node a
    integer(4), intent(in) :: bi             ! idx of leaf node b
    logical :: overlap                       ! if some pair of trianlgles overlapped
    real(8) :: c1(3,3), c2(3,3)
    integer(4) :: i, j
#define a m%n(ai)
#define b m%n(bi)
    do i = 1,2
        if( a%c(i) == 0) cycle
        do j = 1,2
            if( b%c(j) == 0) cycle
            c1 = m%c(:,m%t(:,a%c(i)))
            c2 = m%c(:,m%t(:,b%c(j)))
            c2 = matmul(tf%rot,c2) + spread(tf%trans,2,3)
            tric = tric + 1
            overlap = TriTri(c1,c2)
            if(overlap) return
        end do
    end do
#undef a
#undef b
end function OverlapLeaf

end module

