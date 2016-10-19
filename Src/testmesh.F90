module TestMesh
    use MeshModule
    implicit none
    real(8), parameter :: pi = 3.1415926535897932d0
contains
    function benchmark(lvls,d1,d2,n) result(t)
        type(TriMesh) :: m
        real(8) :: ang(3),phi,th,trlen, tr(3)
        real(8) :: s,f, d1, d2, t
        integer(4) :: lvls,n,i,k
        logical :: q
        call srand(1)
        m = buildSuperBall(1d0,2.5d0,lvls)
        !print*, 'lvls', m%levels
        call cpu_time(s)
        k = 0
        tric = 0
        dopc = 0
        dopac = 0
        do i=1,n
            ang = [rand()*2*pi, rand()*pi, rand()*2*pi]
            phi = rand()*2*pi
            th = rand()*pi
            trlen = d1 + (d2-d1)*rand()
            tr = trlen*[ sin(th)*cos(phi), sin(th)*sin(phi), cos(th)]
            if (testOverlap(m,ang,tr)) k = k + 1
        end do
        call cpu_time(f)
        t = (f-s)/n
        write(*,*) d1, d2, real(k)/n, t
        write(*,*) 'dop', real(dopc)/n, 'tri', real(tric)/n
        !print*, 'dopac', real(dopac)/dopc
    end function

    subroutine err_meas(mlvl,m)
        integer(4) :: mlvl, i
        real(8) :: m, max_err, rms_err
        type(TriMesh) :: mesh
        do i=1,mlvl
            mesh = buildSuperBall(1.0d0,m,i,max_err,rms_err)
            write(*,*) i, max_err, rms_err
        end do
    end subroutine

end module

program TestMeshProgram
    use TestMesh
    use MeshModule
    type(TriMesh) :: m
    integer(4) :: i
    real(8) :: s,f, u(3,3), v(3,3)
    logical :: d, q
    real(8) :: rad

    write(*,*) 'm = 2'
    call err_meas(5,2d0)
    write(*,*) 'm = 2.5'
    call err_meas(5,2.5d0)
    write(*,*) 'm = 3'
    call err_meas(5,3d0)
    write(*,*) 'm = 10'
    call err_meas(5,10d0)

    u(:,1) = [0,0,0]
    u(:,2) = [1,0,0]
    u(:,3) = [0,1,0]
    v(:,1) = [0.1,0.1,0.0]
    v(:,2) = [0.9,0.1,0.0]
    v(:,3) = [0.1,0.9,0.0]
    write(*,*) tritri(u,v)

    m = buildSuperBall(1d0,2.5d0,4)
    do i=-1,13
#define str (0.02d0)
        rad = 2.0+i*str
        s = benchmark(4,rad,rad+str,10000)
    end do
    !print*, testOverlap(m,[0.,0.,0.],[0.05,0.8,0.0])
    !print*, testOverlap(m,[pi/3,pi/2,0.0d0],[1.9d0,0.8d0,0.0d0])
    !print*, testOverlap(m,[pi/3,pi/2,0.0d0],[1.9d0,0.9d0,0.0d0])
    do i=1,5
        write(*,*) 'lvl', i
        write(*,*) benchmark(i,1.8d0,2.0d0,100000)
    end do
end program


