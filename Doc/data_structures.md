Data structures
===============
The following data structures are used for input variables:

**bonds etc**
```fortran
type bond_var
   real(8) :: k ! force constant
   integer(4) :: p ! power
   real(8) :: eq ! equilibrium value
end type bond_var
```

**adsorption conditions**
```fortran
type adscond_var
   character(8) :: txplane !'xy-plane' is yet the only option
   character(3) :: txend ! '-','+','-|+' for negative, positive, both ends
   real(8) :: offset ! threshold distance of adsorption for the com of a particle
end type adscond_var
```

**stochastic data**
```fortran
type sf_var
   integer(4) :: nplow ! lower particle number
   integer(4) :: npupp ! upper particle number
   integer(4) :: ndim ! number of dimension of the stochastic variable
   integer(4) :: fac ! factor describing the separation of stoch. data
   integer(4) :: ngr ! logical flag if for single sf value
end type sf_var
```

**correlation functions (input)**
```fortran
type cf_input_var
   integer(4) :: nmean ! length (number of values) forming a mean
   integer(4) :: nlevel ! length (number of values) forming a level
   integer(4) :: nolevel ! number of levels
   integer(4) :: legendre ! order of Legendre polynomial (1, 2, or 3)
   logical :: lsvalue ! logical flag if for single sf value
   logical :: lsubmean ! logical flag if subtraction of mean of cf
   logical :: lnorm ! logical flag if normalization of cf
end type cf_input_var
```

**correlation functions**
```fortran
type cf_var
   integer(4)     :: nmean                            ! length (number of values) forming a mean
   integer(4)     :: nlevel                           ! length (number of values) forming a level
   integer(4)     :: ratio                            ! nlevel/nmean
   integer(4)     :: nolevel                          ! number of levels
   integer(4)     :: legendre                         ! order of Legendre polynomial (1, 2, or 3)
   logical        :: lsvalue                          ! logical flag if for single sf value
   logical        :: lsubmean                         ! logical flag if subtraction of mean of cf
   logical        :: lnorm                            ! logical flag if normalisation of cf

   real(8), allocatable    :: sf(:,:,:,:)             ! stochastic function
   real(8), allocatable    :: sf_aver(:,:,:)          ! stochastic function, local average (of particle)
   real(8), allocatable    :: sf_mean(:,:)            ! stochastic function, gobal mean
   real(8), allocatable    :: cf(:,:,:)               ! time correlation function
   real(8), allocatable    :: cf2(:,:,:)              ! time correlation function squared
   integer(4), allocatable :: Np(:,:,:)               ! counter of correlation sampling (of particle)
   integer(4), allocatable :: Ngr(:,:,:)              ! counter of correlation sampling (of group)
   integer(8), allocatable :: Nlev(:,:)               ! counter of level (of particle)
end type cf_var
```

**1D static variables**
```fortran
type static1D_var
   logical :: l ! logical flag for engagement
   real(8) :: min ! minimum value
   real(8) :: max ! maximum value
   integer(4) :: nbin ! number of bins
   logical :: lnorm ! logical flag for normalization
   character(10):: label ! title
   real(8) :: nvar ! expanded into nvar variables (not read)
end type static1D_var
```

**2D static variables**
```fortran
type static2D_var
   logical :: l ! logical flag for engeagement
   real(8) :: min(2) ! minimum value
   real(8) :: max(2) ! maximum value
   integer(4) :: nbin(2) ! number of bins
   logical :: lnorm ! logical flag for normalization
   character(10):: label ! title
   real(8) :: nvar ! expanded into nvar variables (not read)
end type static2D_var
```

**potential energy**
```fortran
   type potenergy_var
real(8), allocatable :: twob(:)       ! two-body contribution (excluding ewald contribution)
   real(8), allocatable :: oneb(:)       ! one-body contribution (dielectric discontinuity)
   real(8)              :: tot           ! total
   real(8)              :: rec           ! reciprocal electrostatic space contribution (UEwald)
   real(8)              :: stat          ! static electrostatic contribution (umbodyp)
   real(8)              :: pol           ! polarization electrostatic contribution (umbodyp)
   real(8)              :: bond
   real(8)              :: angle
   real(8)              :: crosslink
   real(8)              :: external      ! external contribution
end type potenergy_var
```

**chain properties**
```fortran
type chainprop_var
   real(8)    :: ro(3)                   ! center of mass
   real(8)    :: rbb2                    ! bead-to-bead distance squared
   real(8)    :: angle                   ! angle between consecutive beads
   real(8)    :: cos                     ! cos(180 - angle between consecutive beads)
   real(8)    :: ree(3)                  ! end-to-end vector
   real(8)    :: ree2                    ! end-to-end distance squared
   real(8)    :: rg2                     ! radius of gyration squared
   real(8)    :: rg2s                    ! square extention along principal axes (smallest)
   real(8)    :: rg2m                    ! square extention along principal axes (middle)
   real(8)    :: rg2l                    ! square extention along principal axes (largest)
   real(8)    :: rg2z                    ! radius of gyration squared projected on the z-axis
   real(8)    :: rg2xy                   ! radius of gyration squared projected on the xy-plane
   real(8)    :: lpree                   ! persistence length based on end-to-end separation
   real(8)    :: lprg                    ! persistence length based on radius of gyration
   real(8)    :: shape                   ! square end-to-end distance / square of radius of gyration
   real(8)    :: asph                    ! asphericity (JCP 100, 636 (1994))
   real(8)    :: torp                    ! toroidicity
end type chainprop_var
```

**network properties**
```fortran
type networkprop_var
   real(8)    :: ro(3)                   ! center of mass
   real(8)    :: rg2                     ! radius of gyration squared
   real(8)    :: rg2x                    ! radius of gyration squared projected on the x-axis
   real(8)    :: rg2y                    ! radius of gyration squared projected on the y-axis
   real(8)    :: rg2z                    ! radius of gyration squared projected on the z-axis
   real(8)    :: rg2s                    ! square extention along principal axes (smallest)
   real(8)    :: rg2m                    ! square extention along principal axes (middle)
   real(8)    :: rg2l                    ! square extention along principal axes (largest)
   real(8)    :: eivr(3,3)               ! normalized eigenvectors of the principal frame
   real(8)    :: theta(3)                ! angles of axes of largest extension and x-, y-, and z-axes of main frame
   real(8)    :: asph                    ! asphericity (JCP 100, 636 (1994))
   real(8)    :: alpha                   ! degree of ionization (for titrating systems)
end type networkprop_var
```

**simple averaging**
```fortran
type aver_var
   real(8)      :: s2                    ! summation/averaging over steps
   real(8)      :: s1                    ! summation/averaging over macrosteps
end type aver_var
```

**complexation**
```fortran
type cluster_var
   integer(4), allocatable :: ip
   integer(4), allocatable :: np
end type cluster_var
```

**blocks in chains**
```fortran
type :: block_type
   integer(4)  :: pt  !particle type
   integer(4)  :: np  !number of particles
end type block_type
```

**scalar quantities**
```fortran
type scalar_var
   character(40) :: label                             ! label
   real(8)       :: norm                              ! normalization factor
   integer(4)    :: nsamp1                            ! number of macrosteps sampled
   integer(4)    :: nsamp2                            ! number of values sampled per macrostep
   real(8)       :: avs1                              ! average of the run
   real(8)       :: avsd                              ! precision of average of the run
   real(8)       :: avs2                              ! average of a macrostep
   real(8)       :: fls1                              ! fluctuation of the run
   real(8)       :: flsd                              ! precision of fluctuation of the run
   real(8)       :: fls2                              ! fluctuation of a macrostep
   real(8)       :: value                             ! value of a configuration
   integer(4)    :: nsamp                             ! number of samplings
   integer(4)    :: nblocklen                         ! for sampling with variable blocklen
   integer(4)    :: nblock(mnblocklen)                ! for sampling with variable blocklen
   real(8)       :: av_s1(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: av_sd(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: av_s2(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: fl_s1(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: fl_sd(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: fl_s2(mnblocklen)                 ! for sampling with variable blocklen
   real(8)       :: av_sd_extrap                      ! for sampling with variable blocklen
   real(8)       :: av_sd_stateff                     ! for sampling with variable blocklen
   real(8)       :: fl_sd_extrap                      ! for sampling with variable blocklen
   real(8)       :: fl_sd_stateff                     ! for sampling with variable blocklen
end type scalar_var
```

**1D distribution functions**
```fortran
type df_var
   character(27) :: label                             ! label
   real(8)       :: norm                              ! normalization factor
   integer(4)    :: nsamp1                            ! number of macrosteps sampled
   integer(4)    :: nsamp2                            ! number of values sampled per macrostep
   real(8)       :: min                               ! minimum value of df
   real(8)       :: max                               ! maximum value of df
   integer(4)    :: nbin                              ! number of grid points
   real(8)       :: bin                               ! grid length of df
   real(8)       :: bini                              ! inverse of bin
   real(8)       :: nsampbin(-1:mnbin_df)             ! number of values sampled in each bin during macrostep
   real(8)       :: avs1(-1:mnbin_df)                 ! average of the run
   real(8)       :: avsd(-1:mnbin_df)                 ! precision of average of the run
   real(8)       :: avs2(-1:mnbin_df)                 ! average of a macrostep
end type df_var
```

**2D distribution function**
```fortran
type df2d_var
   character(27) :: label                             ! label of df
   real(8)       :: norm                              ! normalization factor
   integer(4)    :: nsamp1                            ! number of macrosteps sampled
   integer(4)    :: nsamp2                            ! number of values sampled per macrostep
   real(8)       :: min(2)                            ! minimum value of df
   real(8)       :: max(2)                            ! maximum value of df
   integer(4)    :: nbin(2)                           ! number of grid points
   real(8)       :: bin(2)                            ! grid length of df
   real(8)       :: bini(2)                           ! 1/bini
   real(8)       :: avs1(-1:mnbin_df2d,-1:mnbin_df2d) ! average of the run
   real(8)       :: avsd(-1:mnbin_df2d,-1:mnbin_df2d) ! precision of average of the run
   real(8)       :: avs2(-1:mnbin_df2d,-1:mnbin_df2d) ! average of a macrostep
end type df2d_var
```

**cluster1 trial move**
```fortran
type cluster1_tm_var
   logical      :: l                    ! flag for cluster1 tm
   real(8)      :: p                    ! relative probability of cluster1 move
   real(8)      :: rad                  ! radius of region for cluster members
   real(8)      :: psel                 ! probability to select a particle within rad
   real(8)      :: dtran                ! maximal translational trial displacement
   real(8)      :: drot                 ! maximal rotational trial displacement
end type cluster1_tm_var
```

**cluster 2 trial move**
```fortran
type cluster2_tm_var
   logical      :: l                    ! flag for cluster1 tm
   real(8)      :: p                    ! relative probability of cluster2 move
   real(8)      :: dtran                ! maximal translational trial displacement
   real(8)      :: drot                 ! maximal rotational trial displacement
   integer(4)   :: mode                 ! =0 : search members only of type iptmove
                                        ! =1 : search members across all particle types
end type cluster2_tm_var
```

**trial move**
```fortran
type trialmove_var
   logical      :: l                    ! flag for type of trial move
   real(8)      :: p                    ! relative probability of type of trial move
   real(8)      :: dtran                ! maximal translational trial displacement
   real(8)      :: drot                 ! maximal rotational trial displacement
   integer(4)   :: mode                 ! specific for type of trial displacment
   logical      :: lcl1                 ! flag for cluster1 trial move
   real(8)      :: radcl1               ! radius of region for cluster members
   real(8)      :: pselcl1              ! probability to select a particle within rad
end type trialmove_var
```

**node in the DOP-tree**
```fortran
type Node                             ! node in the DOP-tree
  !private
  real(8) :: dop(6)                 ! bounding box in object coordinate system
  integer(4) :: c(2)                ! id:s of children
  logical :: leaf                   ! leaf node: children is triangles
end type
```

**triangle mesh, with DOP-tree**
```fortran
type TriMesh  ! triangle mesh, with DOP-tree
  !private
  real(8), allocatable :: c(:,:)    ! c(3,np) coordinates of triangle verticies
  integer(4), allocatable :: t(:,:) ! t(3,nt) triangles as index as verticies into c
  type(node), allocatable :: n(:)   ! nodes in tree
  integer(4) :: levels              ! levels of subdivisions of triangles
end type
```

**affine transformation**
```fortran
type AffineTrans                      ! affine transformation
  real(8) :: trans(3)               ! location of object origin in lab system
  real(8) :: rot(3,3)               ! rotation matrix applied to object
  integer(4) :: sel(3,6)            ! for each directed axis in lab system which 3 axises in object system contribute positively
end type
```

**sso step**
```fortran
type :: step
   integer(8)  :: n  !number of steps
   real(8)     :: d2  !squared displacement
   real(8)     :: d4 !displacement**4
end type step
```

**SSOPart**
```fortran
type  :: ssopart_var
   real(8)     :: fac      !increment of part length
   integer(4)  :: nextstep !step at which next part starts
   integer(4)  :: i        !current part
   integer(4)  :: n        !number of parts
end type ssopart_var
type(ssopart_var), save  :: SSOPart
```

**SSOParameters**
```fortran
type  :: ssoparam_var
   real(8)     :: used     ! used dtran
   real(8)     :: opt      ! dtran with the highest mobility
   real(8)     :: err   ! accuracy of opt
end type ssoparam_var
type(ssoparam_var), save, allocatable  :: SSOParameters(:,:)
```

**sso mobility**
```fortran
type  :: mobility_var
   real(8)     :: val         !value
   real(8)     :: error       !error
   real(8)     :: smooth      !smooth
end type mobility_var
type(mobility_var), allocatable, save :: Mobility(:)
```

**celllist cell-pointer-array**
```fortran
type cell_pointer_array
   type(cell_type), pointer              :: p => null()        ! pointer to a cell, usefull to create an array of pointers
end type cell_pointer_array
```

**celllist cell-type**
```fortran
type cell_type
   integer(4)                            :: id                 ! for easy recognition
   integer(4)                            :: npart              ! number of particles per cell
   integer(4)                            :: nneighcell         ! number of neighbouring cells
   type(cell_pointer_array), allocatable :: neighcell(:)       ! pointer to the neighbouring cells
   integer(4)                            :: iphead             ! first particle in the linked list
end type cell_type
```
