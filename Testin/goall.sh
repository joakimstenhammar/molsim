#!/bin/bash
#*************************************************************
#*                   Start of job                            *
#*************************************************************

#set -e

changelink () {
job=$1

for i in pos ori liv anv for tor group utot

do 
   rm -f $job.ana.$i
   ln -s $job.$i $job.ana.$i
done
}

if [ -f "../version.conf" ]
then
        version=".$(cat ../version.conf)"
else
        verson=""
fi
#####################
molsim="molsim_ser$version"
core=1
####################
rm -r out
cp -r in out

cd out
 
#gotest
$molsim test                      $core

changelink test
$molsim test.ana                  $core
$molsim test.mix                  $core

#gofree
$molsim  free.mc.uext1            $core
$molsim  free.mc.uext2            $core
$molsim  free.md1                 $core
$molsim  free.md2                 $core
$molsim  di.mc                    $core
 
#gohs
$molsim  hs.mc                    $core
$molsim  hs.mc.percolation1       $core
$molsim  hs.mc.percolation2       $core
$molsim  hs.mc.dimer              $core
$molsim  hs2d.mc                  $core
$molsim  hs2d.mc.lradatbox        $core
$molsim  ramp.mc                  $core
$molsim  ramp.mc.ads              $core
$molsim  sw.mc                    $core
$molsim  sw.mc.ads                $core

#gocoul
$molsim  holm.ew1                 $core
$molsim  holm.ew2                 $core
$molsim  holm.ew3                 $core
$molsim  salt.mc.mi               $core
$molsim  salt.mc.ew               $core
$molsim  salt.mccl2.ew            $core
$molsim  salt.mc.cyl.mol          $core
$molsim  salt.squarewell          $core
$molsim  ion.mc.umb               $core

#godip
$molsim  qq.ew.surf               $core
$molsim  qq.ew.nosurf             $core
$molsim  d.ew.surf                $core
$molsim  d.ew.nosurf              $core
$molsim  dd.ew.surf               $core
$molsim  dd.ew.nosurf             $core
$molsim  dip.mc.ew                $core
$molsim  dip.mc.ewspm             $core
$molsim  dip.mc.rd.ew             $core
$molsim  dip.mc.to.ew             $core
$molsim  dip.mc.sphimage          $core
cp dip.md.ew.cnf.save dip.md.ew.cnf
$molsim  dip.md.ew                $core
$molsim  dip.md.rd.ew             $core
$molsim  dip.md.to.ew             $core
$molsim  dippol.md.ew1            $core
$molsim  dippol.md.ew2            $core
$molsim  dippol.md.rf             $core

#gomic
$molsim  mic.mc.ocm               $core
$molsim  mic.mc.sph               $core
$molsim  mic.mc.cyl.mf            $core
$molsim  mic.mc.cyl.pmf           $core
$molsim  mic.mc.ew                $core
$molsim  mic.mc.ewspm             $core
$molsim  mic.mccl1.ew             $core
$molsim  mic.mccl2.ew             $core
$molsim  mic2.mc.ew               $core
$molsim  linecharge.mc.cyl        $core

#godieldis
$molsim  dieldis.p11              $core
$molsim  dieldis.p12              $core
$molsim  dieldis.p21              $core
$molsim  dieldis.p22              $core
$molsim  dieldis.p23              $core
$molsim  dieldis.p4               $core
$molsim  dieldis.s11a             $core
$molsim  dieldis.s11b             $core
$molsim  dieldis.s11c             $core
$molsim  dieldis.s12a             $core
$molsim  dieldis.s12b             $core
$molsim  dieldis.s12c             $core
$molsim  dieldis.s21a             $core
$molsim  dieldis.s21b             $core
$molsim  dieldis.s21c             $core
$molsim  dieldis.s22a             $core
$molsim  dieldis.s22b             $core
$molsim  dieldis.s22c             $core
$molsim  dieldis.s23a             $core
$molsim  dieldis.s23b             $core
$molsim  dieldis.s23c             $core
$molsim  dieldis.s21a_bg          $core
$molsim  dieldis.s2_center        $core
$molsim  dieldis.z0.q             $core
$molsim  dieldis.z0.p1            $core
$molsim  dieldis.z0.p2            $core
$molsim  dieldis.z10.q            $core
$molsim  dieldis.z10.p1           $core
$molsim  dieldis.z10.p2           $core

#gomol
$molsim bw.mc                     $core
$molsim b.md                      $core
$molsim pt.sph1                   $core
$molsim pt.sph2                   $core
$molsim lj.mc                     $core
$molsim lj.mc.rd                  $core
$molsim lj.mc.to                  $core
$molsim lj.mc.llist               $core
$molsim lj.mc.vlistllist          $core
$molsim lj.mc.mvt                 $core
$molsim lj.mc.npt                 $core
$molsim lj.md                     $core
$molsim w.mc                      $core
changelink w.mc
$molsim w.mc.ana                  $core
$molsim w.md                      $core

#gonemo
$molsim  wnemo.mcall.ew           $core
$molsim  wnemo.mcall.rf           $core
$molsim  wnemo.mcall.umb          $core
$molsim  wnemo.md                 $core
$molsim  wnemo.md.ew              $core
$molsim  wnemo.md.ew1             $core
$molsim  wnemo.md.ew2             $core
$molsim  wnemo.md.ewspm           $core
$molsim  wnemo.md.rf              $core
$molsim  wnemo.md.int             $core

#gochain
$molsim  chain.mc.002             $core
$molsim  chain.mc.003             $core
$molsim  chain.mc.004             $core
$molsim  chain.mc.trico           $core
$molsim  chain.mc.diblock         $core
$molsim  chain.mc.triblock        $core
$molsim  chain.mc.homo.diblock    $core
$molsim  chain.mc.nohs            $core
$molsim  chain.mc.hs              $core
$molsim  chain.mc.confined        $core
$molsim  chain.mc.line            $core
$molsim  chain.mc.spart           $core
$molsim  chain.mc.pivot1          $core
$molsim  chain.mc.pivot2          $core
$molsim  chain.mc.pivot3          $core
$molsim  chain.mc.chain           $core
$molsim  chain.mc.slither         $core
$molsim  chain.mc.2chains         $core
$molsim  chain.mc.2chaintypes     $core
$molsim  chain.mc.static          $core
$molsim  chain.md.nohs            $core
$molsim  chain.bd.nohs            $core

#gopads
$molsim p.mc.ads                  $core
$molsim p.bd.ads                  $core # run
$molsim p.bd.ads1                 $core # run 1
cp p.bd.ads1.cnf p.bd.ads2.cnf
cp p.bd.ads1.user p.bd.ads2.user
$molsim p.bd.ads2                 $core # run 2
cp ads_20-mer ads_20-mer2
echo "comparison of FUSER and ads_20-mer"
diff p.bd.ads2.user p.bd.ads.user
diff ads_20-mer2 ads_20-mer
echo "comparison done"

#gope
$molsim  pe.mc.sph1               $core
$molsim  pe.mc.sph2               $core
$molsim  pe.mc.sph3               $core
$molsim  pe.mc.sph4               $core
$molsim  pe.mc.ell                $core
$molsim  pe.mc.ew.catanionic      $core
$molsim  pe.mc.mol                $core
$molsim  pe.md.ew                 $core
$molsim  pepe.mc                  $core
$molsim  pem.mc.ew                $core
$molsim  peprot.mc                $core
$molsim  pamf.md.ew               $core
$molsim  prot.mccl1.ew            $core
cp daniel.cnf.save daniel.cnf
$molsim  daniel                   $core

#goweakcharge
$molsim  pe.mc.weakcharge         $core
# $molsim  pnw.mc.weakcharge        $core
$molsim  lysozyme.mc.weakcharge   $core
$molsim  ipma.mc.weakcharge       $core

#gonetwork
$molsim  pnw.mc1                  $core
$molsim  pnw.mc2                  $core
# $molsim  pnw.mc.finite            $core
# $molsim  pnw.mc.finite2           $core
$molsim  penw.mc.ew               $core
$molsim  penw.mc.finite           $core
$molsim  steffi_ion_hollowsphere  $core
$molsim  steffi_ion_coreshell     $core

#gobrush
$molsim  pamf.planbrush           $core
$molsim  pamf.sphbrush.sph        $core
$molsim  pamf.sphbrush.cyl        $core
$molsim  pamf.sphbrush.ew         $core

#gohierarchical
$molsim  bbrush.mc                $core
$molsim  bbrush.mc.ads            $core
$molsim  dendri.mc                $core

#gosuperball
$molsim sb_x                      $core
$molsim sb_xz                     $core
$molsim sb_xyz                    $core
$molsim sb_y                      $core
$molsim sb_z                      $core
$molsim sb_z_y090                 $core
$molsim sb_z_y180                 $core
$molsim sb_z_y030                 $core
$molsim sb_z_y045                 $core
$molsim sb_z_xyz                  $core
$molsim sb_nr                     $core
$molsim sb_mesh                   $core
$molsim sb.mc.ads                 $core
 
$molsim  capsid.mc.charge         $core
$molsim  capsid.mc.uniform        $core
$molsim  jurij1                   $core
$molsim  jurij2                   $core
$molsim  ao.mc                    $core
$molsim  ellipsoidob.mc           $core
$molsim  ellipsoidpro.mc          $core
$molsim  ellipsoidpro.dip.mc      $core

#gomix
$molsim  w.mc.mix                 $core
$molsim  b.md.mix                 $core
$molsim  hs.b2.mix                $core

#bugfixes
$molsim  ltime                    $core

echo
echo "Remove group and dump files"
echo
rm -f *.group
rm -f *.pos;  rm -f *.ori;  rm -f *.liv;  rm -f *.anv;  rm -f *.for;  rm -f *.tor; rm -f *.utot
 
echo
echo "Your job is completed!"
echo

cd ../

diffall.sh
