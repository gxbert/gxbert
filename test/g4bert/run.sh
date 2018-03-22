#!/bin/sh
#This script will submit 
#12(particle type)x5(target type)x3(projectile energy)
#180 jobs
#
#Testing machine 
#
#MemTotal:       65,681,388 kB
#16 x Intel(R) Xeon(R) CPU E5-2640 v3 @ 2.60GHz
#
#
#1k events 
#real    0m31.321s
#user    5m23.352s
#sys     9m23.421s
#
#64k events
#real    11m56.846s
#user    319m8.287s
#sys     9m48.965s

run_() {
   proj=$1
   targ=$2
   energy=$3
   ./g4bert $proj $targ $energy > result.${proj}.${targ}.${energy} &
   #./g4bert $proj $targ $energy | grep RESULT > result.${proj}.${targ}.${energy} &
}

energy_() {
   proj=$1
   targ=$2
   run_ $proj $targ 800MeV
   run_ $proj $targ 1.5GeV
   run_ $proj $targ 3GeV
}

targ_() {
  proj=$1
  energy_ $proj C
  energy_ $proj Al
  energy_ $proj Fe
  energy_ $proj In
  energy_ $proj Pb
}

targ_ proton
targ_ neutron
targ_ gamma
targ_ kaon-
targ_ kaon+
targ_ lambda
targ_ pi-
targ_ pi+
targ_ sigma-
targ_ sigma+
targ_ omega-
targ_ xi-
wait
