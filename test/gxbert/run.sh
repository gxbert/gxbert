#!/bin/sh

#64k events for each jobs
#real    52m34.128s
#user    1405m34.008s
#sys     0m27.701s

run_() {
   proj=$1
   targ=$2
   energy=$3
   ./gxbert $proj $targ $energy > result.${proj}.${targ}.${energy} &
   #( ./gxbert $proj $targ $energy | grep RESULT > result.${proj}.${targ}.${energy}; root -l -b 'compair_only_ratio.C( "'$proj'", "'$targ'" , "'$energy'" )' ) &
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

