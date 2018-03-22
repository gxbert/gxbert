source /afs/slac.stanford.edu/package/geant4/vol52/gcc/setup-gcc.sh

ulimit -c unlimited

alias vi="vim"

G4INSTALL=/home/tkoi/BERT/standalone/rev0/geant4-install

if [ $?LD_LIBRARY_PATH == 1 ]
then
   export LD_LIBRARY_PATH=$G4INSTALL/lib64
else
   export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$G4INSTALL/lib64
fi

export G4ENSDFSTATEDATA=/afs/slac.stanford.edu/package/geant4/vol55/data/G4ENSDFSTATE2.2
export G4LEVELGAMMADATA=/afs/slac.stanford.edu/package/geant4/vol55/data/PhotonEvaporation5.1
