device=0

confFile=confs/N00010000-T0.7-100x020x020.gro
dt=0.005
nstep=3000000
confFeq=2000
thermoFeq=100
rcut=5.0
nlistExten=0.49

refT=0.7
tauT=1.0

NTCell=96
NTAtom=96

# adapt rcut parameters
densityProfileSamplingFeq=40
rcutAssignFeq=40
rcutUpdateFeq=20000
rcutNumRefine=2
refh=1.0
rcmin=03.0
rcmax=10.0
rcstep=0.5
targetPrec=0.010

# force correction parameters
assignForceCorrFeq=20
updateForceCorrFeq=20000
