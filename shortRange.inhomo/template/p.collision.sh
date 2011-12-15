device=0

NTCell=256
NTAtom=96

# ball generation ##################################
ball_init_conf=none
ball_output_conf=ball.init.gro

ball_dt=0.002
ball_nstep=1000000
ball_confFeq=200000
ball_thermoFeq=100
ball_rcut=2.5
ball_nlistExten=0.49
ball_nlistSizeFactor=1000.0
ball_clistDivision=1
# ##################################################
ball_refT=0.33
ball_tauT=1.0


# collision ########################################
collision_init_conf=confs/ball.init.gro

collision_dt=0.002
collision_nstep=100000
collision_confFeq=2000
collision_thermoFeq=100
collision_rcut=5.0
collision_nlistExten=0.49
collision_nlistSizeFactor=100.0
collision_clistDivision=1
collision_prec=0.030

collision_x=0.36
collision_u=3.17

# density profile parameters
refh=1.0
densityProfileSamplingFeq=10

# adapt rcut parameters
rcutAssignFeq=10
rcutUpdateFeq=200
rcutNumRefine=1
rcmin=02.5
rcmax=05.0
rcstep=0.25

assignForceCorrFeq=10
updateForceCorrFeq=200
