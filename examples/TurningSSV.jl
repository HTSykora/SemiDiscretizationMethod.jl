using SemiDiscretizationMethod

using Plots

#--------------------Problem definitnion------------------------
# computing the stablility of the turning process in case of SipindeSpeedVariation (in dimensionless form)
# ζ damping of the oscillator (woth normalized natural frequency = 1)
# kw cuttin gcoefficients
# Ω spindle speed
# ASSV amplitude of the delay (spindle speed) relative SipindeSpeedVariation (ASSV<<1 typically)
# T timeperiod of the spindle speed variation

function createTurningSSVProblem(kw,ζ,Ω,ASSV,T)
    AMx =  ProportionalMX(t->@SMatrix [0. 1.; -1-kw -ζ]);
    τ1=t -> 2π/Ω *(1+ ASSV*sin(t/T*2π)) # if function is needed, the use τ1 = t->foo(t)
    BMx1 = DelayMX(τ1,t->@SMatrix [0. 0.; kw 0.]);
    cVec = Additive(t->@SVector [0.,cos(t*2π)])#τ1(t)*kw
    LDDEProblem(AMx,[BMx1],cVec)
end;

kw=0.2
ζ=0.1
Ω=0.3
ASSV=0.1
NT=10 
T=2π / Ω * NT  #spindle speed variation in integer number of revolution (NT)

τmax=2π/Ω *(1+ ASSV) # the largest τ of the system

turningSSV_lddep=createTurningSSVProblem(kw,ζ,Ω,ASSV,T); # LDDE problem for Turning SSV problem
method=SemiDiscretization(1,0.025) # 3rd order semi discretization with Δt=0.05
### -------- Try to increase the resolutionby decreasing Δt


#<<<<<<<<<<<<<<<<<<<<<<Stability calculation in one point>>>>>>>>>>>>>>>>>>>>>>>

#--------------------------------------------------------------------------------
#-----------------traditional mapping of SemiDiscretization----------------------
#-----------------------faster computation if T>>τmax----------------------------
#--------------------------------------------------------------------------------

@time mapping=DiscreteMapping_1step(turningSSV_lddep,method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true); #The discrete mapping of the system

@time spectralRadiusOfMapping(mapping); # spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
@time fixPointOfMapping(mapping); # stationary solution of the system (equilibrium position)



#----------------------------------------------------------------------------------------------
#--------------------#updated Left\Right mapping matrices of SemiDiscretization----------------
#----------------------- Order of magnitude faster computation if T<=τmax----------------------
#----------------------------------------------------------------------------------------------
@time mappingLR=DiscreteMapping_LR(turningSSV_lddep,method,τmax,
    n_steps=Int((T+100eps(T))÷method.Δt),calculate_additive=true); #The discrete mapping of the system

@time spectralRadiusOfMapping(mappingLR); # spectral radius ρ of the mapping matrix (ρ>1 unstable, ρ<1 stable)
@time fixPointOfMapping(mappingLR); # stationary solution of the system (equilibrium position)


