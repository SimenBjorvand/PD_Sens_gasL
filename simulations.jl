# Date: 25.10.2022
# Modified: 27.09.2023
# Author: Simen Bjorvand
# Description: Run the MPC simulation here for the nominal case

# Libraries:
using JuMP, Ipopt, Statistics, SparseArrays, DifferentialEquations, Sundials, DelimitedFiles, PyPlot, LinearAlgebra

# My files/codes:
include("myParams.jl")
include("plant.jl")
include("/home/simen/PhD/utils/colMatrix.jl")
include("msUtils.jl")
include("decMPCUtil.jl")
include("sensDecMPC.jl")
include("disturbance.jl")


# Simulation parameters:
max_iterations = 60                 # [-]
dt             = 30.                 # [s]
simTime        = max_iterations*dt*10   # [s]

# Model parameters:
par = gasLiftPar()
Nu = par.Nw
Nx = 3*par.Nw + 2
Nz = 13*par.Nw + 7
ΔGOR1,ΔGOR2 = dist3()


###
# Get initial condtions and guesses for mpc
ω0 = readdlm("ss_data/init_cond/noGasLift.txt")[:,1]
u0 = ω0[1:Nu]
x0 = ω0[Nu+1:end][1:Nx]
z0 = ω0[Nu+1:end][Nx+1:end]

u0 = round.(u0,digits=4)
x0 = round.(x0,digits=4)
z0 = round.(z0,digits=4)

ωSS = readdlm("ss_data/nominalOpt.txt")[:,1]
uSS = ωSS[1:Nu]
xSS = ωSS[Nu+1:end][1:Nx]
zSS = ωSS[Nu+1:end][Nx+1:end]

# Set up parameters to store the solution
uL = Array{Array{Float64,1},1}(undef,Nu)
xL = Array{Array{Float64,1},1}(undef,Nx)
zL = Array{Array{Float64,1},1}(undef,Nz)

for i in 1:Nu
    uL[i] = [u0[i]]
end
for i in 1:Nx
    xL[i] = [x0[i]]
end
for i in 1:Nz
    zL[i] = [z0[i]]
end

solveTimes = Array{Float64,1}(undef,0)
limSolveTimes = Array{Float64,1}(undef,0)
iterList = Array{Float64,1}(undef,0)
PCIterList = Array{Float64,1}(undef,0)
dirList = Array{Float64,1}(undef,0)
gradList = Array{Float64,1}(undef,0)


# Set up the controller:
lprev = [uSS]


# Set up the iterative simulation parameters:
uk = copy(u0)
xk = copy(x0)
zk = copy(z0)

for i in 1:max_iterations
    println("_______________________")
    println("Iteration: ", i)

    # Call on the distributed MPC
    uk,iter,solveT,lprev,npk,ngk,nIL,nFTL,nBTL,nLT = sensDistMPCN(xk,uk,lprev,par,uSS,xSS,zSS,80)

    st = sum(sum(solveT))

    snIL = sum(sum(nIL))

    snLT = sum(nLT)

    uk = round.(uk,digits=16)

    for j in 1:10
        # Simulate plant for next state:
        xk,zk = Plant(xk,zk,uk,[ΔGOR1[10*(i-1)+j],ΔGOR2[10*(i-1)+j]],par,dt)
        xk = round.(xk,digits=16)
        zk = round.(zk,digits=16)

        # Save the data from this timestep
        for i in 1:Nu
            push!(uL[i],uk[i])
        end
        for i in 1:Nx
            push!(xL[i],xk[i])
        end
        for i in 1:Nz
            push!(zL[i],zk[i])
        end
    end

    # Prepare the  mpc for the next iteration
    push!(solveTimes,st)
    push!(limSolveTimes,snLT)
    push!(iterList,iter)
    push!(PCIterList,snIL)
    push!(dirList,npk)
    push!(gradList,ngk)
end


fig3,ax3 = subplots(1,figsize=[8,8])
#ax3.plot(1:60,itSD)
ax3.step(1:60,iterList)
#ax3.plot(1:60,solveTimes*0.2)
fig3


####
fig4,ax4 = subplots(1,figsize=[8,8])
#ax3.plot(1:60,itSD)
#ax4.step(1:60,solveTimes)
ax4.step(1:60,msTime)
#ax3.plot(1:60,solveTimes*0.2)
fig4

nsTime = readdlm("sim_data/newton_sens/solTime.txt")
nTime  = readdlm("sim_data/newton/solTime.txt")
msTime = readdlm("sim_data/msSim/solTime.txt")
ssTime = readdlm("sim_data/SD_sens/solTime.txt")
sTime = readdlm("sim_data/SD/solTime.txt")

fig4,ax4 = subplots(1,figsize=[8,8])
#ax4.step(1:60,msTime)
ax4.step(1:60,nTime)
ax4.step(1:60,nsTime)
fig4

fig4,ax4 = subplots(1,figsize=[8,8])
#ax4.step(1:60,msTime)
ax4.step(1:60,nTime)
ax4.step(1:60,sTime)
fig4

fig4,ax4 = subplots(1,figsize=[8,8])
#ax4.step(1:60,msTime)
ax4.step(1:60,nTime)
ax4.step(1:60,ssTime)
fig4

fig4,ax4 = subplots(1,figsize=[8,8])
ax4.step(1:60,msTime)
ax4.step(1:60,sTime)
ax4.step(1:60,nsTime)
ax4.step(1:60,ssTime)
fig4

timeS = 0:dt/3600:simTime/3600


#=
fig1,ax1 = subplots(5,2,figsize=[12,8])
ax1[1,1].plot(timeS,xL[1])
ax1[1,2].plot(timeS,xL[2])
ax1[2,1].plot(timeS,xL[3])
ax1[2,2].plot(timeS,xL[4])
ax1[3,1].plot(timeS,xL[5])
ax1[3,2].plot(timeS,xL[6])
ax1[4,1].step(timeS,uL[1])
ax1[4,2].step(timeS,uL[2])
ax1[5,1].plot(timeS,xL[7])
ax1[5,2].plot(timeS,xL[8])

fig1


fig2,ax2 = subplots(4,2,figsize=[12,8])
ax2[1,1].plot(timeS,zL[19])
ax2[1,2].plot(timeS,zL[20])
ax2[2,1].plot(timeS,zL[21])
ax2[2,2].plot(timeS,zL[22])
ax2[3,1].step(timeS,uL[1])
ax2[3,2].step(timeS,uL[2])
ax2[4,1].plot(timeS,zL[32])
ax2[4,2].plot(timeS,zL[33])
fig2
=#


# Save the simulation data
#=
open("sim_data/SD_sens/x.txt","w") do io
    writedlm(io,xL)
end

open("sim_data/SD_sens/u.txt","w") do io
    writedlm(io,uL)
end

open("sim_data/SD_sens/z.txt","w") do io
    writedlm(io,zL)
end


open("sim_data/SD_sens/time.txt","w") do io
    writedlm(io,timeS)
end

open("sim_data/SD_sens/solTime.txt","w") do io
    writedlm(io,solveTimes)
end

open("sim_data/SD_sens/limitSolTime.txt","w") do io
    writedlm(io,limSolveTimes)
end

open("sim_data/SD_sens/iterations.txt","w") do io
    writedlm(io,iterList)
end

open("sim_data/SD_sens/predCorIters.txt","w") do io
    writedlm(io,PCIterList)
end

open("sim_data/SD_sens/dir.txt","w") do io
    writedlm(io,dirList)
end

open("sim_data/SD_sens/grad.txt","w") do io
    writedlm(io,gradList)
end


=#


#=
open("sim_data/stepest_decent_2step/GOR.txt","w") do io
    writedlm(io,GOR)
end
=#