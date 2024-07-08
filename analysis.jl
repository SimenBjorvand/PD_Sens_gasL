# function to make plots of my simulation results so I can plot them
using DelimitedFiles, PyPlot, Statistics

# read simulation data:
pathC  = "sim_data/msSim/"
pathSD = "sim_data/SD/"
pathN  = "sim_data/newton/"
pathSDsens = "sim_data/SD_sens/"
pathNsens  = "sim_data/newton_sens/"

solTC  = readdlm(pathC*"solTime.txt")
solTSD = readdlm(pathSD*"solTime.txt")
solTN  = readdlm(pathN*"solTime.txt")
solTSDsens = readdlm(pathSDsens*"solTime.txt")
solTNsens  = readdlm(pathNsens*"solTime.txt")


####
times = readdlm(pathN*"time.txt")[:,1]
timesShort = 0:times[end]/(length(solTN)):times[end]


fig,ax = subplots(1,figsize = [14*0.7,10*0.7])
ax.step(timesShort,[solTC[1];solTC],label="Central")
ax.step(timesShort,[solTSD[1];solTSD]/5,label="SD")
ax.step(timesShort,[solTN[1];solTN]/5,label="N")
ax.step(timesShort,[solTSDsens[1];solTSDsens]/5,label="SDsens")
ax.step(timesShort,[solTNsens[1];solTNsens]/5,label="Nsens")
ax.set_ylabel("Solution time")
ax.set_xlabel("Time [h]")
ax.set_ylim([0,3])
ax.legend()
fig


mC   = mean(solTC)
mSD  = mean(solTSD)/5
mN   = mean(solTN)/5
mSDs = mean(solTSDsens)/5
mNs  = mean(solTNsens)/5

#####
mSD/mC
mN/mC
mSDs/mC
mNs/mC

####
