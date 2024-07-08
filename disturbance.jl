# File to create disturbance trajectories for my uncertain variable ΔGOR

# ΔGOR1 ∈ [-0.04 0.04]
# ΔGOR1 ∈ [-0.01 0.01]

function dist()
    ΔGOR11 = -0.03*ones(12)
    ΔGOR21 = zeros(12)

    ΔGOR12 = [-0.03:(0.05)/11:0.02;0.02]
    ΔGOR22 = zeros(12)

    ΔGOR13 = 0.02*ones(12)
    ΔGOR23 = zeros(12)

    ΔGOR14 = 0.02*ones(12)
    ΔGOR24 = Array(0:(0.01)/11:0.01)

    ΔGOR15 = 0.02*ones(12)
    ΔGOR25 = 0.01*ones(12)

    ΔGOR1 = [ΔGOR11;ΔGOR12;ΔGOR13;ΔGOR14;ΔGOR15]
    ΔGOR2 = [ΔGOR21;ΔGOR22;ΔGOR23;ΔGOR24;ΔGOR25]
    return ΔGOR1,ΔGOR2

    #=
    fig,ax = subplots(figsize = [8,8])
    ax.plot(1:60,ΔGOR1.+par.GOR[1])
    ax.plot(1:60,ΔGOR2.+par.GOR[2])
    fig=#
end


function dist2()
    ΔGOR11 = -0.03*ones(120)
    ΔGOR21 = zeros(120)

    ΔGOR12 = -0.03:(0.05)/119:0.02
    ΔGOR22 = zeros(120)

    ΔGOR13 = 0.02*ones(120)
    ΔGOR23 = zeros(120)

    ΔGOR14 = 0.02*ones(120)
    ΔGOR24 = Array(0:(0.01)/119:0.01)

    ΔGOR15 = 0.02*ones(120)
    ΔGOR25 = 0.01*ones(120)

    ΔGOR1 = [ΔGOR11;ΔGOR12;ΔGOR13;ΔGOR14;ΔGOR15]
    ΔGOR2 = [ΔGOR21;ΔGOR22;ΔGOR23;ΔGOR24;ΔGOR25]
    return ΔGOR1,ΔGOR2

end



function dist3()
    ΔGOR11 = -0.04*ones(100)
    ΔGOR21 = -0.01*ones(100)

    ΔGOR12 = -0.04:(0.08)/199:0.04
    ΔGOR22 = -0.01*ones(200)

    ΔGOR13 = 0.04*ones(200)
    ΔGOR23 = -0.01:(0.02)/199:0.01

    ΔGOR14 = 0.04*ones(100)
    ΔGOR24 = 0.01*ones(100)


    ΔGOR1 = [ΔGOR11;ΔGOR12;ΔGOR13;ΔGOR14]
    ΔGOR2 = [ΔGOR21;ΔGOR22;ΔGOR23;ΔGOR24]
    
    return ΔGOR1,ΔGOR2
end