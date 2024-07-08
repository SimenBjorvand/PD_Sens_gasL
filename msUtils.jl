# Author: Simen Bjorvand
# Date: 07.12.21  Original(28.09.21)
# Description:
# This file will contain help functions for msMPC

# Required packages:
# ~

#Required files:
# ~ 


# Function definition:

function scenarioGeneration(pBounds,NFE,NRH)
    Np = length(pBounds)
    NPD = length(pBounds[1])
    NPR = NPD^Np
    NSC = NPR^NRH
    ty = typeof(pBounds[1][1])

    pSet = Array{Array{ty,1},1}(undef,NPR)
    for i in 1:NPR
        pSet[i] = Array{ty,1}(undef,Np)
    end
    for np in 1:Np
        for i in 1:NPD^(np-1)
            for j in 1:NPD
                for k in 1:NPD^(Np-np)
                    pSet[(i-1)*NPD^(Np-np+1) + (j-1)*NPD^(Np-np) + k][np] = pBounds[np][j]
                end
            end
        end
    end

    scenarios = Array{Array{Array{ty,1},1},1}(undef,NSC)
    for i in 1:NSC
        scenarios[i] = Array{Array{ty,1},1}(undef,NFE)
    end
    for i in 1:NFE
        if i <= NRH
            for j = 1:NPR^i
                for k = 1:NPR^(NRH-i)
                    scenarios[(j-1)*(NPR^(NRH-i))+k][i] = pSet[mod(j-1,NPR)+1]
                end
            end   
        else
            for j = 1:NSC
                scenarios[j][i] = pSet[mod(j-1,NPR)+1]     #pDist[np,(NPD+1)÷2]      # pSet[(length(pSet)+1)÷2] 
            end
        end
    end

    return scenarios, pSet
end


function scenarioGenerationBox(pBounds,NFE,NRH)
    Np = length(pBounds)
    NPR = 2^Np + 1
    NSC = NPR^NRH
    ty = typeof(pBounds[1][1])

    pSet = Array{Array{ty,1},1}(undef,NPR)
    for i in 1:NPR
        pSet[i] = Array{ty,1}(undef,Np)
    end

    function multVec(v,w)
        Nv = length(v)
        Nw = length(w)
        d = Array{Array{Float64,1},1}(undef,Nv*Nw)
        for nv in 1:Nv
            for nw in 1:Nw
                i = (nv-1)*Nw+nw
                d[i] = [v[nv];w[nw]]
            end
        end
        return d
    end
    hold = pBounds[Np][[1,end]]
    if Np > 1
        for i in 1:Np-1
            hold = multVec(pBounds[Np-i][[1,end]],hold)
        end
    end
    pSet[1:2^(Np-1)] = hold[1:2^(Np-1)]
    pSet[2^(Np-1)+1] = [i[2] for i in pBounds]
    pSet[2^(Np-1)+2:end] = hold[2^(Np-1)+1:end]


    scenarios = Array{Array{Array{ty,1},1},1}(undef,NSC)
    for i in 1:NSC
        scenarios[i] = Array{Array{ty,1},1}(undef,NFE)
    end
    for i in 1:NFE
        if i <= NRH
            for j = 1:NPR^i
                for k = 1:NPR^(NRH-i)
                    scenarios[(j-1)*(NPR^(NRH-i))+k][i] = pSet[mod(j-1,NPR)+1]
                end
            end   
        else
            for j = 1:NSC
                scenarios[j][i] = pSet[mod(j-1,NPR)+1]     #pDist[np,(NPD+1)÷2]      # pSet[(length(pSet)+1)÷2] 
            end
        end
    end

    return scenarios, pSet
end



function weightPerScenario1(wpSet,NRH)
    if NRH == 1
        return wpSet
    else
        bSet = weightPerScenario1(wpSet,NRH-1)
        return vcat([bSet*i for i in wpSet]...)
    end
end