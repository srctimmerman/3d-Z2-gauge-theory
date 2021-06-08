using DelimitedFiles
using Plots
using Random; Random.seed!(0);
rng = MersenneTwister(1234);

include("CubicLattice.jl")

const EqSteps = 10000   #number of burn in steps
const NumSteps = 100000 #number of update steps counted in the simulation

#parameters of the 3D Z2 gauge theory
J = 1
B = 0

#size of the lattice
L = 8
M = L^3
N = 3*L^3
Spin = fill(1,N);

#generating data structures to enable the simulation to efficiently
#navigate the 3D cubic lattice
UpCube = make_up_cube(L);
println("Done w UpCube!")
AssociatedCube = make_associated_cube(L);
println("Done w Associated Cube!")
CubeStars = make_cube_stars(L);
println("Done w Cube Stars!")

#here is the brute force calculation of the energy
function Energy_Total(Spin, UpCube)

    H_trans = sum(Spin)
    H_plaq = 0
    for i = 1:M     #iterate over up-cubes
        for j = 1:3 #iterate over all of the plaquettes associated with a single lattice site
            H_plaq += Spin[UpCube[i,j,1]] * Spin[UpCube[i,j,2]] * Spin[UpCube[i,j,3]] * Spin[UpCube[i,j,4]]
        end # limit all arithmetic to integers for as long as we can
    end
    return -J*H_plaq - B*H_trans
end #Energy_Total

#here is the energy DIFFERENCE when the spin at spin_index is flipped,
#calculated from the local square plaquettes
function Energy_Diff(Spin,spin_index, AssociatedCube)

    plaquette_factor = 0
    #iterate over plaquettes involving this spin
    for i = 1:4
        Plaq = AssociatedCube[spin_index, i,:]
        plaquette_factor -= Spin[Plaq[1]] * Spin[Plaq[2]] * Spin[Plaq[3]] * Spin[Plaq[4]]
    end

    return -2*J*plaquette_factor + 2*B*Spin[spin_index]
end

#computes whether an update with energy cost DeltaE should be accepted
function MetropolisAccept(DeltaE, beta)::Bool
    if DeltaE <= 0
        return true
    else
        rnum = rand(rng)  #random number for Metropolis
        if (exp(-beta*DeltaE) > rnum)
            return true
        end
    end
    return false
end

#data structures to hold thermodynamic output
Tarr = 2:-0.02:1
EMC = zeros(Float64,0)
SpecHeat = zeros(Float64,0)
Eexact = zeros(Float64,0)
Acceptance = zeros(Float64,0)

#loops over different temperatures
for T in Tarr  #count down

    beta = 1.0/T

    #initialize the energy
    Energy = 0
    Energy = Energy_Total(Spin, UpCube)

    #burning in to an initial config
    for step = 1:EqSteps

        #multiple single spin flips
        for i = 1:N
            spin_i = rand(1:N)
            DeltaE = Energy_Diff(Spin,spin_i, AssociatedCube)
            if MetropolisAccept(DeltaE, beta)
                Energy += DeltaE
                Spin[spin_i] *= -1
            end
        end #i

        #star flips
        for i = 1:L^3
            site_i = rand(1:L^3)
            for j = 1:6
                Spin[CubeStars[site_i,j]] *= -1
            end
        end
    end #Equilibration

    E_avg = 0
    E2 = 0
    A_rate = 0
    #actual Monte Carlo prouction steps
    for step = 1:NumSteps
        Accept = 0
        for i = 1:N  #multiple single spin flips
            spin_i = rand(1:N)
            DeltaE = Energy_Diff(Spin,spin_i, AssociatedCube)
            if MetropolisAccept(DeltaE, beta)
                Energy += DeltaE
                Spin[spin_i] *= -1
                Accept += 1
            end
        end #i

        #star flips
        for i = 1:L^3
            site_i = rand(1:L^3)
            for j = 1:6
                Spin[CubeStars[site_i,j]] *= -1
            end
        end

        E_avg += Energy
        E2 += Energy*Energy
        A_rate += Accept/N
    end #Monte Carlo production step

    #compute, print out, and store observables for this temperature
    Cv = E2/NumSteps - (E_avg/NumSteps)^2
    println(T," ",E_avg/NumSteps," ",E2/NumSteps," ",Cv/N/T/T," ",A_rate/NumSteps)

    push!(EMC,E_avg/NumSteps/N)
    push!(SpecHeat,Cv/(T*T*N))
    push!(Acceptance,A_rate/NumSteps)

end #T

open(string(L, "_.txt"), "w") do io
    writedlm(io, [Tarr EMC SpecHeat Acceptance])
end
