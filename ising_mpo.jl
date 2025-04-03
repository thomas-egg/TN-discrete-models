using ITensors, ITensorMPS

# Time evolution using a simple Trotter step (this is schematic)
function time_evolve_master(psi, W)
    psi = tdvp(W, -1.5im, psi; nsweeps=10)
    println("Time evolution step completed.")
    psi = normalize(psi)
    return psi
end

let
    L = 64            # Number of sites
    β = 1.0           # Inverse temperature
    τ = 1.0           # Time step
    
    # Define sites with spin-1/2 for the Ising model
    sites = siteinds(n -> "S=1/2", L)
    
    # Create the AutoMPO for the Hamiltonian
    ampo = AutoMPO()
    
    # Adding terms to the Hamiltonian
    for j in 1:L
        add!(ampo, 1/τ, "Id", j)           # Identity term for each site
        add!(ampo, -1/τ, "Sx", j)          # Spin interaction in x-direction
        add!(ampo, -1/τ * exp(-β), "Sz", j) # Spin interaction in z-direction (with Boltzmann factor)
    end
    
    # Construct MPO from the AutoMPO
    W = MPO(ampo, sites)
    
    # Random initial state in MPS form
    psi0 = randomMPS(sites; linkdims=4)
    
    # Time evolution using the master equation
    psi_final = time_evolve_master(psi0, W)
    
    # Sample the final state
    result = sample(psi_final)
    exp_en = inner(psi_final', W, psi_final)
    
    # Print and write the final state
    println("Expected Energy: ", exp_en)
    println("Final state: ", result)
    write("result.txt", result)
end