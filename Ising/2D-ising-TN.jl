#=
Tensor Network  Julia implementation of a 2D square Ising model

Tensor Renormalization Group (TRG) algorithm can be used to com-
pute the partition function of a 2-dimensional lattice model. In
this case, we are looking at a 2D lattice model
=#

# Function to return Ising model values
using ITensors
function Sig(s)
    return 1.0 - 2.0 * (s-1)
end

function Energy(si, sj, sk, sl)
    E = (Sig(si) * Sig(sj)) + (Sig(sj) * Sig(sk)) + (Sig(sk) * Sig(sl)) + (Sig(sl) * Sig(si))
    return E
end

# MAIN
function main(kbT)
    
    let
        # Parameters
        dim = 2
        #kbT = parse(Float64, arg["kbT"])
        maxdim = 20
        topscale = 7

        # Set indices
        s = Index(dim, "scale=0")
        l = addtags(s, "left")
        r = addtags(s, "right")
        u = addtags(s, "up")
        d = addtags(s, "down")

        # Construct tensor
        A = ITensor(l, r, u, d)

        # Loop over all tensors
        for sl in 1:dim
            for sd in 1:dim
                for sr in 1:dim
                    for su in 1:dim

                        # Build up A
                        E = Energy(sl, sd, sr, su)
                        P = exp(-E / kbT)
                        A[l=>sl, d=>sd, r=>sr, u=>su] = P

                    end
                end
            end
        end

        # Partition function per-site
        z = 1.0

        # TRG algorithm
        for scale in 1:topscale

            println("\n---------- Scale $(scale) -> $(1 + scale)  ----------")

            # Factorize tensor A into Fr, Fl
            Fl, Fr = factorize(A, (r, d), which_decomp="svd", maxdim=maxdim)

            # Set tags of Fd and Fu
            Fu, Fd = factorize(A, (l, d), which_decomp="svd", maxdim=maxdim)

            # Set tags
            Fl = replacetags(Fl, "Link,fact", "left")
            Fr = replacetags(Fr, "Link,fact", "left")
            Fu = replacetags(Fu, "Link,fact", "up")
            Fd = replacetags(Fd, "Link,fact", "up")

            # Define shared indices for Fr
            l_new = commoninds(Fl, Fr)
            r_new = replacetags(l_new, "left", "right")
            Fr *= delta(l_new, r_new)

            # " " for Fd
            u_new = commoninds(Fu, Fd)
            d_new = replacetags(u_new, "up", "down")
            Fd *= delta(u_new, d_new)

            # Renormalize A
            Fl *= delta(r, l)
            Fu *= delta(d, u)
            Fr *= delta(l, r)
            Fd *= delta(u, d)
            A = Fl * Fu * Fr * Fd

            # Update variables
            l = l_new
            r = r_new
            u = u_new
            d = d_new

            # Normalize tensor and update Z
            trace_A = (A * delta(l, r) * delta(u, d))[1, 1, 1, 1]
            A /= trace_A
            z *= trace_A^(1.0 / (2^(1+scale)))

        end

        # Print log of Z
        println(log(z))

        # Final value
        return log(z)
    end
end