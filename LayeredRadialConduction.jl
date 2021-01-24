using Plots

mutable struct Grid
    R   # radius array
    nR  # number of grid points
    dR  # grid spacing
    D   # diffusivity array
    K   # conductivity array
    T   # temperature array
    j   # boundary index
end

function setup(diffProbe, condProbe, diffSoil, KM)
    #### USER INPUTS ##########################
    rProbe = 0.00476; # radius of probe in cm
    nR = 101;     # number of grid points
    modelSize = 5;   # model size vs probe radius
    ###########################################


    # grid spacing
    RModel = rProbe * modelSize;
    dR = RModel/(nR-1);
    R = LinRange(0, RModel, nR); # set the radius coordinates

    ### Make diffusivity and conductivity grid
    K = KM * ones(nR);
    #K[1:1+(nR-1)/modelSize] .= condProbe;
    K[1:trunc(Int, 1+(nR-1)/modelSize)] .= condProbe;
    D = diffSoil * ones(nR);
    D[1:trunc(Int, 1+(nR-1)/modelSize)] .= diffProbe;

    ### initial temperature & boundary
    T = zeros(nR);         # temperature is 0 in the medium
    j = trunc(Int, 1+(nR-1)/modelSize); # j is last point inside the probe
    T[1:j] .= 1; # temperature is 1 in the probe

    return Grid(R, nR, dR, D, K, T, j)
end

function rDiffusion(G, t) #diffusivity of medium, conductivity of medium
    # time step
    Sigma = 0.2; # timestep stability size
    dt = (Sigma*(g.dR^2))/max(g.D[1],g.D[end]);
    nt = t/dt;     # number of time steps
    n = 0;

    un = zeros(g.nR); #initialize a temporary array
    while n < nt   #loop for values of n from 1 to nt
        copy!(un,g.T)  #copy existing values of g.T into un

        # r = 0
        g.T[1] = un[1] + 4*dt*diffSoil*(un[2] - un[1])*(1/(g.dR^2));

        # general interior
        for i in 2:g.nR-1
            delta = dt * ((1/g.R[i])*(1/g.dR^2) *
                     (g.R[i+1]*g.D[i+1]*(un[i+1] - un[i]) -
                     g.R[i-1]*g.D[i-1]*(un[i] - un[i-1])));
            g.T[i] = un[i] + delta;
        end

        # heat flux at boundary
        j = g.j; # find the index of the boundary
        k1 = g.K[1];
        k2 = g.K[end];
        d1 = g.D[1];
        d2 = g.D[end];

        # derivatives at boundary - calculated using linear flux equations
        # dun[j+1] = k1 * (un[j-2] - 4*un[j-1] - un[j+1] + 4*un[j+2])/(2*g.dR*(k1+k2));
        # dTm = k2 * (un[j-2] - 4*un[j-1] + 4*un[j+1] - un[j+2])/(2*g.dR*(k1+k2));
        # d2tp = (un[j-2]*(-k1) + un[j-1]*(4*k1) - un[j+1]*(6*k1+2*k2) + un[j+2]*(3*k1+2*k2)) /
        #         (3 * g.dR^2 * (k1 + k2));
        # d2tm = (un[j-2]*(2*k1+3*k2) - un[j-1]*(2*k1+6*k2) + un[j+1]*(4*k2) - un[j+2]*(k2)) /
        #         (3 * g.dR^2 * (k1 + k2));

        #derivatives at boundary v2 - calculated using radial flux equations
        dTp = (3*un[j-2]*g.dR*k1 - 3*un[j+2]*g.dR*k1 - 12*un[j-1]*g.dR*k1 + 12*un[j+1]*g.dR*k1 + 6*un[j-2]*k1*g.R[j] - 2*un[j+2]*k1*g.R[j] - 4*un[j+2]*k2*g.R[j] - 12*un[j-1]*k1*g.R[j] + 8*un[j+1]*k1*g.R[j] + 4*un[j+1]*k2*g.R[j])/(2*g.dR*(3*g.dR*k1 + 3*g.dR*k2 + 2*k1*g.R[j] - 2*k2*g.R[j]))
        dTm = (3*un[j-2]*g.dR*k2 - 3*un[j+2]*g.dR*k2 - 12*un[j-1]*g.dR*k2 + 12*un[j+1]*g.dR*k2 - 4*un[j-2]*k1*g.R[j] - 2*un[j-2]*k2*g.R[j] + 6*un[j+2]*k2*g.R[j] + 4*un[j-1]*k1*g.R[j] + 8*un[j-1]*k2*g.R[j] - 12*un[j+1]*k2*g.R[j])/(2*g.dR*(3*g.dR*k1 + 3*g.dR*k2 + 2*k1*g.R[j] - 2*k2*g.R[j]))
        d2tp = -(un[j-2]*g.dR*k1 - 3*un[j+2]*g.dR*k1 - 2*un[j+2]*g.dR*k2 - 4*un[j-1]*g.dR*k1 + 6*un[j+1]*g.dR*k1 + 2*un[j+1]*g.dR*k2 + 2*un[j-2]*k1*g.R[j] - 2*un[j+2]*k1*g.R[j] - 4*un[j-1]*k1*g.R[j] + 4*un[j+1]*k1*g.R[j])/(g.dR^2*(3*g.dR*k1 + 3*g.dR*k2 + 2*k1*g.R[j] - 2*k2*g.R[j]))
        d2tm = (2*un[j-2]*g.dR*k1 + 3*un[j-2]*g.dR*k2 - un[j+2]*g.dR*k2 - 2*un[j-1]*g.dR*k1 - 6*un[j-1]*g.dR*k2 + 4*un[j+1]*g.dR*k2 - 2*un[j-2]*k2*g.R[j] + 2*un[j+2]*k2*g.R[j] + 4*un[j-1]*k2*g.R[j] - 4*un[j+1]*k2*g.R[j])/(g.dR^2*(3*g.dR*k1 + 3*g.dR*k2 + 2*k1*g.R[j] - 2*k2*g.R[j]))



        # temperature change rate at boundary
        dt1 = (d1/g.R[j])*(dTp + g.R[j]*d2tp);
        dt2 = (d2/g.R[j])*(dTm + g.R[j]*d2tm);

        g.T[j] = un[j] + dt*(0.5*dt1 + 0.5*dt2);

        # outer boundary of the model - make it same as next to the edge
        ghost = 2*un[end-1] - un[end-2];
        g.T[end] = un[end] + dt * ((1/g.R[end])*(1/g.dR^2) *
                 ((g.dR + g.R[end])*g.D[end]*(ghost - un[end]) -
                 g.R[end-1]*g.D[end-1]*(un[end] - un[end-1])));

        n += 1;
    end
    return G;
end

function rDiffusionOld(gNaive, t) #diffusivity of medium, conductivity of medium
    # time step
    Sigma = 0.2; # timestep stability size
    dt = (Sigma*(gNaive.dR^2))/max(diffProbe, diffSoil)
    nt = t/dt;     # number of time steps
    n = 0;


    un = zeros(gNaive.nR); #initialize a temporary array
    while n < nt   #loop for values of n from 1 to nt
        copy!(un,gNaive.T)  #copy existing values of gNaive.T into un

        # r = 0
        gNaive.T[1] = un[1] + 4*dt*diffSoil*(un[2] - un[1])*(1/(gNaive.dR^2));

        # general interior
        for i in 2:gNaive.nR-1
            delta = dt * ((1/gNaive.R[i])*(1/gNaive.dR^2) *
                     (gNaive.R[i+1]*gNaive.D[i+1]*(un[i+1] - un[i]) -
                     gNaive.R[i-1]*gNaive.D[i-1]*(un[i] - un[i-1])));
            gNaive.T[i] = un[i] + delta;
        end

        # outer boundary of the model - make it same as next to the edge
        ghost = 2*un[end-1] - un[end-2];
        gNaive.T[end] = un[end] + dt * ((1/gNaive.R[end])*(1/gNaive.dR^2) *
                 ((gNaive.dR + gNaive.R[end])*gNaive.D[end]*(ghost - un[end]) -
                 gNaive.R[end-1]*gNaive.D[end-1]*(un[end] - un[end-1])));

        n += 1;
    end
    return gNaive;
end

####### MATERIAL PROPERTIES #########
diffProbe = 9e-8; #M^2/s
condProbe = 0.185; # W/m-K
diffSoil = 3e-7;
condSoil = 0.4;
#####################################

g = setup(diffProbe, condProbe, diffSoil, condSoil)  # get grid of soil properties and coordinates
gNaive = setup(diffProbe, condProbe, diffSoil, condSoil);

nSteps = 300; # number of timesteps for the gif
@gif for t ∈ 1:nSteps;
    println(g.T[21])
    g != rDiffusion(g, 0.5); # medium D, medium k, time 10s
    gNaive != rDiffusionOld(gNaive, 0.5);
    plot(gNaive.R, gNaive.T, color = :black, label = "Naive Boundary")
    plot!(g.R, g.T, color = :red, label = "Heat Conserving Boundary")
    plot!([0.00476], seriestype="vline", color = :gray, label="")
    xlabel!("Radius [m]")
    ylabel!("Temperature [c]")
    ylims!(-0.5,1)
    xlims!(0,0.01)
    title!(string("Time = ", 0.5*t, "s"))
end
