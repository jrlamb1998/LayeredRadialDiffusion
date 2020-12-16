using Plots
mutable struct Grid
    R   # radius array
    NR  # number of grid points
    DR  # grid spacing
    D   # diffusivity array
    K   # conductivity array
    T   # temperature array
    j   # boundary index
end

function setup(D1, K1, D2, KM)
    RProbe = 0.00476; # radius of probe in cm
    NR = 101;     # number of grid points

    ModelSize = 5;   # model size vs probe radius

    # grid spacing
    RModel = RProbe * ModelSize;
    DR = RModel/(NR-1);
    R = LinRange(0, RModel, NR); # set the radius coordinates


    # D and K
    K = KM * ones(NR);
    #K[1:1+(NR-1)/ModelSize] .= K1;
    K[1:trunc(Int, 1+(NR-1)/ModelSize)] .= K1;
    D = D2 * ones(NR);
    D[1:trunc(Int, 1+(NR-1)/ModelSize)] .= D1;

    # initial temperature
    T = zeros(NR);         # temperature is 0 in the medium
    T[1:trunc(Int, 1+(NR-1)/ModelSize)] .= 1; # temperature is 1 in the probe

    j = trunc(Int, 1+(NR-1)/ModelSize);

    return Grid(R, NR, DR, D, K, T, j)
end

function rDiffusion(G, t) #diffusivity of medium, conductivity of medium
    # time step
    Sigma = 0.2; # timestep stability size
    dt = (Sigma*(G.DR^2))/max(D1, D2)
    nt = t/dt;     # number of time steps
    n = 0;


    un = zeros(G.NR); #initialize a temporary array
    while n < nt   #loop for values of n from 1 to nt
        copy!(un,G.T)  #copy existing values of G.T into un

        # r = 0
        G.T[1] = un[1] + 4*dt*D2*(un[2] - un[1])*(1/(G.DR^2));

        # general interior
        for i in 2:G.NR-1
            delta = dt * ((1/G.R[i])*(1/G.DR^2) *
                     (G.R[i+1]*G.D[i+1]*(un[i+1] - un[i]) -
                     G.R[i-1]*G.D[i-1]*(un[i] - un[i-1])));
            G.T[i] = un[i] + delta;
        end

        # heat flux at boundary
        j = G.j; # find the index of the boundary
        k1 = G.K[1];
        k2 = G.K[end];
        d1 = G.D[1];
        d2 = G.D[end];

        # derivatives at boundary - calculated using linear flux equations
        # dun[j+1] = k1 * (un[j-2] - 4*un[j-1] - un[j+1] + 4*un[j+2])/(2*G.DR*(k1+k2));
        # dTm = k2 * (un[j-2] - 4*un[j-1] + 4*un[j+1] - un[j+2])/(2*G.DR*(k1+k2));
        # d2tp = (un[j-2]*(-k1) + un[j-1]*(4*k1) - un[j+1]*(6*k1+2*k2) + un[j+2]*(3*k1+2*k2)) /
        #         (3 * G.DR^2 * (k1 + k2));
        # d2tm = (un[j-2]*(2*k1+3*k2) - un[j-1]*(2*k1+6*k2) + un[j+1]*(4*k2) - un[j+2]*(k2)) /
        #         (3 * G.DR^2 * (k1 + k2));

        #derivatives at boundary v2 - calculated using radial flux equations
        dTp = (3*un[j-2]*G.DR*k1 - 3*un[j+2]*G.DR*k1 - 12*un[j-1]*G.DR*k1 + 12*un[j+1]*G.DR*k1 + 6*un[j-2]*k1*G.R[j] - 2*un[j+2]*k1*G.R[j] - 4*un[j+2]*k2*G.R[j] - 12*un[j-1]*k1*G.R[j] + 8*un[j+1]*k1*G.R[j] + 4*un[j+1]*k2*G.R[j])/(2*G.DR*(3*G.DR*k1 + 3*G.DR*k2 + 2*k1*G.R[j] - 2*k2*G.R[j]))
        dTm = (3*un[j-2]*G.DR*k2 - 3*un[j+2]*G.DR*k2 - 12*un[j-1]*G.DR*k2 + 12*un[j+1]*G.DR*k2 - 4*un[j-2]*k1*G.R[j] - 2*un[j-2]*k2*G.R[j] + 6*un[j+2]*k2*G.R[j] + 4*un[j-1]*k1*G.R[j] + 8*un[j-1]*k2*G.R[j] - 12*un[j+1]*k2*G.R[j])/(2*G.DR*(3*G.DR*k1 + 3*G.DR*k2 + 2*k1*G.R[j] - 2*k2*G.R[j]))
        d2tp = -(un[j-2]*G.DR*k1 - 3*un[j+2]*G.DR*k1 - 2*un[j+2]*G.DR*k2 - 4*un[j-1]*G.DR*k1 + 6*un[j+1]*G.DR*k1 + 2*un[j+1]*G.DR*k2 + 2*un[j-2]*k1*G.R[j] - 2*un[j+2]*k1*G.R[j] - 4*un[j-1]*k1*G.R[j] + 4*un[j+1]*k1*G.R[j])/(G.DR^2*(3*G.DR*k1 + 3*G.DR*k2 + 2*k1*G.R[j] - 2*k2*G.R[j]))
        d2tm = (2*un[j-2]*G.DR*k1 + 3*un[j-2]*G.DR*k2 - un[j+2]*G.DR*k2 - 2*un[j-1]*G.DR*k1 - 6*un[j-1]*G.DR*k2 + 4*un[j+1]*G.DR*k2 - 2*un[j-2]*k2*G.R[j] + 2*un[j+2]*k2*G.R[j] + 4*un[j-1]*k2*G.R[j] - 4*un[j+1]*k2*G.R[j])/(G.DR^2*(3*G.DR*k1 + 3*G.DR*k2 + 2*k1*G.R[j] - 2*k2*G.R[j]))



        # temperature change rate at boundary
        dt1 = (d1/G.R[j])*(dTp + G.R[j]*d2tp);
        dt2 = (d2/G.R[j])*(dTm + G.R[j]*d2tm);

        G.T[j] = un[j] + dt*(0.5*dt1 + 0.5*dt2);

        # outer boundary of the model - make it same as next to the edge
        ghost = 2*un[end-1] - un[end-2];
        G.T[end] = un[end] + dt * ((1/G.R[end])*(1/G.DR^2) *
                 ((G.DR + G.R[end])*G.D[end]*(ghost - un[end]) -
                 G.R[end-1]*G.D[end-1]*(un[end] - un[end-1])));

        n += 1;
    end
    return G;
end

function rDiffusionOld(GOld, t) #diffusivity of medium, conductivity of medium
    # time step
    Sigma = 0.2; # timestep stability size
    dt = (Sigma*(GOld.DR^2))/max(D1, D2)
    nt = t/dt;     # number of time steps
    n = 0;


    un = zeros(GOld.NR); #initialize a temporary array
    while n < nt   #loop for values of n from 1 to nt
        copy!(un,GOld.T)  #copy existing values of GOld.T into un

        # r = 0
        GOld.T[1] = un[1] + 4*dt*D2*(un[2] - un[1])*(1/(GOld.DR^2));

        # general interior
        for i in 2:GOld.NR-1
            delta = dt * ((1/GOld.R[i])*(1/GOld.DR^2) *
                     (GOld.R[i+1]*GOld.D[i+1]*(un[i+1] - un[i]) -
                     GOld.R[i-1]*GOld.D[i-1]*(un[i] - un[i-1])));
            GOld.T[i] = un[i] + delta;
        end

        # heat flux at boundary
        j = GOld.j; # find the index of the boundary
        k1 = GOld.K[1];
        k2 = GOld.K[end];
        d1 = GOld.D[1];
        d2 = GOld.D[end];

        # outer boundary of the model - make it same as next to the edge
        ghost = 2*un[end-1] - un[end-2];
        GOld.T[end] = un[end] + dt * ((1/GOld.R[end])*(1/GOld.DR^2) *
                 ((GOld.DR + GOld.R[end])*GOld.D[end]*(ghost - un[end]) -
                 GOld.R[end-1]*GOld.D[end-1]*(un[end] - un[end-1])));

        n += 1;
    end
    return GOld;
end

D1 = 9e-8; #M^2/s
K1 = 0.185; # W/m-K
D2 = 3e-7;
K2 = 0.4;
G = setup(D1, K1, D2, K2)  # get grid of soil properties and coordinates
GOld = setup(D1, K1, D2, K2);

@gif for t ∈ 1:300; # do 100 time steps
    println(G.T[21])
    G != rDiffusion(G, 0.5); # medium D, medium k, time 10s
    GOld != rDiffusionOld(GOld, 0.5);
    plot(GOld.R, GOld.T, color = :black, label = "Naive Boundary")
    plot!(G.R, G.T, color = :red, label = "Heat Conserving Boundary")
    plot!([0.00476], seriestype="vline", color = :gray, label="")
    xlabel!("Radius [m]")
    ylabel!("Temperature [c]")
    ylims!(-0.5,1)
    xlims!(0,0.01)
    title!(string("Time = ", 0.5*t, "s"))
end
