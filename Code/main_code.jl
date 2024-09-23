using Plots
pyplot()

include("Functions/lagoon_bed_shear_stress.jl")
include("Functions/marsh_bed_shear_stress.jl")
include("Functions/wave_height.jl")
include("Functions/wave_number.jl")
include("Functions/wave_period.jl")
include("Functions/wave_power.jl")

# Barrier Dynamics Parameters

αe = 0.02
bbmc = 1000
Dt = 10
He = 2
K = 2000
Qow_max = 100
We = 800
zdot = 0.005
Vd_max = He*We

# Marsh-Lagoon Dynamics Parameters

β = 10^10
Bpeak = 2.5
χref = 0.158
Co = 0.03
Dmin = 0 
g = 9.80171
ka = 2 
ke = 0.15
ko = 0.001
λ = 0.0001
νGp = 0.0138
P = 12.5/(24*365)
por = 1000/2650
r = 1.4
ρ = 1000
ρo = 1000
τcr = 0.1
U = 10
ws = 0.5*10^-3*(60*60*24*365)
x = 10
Dmax = 0.7167*r-0.0483
AA = 0.25*(Dmax-Dmin)^2

# Carbon Storage Parameters



# Barrier Intital Conditions

α0 = αe
W0 = We
H0 = He
xt0 = 0
xs0 = Dt/αe
xb0 = Dt/αe+We
Z0 = Dt

# Marsh-Lagoon Initial Conditions

bbm0 = 1000
bL0 = 10000
bim0 = 2000
zm0 = Dmax/2
zL0 = 2
xbm0 = xb0+bbm0
xim0 = xb0+bbm0+bL0
xmm0 = xb0+bbm0+bL0+bim0-(zm0-r/2)/β

# Computational Parameters

t0 = 0
tmax = 1000
n = 100000
dt = (tmax-t0)/n
t = t0:dt:tmax

# Array Preallocation

α = zeros(n+1); α[1] = α0; 
W = zeros(n+1); W[1] = W0; 
H = zeros(n+1); H[1] = H0; Hdot = zeros(n)
xt = zeros(n+1); xt[1] = xt0; xtdot = zeros(n)
xs = zeros(n+1); xs[1] = xs0; xsdot = zeros(n)
xb = zeros(n+1); xb[1] = xb0; xbdot = zeros(n)
Z = zeros(n+1); Z[1] = Z0
bbm = zeros(n+1); bbm[1] = bbm0
bL = zeros(n+1); bL[1] = bL0
bim = zeros(n+1); bim[1] = bim0
zm = zeros(n+1); zm[1] = zm0; zmdot = zeros(n)
zL = zeros(n+1); zL[1] = zL0; zLdot = zeros(n)
xbm = zeros(n+1); xbm[1] = xbm0; xbmdot = zeros(n)
xim = zeros(n+1); xim[1] = xim0; ximdot = zeros(n)
xmm = zeros(n+1); xmm[1] = xmm0; xmmdot = zeros(n)
Qsf = zeros(n)
Qow_H = zeros(n)
Qow_Bl = zeros(n)
Qow_Bm = zeros(n)
Fc = zeros(n)
Fm = zeros(n)
O = zeros(n)

# Main Code

for i in 1:n

    # Deficit Volume Calculations 

    ϕ = min(1,bbm[i]/bbmc)
    Vd_B = max(0,(We-W[i])*(H[i]+ϕ*(zm[i]-r/2)+(1-ϕ)*(zL[i]-r/2)))
    Vd_H = max(0,(He-H[i])*(W[i]))
    Vd = Vd_B+Vd_H

    # Overwash Calculations

    Qsf[i] = K*(αe-α[i])

    if Vd<Vd_max

        Qow_H[i] = Qow_max*Vd_H/Vd_max
        Qow_B = Qow_max*Vd_B/Vd_max

    else 

        Qow_H[i] = Qow_max*Vd_H/Vd
        Qow_B = Qow_max*Vd_B/Vd

    end

    Qow = Qow_H[i]+Qow_B
    Qow_Bl[i] = (1-ϕ)*Qow_B
    Qow_Bm[i] = ϕ*Qow_B

    # Barrier Equations

    xtdot[i] = 4*Qsf[i]*((H[i]+Dt)/(Dt*(2*H[i]+Dt)))+2*zdot/α[i]
    xsdot[i] = 2*Qow/(2*H[i]+Dt)-4*Qsf[i]*((H[i]+Dt)/(2*H[i]+Dt)^2)
    xbdot[i] = Qow_Bm[i] /(H[i]+zm[i]-r/2)
    Hdot[i] = Qow_H[i]/W[i]-zdot

    # Marsh-Lagoon Dynamics

    # Lagoon 

    ZL = (zL[i]+(zL[i]-min(r,zL[i])))/2
    τL = lagoon_bed_shear_stress(bL[i],ZL,g,ko,U)
    SL = max((τL-τcr)/τcr,0)
    Cr = ρ*λ*SL/(1+SL*λ)
    Fc[i] = (Cr-Co)*min(r,zL[i])/P/ρ

    # Marshes

    Zm = (zm[i]+(zm[i]-min(r,zm[i])))/2
    B = Bpeak*(Dmax-zm[i])*(zm[i]-Dmin)/AA

    if B <= 1*ℯ^-3

        B = 0

    end

    Bfrac = B/Bpeak
    AMC = 180*νGp*B
    Rref = AMC*χref
    O[i] = 1/por*(Rref/ρo)

    if Zm > 1*ℯ^-4

        τm = marsh_bed_shear_stress(bL,Zm,Bfrac,g,ko,U)

    else

        τm = 0 

    end

    Sm = max((τm-τcr)/τcr,0)
    Cm = ρ*λ*Sm/(1+Sm*λ)
    Fm[i] = (Cr-Cm)*min(r,zm[i])/P/ρ

    # Marsh-Lagoon Edges

    zs = zm[i]+(zL[i]-zm[i])*(1-exp(-x*0.1/zL[i]))
    Zs = (zs+(zs-min(r,zs)))/2
    WP = wave_power(bL[i],Zs,g,U)
    E = ke*WP/(zs-zm[i])-ka*Cr*ws/ρ
   
    # Moving Boundaries
    
    zmdot[i] = -Fm[i]-O[i]+zdot
    zLdot[i] = -2*E*(zL[i]-zm[i])/bL[i]+Fm[i]*(bbm[i]+bim[i])/bL[i]+Fc[i]+zdot
    xbmdot[i] = -E+Qow_Bl[i] /(zL[i]-zm[i])
    ximdot[i] = E
    xmmdot[i] = (zdot-zmdot[i])/β

    # Foward Euler Method Implementation

    xt[i+1] = xt[i]+dt*xtdot[i]
    xs[i+1] = xs[i]+dt*xsdot[i]
    xb[i+1] = xb[i]+dt*xbdot[i]
    H[i+1] = H[i]+dt*Hdot[i]
    Z[i+1] = Z[i]+dt*zdot
    zm[i+1] = zm[i]+dt*zmdot[i]
    zL[i+1] = zL[i]+dt*zLdot[i]
    xbm[i+1] = xbm[i]+dt*xbmdot[i]
    xim[i+1] = xim[i]+dt*ximdot[i]
    xmm[i+1] = xmm[i]+dt*xmmdot[i]
    
    # Check for Complete Marsh Erosion

    if xbm[i+1] <= xb[i+1]

        xbm[i+1] = xb[i+1]

    end

    if xmm[i+1] <= xim[i+1]

        xim[i+1] = xmm[i+1]

    end

    # Check for Drowning or Marsh Filling

    if zm[i+1] > Dmax

        break

    end

    if zL[i+1] <= zm[i+1]

        xbm[i+1] = xb[i+1]
        xim[i+1] = xmm[i+1]

    end

    # Updating Additional Values
    
    α[i+1] = Dt/(xs[i+1]-xt[i+1])
    W[i+1] = xb[i+1]-xs[i+1]
    bbm[i+1] = xbm[i+1]-xb[i+1]
    bL[i+1] = xim[i+1]-xbm[i+1]
    bim[i+1] = xmm[i+1]-xim[i+1]

end

plot(t,[bbm[:] bL[:]/10^3 W[:] zL[:]], 

    layout = (2,2),
    size = [1200 800],
    grid = false,
    legend = false,
    title = ["Backbarrier Marsh Width (bbm)" "Lagoon Width (bL)" "Barrier Width (W)" "Lagoon Depth (zL)"],
    xlabel = "Time (yrs)",
    ylabel = ["bbm (m)" "bl (km)" "W (m)" "zL (m)"],
    lw = 3,
    lc = :black
    
)