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

α = zeros(1,n+1); α[1,1] = α0; 
W = zeros(1,n+1); W[1,1] = W0; 
H = zeros(1,n+1); H[1,1] = H0; Hdot = zeros(1,n)
xt = zeros(1,n+1); xt[1,1] = xt0; xtdot = zeros(1,n)
xs = zeros(1,n+1); xs[1,1] = xs0; xsdot = zeros(1,n)
xb = zeros(1,n+1); xb[1,1] = xb0; xbdot = zeros(1,n)
Z = zeros(1,n+1); Z[1,1] = Z0
bbm = zeros(1,n+1); bbm[1,1] = bbm0
bL = zeros(1,n+1); bL[1,1] = bL0
bim = zeros(1,n+1); bim[1,1] = bim0
zm = zeros(1,n+1); zm[1,1] = zm0; zmdot = zeros(1,n)
zL = zeros(1,n+1); zL[1,1] = zL0; zLdot = zeros(1,n)
xbm = zeros(1,n+1); xbm[1,1] = xbm0; xbmdot = zeros(1,n)
xim = zeros(1,n+1); xim[1,1] = xim0; ximdot = zeros(1,n)
xmm = zeros(1,n+1); xmm[1,1] = xmm0; xmmdot = zeros(1,n)
Qsf = zeros(1,n)
Qow_H = zeros(1,n)
Qow_Bl = zeros(1,n)
Qow_Bm = zeros(1,n)
Fc = zeros(1,n)
Fm = zeros(1,n)
O = zeros(1,n)

# Main Code

for i in 1:n

    # Deficit Volume Calculations 

    ϕ = min(1,bbm[1,i]/bbmc)
    Vd_B = max(0,(We-W[1,i])*(H[1,i]+ϕ*(zm[1,i]-r/2)+(1-ϕ)*(zL[1,i]-r/2)))
    Vd_H = max(0,(He-H[1,i])*(W[1,i]))
    Vd = Vd_B+Vd_H

    # Overwash Calculations

    Qsf[1,i] = K*(αe-α[1,i])

    if Vd<Vd_max

        Qow_H[1,i] = Qow_max*Vd_H/Vd_max
        Qow_B = Qow_max*Vd_B/Vd_max

    else 

        Qow_H[1,i] = Qow_max*Vd_H/Vd
        Qow_B = Qow_max*Vd_B/Vd

    end

    Qow = Qow_H[1,i]+Qow_B
    Qow_Bl[1,i] = (1-ϕ)*Qow_B
    Qow_Bm[1,i] = ϕ*Qow_B

    # Barrier Equations

    xtdot[1,i] = 4*Qsf[1,i]*((H[1,i]+Dt)/(Dt*(2*H[1,i]+Dt)))+2*zdot/α[1,i]
    xsdot[1,i] = 2*Qow/(2*H[1,i]+Dt)-4*Qsf[1,i]*((H[1,i]+Dt)/(2*H[1,i]+Dt)^2)
    xbdot[1,i] = Qow_Bm[1,i] /(H[1,i]+zm[1,i]-r/2)
    Hdot[1,i] = Qow_H[1,i]/W[1,i]-zdot

    # Marsh-Lagoon Dynamics

    # Lagoon 

    ZL = (zL[1,i]+(zL[1,i]-min(r,zL[1,i])))/2
    τL = lagoon_bed_shear_stress(bL[1,i],ZL,g,ko,U)
    SL = max((τL-τcr)/τcr,0)
    Cr = ρ*λ*SL/(1+SL*λ)
    Fc[1,i] = (Cr-Co)*min(r,zL[1,i])/P/ρ

    # Marshes

    Zm = (zm[1,i]+(zm[1,i]-min(r,zm[1,i])))/2
    B = Bpeak*(Dmax-zm[1,i])*(zm[1,i]-Dmin)/AA

    if B <= 1*ℯ^-3

        B = 0

    end

    Bfrac = B/Bpeak
    AMC = 180*νGp*B
    Rref = AMC*χref
    O[1,i] = 1/por*(Rref/ρo)

    if Zm > 1*ℯ^-4

        τm = marsh_bed_shear_stress(bL,Zm,Bfrac,g,ko,U)

    else

        τm = 0 

    end

    Sm = max((τm-τcr)/τcr,0)
    Cm = ρ*λ*Sm/(1+Sm*λ)
    Fm[1,i] = (Cr-Cm)*min(r,zm[1,i])/P/ρ

    # Marsh-Lagoon Edges

    zs = zm[1,i]+(zL[1,i]-zm[1,i])*(1-exp(-x*0.1/zL[1,i]))
    Zs = (zs+(zs-min(r,zs)))/2
    WP = wave_power(bL[1,i],Zs,g,U)
    E = ke*WP/(zs-zm[1,i])-ka*Cr*ws/ρ
   
    # Moving Boundaries
    
    zmdot[1,i] = -Fm[1,i]-O[1,i]+zdot
    zLdot[1,i] = -2*E*(zL[1,i]-zm[1,i])/bL[1,i]+Fm[1,i]*(bbm[1,i]+bim[1,i])/bL[1,i]+Fc[1,i]+zdot
    xbmdot[1,i] = -E+Qow_Bl[1,i] /(zL[1,i]-zm[1,i])
    ximdot[1,i] = E
    xmmdot[1,i] = (zdot-zmdot[1,i])/β

    # Foward Euler Method Implementation

    xt[1,i+1] = xt[1,i]+dt*xtdot[1,i]
    xs[1,i+1] = xs[1,i]+dt*xsdot[1,i]
    xb[1,i+1] = xb[1,i]+dt*xbdot[1,i]
    H[1,i+1] = H[1,i]+dt*Hdot[1,i]
    Z[1,i+1] = Z[1,i]+dt*zdot
    zm[1,i+1] = zm[1,i]+dt*zmdot[1,i]
    zL[1,i+1] = zL[1,i]+dt*zLdot[1,i]
    xbm[1,i+1] = xbm[1,i]+dt*xbmdot[1,i]
    xim[1,i+1] = xim[1,i]+dt*ximdot[1,i]
    xmm[1,i+1] = xmm[1,i]+dt*xmmdot[1,i]
    
    # Check for Complete Marsh Erosion

    if xbm[1,i+1] <= xb[1,i+1]

        xbm[1,i+1] = xb[1,i+1]

    end

    if xmm[1,i+1] <= xim[1,i+1]

        xim[1,i+1] = xmm[1,i+1]

    end

    # Check for Drowning or Marsh Filling

    if zm[1,i+1] > Dmax

        break

    end

    if zL[1,i+1] <= zm[1,i+1]

        xbm[1,i+1] = xb[1,i+1]
        xim[1,i+1] = xmm[1,i+1]

    end

    # Updating Additional Values
    
    α[1,i+1] = Dt/(xs[1,i+1]-xt[1,i+1])
    W[1,i+1] = xb[1,i+1]-xs[1,i+1]
    bbm[1,i+1] = xbm[1,i+1]-xb[1,i+1]
    bL[1,i+1] = xim[1,i+1]-xbm[1,i+1]
    bim[1,i+1] = xmm[1,i+1]-xim[1,i+1]

end

plot(t,[bbm[1,:] bL[1,:]/10^3 W[1,:] zL[1,:]], 

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