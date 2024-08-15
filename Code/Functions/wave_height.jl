function wave_height(bL,z,g,U)

    delta = z*g/U^2
    chi = bL*g./U^2
    epsilon = 3.64*10^-3*(tanh(0.493*delta^0.75)*tanh(3.13*10^-3*chi^0.57/tanh(0.493*delta^0.75)))^1.74
    Hs = 4*sqrt(U^4*epsilon/g^2)
    
    return Hs

end