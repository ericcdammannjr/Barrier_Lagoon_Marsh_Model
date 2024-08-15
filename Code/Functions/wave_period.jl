function wave_period(bL,z,g,U)

    delta = z*g/U^2
    chi = bL*g./U^2
    ni = 0.133*(tanh(0.331*delta^1.01)*tanh(5.215*10^-4*chi^0.73/tanh(0.331*delta^1.01)))^-0.37
    Tp = U/ni/g

    return Tp

end