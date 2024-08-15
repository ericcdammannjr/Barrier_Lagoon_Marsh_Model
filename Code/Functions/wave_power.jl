function wave_power(bL,Zs,g,U)

    Hs = wave_height(bL,Zs,g,U)
    Tp = wave_period(bL,Zs,g,U)
    k = wave_number(1/Tp,Zs,g)
    cg = 2*pi/k/Tp*0.5*(1+2*k*Zs/(sinh(2*k*Zs)))
    WP = cg*9800/16*abs(Hs).^2

    return WP

end