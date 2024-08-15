function lagoon_bed_shear_stress(bL,ZL,g,ko,U)

    Hs = wave_height(bL,ZL,g,U)
    Tp = wave_period(bL,ZL,g,U)
    k = wave_number(1/Tp,ZL,g)
    Um = (pi*Hs/Tp/sinh(k*ZL))
    aw = Tp*Um/pi
    fw = 0.4*(aw/ko)^-0.75
    τ = 1/2*1020*fw*Um^2

    return τ

end