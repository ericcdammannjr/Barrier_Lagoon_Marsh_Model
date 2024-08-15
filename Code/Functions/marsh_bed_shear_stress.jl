function marsh_bed_shear_stress(bL,Zm,Bfrac,g,ko,U)

    if Bfrac == 0

        Hs = wave_height(bL,Zm,g,U)
        Tp = wave_period(bL,Zm,g,U)
        k = wave_number(1/Tp,Zm,g)
        Um = (pi*Hs/Tp/sinh(k*Zm))
        aw = Tp*Um/(2*pi)
        fw = 0.4*(aw/ko)^-0.75
        τ = 1/2*1020*fw*Um^2

    else

        τ = 0

    end 

    return τ

end