function brent_min(f_and_df::Function, x0::Real; xatol::Real=√(eps(typeof(x0))), xrtol::Real=0, maxiters::Int=25, step_max::Real=floatmax(typeof(x0)), verbose::Bool=false)
    # 1: find a bracket using Newton's method
    # 2: apply Brent's method to the derivative to find the minimum

    x = x0
    f, df = f_and_df(x)
    if df == 0 return x end

    x_prev, df_prev = x, df
    it = 1

    if verbose
        @printf("|       Starting iterations brent_min     |\n")
        @printf("| %3s | %9s | %9s | %9s |\n", "it", "x", "f(x)", "df(x)")
    end

    while df * df_prev > 0 && it < maxiters
        if verbose
            @printf("| %3d | %9.6f | %9.2e | %9.2e |\n", it, x, f, df)
        end

        x_prev, df_prev = x, df
        step = - 2f / df

        if abs(step) > step_max
            step = sign(step) * step_max
        end

        x += step

        f, df = f_and_df(x)

        it += 1
    end

    if df * df_prev > 0
        if abs(df) > 100eps(typeof(df)) @warn "brent_min: failed to find a bracket, final x = $x, f = $f, df = $df" end
        return x
    elseif verbose
        @printf("| %3d | %9.6f | %9.2e | %9.2e |\n", it, x, f, df)
        @printf("|        Finished Newton iterations       |\n")
    end

    find_zero(x -> f_and_df(x)[2], (x_prev, x), Roots.Brent(), maxiters=maxiters-it, verbose=verbose, xatol=xatol, xrtol=xrtol)
end

function parabola_roots(A::T, B::T, C::T) where T

    imag = zero(T)
    if A == 0
        if B == 0
            x1 = NaN
            x2 = NaN
        else
            x1 = -C / B
            x2 = x1
        end
    elseif C == 0 then
        x1 = -B / A
        x2 = 0.0
    else

        # Citarduaq's formula is used to allow for small c(1) in a numerically stable way
        Bscaled = B / A
        Cscaled = C / A

        sqrt_max = √floatmax(T)

        if Bscaled > sqrt_max || Bscaled < -sqrt_max
            # Bscaled^2 would overflow, so we let √D ≈ |Bscaled| instead
            x1 = -Cscaled / Bscaled
            x2 = Cscaled / x1
        else
            D = Bscaled^2 - 4 * Cscaled

            if D ≤ 0
                if D < 0
                    imag = sqrt(-D) / 2
                end
                x1 = -Bscaled / 2
                x2 = x1
            else
                x1 = -2 * Cscaled / (Bscaled + sqrt(D) * sign(Bscaled))
                x2 = Cscaled / x1
            end
        end
    end

    return x1, x2, imag
end
