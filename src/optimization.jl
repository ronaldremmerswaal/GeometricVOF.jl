function brent_min(f_and_df::Function, x0::Real; xtol::Real=√(eps(typeof(x0))), maxiters::Int=25, step_max::Real=floatmax(typeof(x0)))
    # 1: find a bracket using Newton's method
    # 2: apply Brent's method to the derivative to find the minimum

    x = x0
    f, df = f_and_df(x)
    if df == 0 return x end

    x_prev, df_prev = x, df
    it = 1

    converged = false

    while df * df_prev > 0 && it < maxiters

        x_prev, df_prev = x, df
        step = - 2f / df

        if abs(step) > step_max
            step = sign(step) * step_max
        end

        x += step

        converged = abs(step) < xtol
        if converged break end

        f, df = f_and_df(x)

        # println("it: $it, x: $x, f: $f, df: $df")
        it += 1
    end

    if converged
        return x
    elseif df * df_prev > 0
        @warn "Failed to find a bracket"
        return x
    elseif it >= maxiters
        @warn "Maximum iterations reached"
        return x
    end

    find_zero(x -> f_and_df(x)[2], (x_prev, x), Roots.Brent(), maxiters=maxiters-it)
end
