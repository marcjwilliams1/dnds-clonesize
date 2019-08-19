using Distances
using Optim
using DataFrames
using RCall
using Statistics
using Distributions
using Optim, NLSolversBase, Random
using LinearAlgebra: diag

function myoptimfunc(dnds, A, Δ, rλ; t = 10.0, λ = 1.0, Amin = 0.01, ρ = 100.0)
    euclidean(cumulativednds(Amin, A, ρ; t = t, r = rλ, λ = λ, Δ = Δ), dnds)
end
mutable struct dndsfit
    Δ
    rλ
    DF
    plot
end

using Bridge
#myexpint(x) = map(val -> -Bridge.expint(-val), x)
function myexpint(x)
    if length(x) == 1
        if x > 0.0
            x = -1.0
        end
        x = -Bridge.expint(-x)
        return x
    else
        ei = zeros(Float64, length(x))
        for i in 1:length(x)
            if x[i] > 0.0
                x[i] = -1.0
            end
            ei[i] = -Bridge.expint(-x[i])
        end
        return ei
    end
end

expintdiff(N, Amax, Amin, ρ) = (myexpint(-(Amax .* ρ) ./ N) .- myexpint(-(Amin .* ρ) / N))
EulerMaclaurin(N, Amax, Amin, ρ) = ((exp.(-(Amax .* ρ) ./ N) ./ (Amax .* ρ))
                .+ (exp.(-Amin .* ρ / N) ./ (Amin .* ρ))) ./ 2


function Cnintervaltheory(Amax, Amin, ρ; t = 10.0, r = 0.5, λ = 1.0, Δ = 0.0, μ = 0.001, n0 = 1)
    if Δ > 0.0
        Nave = n0 / (1 - ((r * λ * t) / (1 + r * λ * t))^n0)
        Nt = ((1 + Δ) * exp.(2 * r * λ * t * Δ) - (1 - Δ)) / (2 * Δ)
    else
        Nave = n0 / (1 - ((r * λ * t) / (1 + r * λ * t))^n0)
        Nt = 1 + r * λ * t
    end
    corr = EulerMaclaurin(Nt, Amax, Amin, ρ)
    Cnint = μ * Nave .* (1 / (1 + Δ)) .* ((myexpint(-(Amax .* ρ) ./ Nt) .- myexpint(-(Amin .* ρ) / Nt)) .+ corr)
    return Cnint
end

function cumulativednds(Amin, A, ρ; t = 10.0, r = 0.5, λ = 1.0, Δ = 0.00001)
    if Δ == 0.0
        Δ = 0.00001
    end
    Nt1 = ((1 + Δ) * exp.(2 * r * λ * t * Δ) - (1 - Δ)) / (2 * Δ)
    Nt2 = 1 + r * λ * t
    cdnds = (1 / (1 + Δ)) .* (expintdiff(Nt1, A, Amin, ρ) .+ EulerMaclaurin(Nt1, A, Amin, ρ) ) ./
    (expintdiff(Nt2, A, Amin, ρ) .+ EulerMaclaurin(Nt2, A, Amin, ρ))
    return cdnds
end

function optimizationresults(dnds, A; t = 10.0, λ = 1.0, Amin = 0.01, ρ = 100.0)
    res = optimize(x -> myoptimfunc(dnds, A, x[1], x[2]; t = t, λ = λ, Amin = Amin, ρ = ρ), [-1.0, 0.0001], [1.0, 100.0], [0.1, 0.1])
    deltares, lambdares = Optim.minimizer(res)
    DF = DataFrame(dnds = dnds, A = A,
    dndsfit = cumulativednds(Amin, A, ρ; t = t, r = lambdares, λ = λ, Δ = deltares),
    deltafit = deltares,
    lambdarfit = lambdares)

    residuals = DF[:dnds] .- DF[:dndsfit]
    SSres = sum(residuals.^2)
    SStot = sum( (DF[:dnds] .- mean(DF[:dnds])).^2)
    rsq = 1 - (SSres / SStot)
    DF[:rsq] = rsq
    @rput DF
    @rput deltares
    @rput lambdares
    x = R"""
    library(ggplot2)
    library(cowplot)
    g <- ggplot(DF, aes(x = A, y = dnds)) +
    geom_point(alpha = 0.5, col = "plum4") +
    geom_line(aes(y = dndsfit), alpha = 0.7, col = "firebrick") +
    xlab("Clone area") +
    ylab("Interval dN/dS") +
    ggtitle(paste0("Δ = ", round(deltares, 2), ", λr = ", round(lambdares, 2)))
    """
    return dndsfit(deltares, lambdares, DF, x)
end

function LLfunc(dnds, A, Δ, rλ, sigma; t = 10.0, λ = 1.0, Amin = 0.01, ρ = 100.0)
    #print(c(A, s))
    R = dnds - cumulativednds(Amin, A, ρ; t = t, r = rλ, λ = λ, Δ = Δ)
    R = sum(logpdf.(Normal(0, sigma), R))
    -sum(R)
end

function LLoptimizationresults(dnds, A; t = 10.0, λ = 1.0, Amin = 0.01, ρ = 100.0)
    function myLLfunc(x)
        LLfunc(dnds, A, x[1], x[2], x[3]; t = t, λ = λ, Amin = Amin, ρ = ρ)
    end
    res = optimize(myLLfunc, [-1.0, 0.01, 0.001], [1.0, 100.0, 1.0], [0.1, 0.1, 0.1])
    #deltares, lambdares = Optim.minimizer(res)
    params = Optim.minimizer(res)
    func = TwiceDifferentiable(myLLfunc, ones(3))
    numerical_hessian = hessian!(func, params)
    if (sum(numerical_hessian.==0.0) > 5 )
        CI = zeros(Float64, 2)
    else
        var_cov_matrix = inv(numerical_hessian)
        CI = 1.96 .* sqrt.(abs.(diag(var_cov_matrix)))
    end

    deltares = params[1]
    lambdares = params[2]

    dndsmodel = hcat(cumulativednds(Amin, A, ρ; t = t, r = lambdares - CI[2], λ = λ, Δ = deltares - CI[1]),
    cumulativednds(Amin, A, ρ; t = t, r = lambdares + CI[2], λ = λ, Δ = deltares - CI[1]),
    cumulativednds(Amin, A, ρ; t = t, r = lambdares - CI[2], λ = λ, Δ = deltares + CI[1]),
    cumulativednds(Amin, A, ρ; t = t, r = lambdares + CI[2], λ = λ, Δ = deltares + CI[1]))

    mindnds = minimum(dndsmodel, dims = 2)[:]
    maxdnds = maximum(dndsmodel, dims = 2)[:]

    #sometimes the CI calculation fails due to inability to calculate the hessian
    # this is often due to a poor fit, if this is the case we'll just set the CI to the dnds

    if sum(mindnds .< 10^-5) > 2
        @warn "Cannot compute CI intervals, returning MLE for CIs"
        mindnds = cumulativednds(Amin, A, ρ; t = t, r = lambdares, λ = λ, Δ = deltares)
    end

    if sum(maxdnds .> 10^6) > 2
        maxdnds = cumulativednds(Amin, A, ρ; t = t, r = lambdares, λ = λ, Δ = deltares)
    end

    DF = DataFrame(dnds = dnds, A = A,
    dndsfit = cumulativednds(Amin, A, ρ; t = t, r = lambdares, λ = λ, Δ = deltares),
    dndsfitlq = mindnds,
    dndsfituq = maxdnds,
    deltafit = deltares,
    lambdarfit = lambdares,
    deltafitlq = deltares - CI[1],
    lambdarfitlq = lambdares - CI[2],
    deltafituq = deltares + CI[1],
    lambdarfituq = lambdares + CI[2],
    sedelta = CI[1],
    selambda = CI[2])

    residuals = DF[:dnds] .- DF[:dndsfit]
    SSres = sum(residuals.^2)
    SStot = sum( (DF[:dnds] .- mean(DF[:dnds])).^2)
    rsq = 1 - (SSres / SStot)
    DF[:rsq] = rsq
    @rput DF
    @rput deltares
    @rput lambdares
    x = R"""
    library(ggplot2)
    library(cowplot)
    g <- ggplot(DF, aes(x = A, y = dnds)) +
    geom_point(alpha = 0.5, col = "plum4") +
    geom_line(aes(y = dndsfit), alpha = 0.8, col = "firebrick") +
    geom_ribbon(aes(ymin = dndsfitlq, ymax = dndsfituq),alpha = 0.3, fill = "firebrick") +
    xlab("Clone area") +
    ylab("Interval dN/dS") +
    ggtitle(paste0("Δ = ", round(deltares, 2), ", λr = ", round(lambdares, 2)))
    """
    return dndsfit(deltares, lambdares, DF, x)
end
