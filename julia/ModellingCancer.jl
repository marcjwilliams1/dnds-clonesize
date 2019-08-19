using CancerSeqSim
using DataFrames
using RCall
using ProgressMeter
using Distributions
using Random
using ArgParse
R"""
library(viridis)
library(cowplot)
library(tidyverse)
library(popgendnds) #this is a package written for the project that does a ML fit
""";

s = ArgParseSettings()
@add_arg_table s begin
    "--syntheticcohort"
        help = "Synthetic cohort"
    "--syntheticcohort_its"
        help = "# Iterations synthetic cohort"
        arg_type = Int
        default = 500
    "--syntheticcohort_diffmu"
        help = "Synthetic cohort different mu"
    "--syntheticcohort_diffmu_its"
        help = "# Iterations synthetic cohort"
        arg_type = Int
        default = 2000
    "--syntheticcohort_power"
        help = "Synthetic cohort power"
    "--syntheticcohort_power_its"
        help = "# Iterations synthetic cohort"
        arg_type = Int
        default = 2000
    "--syntheticcohort_fmin"
        help = "Synthetic cohort fmin"
    "--syntheticcohort_fmin_its"
        help = "# Iterations synthetic cohort"
        arg_type = Int
        default = 2000
    "--syntheticcohort_inferreds"
        help = "Inferred s"
end

parsed_args = parse_args(ARGS, s)

println("Generate random cohort with different selection coefficients")
Random.seed!(123)
mup = 0.1
Nmax = 10^5
mud = mup/(Nmax)
mud = 0.005
println(Nmax * mud)
its = parsed_args["syntheticcohort_its"]
svec = [0.0, 0.25, 0.5, -0.1]
Nd = zeros(Float64, its * length(svec))
Np = zeros(Float64, its * length(svec))
clonesized = zeros(Float64, its * length(svec))
clonesizep = zeros(Float64, its * length(svec))
sel = String[]
j = 1

driverVAF = Float64[]
passengerVAF = Float64[]
simnumdriver = Int64[]
simnumpass = Int64[]
seld = Float64[]
selp = Float64[]

for s in svec
    println(s)
    sfunc() = s
    @showprogress for i in 1:its
        x = simulatedifferentmutations(Nmax = Nmax, μp = mup, μd = mud, clonalmutations = 100, s = sfunc,
        fitnessfunc = CancerSeqSim.nonmultiplicativefitness);
        append!(driverVAF,
            map((depth, freq) -> rand(Binomial(depth, freq))/depth, rand(Poisson(100), length(x.output.trueVAFd)),
            x.output.trueVAFd./Nmax))
        append!(passengerVAF,
            map((depth, freq) -> rand(Binomial(depth, freq))/depth, rand(Poisson(100), length(x.output.trueVAFp)),
            x.output.trueVAFp./Nmax))
        #append!(driverVAF, x.output.trueVAFd./Nmax)
        #append!(passengerVAF, x.output.trueVAFp./Nmax)
        append!(simnumdriver, fill(i, length(x.output.trueVAFd)))
        append!(simnumpass, fill(i, length(x.output.trueVAFp)))
        append!(seld, fill(s, length(x.output.trueVAFd)))
        append!(selp, fill(s, length(x.output.trueVAFp)))
    end
end

DFd = DataFrame(VAF = driverVAF, s = seld, simnum = simnumdriver)
DFp = DataFrame(VAF = passengerVAF, s = selp, simnum = simnumpass)

dfd = DFd[DFd[:VAF] .> 0.01, :]
dfp = DFp[DFp[:VAF] .> 0.01, :]

@rput dfp
@rput dfd

R"""

dnds <- data.frame()
for (selection in unique(dfd$s)){
    drivers <- dfd %>%
    filter(s == selection)
    passengers <- dfp %>%
    filter(s == selection)

out <- get_cumulative(drivers$VAF, passengers$VAF, ccfmin = 0.01, mud = 0.005, mup = 0.1, step = 0.025) %>%
    mutate(w = nmuts/max(nmuts)) %>%
    filter(VAF < 0.5) %>%
    select(VAF, dnds, w) %>%
    fitintervaldNdSlsq(., trues = selection, ccfmin = 0.01)
    print(out)
    dnds <- rbind(dnds, out$data)
}

write_csv(dnds, $(parsed_args["syntheticcohort"]))

"""

println("Generate cohort with different mutation rates")

using Random
Random.seed!(123)
mup = 0.5
Nmax =  10^4
mudvec = [0.2, 0.1, 0.01, 0.005, 0.001]
println(Nmax.*mudvec)
its = parsed_args["syntheticcohort_diffmu_its"]
Nd = zeros(Float64, its * length(mudvec))
Np = zeros(Float64, its * length(mudvec))
clonesized = zeros(Float64, its * length(mudvec))
clonesizep = zeros(Float64, its * length(mudvec))
sel = String[]
j = 1

driverVAF = Float64[]
passengerVAF = Float64[]
simnumdriver = Int64[]
simnumpass = Int64[]
mudout = Float64[]
mupout = Float64[]
soutd = Float64[]
soutp = Float64[]
svec = [0.0, 0.2, 0.5]

for sel in svec
    println(sel)
    sfunc() = sel
    println(sfunc())
    for mud in mudvec
        println(mud)
        @showprogress for i in 1:its
            x = simulatedifferentmutations(Nmax = Nmax, μp = mup, μd = mud, clonalmutations = 0, s = sfunc,
            fitnessfunc = CancerSeqSim.nonmultiplicativefitness);
            append!(driverVAF, x.output.trueVAFd./Nmax)
            append!(passengerVAF, x.output.trueVAFp./Nmax)
            append!(simnumdriver, fill(i, length(x.output.trueVAFd)))
            append!(simnumpass, fill(i, length(x.output.trueVAFp)))
            append!(mudout, fill(mud, length(x.output.trueVAFd)))
            append!(mupout, fill(mud, length(x.output.trueVAFp)))
            append!(soutp, fill(sel, length(x.output.trueVAFp)))
            append!(soutd, fill(sel, length(x.output.trueVAFd)))
        end
    end
end

DFd = DataFrame(VAF = driverVAF, mud = mudout, trues = soutd, simnum = simnumdriver)
DFp = DataFrame(VAF = passengerVAF, mud = mupout, trues = soutp, simnum = simnumpass)

dfd = DFd[DFd[:VAF] .> 0.01, :]
dfp = DFp[DFp[:VAF] .> 0.01, :]

@rput dfp
@rput dfd

# we'll move over to R to fit the data and make the plot using popgendnds
R"""
ccfmin <- 0.01
dnds <- data.frame()
for (sel in rev(unique(dfp$trues))){
print(paste0("selection = ", sel))
    for (mu in unique(dfp$mud)){
        print(mu)
        drivers <- dfd %>%
            filter(VAF > ccfmin) %>%
            filter(mud == mu, trues == sel)
        passengers <- dfp %>%
            filter(VAF > ccfmin) %>%
            filter(mud == mu, trues == sel)

        cum1 <- get_cumulative(drivers$VAF, passengers$VAF,
                ccfmin = ccfmin, mud = mu, mup = 0.5, step = 0.01) %>%
                mutate(w = nmuts/max(nmuts)) %>%
                select(VAF, dnds, w)
        out <- cum1 %>%
            fitintervaldNdSlsq(., trues = sel, ccfmin = ccfmin)
        print(out)
        dnds <- rbind(dnds, mutate(out$data, mud = mu))
    }
}

write_csv(dnds, $(parsed_args["syntheticcohort_diffmu"]))
"""


println("Power calculations")

Random.seed!(123)
mup = 0.005
Nmax =  10^5
mudvec = [2.7 * 0.005]
println(Nmax.*mudvec)
its = parsed_args["syntheticcohort_power_its"]
Nd = zeros(Float64, its * length(mudvec))
Np = zeros(Float64, its * length(mudvec))
clonesized = zeros(Float64, its * length(mudvec))
clonesizep = zeros(Float64, its * length(mudvec))
sel = String[]
j = 1

driverVAF = Float64[]
passengerVAF = Float64[]
simnumdriver = Int64[]
simnumpass = Int64[]
mudout = Float64[]
mupout = Float64[]
soutd = Float64[]
soutp = Float64[]
svec = [0.25]

for sel in svec
    println(sel)
    sfunc() = sel
    for mud in mudvec
        println(mud)
        @showprogress for i in 1:its
            x = simulatedifferentmutations(Nmax = Nmax, μp = mup, μd = mud, clonalmutations = 0, s = sfunc,
            fitnessfunc = CancerSeqSim.nonmultiplicativefitness);
            append!(driverVAF, x.output.trueVAFd./Nmax)
            append!(passengerVAF, x.output.trueVAFp./Nmax)
            append!(simnumdriver, fill(i, length(x.output.trueVAFd)))
            append!(simnumpass, fill(i, length(x.output.trueVAFp)))
            append!(mudout, fill(mud, length(x.output.trueVAFd)))
            append!(mupout, fill(mud, length(x.output.trueVAFp)))
            append!(soutp, fill(sel, length(x.output.trueVAFp)))
            append!(soutd, fill(sel, length(x.output.trueVAFd)))
        end
    end
end

DFd = DataFrame(VAF = driverVAF, mud = mudout, trues = soutd, simnum = simnumdriver)
DFp = DataFrame(VAF = passengerVAF, mud = mupout, trues = soutp, simnum = simnumpass)

dfd = DFd[DFd[:VAF] .> 0.01, :]
dfp = DFp[DFp[:VAF] .> 0.01, :]

@rput dfp
@rput dfd

# we'll move over to R to fit the data and make the plot using popgendnds
R"""

dnds <- data.frame()
j <- 1

samplevec <- c(10, 20, 30, 50, 75, 100, 200, 300, 500, 750, 1000) / length(dfp$VAF)

plots <- list()
i <- 1

for (samplefrac in samplevec){
    for (j in 1:30){
        print(samplefrac)
            drivers <- dfd %>%
        sample_frac(samplefrac)
            passengers <- dfp %>%
        sample_frac(samplefrac)
        print(length(drivers$VAF))
        print(length(passengers$VAF))
        ndrivers = length(drivers$VAF)
        npassengers = length(passengers$VAF)
        out <- get_cumulative(drivers$VAF, passengers$VAF,
        ccfmin = 0.01, mud = 2.7 * 0.005, mup = 0.005, step = 0.01) %>%
                mutate(w = nmuts/max(nmuts)) %>%
                select(VAF, dnds, w) %>%
                fitintervaldNdSlsq(., trues = 0.25, ccfmin = 0.01)
                print(out)
        dnds <- rbind(dnds, mutate(out$data, sample = samplefrac, simnum = j,
                    ndrivers = ndrivers, npassengers = npassengers))
        j <- j + 1
    }
plots[[i]] <- out$plot
i <- i + 1
}

write_csv(dnds, $(parsed_args["syntheticcohort_power"]))
"""

println("dN/dS as a function of fmin")
Random.seed!(1)
mup = 0.001
Nmax =  10^3 #set small Nmax otherwise there are too many mutations
mud = 2.7 * 0.001
println(Nmax.*mud)
j = 1

sel = 0.25
sfunc() = sel
nsets = 50
its = parsed_args["syntheticcohort_fmin_its"]

driverVAF = Float64[]
passengerVAF = Float64[]
simnumdriver = Int64[]
simnumpass = Int64[]
mudout = Float64[]
mupout = Float64[]
soutd = Float64[]
soutp = Float64[]
setnump = Int64[]
setnumd = Int64[]
Nd = zeros(Float64, its * length(mudvec))
Np = zeros(Float64, its * length(mudvec))
clonesized = zeros(Float64, its * length(mudvec))
clonesizep = zeros(Float64, its * length(mudvec))
#sel = String[]

for n in 1:nsets
    @showprogress for i in 1:its
        x = simulatedifferentmutations(Nmax = Nmax, μp = mup, μd = mud, clonalmutations = 0, s = sfunc,
        fitnessfunc = CancerSeqSim.nonmultiplicativefitness);
        append!(driverVAF, x.output.trueVAFd./Nmax)
        append!(passengerVAF, x.output.trueVAFp./Nmax)
        append!(simnumdriver, fill(i, length(x.output.trueVAFd)))
        append!(simnumpass, fill(i, length(x.output.trueVAFp)))
        append!(mudout, fill(mud, length(x.output.trueVAFd)))
        append!(mupout, fill(mud, length(x.output.trueVAFp)))
        append!(soutp, fill(sel, length(x.output.trueVAFp)))
        append!(soutd, fill(sel, length(x.output.trueVAFd)))
        append!(setnump, fill(n, length(x.output.trueVAFp)))
        append!(setnumd, fill(n, length(x.output.trueVAFd)))
    end
end

DFd = DataFrame(VAF = driverVAF, mud = mudout, trues = soutd, simnum = simnumdriver, simulationset = setnumd)
DFp = DataFrame(VAF = passengerVAF, mud = mupout, trues = soutp, simnum = simnumpass, simulationset = setnump)

dfd = DFd
dfp = DFp

@rput dfp
@rput dfd

# we'll move over to R to fit the data and make the plot using popgendnds
# dN/dS calcuations are weighted
R"""

dnds <- data.frame()
sout <- data.frame()
fminvec <- c(0.0, 0.01, 0.05)
for (set in unique(dfd$simulationset)){
    for (fmin in fminvec){
        drivers <- dfd %>%
            filter(simulationset == set) %>%
            filter(VAF > fmin)
        passengers <- dfp %>%
            filter(simulationset == set) %>%
            filter(VAF > fmin)

        out <- get_cumulative(drivers$VAF, passengers$VAF,
            ccfmin = fmin + 0.000001, mud = 2.7 * 0.05, mup = 0.05, step = 0.025) %>%
            mutate(w = nmuts/max(nmuts)) %>%
            #filter(VAF < 0.5) %>%
            select(VAF, dnds, w) %>%
            fitintervaldNdSlsq(., trues = 0.25, ccfmin = fmin + 0.000001)
        dnds <- rbind(dnds, mutate(out$data, fmin = fmin, simulationset = set))
        sout <- rbind(sout, data.frame(s = out$parameters$s[1], fmin = fmin))
    }
}

write_csv(dnds, $(parsed_args["syntheticcohort_fmin"]))
write_csv(sout, $(parsed_args["syntheticcohort_inferreds"]))
"""
