using LinearAlgebra
using Statistics
using StatsBase
using Distributions
using ProgressMeter

"""
	metacommunity(X, nsim=200)

	Compute coherence, turnover, clumping of a metacommunity. Equivilant to `metacommunity(X,ord=1,model=5,Comm=0,iter=nsim,ax=1)` in EMS by CHRIS L. HIGGINS (Presley et al. 2010)

	`X` is assumed to be site-by-species incidence matrix. Site as row, species as col.
	`nsim` specifies numbers of simulation (null model generation).
	return values are identically named and ordered to EMS.
"""
function metacommunity(X, nsim=200)
	Abs, Apr, MA, SA = coherence(X, nsim)
	Re, Rpr, MR, SR = turnover(X, nsim)
	M, Mpr = clumping(X)
	ret = [Abs, Apr, MA, SA, Re, Rpr, MR, SR, M, Mpr]
	return ret
end

function summarize(ret)
	Abs, Apr, MA, SA, Re, Rpr, MR, SR, M, Mpr = ret
	coh = MA - Abs
	clu = M - 1
	tur = Re - MR
	if coh > 0; print("POSITIVE coherence");else;print("NEGATIVE coherence");end
	print(", p = ");println(Apr);
	print("\tEmbedded absences: "); println(Abs)
	print("\tNull model expected: "); println(MA)
	if tur > 0; print("POSITIVE turnover");else;print("NEGATIVE turnover");end
	print(", p = ");println(Rpr);
	print("\tSpecies replacement: "); println(Re)
	print("\tNull model expected: "); println(MR)
	if clu > 0; print("POSITIVE clumping");else;print("NEGATIVE clumping");end
	print(", p = ");println(Mpr);
	print("\tMorista's I: "); print(M)
end

function coherence(X, nsim=200)
	X = makebinary(X)
	X = allzeros(X)
	Y = ordinate(X)
	Abs = embedabs(Y)
	null = randomnull_model5(X, nsim)
	EAbs = zeros(Int, nsim)
	EX = zeros(Int64, size(X))
	@showprogress 1 "examing coherence of null..." for i = 1:nsim
		EX = allzeros(@view null[:,:,i])
		EY = ordinate(EX)
		EAbs[i] = embedabs(EY)
	end
	MA = mean(EAbs)
	SA = std(EAbs)
	Apr = ztest(Abs, MA, SA, nsim)
	return Abs, Apr, MA, SA
end

function turnover(X, nsim=200)
	X = makebinary(X)
	X = allzeros(X)
	Y = ordinate(X)
	Y = boundfill_comm0(Y)
	Re = checkstat(Y)
	null = replacenull_comm0(Y, nsim)
	ERe = zeros(Int, nsim)
	@showprogress 1 "examing turnover of null... " for i = 1:nsim
		ERe[i] = checkstat(@view null[:,:,i])
	end
	MR = mean(ERe)
	SR = std(ERe)
	Rpr = ztest(Re, MR, SR, nsim)
	return Re, Rpr, MR, SR
end

function clumping(X)
	X = transpose(X)
	X = makebinary(X)
	X = allzeros(X)
	Y = ordinate(X)
	r, c = size(Y)
	ComBnd = zeros(c, 1)
	M = 0
	ComBndChi = 0
	@inbounds @simd for i = 1:r
		ind = findall(Y[i,:] .== 1)
		First = ind[1];
		Last = ind[length(ind)]
		for j = 1:c
			if First == j
				ComBnd[j] += 1
			end
			if Last == j
				ComBnd[j] += 1
			end
		end
	end
	TotComBnd = r * 2 - ComBnd[1] - ComBnd[c]
	ExpComBnd = TotComBnd / (c - 2)
	df = -1
	@inbounds @simd for i = 2:(c-1)
		M += ((ComBnd[i]/TotComBnd)*((ComBnd[i]-1)/(TotComBnd-1)))
		ComBndChi += ((ComBnd[i] - ExpComBnd)^2/ExpComBnd)
		df += 1
	end
	M *= c-2
	if M < 1
		Mpr = chi2test(ComBndChi, df, "left")
	else
		Mpr = chi2test(ComBndChi, df, "right")
	end
	return M, Mpr
end

# two-tailed z-test
function ztest(mu, xbar, stddev, n)
    z = (xbar-mu)/(stddev/sqrt(n))
    p = 2 * cdf(Normal(), -abs(z))
    return p
end

function chi2test(x, v, tail="left")
	if tail == "left"
		p = cdf(Chisq(v), x)
	elseif tail == "right"
		p = 1 - cdf(Chisq(v), x)
	end
	return p
end

function boundfill_comm0(X)
	sites, species = size(X)
	Y = copy(X)
	for i = 1:species
		ind = findall(X[:,i] .== 1)
		First = ind[1]
		Last = ind[length(ind)]
		Y[First:Last,i] .= 1
	end
	return Y
end

function checkstat(X)
	r, c = size(X)
	C = 0
	Xij = zeros(Int64, r)
	@inbounds @simd for i = 1:(c-1)
		Xi = @view X[:,i]
		R1 = sum(Xi)
		for j = (i+1):c
			Xj = @view X[:,j]
			S = sum(Xi .& Xj)
			R2 = sum(Xj)
			C += (R1-S)*(R2-S)
		end
	end
	return C
end

function embedabs(X)
	sites, species = size(X)
	FirstSpecies, LastSpecies, FirstSite, LastSite = boundfind(X)
	Abs = 0
	for i = 1:sites
		for j = 1:species
			if X[i,j] == 0
				flag = 0
				if j > FirstSpecies[i] && j < LastSpecies[i]
					Abs += 1
					flag = 1
				end
				if flag == 0 && i > FirstSite[j] && i < LastSite[j]
					Abs += 1
				end
			end
		end
	end
	return Abs
end

function boundfind(X)
	r,c = size(X)
	FirstSpecies = zeros(r,1)
    LastSpecies = zeros(r,1)
    FirstSite = zeros(1,c)
	LastSite = zeros(1,c)
	for i = 1:r
        ind = findall(X[i,:] .== 1)
        FirstSpecies[i] = ind[1];
        LastSpecies[i] = ind[length(ind)];
	end
	for i = 1:c
        ind = findall(X[:,i] .== 1);
        FirstSite[i] = ind[1];
        LastSite[i] = ind[length(ind)];
	end
	return FirstSpecies, LastSpecies, FirstSite, LastSite
end

# EMS randomnull model=5
function randomnull_model5(X, nsim)
    nrow, ncol = size(X)
    Incidence = sum(X, dims=1)
    Richness = sum(X, dims=2)
    FixRichness = ones(nrow, 1)
    PropIncidence = Incidence / sum(Incidence)
    P = FixRichness * PropIncidence
    null = zeros(UInt8, nrow, ncol, nsim)
    stop = 0
	wv = FrequencyWeights(@view P[1,:])
	pbar = Progress(nsim, 1, "generating random null...   ")
    @inbounds @simd for k = 1:nsim
        for i = 1:nrow
            wv = FrequencyWeights(@view P[i,:])
            for j = 1:Richness[i]
                stop = 0
                while stop == 0
                    ind = sample(1:ncol, wv)
                    if null[i, ind, k] != 1
                        null[i, ind, k] = 1
                        stop = 1
                    end
				end
            end
		end
		next!(pbar)
    end
    return null
end

function replacenull_comm0(X, nsim)
	sites, species = size(X)
	null = zeros(Int, sites, species, nsim)
	pbar = Progress(nsim, 1, "generating replace null...  ")
	@inbounds @simd for k = 1:nsim
		for i = 1:species
			len = sum(@view X[:,i])
			rd = Int(ceil(rand()*(sites-len+1)))
			null[rd:(rd+len-1),i,k] .= 1
		end
		next!(pbar)
	end
	return null
end

function makebinary(X)
	return Int.(X .> 0)
end

function allzeros(X)
	Rsum = sum(X, dims=2)
	Rind = findall(Rsum[:] .> 0)
	Y = @view X[Rind,:]
	Csum = sum(X, dims=1)
	Cind = findall(Csum[:] .> 0)
	Y = @view Y[:,Cind]
	return Y
end

function ordinate(X, ax=1)
	V, W = ca(X)
	ir = sortperm(@view V[ax+1,:])
	ic = sortperm(@view W[ax+1,:])
	Y = @view X[ir,ic]
	Y = @view Y[size(Y,1):-1:1,:]
	return Y
end

function ca(X)
	tot = sum(X)
	Y = X / tot
	rw = sum(Y, dims=2)
	cw = sum(Y, dims=1)
	rc = rw * cw
	Y = (Y-rc) ./ sqrt.(rc)
	F = svd(Y, full=true)
	V = F.U ./ sqrt.(rw)
	W = F.V ./ sqrt.(cw)
	return V, W
end

function eig(X)
	vals, vecs = eigen(X)
	if !isreal(vals)
		mask = isreal.(vals)
		vals = Float64.(mask .* vals + abs.(.!mask .* vals))
	end
	if !isreal(vecs)
		mask = isreal.(vecs)
		vecs = Float64.(mask .* vecs + abs.(.!mask .* vecs))
	end
	ind = sortperm(-vals)
	vecs = @view vecs[:,ind]
	for j = 1:size(vecs, 2)
		if sum(@view vecs[:,j]) < 0
			vecs[:,j] = - @view vecs[:,j]
		end
	end
	return vecs
end

println()
println("########### fastEMS.jl ###########")
println("#       AUTHOR: YUANFEI PAN      #")
println("#         DATE: 2020-02-27       #")
println("########### fastEMS.jl ###########")
println()
println("INFO: Precompiling fastEMS library.")
println("INFO: This may take a few seconds.")
metacommunity(rand(0:1, 10, 10))
println("INFO: fastEMS library loaded.")
println("INFO: See ?metacommunity for help.")
