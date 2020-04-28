using LinearAlgebra
using Statistics
using Random

function coherence(X, ord=1, model=5, iter=200, ax=1)
	X = makebinary(X)

	if ord == 1
		X = allzeros(X)
		Y = ordinate(X, ax)
	end

	Abs = embedabs(Y)

	EAbs = zeros(Int, iter)
	for i = 1:iter
		null = randomnull(X, model)
		if ord == 1
			EX = allzeros(null)
			EY = ordinate(EX, ax)
		end
		EAbs[i] = embedabs(EY)
	end

	MA = mean(EAbs)
	SA = std(EAbs)

	# z-test
	return Abs, MA, SA
end

function turnover(X, ord=1, Comm=0, iter=200, ax=1)
	X = makebinary(X)

	if ord == 1
		X = allzeros(X)
		Y = ordinate(X, ax)
	end

	Y = boundfill(Y, Comm)
	Re = checkstat(Y)

	ERe = zeros(iter)
	for i = 1:iter
		null = replacenull(Y, Comm)
		ERe[i] = checkstat(null)
	end
	MR = mean(ERe)
	SR = std(ERe)

	return Re, MR, SR
end

function clumping(X, ord=1, Comm=0, ax=1)
	if Comm == 0
		X = transpose(X)
	end

	X = makebinary(X)

	if ord == 1
		X = allzeros(X)
		Y = ordinate(X)
	end

	r,c = size(Y)
	ComBnd = zeros(c,1)
	M = 0
	ComBndChi = 0

	for i = 1:r
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

	for i = 2:(c-1)
		M += ((ComBnd[i]/TotComBnd)*((ComBnd[i]-1)/(TotComBnd-1)))
		ComBndChi += ((ComBnd[i] - ExpComBnd)^2/ExpComBnd)
		df += 1
	end

	M *= c-2

	return M, ComBndChi, df
end

function replacenull(X, Comm)
	sites, species = size(X)
	null = zeros(Int, sites, species)
	if Comm == 0
		for i = 1:species
			ind = findall(X[:,i] .== 1)
			len = length(ind)
			rd = Int(ceil(rand()*(sites-len+1)))
			null[rd:(rd+len-1),i] .= 1
		end
	elseif Comm == 1
		for i = 1:sites
			ind = findall(X[i,:] .== 1)
			len = length(ind)
			rd = Int(ceil(rand()*(species-len+1)))
			null[i,rd:(rd+len-1)] .= 1
		end
	end
	return null
end

function checkstat(X)
	r,c = size(X)
	C = 0
	for i = 1:(c-1)
		Xi = X[:,i]
		R1 = sum(Xi)
		for j = (i+1):c
			Xj = X[:,j]
			S = sum(Xi.&Xj)
			R2 = sum(Xj)
			C += (R1-S) * (R2-S)
		end
	end
	return C
end

function boundfill(X, Comm)
	sites, species = size(X)
	Y = copy(X)
	if Comm == 0
		for i = 1:species
			ind = findall(X[:,i] .== 1)
			First = ind[1]
			Last = ind[length(ind)]
			Y[First:Last,i] .= 1
		end
	elseif Comm == 1
		for i = 1:sites
			ind = findall(X[i,:] .== 1)
			First = ind[1]
			Last = ind[length(ind)]
			Y[i,First:Last] .= 1
		end
	end
	return Y
end

function makebinary(X)
	Y = Int.(X .> 0)
	return Y
end

function allzeros(X, rc=3)
	Y = X
	if rc == 3 || rc == 1
		Rsum = sum(X, dims=2)
		Rind = findall(Rsum[:] .> 0)
		Y = X[Rind,:]
	end
	if rc == 3 || rc == 2
		Csum = sum(X, dims=1)
		Cind = findall(Csum[:] .> 0)
		Y = Y[:,Cind]
	end
	return Y
end

function ordinate(X, ax=1)
	V, W = ca(X)
	ir = sortperm(V[ax+1,:])
	ic = sortperm(W[ax+1,:])
	Y = X[ir,ic]
	Y = Y[size(Y,1):-1:1,:]
	return Y
end

function ca(X)
	r = sum(X, dims=2)
	c = sum(X, dims=1)
	N = sum(r)
	R = diagm(0 => r[:])
	C = diagm(0 => c[:])

	Risq = sqrt(inv(R))
	Cisq = sqrt(inv(C))

	M = Risq * X * Cisq
	Mt = transpose(M)
	P = M * Mt
	Q = Mt * M

	vals, vecs = eig(P)
	UP = transpose(vecs)
	vals, vecs = eig(Q)
	UQ = transpose(vecs)

	V = sqrt(N) * UP * Risq
	W = sqrt(N) * UQ * Cisq

	return V, W
end

function eig(X, n=size(X,1))
	vals, vecs = eigen(X)
	
	if !isreal(vals)
		mask = isreal.(vals)
		image_ = abs.(.!mask .* vals)
		real_ = mask .* vals
		vals = Float64.(real_ + image_)
	end
	if !isreal(vecs)
		mask = isreal.(vecs)
		image_ = abs.(.!mask .* vecs)
		real_ = mask .* vecs
		vecs = Float64.(real_ + image_)
	end
	
	ind = sortperm(-vals)
	vals = vals[ind]
	vecs = vecs[:,ind]

	for j = 1:size(vecs, 2)
		if sum(vecs[:,j]) < 0
			vecs[:,j] = -vecs[:,j];
		end
	end

	return vals[1:n], vecs[:,1:n]
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

# only binary mode is implemented
function randomnull(X, model=1)
	r,c = size(X)
	
	Richness = sum(X, dims=2)
	Incidence = sum(X, dims=1)
	N = sum(Richness)

    FixRichness = ones(r,1)
    EqRichness = ones(r,1) * (1/r)
    PropRichness = Richness / N
    FixIncidence = ones(1,c)
    EqIncidence = ones(1,c) * (1/c)
	PropIncidence = Incidence / N
	
	null = zeros(r,c)
	if model == 1
		P = EqRichness * EqIncidence
		null = assignsp(P, N)
	elseif model == 2
		P = EqRichness * FixIncidence
		for i = 1:c
			Pi = reshape(P[:,i], (:,1))
			null[:,i] = assignsp(Pi, Incidence[i])
		end
	elseif model == 3
		P = FixRichness * EqIncidence
		for i = 1:r
			Pi = reshape(P[i,:], (1,:))
			null[i,:] = assignsp(Pi, Richness[i])
		end;
	elseif model == 4
		P = PropRichness * FixIncidence
		for i = 1:c
			Pi = reshape(P[:,i], (:,1))
			null[:,i] = assignsp(Pi, Incidence[i])
		end;
	elseif model == 5
		P = FixRichness * PropIncidence
		for i = 1:r
			Pi = reshape(P[i,:], (1,:))
			null[i,:] = assignsp(Pi, Richness[i])
		end
	elseif model == 6
		P = PropRichness * EqIncidence
        null = assignsp(P, N)
	elseif model == 7
		P = EqRichness * PropIncidence
		null = assignsp(P, N)
	elseif model == 8
		P = PropRichness * PropIncidence
		null = assignsp(P, N)
	elseif model == 9
		null = swap(X)
	elseif model == 10
		null = randnest(X)
	end

	return null
end

# not functionally idential to Swap.m, but similar and quicker
function swap(X)
	r,c = size(X)
	ir = randperm(r)
	return X[ir,:]
end

function randnest(X)
	r,c = size(X)
	Incidence = sum(X, dims=1)
	SpecFreq = Incidence / r
	null = ones(r, c)
	mask = rand(r, c) .< repeat(SpecFreq, 5)
	return mask .* null
end

function assignsp(P, N, con=1)
	r,c = size(P)
	null = zeros(r,c)
	X = P[:]
	
	X = X / sum(X)
	C = cumsum(X)

	while true
		for i = 1:N
			while true
				ind = findall(C .>= rand())
				if null[ind[1]] != 1
					null[ind[1]] = 1
					break
				end
			end
		end
		if con == 1
			break
		end
	end
	
	return null
end