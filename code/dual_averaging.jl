using Convex, PyPlot

#srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 

for _ in 1:6
	# problem data 
	n = 50
	A = randn(n, n) #[-7 1; 1 -7]
	b = randn(n)
	R = 1 
	eta = 1 / norm(A)
	
	P(x) = proj(x, R)
	f(x) = dot(x, A*x) + 2*dot(b, x)
	grad(x) = 2 * (A * x + b)
	alpha(k) = eta / (k + 1)

	# dual averaging 

	x = P(randn(n))
	xs = [x]
	fs = [f(x)]
	z = zeros(n)
	for k in 1:50000
		z += grad(x)
		x = -alpha(k - 1) * z
		if norm(x) > R   x = P(x)   end
		push!(fs, f(x))
	end
	diffs = sign.(fs[1:end-1] - fs[2:end])
	diffs = sum(1 + diffs)/2
	plot(fs, label="$diffs / 50000 descent, Op norm = $(norm(A))")
end
legend()
savefig("output.png", dpi=1000)
