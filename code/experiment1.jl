using Convex, PyPlot

#srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
x_opt = [-1, 0]
for _ in 1:1
	# problem data 
	n = 50
	eigs = randn(n)
	A = diagm(eigs) #[-7 1; 1 -7]
	b = randn(n)
	R = 1 
	eta = 1 / (4 * norm(A))
	println("OP NORM = $(norm(A))")
	P(x) = proj(x, R)
	f(x) = dot(x, A*x) + 2*dot(b, x)
	grad(x) = 2 * (eigs .* x + b)
	x = P(randn(n)); xs = [x]
	fs = [f(x)]; grad_norms = []
        for k in 1:30
		push!(grad_norms, norm(grad(x)) ^ 2)
		x -= eta * (grad(x))
		if norm(x) > R   x = P(x)   end
		push!(fs, f(x))
		push!(xs, x)
	end
	diff_norms = [norm(xs[i + 1] - xs[i]) ^ 2 for i in  1:(length(xs) -1)]
	# plot(diff_norms, label="norms of points")
	# plot(eta * grad_norms, label="upper bound")
	plot(fs[2:end]-fs[1:end - 1], label = "fn value decrease")
	plot(-norm(A)* diff_norms, label="upper bound 1")
	plot(-1/(16 * norm(A)) * grad_norms, label="upper bound 2")

end
legend()
savefig("output.png", dpi=1000)
