
using Convex
using PyPlot

n = 100
X = randn(n, n)
eigs = randn(n)
A = diagm(eigs)
b = randn(n);
x = Variable(n)
X = Semidefinite(n)
constraints = [ trace(X) <= 1, ([X x; x' 1] ⪰ 0)]
problem = minimize(trace(A*X) + 2*b'*x, constraints)
solve!(problem);
f_star = problem.optval/2 # your value from CVX goes here 

x = randn(n)
x = x / norm(x)
f(x) = 0.5 * dot(x, A*x) + dot(b, x)
grad(x) = A*x + b
fs = [f(x)]
norms = [norm(x)]
for k in 1:3000
	x -= 0.1 * grad(x)
	if norm(x) > 1
		x /= norm(x)
	end
	push!(fs, f(x) - f_star)
	push!(norms, norm(x))
end

semilogy(fs[2:300])
xlabel(L"$k$")
ylabel(L"$f(x^k) - f^\star$")
savefig("output.png", bbox_inches="tight", dpi=1000)

semilogy(norms[2:300])
xlabel(L"$k$")
ylabel(L"$\|x^k\|_2$")
savefig("output.png", bbox_inches="tight", dpi=1000)
