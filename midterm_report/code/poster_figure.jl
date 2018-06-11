using Convex, PSPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
ctrs = 1
for i in 1:1
    n = 2; Q = qr(randn(2, 2))[1];
    eig = rand(1:10, 2)
    eig[1] *= -1
    A = Q * diagm(eig) * Q'; b = [-1, 0]; R = 1; eta = 1/(2 * norm(A))
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(2); xs = [x]
    z = zeros(2); zs = [z]
    println("Figure 1, eigenvalues are $eig")
    println("b^T u_i, $(Q'*b)")
    fs = [f(x)]; grad_norms = []
    kmax = 5000
    for k in 1:kmax
        x -= eta * (grad(x))
        z -= eta * (grad(z) - minimum(eig) * z)
        if norm(x) > R x = P(x) end
        if norm(z) > R z = P(z) end
	push!(fs, f(x))
	push!(xs, x)
        push!(zs, z)
        # println("x(k)^T \nabla f(x^(k): $(dot(x, grad(x)))")
        # println("x(k)^T A x(k), $(dot(x, A*x))")
    end
    xs = hcat(xs...)
    zs = hcat(zs...)
    X = linspace(-1.1, 1.1, 2000); Y = linspace(-1, 1, 2000)
    ctrs = contour(X, Y, [f([y, x]) for x in X, y in Y], levels=[-20:0.15:20],linewidths=0.5)
    rands = [P(randn(2)) for _ in 1:2000]
    rands = hcat(rands...)
    plot(x[1], x[2], marker="o", markersize=10, color="r", linewidth=2, label="xopt")
    scatter(rands[1,:], rands[2,:], linewidth=1, s=3, label="BBBB")
    plot(xs[1, 1:500], xs[2, 1:500], marker="d", color="k", markersize=6, linewidth=1.25,label="xxxx")
    colorbar(ctrs)
    println("final norm: $(norm(x))")
    println("final value $(f(x))")
end
printfig("../figs/convergence_poster.eps", bbox_inches="tight")# dpi=2000)
