using Convex, PSPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
fig = figure(figsize=(10, 3))
ctrs = 1
for i in 1:3
    n = 2; Q = qr(randn(2, 2))[1];
    eig = rand(1:10, 2)
    eig[1] *= -1
    A = Q * diagm(eig) * Q'; b = [-1, 0]; R = 1; eta = 1/(2 * norm(A))
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(2); xs = [x]
    println("Figure 1, eigenvalues are $eig")
    println("b^T u_i, $(Q'*b)")
    fs = [f(x)]; grad_norms = []
    kmax = 500
    for k in 1:kmax
        x -= eta * (grad(x))
        if norm(x) > R x = P(x) end
	push!(fs, f(x))
	push!(xs, x)
    end
    xs = hcat(xs...)
    subplot(130 + i)
    X = linspace(-1.1, 1.1, 2000); Y = linspace(-1, 1, 2000)
    ctrs = contour(X, Y, [f([y, x]) for x in X, y in Y], levels=[-15:0.15:15],linewidths=0.5)
    plot(xs[1, :], xs[2, :], marker="d", color="k", linewidth=0.75)
    #plot(x[1], x[2], marker="o", linewidth=2)
    print(x)
    colorbar(ctrs)
    println("final norm: $(norm(x))")
    println("final value $(f(x))")
end
subplots_adjust(wspace=0.35)
# cax = axes([0.91, 0.1, 0.015, 0.8])
# colorbar(ctrs, cax=cax)
tight_layout()
printfig("../figs/convergence.eps", bbox_inches="tight")# dpi=2000)
