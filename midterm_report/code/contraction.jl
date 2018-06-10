using Convex, PSPlot, ProgressMeter

#srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
fig = figure(figsize=(10, 3))
ctrs = 1
for i in 1:1
    n = 2; Q = qr(randn(2, 2))[1];
    eig = rand(1:10, 2)
    eig[1] *= -1
    A = Q * diagm(eig) * Q'; b = [-1, 0]; R = 1; eta = 1/(2 * norm(A))
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(2); xs = [x]; ys = [x];
    println("Figure 1, eigenvalues are $eig")
    println("b^T u_i, $(Q'*b)")
    fs = [f(x)]; grad_norms = []
    kmax = 3000
    @showprogress "Doing PGD..." for k in 1:kmax
        x -= eta * (grad(x))
        push!(ys, x)
        if norm(x) > R x = P(x) end
	push!(fs, f(x))
	push!(xs, x)
    end
    distsx = [norm(x - xs[end]) for x in xs[1:200]]
    distsy = [norm(y - xs[end]) for y in ys[1:200]]
    plot(distsx[2:end] ./ distsx[1:end-1], color="b")
    plot(distsy[2:end] ./ distsy[1:end-1], color="r")
    plot([norm(x) for x in xs[2:end]], color="k")
    println("max: $(maximum(distsx[2:25] ./ distsx[1:24]))")
    #dists = [norm(x - xs[end]) for x in xs]
    #factor = exp.(log.(dists[2:end]) - log.(dists[1:end-1]))
    #plot(factor, label="Op norm = $(norm(A))")
    #println("final norm: $(norm(x))")
    #println("final value $(f(x))")
end
legend()
tight_layout()
printfig("../figs/contraction.eps", bbox_inches="tight")# dpi=2000)
