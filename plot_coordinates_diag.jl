using Convex, PSPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
fig = figure(figsize=(10, 3))
ctrs = 1
for i in 1:10
    n = 6
    eig = sort(rand(1:10, n) .* rand([-1, 1], n))
    A = diagm(eig); R = 1; eta = 1/(2 * norm(A)); b = randn(n)
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(n); xs = [x]
    println("Evals: $eig")
    fs = [f(x)]; grad_norms = []
    kmax = 50
    khit = -1
    for k in 1:kmax
        x -= eta * (grad(x))
        if norm(x) > R
            x = P(x)
            if khit == -1  khit = k end 
        end
	push!(fs, f(x))
	push!(xs, x)
    end
    xs = hcat(xs...)
    figure()
    for j in 1:n
        plot(0:kmax, xs[j, :], label="x_$j  (lambda = $(eig[j]))")
    end
    ymin, ymax = minimum(xs), maximum(xs)
    plot((khit, khit), (ymin, ymax), "k-")
    legend()
    printfig("no_abs/figure$(i).eps", bbox_inches="tight")
end
#subplots_adjust(wspace=0.35)
# dpi=2000)
