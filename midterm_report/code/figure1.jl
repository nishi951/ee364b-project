using Convex, PSPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 

# problem data
n = 2
Q = qr(randn(2, 2))[1]
A = Q * diagm([-5, 5]) * Q'; b = [-1, 0]; R = 1; eta = 0.09
P(x) = proj(x, R)
f(x) = dot(x, A*x) + 2*dot(b, x) 
grad(x) = 2 * (A * x + b)
x_opt = [-1, 0]

# contour plot
X = linspace(-1.1, 1.1, 2000); Y = linspace(-1, 1, 2000)
contour(X, Y, [f([x, y]) for x in X, y in Y], levels=[-200:0.15:200],linewidths=0.5)

for _ in 1:3
    c = rand()
    x = -c * b; xs = [x]
    fs = [f(x)]; grad_norms = []
    for k in 1:500
        println("x_t^Tb = $(dot(x, b))")
        x -= eta * (grad(x))
        if norm(x) > R x = P(x) end
	push!(fs, f(x))
	push!(xs, x)
    end
    xs = hcat(xs...)
    println(fs[end])
    plot(xs[1,:], xs[2,:], marker="o")
end
xlabel("x1")
ylabel("x2")
printfig("../figs/convergence.eps")#, bbox_inches="tight", dpi=2000)
