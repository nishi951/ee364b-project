using Convex, PyPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 

# problem data
n = 7

for _ in 1:4
    Q = qr(randn(n, n))[1]
    A = diagm([-1, 4, 3, 2, 1, 4, 5]);
    b = [1, 0, 0, 0, 0, 0, 0]; R = 1; eta = 1 / (3 * norm(A))
    P(x) = proj(x, R)
    f(x) = dot(x, A*x) + 2*dot(b, x) 
    grad(x) = 2 * (A * x + b)
    x = zeros(n); xs = [x]
    fs = [f(x)]; grad_norms = []
    grads = [grad(x)]
    for k in 1:500
        x -= eta * (grad(x))
        if norm(x) > R x = P(x) end
	push!(fs, f(x))
	push!(xs, x)
        push!(grads, grad(x))
        println("$(dot(x, grad(x)))")
        println("A term: $(dot(x, A * x)), \t b term: $(dot(b, x))")
    end
    bounds = []
    for i in 2:500
        b = dot(grads[i] - grads[i - 1], xs[i-1]) - eta * dot(grads[i], grads[i-1])
        push!(bounds, b)
    end
    for b in bounds
        #println("$b")
    end
    #xs = hcat(xs...)
    println(fs[end])
end
