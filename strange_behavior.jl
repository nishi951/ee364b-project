using Convex, PSPlot

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
fig = figure(figsize=(10, 3))
ctrs = 1
include("/Users/rpathak/Documents/Classes/EE364B/project/proposal/trlib/bindings/julia/trlib.jl")
using trlib
import trlib: trlib_solve

function get_x(A, b, R)
    TR = trlib.trlib_data(A, b)
    trlib_solve(TR, R)
    return TR.sol
end

for i in 1:10
    n = 6
    eig = sort(rand(1:10, n) .* rand([-1, 1], n))
    A = diagm(eig); R = 1.0; eta = 1/(2 * norm(A)); b = randn(n)
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(n); xs = [x]
    print("Evals: $eig")
    fs = [f(x)]; grad_norms = []
    kmax = 500
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
    xtrue = get_x(A, b, R)
    print(" norm suboptimality: $(norm(x - xtrue))")
    print(" norm of xtrue: $(norm(xtrue))")
    println(" f(\hat x) - f*: $(f(x) - f(xtrue))")
    xs = hcat(xs...)
    if i != 8 continue end
    figure()
    for j in 1:n
        plot(0:kmax, xs[j, :], label="x_$j  (lambda = $(eig[j]))")
    end
    ymin, ymax = minimum(xs), maximum(xs)
    plot((khit, khit), (ymin, ymax), "k-")
    legend()
    printfig("eps/figure$(i).eps", bbox_inches="tight")
    figure(); plot(fs); printfig("f.eps");
    println("norm of final iterate: $(norm(x))")
    xtrue = get_x(A, b, R)
    println("suboptimality: $(norm(x - xtrue))")
    break
end
#subplots_adjust(wspace=0.35)
# dpi=2000)
