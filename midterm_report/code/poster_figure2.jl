using Convex, PSPlot

include("/Users/rpathak/Documents/Classes/EE364B/project/proposal/trlib/bindings/julia/trlib.jl")
using trlib
import trlib: trlib_solve

function get_x(A, b, R)
    TR = trlib.trlib_data(A, b)
    trlib_solve(TR, 1.0)
    return TR.sol
end

srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
ctrs = 1
for i in 1:8
    n = 2; Q = qr(randn(2, 2))[1];
    eig = rand(1:10, 2)
    eig[1] *= -1
    A = Q * diagm(eig) * Q'; b = [-1, 0]; R = 1; eta = 1/(2 * norm(A))
    P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
    x = zeros(2); xs = [x]
    println("Figure 1, eigenvalues are $eig")
    println("b^T u_i, $(Q'*b)")
    fs = [f(x)]; grad_norms = []
    kmax = 5000
    bd_x = []
    bd = false
    for k in 1:kmax
        x -= eta * (grad(x))
        if norm(x) > R
            x = P(x)
            bd = true
        end
	push!(fs, f(x))
	if !bd  push!(xs, x) end
        if bd push!(bd_x, x) end
        # println("x(k)^T \nabla f(x^(k): $(dot(x, grad(x)))")
        # println("x(k)^T A x(k), $(dot(x, A*x))")
    end
    k = length(xs)
    x = get_x(A, b, R)
    push!(xs, bd_x[1])
    if i == 3 continue end
    if i == 4 continue end
    if i == 1
        semilogy(1:(k+1), [norm(xvar - x) for xvar in xs], color="r", label="iiii")
        semilogy((k+1):25, [norm(bx_x - x) for bx_x in bd_x[1:(25 - k)]], color="b", label="bbbb")
    else
        semilogy(1:(k+1), [norm(xvar - x) for xvar in xs], color="r")
        semilogy((k+1):25, [norm(bx_x - x) for bx_x in bd_x[1:(25 - k)]], color="b")
    end
    legend()
    println("final norm: $(norm(x))")
    println("final value $(f(x))")
end
xlabel("kk")
printfig("../figs/boundary_vs_interior.eps", bbox_inches="tight")# dpi=2000)
