using Convex, PSPlot

#srand(15)
"Project onto B(R), assuming x is not in B(R)"
proj(x, R) = (R / norm(x)) * x 
fig = figure(figsize=(10, 3))
ctrs = 1
for i in 1:10
    for j in 1:10
        n = 2; Q = qr(randn(2, 2))[1];
        eig = [i, j]
        eig[1] *= -1 
        A = diagm(eig); b = [-1, 0]; R = 1; eta = 1/(2 * norm(A))
        P(x) = proj(x, R); f(x) = 0.5 * dot(x, A*x) + dot(b, x); grad(x) = (A * x + b)
        x = zeros(2); xs = [x]
        u = b;
        z = (R/norm(u)) * u; zs = [z]
        fs = [f(x)]; grad_norms = []
        kmax = 1000
        for k in 1:kmax
            x -= eta * (grad(x))
            z -= eta * (grad(z))
            if norm(x) > R
                y = deepcopy(x)
                x = P(x)
                alpha = 1/minimum(eig) * (1 - ((1 - eta * minimum(eig)) ^ (k + 1))) - eta
                term1 = (eta + (R/norm(y)) * alpha)
                term2 = b[findmin(eig)[2]]
                bound = term1^2 * term2^2 
                println("Eigenvalues: $(eig). Bound: $(bound), Term1: $(sign(term1)), Alpha: $(sign(alpha))")
                break
            end
	    push!(fs, f(x))
	    push!(xs, x)
        end
    end
end
