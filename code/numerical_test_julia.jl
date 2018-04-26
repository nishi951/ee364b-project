
# coding: utf-8

# In[41]:


using Convex
using PyPlot


# In[48]:


# Generate random indefinite QP
n = 10
X = randn(n, n)
A = Symmetric(X)
lambda= eigvals(A)
# plt[:hist](lambda)
b = randn(n)


# In[63]:


x = Variable(n)
X = Variable(n, n)
constraints = [
    trace(X) <= 3,
    X == X',
    lambdamin([X x; x' 1]) > 0
]
problem = minimize(trace(A*X) + 2*b'*x, constraints)
solve!(problem)


# In[64]:


1/2

