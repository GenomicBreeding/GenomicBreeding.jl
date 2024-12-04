using StatsBase

n::Int = 10
p::Int = 5
X = rand(n, p)
b = rand(p)
e = rand(n)
y = (X * b) + e

A = hcat(ones(n,1), X)
b_hat = inv(A'*A)*A'*y
b
cor(b_hat[2:end], b)