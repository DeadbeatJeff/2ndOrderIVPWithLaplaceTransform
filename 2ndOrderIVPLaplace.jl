using Symbolics
using SymPy

# tReal = symbols("tReal", real=true)
t = symbols("t", real=true, positive=true)
Y = symbols("Y", real=true)
s = symbols("s", real=true, positive=true)

# Parameters and initial conditions
a1 = -2
a0 = -3
y0 = 4
ydot0 = 5
function u(t)
    t^2
end

println("u(t) = ",u(t))
println("y0 = ", round(N.(y0),digits=4))
println("Dy0 = ", round(N.(ydot0),digits=4))
println(" ")

# Laplace transform
Ly = Y
LDy = s*Y - y0
LD2y = s^2*Y - y0*s - ydot0
LHS = LD2y + a1*LDy + a0*Ly
U = sympy.integrate(u(t) * sympy.exp(-s * t), (t, 0, sympy.∞))  # Laplace transform of the right-hand side 
IVP = LHS - U(s)

# Solve for Y(s)
Y = solve(IVP, Y)[1]
Y = Y.apart(s) # Partial fraction decomposition

# Inverse Laplace transform
function yMatrix(tReal)
    global s
    sympy.inverse_laplace_transform(Y, s, tReal)
end
y = yMatrix(t)[1]
println("y(t) = ",y(t))

# # Check solution
Dy = sympy.diff(y(t), t)[1]
D2y = sympy.diff(y(t), t, 2)[1]
check1 = sympy.simplify(D2y(t) + a1*Dy(t) + a0*y(t))
println("Check 1: u(t) = ",check1)
check2 = subs(y(t), t => 0)
println("Check 2: y0 = ", round(N.(check2),digits=4))
check3 = subs(Dy(t), t => 0)
println("Check 3: Dy0 = ", round(N.(check3),digits=4))
