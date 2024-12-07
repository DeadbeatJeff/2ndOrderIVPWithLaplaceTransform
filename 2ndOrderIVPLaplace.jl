#= #######################################################################
 Author: Dr. Jeff Rolland
 Creation Date: 12/06/2024
 
 This Julia program takes a 2nd order, linear, time-independent (LTI) 
 ordinary differential equation (ODE) dy^2/dt^2 + a1 dy/dt + a0 y = u(t) 
 with initial conditions y(0) = y0 and v(0) = ydot0, takes the Laplace 
 transform of the ODE using the initial conditions, performs a partial 
 fraction decomposition on the Laplace transform, computes the inverse 
 Laplace transform to solve the ODE, then checks the solution against the 
 input function u(t), y(0), and v(0).
####################################################################### =#

# Include modules
using Symbolics
using SymPy

# Define symbolic variables
t = symbols("t", real=true, positive=true)
Y = symbols("Y", real=true)
s = symbols("s", real=true, positive=true)

###########################################################################
# Enter parameters, initial conditions, and input function ################
###########################################################################
a1 = 4                # Coefficient for dy/dt term in ODE (ODE should be "monic": coefficient of d^2y/dt^2 term should be 1)
a0 = 3                # Coefficient for y term in ODE
y0 = 4                # Initial position
ydot0 = 5             # Initial velocity
function u(t)         # Input function
    t^2 
end

###########################################################################
# You Shouldn't Need to Change Anything Below This Line ###################
###########################################################################
println("ODE: dy^2/dt^2 + ", a1, "dy/dt + ", a0, " = ",u(t))
println("Initial Conditions: y0 = ", round(N.(y0),digits=4), ", v0 = ", round(N.(ydot0),digits=4))
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
println("Laplace transform:")
println("Y(s) = ",Y(s))
println(" ")

# Inverse Laplace transform
function yMatrix(t)
    global s
    sympy.inverse_laplace_transform(Y, s, t)
end
y = yMatrix(t)[1]
println("Solution: y(t) = ",y(t))
println(" ")

# # Check solution
Dy = sympy.diff(y(t), t)[1]
D2y = sympy.diff(y(t), t, 2)[1]
check1 = sympy.simplify(D2y(t) + a1*Dy(t) + a0*y(t))
println("Check 1: u(t) = ",check1)
check2 = subs(y(t), t => 0)
println("Check 2: y(0) = ", round(N.(check2),digits=4))
check3 = subs(Dy(t), t => 0)
println("Check 3: Dy(0) = ", round(N.(check3),digits=4))