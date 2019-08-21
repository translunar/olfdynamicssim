using LinearAlgebra
using Plots

function hat(x)
      return [0  -x[3] x[2];
             x[3]  0  -x[1];
            -x[2] x[1]  0];
end

function qtoR(q)
      s = q[1]
      v = q[2:4]
      R = I + 2.0*hat(v)*(s*I + hat(v))
end

function qrot(q,r)
      s = q[1]
      v = q[2:4]
      return r + 2.0*cross(v,(s*r + cross(v,r)))
end

function qmult(q1,q2)
      s1 = q1[1]
      v1 = q1[2:4]
      s2 = q2[1]
      v2 = q2[2:4]

      return [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)]
end

function L(q)
      s = q[1]
      v = q[2:4]
      return [s -v'; v s*I + hat(v)]
end

function R(q)
      s = q[1]
      v = q[2:4]
      return [s -v'; v s*I - hat(v)]
end

function qconj(q)
      return [q[1]; -q[2:4]]
end

function qexp(ϕ)
      mag = norm(ϕ)
      return [cos(mag/2); 0.5*ϕ*sinc(mag/(2*pi))]
end

function gravity(r,t)
      μ_e = 398600.435436 #Earth GM (km^3/s^2)
      μ_m = 4902.800066 #Moon GM (km^3/s^2)

      g_e = -μ_e*r/(norm(r)^3) #Just spherical Earth for now

      #TODO: Add lunar gravity and higher-order gravity terms
      #TODO: Use JPL DE files?
end

function plot_position(soln)
      plot(soln, vars = [(0,1), (0,2), (0,3)], title="Position")
end
function plot_quaternion(soln)
      plot(soln, vars = [(0,4), (0,5), (0,6), (0,7)], title="Quaternion")
end
function plot_velocity(soln)
      plot(soln, vars = [(0,8), (0,9), (0,10)], title="Velocity")
end
function plot_omega(soln)
      plot(soln, vars = [(0,11), (0,12), (0,13)], title="Angular Rate")
end
function plot_mass(soln)
      plot(soln, vars = [(0,14)], title="Mass")
end
