using LinearAlgebra
using Convex, ECOS

function wrench_to_u(w)
    #Converts a 3D wrench [F; τ] to a thruster command

    #TODO: Pull parameters out and put in a common location for dynamcis/control/estimation
    Nu = 16
    # Thruster force Jacobian
    Bf = [0 0 1.0; # +x quad
          0 0 -1.0;
          0 1.0 0;
          0 -1.0 0;
          0 0 1.0; # +y quad
          0 0 -1.0;
          1.0 0 0;
          -1.0 0 0;
          0 0 1.0; # -x quad
          0 0 -1.0;
          0 1.0 0;
          0 -1.0 0;
          0 0 1.0; # -y quad
          0 0 -1.0;
          1.0 0 0;
          -1.0 0 0]';
      # Thruster torque Jacobian
      rxp = [1.0; 0; 0]; #Vector from CoM to +x quad
      rxm = [-1.0; 0; 0]; #Vector from CoM to -x quad
      ryp = [0; 1.0; 0]; #Vector from CoM to +y quad
      rym = [0; -1.0; 0]; #Vector from CoM to -y quad
      Bτ = [hat(rxp)*Bf[:,1:4] hat(ryp)*Bf[:,5:8] hat(rxm)*Bf[:,9:12] hat(rym)*Bf[:,13:16]];

      B = [Bf; Bτ]

      #Set up and solve optimization problem
      u = Variable(Nu)
      e = ones(Nu)
      prob = minimize(e'*u, [u >= 0; B*u == w])
      solve!(prob, ECOSSolver(verbose=0))

      return round.(u.value, digits=5)
end
