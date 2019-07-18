using LinearAlgebra
using Convex, ECOS

function force_to_u(F)
    #Converts a torque to a thruster command
    #Assuming 4-thruster setup

    # Thruster force Jacobian
    Nu = 4
    θ_t = 5.0*pi/180 # Thruster cant angle
    Bf = [0 sin(θ_t) cos(θ_t);
          0 sin(θ_t) cos(θ_t);
          0 -sin(θ_t) cos(θ_t);
          0 -sin(θ_t) cos(θ_t)]'

    #Set up and solve optimization problem
    u = Variable(Nu)
    e = ones(Nu)
    prob = minimize(e'*u, [u >= 0; Bf*u == F])
    solve!(prob, ECOSSolver(verbose=0))

    return round.(u.value, digits=5)
end

function torque_to_u(τ)
    #Converts a torque to a thruster command
    #Assuming 4-thruster setup

    # Thruster force Jacobian
    Nu = 4
    θ_t = 5.0*pi/180 # Thruster cant angle
    Bf = [0 sin(θ_t) cos(θ_t);
          0 sin(θ_t) cos(θ_t);
          0 -sin(θ_t) cos(θ_t);
          0 -sin(θ_t) cos(θ_t)]'

    # Thruster torque Jacobian
    #TODO: Get thruster locations from CAD model
    rul = [-0.5; 0.5; -1.5]; #Vector from CoM to upper-left thruster
    rur = [0.5; 0.5; -1.5]; #Vector from CoM to upper-right thruster
    rlr = [0.5; -0.5; -1.5]; #Vector from CoM to lower-right thruster
    rll = [-0.5; -0.5; -1.5]; #Vector from CoM to lower-left thruster
    Bτ = [cross(rul,Bf[:,1]) cross(rur,Bf[:,2]) cross(rlr,Bf[:,3]) cross(rll,Bf[:,4])];

    #Set up and solve optimization problem
    u = Variable(Nu)
    e = ones(Nu)
    prob = minimize(e'*u, [u >= 0; Bτ*u == τ])
    solve!(prob, ECOSSolver(verbose=0))

    return round.(u.value, digits=5)
end

function wrench_to_u(w)
    #Converts a 3D wrench [F; τ] to a thruster command
    #Assuming 4 RCS thruster quads

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

function tracking(vref, v)

end
