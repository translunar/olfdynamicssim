using LinearAlgebra
using SparseArrays
using OSQP
using ControlSystems

function attitude_tracking(qref, x, params)

      J_dry = params[:J_dry]
      r_tank = params[:r_tank]
      Fmax = params[:Fmax]

      #Attitude tracking while staying as close to full thrust as possible
      Nt = 10 #number of time steps for QP
      h = params[:τ_thrust] #length of time step for QP
      Nx = 6
      Nu = 4

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      ṙ = x[8:10] #vehicle velocity (inertial frame)
      ω = x[11:13] #vehicle angular velocity (body frame)
      m_fuel = x[14] #current fuel mass

      #Total vehicle inertia assuming spherical fuel mass at COM
      J = J_dry + (2*m_fuel*r_tank*r_tank/5)*I

      #attitude error
      qe = qmult(qconj(qref), q)
      ϕ = 2*qe[2:4]

      #Controller state
      x0 = [ϕ; ω]

      # Thruster force Jacobian
      θ_t = 5.0*pi/180 # Thruster angle
      Bf = [cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 sin(θ_t);
            cos(θ_t) 0 sin(θ_t)]'

      # Thruster torque Jacobian
      #TODO: Get thruster locations from CAD model
      rul = [-1.5; -0.5; -0.5]; #Vector from CoM to upper-left thruster
      rur = [-1.5; 0.5; -0.5]; #Vector from CoM to upper-right thruster
      rlr = [-1.5; 0.5; 0.5]; #Vector from CoM to lower-right thruster
      rll = [-1.5; -0.5; 0.5]; #Vector from CoM to lower-left thruster
      Bτ = [cross(rul,Bf[:,1]) cross(rur,Bf[:,2]) cross(rlr,Bf[:,3]) cross(rll,Bf[:,4])];

      #Double-integrator dynamics
      Ac = [zeros(3,3) I; zeros(3,6)]
      Bc = [zeros(3,4); -J\Bτ] #Subtract thrust off of u0
      ABd = exp(h*[Ac Bc; zeros(Nu, Nx+Nu)])
      A = ABd[1:Nx, 1:Nx]
      B = ABd[1:Nx, Nx.+(1:Nu)]

      #Cost weights
      Q = Diagonal([(1/.1)^2*ones(3); (1/.3)^2*ones(3)])
      R = Diagonal((1/10)^2*ones(4))

      #Calculate Infinite-Horizon Cost-to-go
      S = sparse(dare(A,B,Q,R))
      droptol!(S, 1e-6)

      #Build sparse QP
      H = sparse([kron(Diagonal(I,Nt-1),[R zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [R zeros(Nu,Nx); zeros(Nx,Nu) S]])
      g = zeros(Nt*(Nx+Nu))
      C = sparse([[B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]])
      d =[-A*x0; zeros(Nx*(Nt-1))]
      U = kron(Diagonal(I,Nt), [I zeros(Nu,Nx)])
      D = [C; U]
      lb = [d; zeros(Nt*Nu)]
      ub = [d; Fmax*ones(Nt*Nu)]

      prob = OSQP.Model()
      OSQP.setup!(prob; P=H, q=g, A=D, l=lb, u=ub)
      results = OSQP.solve!(prob)

      return U*results.x
end

function attitude_tracking_condensed(qref, x, params)

      J_dry = params[:J_dry]
      r_tank = params[:r_tank]
      Fmax = params[:Fmax]

      #Attitude tracking while staying as close to full thrust as possible
      Nt = 10 #number of time steps for QP
      h = params[:τ_thrust] #length of time step for QP
      Nx = 6
      Nu = 4

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      ṙ = x[8:10] #vehicle velocity (inertial frame)
      ω = x[11:13] #vehicle angular velocity (body frame)
      m_fuel = x[14] #current fuel mass

      #Total vehicle inertia assuming spherical fuel mass at COM
      J = J_dry + (2*m_fuel*r_tank*r_tank/5)*I

      #attitude error
      qe = qmult(qconj(qref), q)
      ϕ = 2*qe[2:4]

      #Controller state
      x0 = [ϕ; ω]

      # Thruster force Jacobian
      θ_t = 5.0*pi/180 # Thruster angle
      Bf = [cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 -sin(θ_t);
            cos(θ_t) 0 sin(θ_t);
            cos(θ_t) 0 sin(θ_t)]'

      # Thruster torque Jacobian
      #TODO: Get thruster locations from CAD model
      rul = [-1.5; -0.5; -0.5]; #Vector from CoM to upper-left thruster
      rur = [-1.5; 0.5; -0.5]; #Vector from CoM to upper-right thruster
      rlr = [-1.5; 0.5; 0.5]; #Vector from CoM to lower-right thruster
      rll = [-1.5; -0.5; 0.5]; #Vector from CoM to lower-left thruster
      Bτ = [cross(rul,Bf[:,1]) cross(rur,Bf[:,2]) cross(rlr,Bf[:,3]) cross(rll,Bf[:,4])];

      #Double-integrator dynamics
      A = I + h*[zeros(3,3) I; zeros(3,6)]
      B = h*[zeros(3,4); -J\Bτ] #Subtract thrust off of u0

      #Cost weights
      Q = Diagonal([(1/.1)^2*ones(3); (1/.3)^2*ones(3)])
      R = Diagonal((1/10)^2*ones(4))

      #Build condensed QP
      Id = Diagonal(I,Nt)
      B̄ = zeros(Nx*Nt,Nu*Nt)
      Ad = Diagonal(I,Nx)
      ā = kron(ones(Nt), x0)
      for k = 1:Nt
            B̄ += kron(Id,Ad*B)
            ā[((k-1)*Nx+1):((k-1)*Nx+Nx)] = Ad*ā[((k-1)*Nx+1):((k-1)*Nx+Nx)]
            Id = [zeros(1,Nt); Id[1:end-1,:]]
            Ad = A*Ad
      end

      Q̄ = kron(Diagonal(I,Nt),Q)
      R̄ = kron(Diagonal(I,Nt),R)

      H = Symmetric(B̄'*Q̄*B̄ + R̄)
      g = B̄'*Q̄*ā

      ū = Variable(Nt*Nu)
      prob = minimize(quadform(ū,H) + 2*g'*ū, [ū >= 0; ū <= Fmax])
      #Convex.solve!(prob, SCSSolver(verbose=0, eps=1e-4, linear_solver=SCS.Direct))
      Convex.solve!(prob, ECOSSolver(verbose=0, abstol=1e-4, reltol=1e-4, feastol=1e-4))
      return ū.value[1:4]
end

function qmult(q1,q2)
      s1 = q1[1]
      v1 = q1[2:4]
      s2 = q2[1]
      v2 = q2[2:4]

      return [s1*s2 - v1'*v2; s1*v2 + s2*v1 + cross(v1,v2)]
end

function qconj(q)
      return [q[1]; -q[2:4]]
end

function force_to_u(F)
    #Converts a torque to a thruster command
    #Assuming 4-thruster setup

    # Thruster force Jacobian
    Nu = 4
    θ_t = 5.0*pi/180 # Thruster angle
    Bf = [cos(θ_t) 0 -sin(θ_t);
          cos(θ_t) 0 -sin(θ_t);
          cos(θ_t) 0 sin(θ_t);
          cos(θ_t) 0 sin(θ_t)]'

    #Set up and solve optimization problem
    u = Variable(Nu)
    e = ones(Nu)
    prob = minimize(e'*u, [u >= 0; Bf*u == F])
    Convex.solve!(prob, ECOSSolver(verbose=0))

    return round.(u.value, digits=5)
end

function torque_to_u(τ)
    #Converts a torque to a thruster command
    #Assuming 4-thruster setup

    # Thruster force Jacobian
    Nu = 4
    θ_t = 5.0*pi/180 # Thruster angle
    Bf = [cos(θ_t) 0 -sin(θ_t);
          cos(θ_t) 0 -sin(θ_t);
          cos(θ_t) 0 sin(θ_t);
          cos(θ_t) 0 sin(θ_t)]'

    # Thruster torque Jacobian
    #TODO: Get thruster locations from CAD model
    rul = [-1.5; -0.5; -0.5]; #Vector from CoM to upper-left thruster
    rur = [-1.5; 0.5; -0.5]; #Vector from CoM to upper-right thruster
    rlr = [-1.5; 0.5; 0.5]; #Vector from CoM to lower-right thruster
    rll = [-1.5; -0.5; 0.5]; #Vector from CoM to lower-left thruster
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
