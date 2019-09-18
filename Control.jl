using LinearAlgebra
using SparseArrays
using OSQP
using ControlSystems

function attitude_tracking_setup(x, params)
      #TODO: Try L1 cost on u and x w/quadratic cost-to-go to minimize thruster-off pulsing
      #TODO: Calculate "effective ISP" based on reasonable errors
      J_dry = params[:J_dry]
      r_tank = params[:r_tank]
      Fmax = params[:Fmax]
      Fmin = params[:Fmin]
      if Fmax != Fmin
            ΔF = Fmax-Fmin
      else
            ΔF = Fmax
      end

      #Attitude tracking while staying as close to full thrust as possible
      Nt = 25 #number of time steps for QP
      τ_thrust = params[:τ_thrust]
      q_thrust = params[:q_thrust]
      h = τ_thrust*(2^q_thrust)
      Nx = 6
      Nu = 4

      #Total vehicle inertia assuming spherical fuel mass at COM
      m_fuel = x[14] #current fuel mass
      J = J_dry + (2*m_fuel*r_tank*r_tank/5)*I

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
      Q = Diagonal([(1/(5*pi/180))^2*ones(3); (1/(10*pi/180))^2*ones(3)])
      Rlqr = Diagonal((1/Fmax)^2*ones(4))
      Rqp = Diagonal((0.1/Fmax)^2*ones(4))
      α = 0.01/Fmax

      #Calculate Infinite-Horizon Cost-to-go
      S = sparse(dare(A,B,Q,Rlqr))
      droptol!(S, 1e-6)

      #Build sparse QP
      U = kron(Diagonal(I,Nt), [I zeros(Nu,Nx)]) #Matrix that picks out all u
      H = sparse([kron(Diagonal(I,Nt-1),[Rqp zeros(Nu,Nx); zeros(Nx,Nu) Q]) zeros((Nx+Nu)*(Nt-1), Nx+Nu); zeros(Nx+Nu,(Nx+Nu)*(Nt-1)) [Rqp zeros(Nu,Nx); zeros(Nx,Nu) S]])
      g = α*(U'*ones(Nt*Nu)) #zeros(Nt*(Nx+Nu))
      C = sparse([[B -I zeros(Nx,(Nt-1)*(Nu+Nx))]; zeros(Nx*(Nt-1),Nu) [kron(Diagonal(I,Nt-1), [A B]) zeros((Nt-1)*Nx,Nx)] + [zeros((Nt-1)*Nx,Nx) kron(Diagonal(I,Nt-1),[zeros(Nx,Nu) Diagonal(-I,Nx)])]])
      d =zeros(Nx*Nt)
      D = [C; U]
      lb = [d; zeros(Nt*Nu)]
      ub = [d; Fmax*ones(Nt*Nu)]

      prob = OSQP.Model()
      OSQP.setup!(prob; P=H, q=g, A=D, l=lb, u=ub, verbose=false)

      return A, prob
end

function attitude_tracking(qref, x, A, qpprob, params)

      #Attitude tracking while staying as close to full thrust as possible
      Nt = 25 #number of time steps for QP
      τ_thrust = params[:τ_thrust]
      q_thrust = params[:q_thrust]
      h = τ_thrust*(2^q_thrust)
      Nx = 6
      Nu = 4
      Fmax = params[:Fmax]
      Fmin = params[:Fmin]
      if Fmax != Fmin
            ΔF = Fmax-Fmin
      else
            ΔF = Fmax
      end

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      ṙ = x[8:10] #vehicle velocity (inertial frame)
      ω = x[11:13] #vehicle angular velocity (body frame)
      m_fuel = x[14] #current fuel mass

      #attitude error
      qe = qmult(qconj(qref), q)
      ϕ = 2*qe[2:4]

      #Controller state
      x0 = [ϕ; ω]

      #Update QP parameters
      d =[-A*x0; zeros(Nx*(Nt-1))]
      lb = [d; zeros(Nt*Nu)]
      ub = [d; Fmax*ones(Nt*Nu)]
      OSQP.update!(qpprob, l=lb, u=ub)

      results = OSQP.solve!(qpprob)
      δu = round.(results.x[1:4], digits=2)

      return δu
end

function effective_Isp(thist, uhist, params)
      #Calculates effective ISP for a pulsed thruster trajectory
      Fmax = params[:Fmax]
      Isp = params[:Isp]
      τ_thrust = params[:τ_thrust]*1000 #put time in ms
      q_thrust = params[:q_thrust]
      h = τ_thrust*(2^q_thrust)
      e = zeros(size(u,1))
      for i = 1:size(uhist,1)
            t_last = 0.0
            t_on = 0.0
            for k = 1:length(thist)
                  on = uhist[i,k]/Fmax
                  if on == 1.0
                        t_on += h
                  elseif on == 0.0
                        #do nothing
                  else
                        t_total = t_on + h
                        t_on += on*h
                        e[i] = (t_total*thruster_efficiency(t_on) + t_last*e[i])/(t_last+t_total)
                        t_on = 0.0
                        t_last += t_total
                  end
            end
      end
      e_avg = sum(e)/length(e)
      return Isp*e_avg, e_avg
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
