using LinearAlgebra
using DifferentialEquations
using Plots

vehicle_params = (m=100.0,
             J=Diagonal([100.0, 100.0, 30.0]),
             m_slosh = 10.0,
             l_slosh = 0.5)

function test()
      r0 = [0; 0; 0]
      q0 = [1.0; 0; 0; 0]
      s0 = [0; 0; 0]
      v0 = [0; 0; 0]
      ω0 = [0; 0; 0]
      ṡ0 = [0; 0; 0]
      x0 = [r0; q0; s0; v0; ω0; ṡ0]
      prob = ODEProblem(open_loop!, x0, (0.0, 60.0), vehicle_params)
      soln = solve(prob)
      plot(soln)
end

function open_loop!(ẋ::AbstractVector,x::AbstractVector,params,t)
      #Just a function to test open-loop inputs

      #No inputs
      u0 = zeros(19)

      #Test inputs for main engine
      u_m0 = [0;0; #main engine gimbal angles (x,y)
              50.0; #main engine thrust
              0;0;0;0; #+x-quad, +y force`
              0;0;0;0; #+y-quad
              0;0;0;0; #-x-quad, -y force
              0;0;0;0] #-y-quad
      u_mx = [pi/6;0; #main engine gimbal angles (x,y)
              50.0; #main engine thrust
              0;0;0;0; #+x-quad, +y force`
              0;0;0;0; #+y-quad
              0;0;0;0; #-x-quad, -y force
              0;0;0;0] #-y-quad
      u_my = [0;pi/6; #main engine gimbal angles (x,y)
              50.0; #main engine thrust
              0;0;0;0; #+x-quad, +y force`
              0;0;0;0; #+y-quad
              0;0;0;0; #-x-quad, -y force
              0;0;0;0] #-y-quad

      #Test inputs for thrusters
      u_fx = [0;0;0; #main engine
              0;0;0;0; #+x-quad, +y force`
              0;0;1.0;0; #+y-quad
              0;0;0;0; #-x-quad, -y force
              0;0;1.0;0] #-y-quad
      u_fy = [0;0;0; #main engine
              0;0;1.0;0; #+x-quad, +y force`
              0;0;0;0; #+y-quad
              0;0;1.0;0; #-x-quad, -y force
              0;0;0;0] #-y-quad
      u_fz = [0;0;0; #main engine
              1.0;0;0;0; #+x-quad, +y force`
              1.0;0;0;0; #+y-quad
              1.0;0;0;0; #-x-quad, -y force
              1.0;0;0;0] #-y-quad
      u_τx = [0;0;0; #main engine
              0;0;0;0; #+x-quad
              1.0;0;0;0; #+y-quad
              0;0;0;0; #-x-quad
              0;1.0;0;0] #-y-quad
      u_τy = [0;0;0; #main engine
              0;1.0;0;0; #+x-quad
              0;0;0;0; #+y-quad
              1.0;0;0;0; #-x-quad
              0;0;0;0] #-y-quad
      u_τz = [0;0;0; #main engine
              0;0;1.0;0; #+x-quad, +y force`
              0;0;0;1.0; #+y-quad, -x force
              0;0;0;1.0; #-x-quad, -y force
              0;0;1.0;0] #-y-quad, +x force

      vehicle_dynamics!(ẋ, x, u_fx, params,t)
end

function vehicle_dynamics!(ẋ::AbstractVector,x::AbstractVector,u::AbstractVector,params,t)

      # Parameters:
      m = params[:m] # vehicle mass
      J = params[:J] # vehicle inertia matrix
      m_slosh = params[:m_slosh] # slosh mass

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      s = x[8:10]/norm(x[8:10]) #slosh mass position vector (body frame)
      ṙ = x[11:13] #vehicle velocity (inertial frame)
      ω = x[14:16] #vehicle angular velocity (body frame)
      ṡ = x[17:19] #slosh mass velocity vector (body frame)

      # Inputs:
      θ = u[1] #main engine gimbal angle about body x-axis
      ϕ = u[2] #main engine gimbal angle about body y-axis
      t = u[3:19] #thruster forces (main engine first)

      # Thruster force Jacobian
      Bf = [sin(ϕ)*cos(θ) -sin(θ) cos(ϕ)*cos(θ); #main engine, assuming gimbal rotation is about x followed by y
            0 0 1.0; # +x quad
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
      #TODO: Get thruster locations from CAD model
      rme = [0; 0; -1.5]; #Vector from CoM to main engine
      rxp = [1.0; 0; 0]; #Vector from CoM to +x quad
      rxm = [-1.0; 0; 0]; #Vector from CoM to -x quad
      ryp = [0; 1.0; 0]; #Vector from CoM to +y quad
      rym = [0; -1.0; 0]; #Vector from CoM to -y quad
      Bτ = [hat(rme)*Bf[:,1] hat(rxp)*Bf[:,2:5] hat(ryp)*Bf[:,6:9] hat(rxm)*Bf[:,10:13] hat(rym)*Bf[:,14:17]];

      #Slosh mass stuff
      #TODO: Add spherical pendulum dynamics
      #TODO: Fit pendulum parameters from CFD model of fuel tank
      Fs = [0.0; 0.0; 0.0]

      #Torques (body frame)
      τ = Bτ*t #Torques from all thrusters

      #Forces (body frame)
      Fb = Bf*t #Forces from all thrusters

      #Rotate forces into inertial frame
      F = qrot(q,Fb)

      # Output:
      ẋ[1:3] = ṙ #Vehicle velocity
      ẋ[4] = -0.5*q[2:4]'*ω #Quaternion kinematics (scalar part)
      ẋ[5:7] = 0.5*(q[1]*ω + cross(q[2:4], ω)) #Quaternion kinematics (vector part)
      ẋ[8:10] = ṡ #Slosh mass velocity
      ẋ[11:13] = F/m + gravity(r,t) #Vehicle acceleration
      ẋ[14:16] = J\(τ - cross(ω,J*ω)) #Vehicle angular acceleration
      ẋ[17:19] = Fs/m_slosh #Slosh mass acceleration
end

function hat(x)
      return [0  -x[3] x[2];
             x[3]  0  -x[1];
            -x[2] x[1]  0];
end

function qrot(q, x)
      s = q[1]
      v = q[2:4]
      return x + 2.0*cross(v, cross(v,x) + s*x);
end

function gravity(r,t)
      #μ_e = 398600.0
      #Fg = -μ_e*r/(norm(r)^3) #Just spherical Earth for now
      Fg = [0.0; 0.0; 0.0]
      #J2 =
      #μ_m =
      #r_m =
      #Fg += -μ_m*(r-r_m)/(norm(r-r_m)^3)
end
