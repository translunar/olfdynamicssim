
using LinearAlgebra
using DifferentialEquations
using Plots

vehicle_params = (m_dry=60.0, #dry mass (kg)
                  m_fuel=48.0, #wet mass (kg)
                  J_dry=[10.0 0 0; 0 30.0 0; 0 0 30.0], #dry inertia (kg*m^2)
                  Isp=285.0, #Isp (s)
                  g0=9.80665, #Standard gravity used to define Isp (m/s^2)
                  Fmax=10.0, #Max thruster force (N)
                  r_tank=0.5 #fuel tank radius (m)
                  )

function test(params)
      m_fuel = params[:m_fuel]
      m0 = m_fuel
      r0 = [0; 0; 0] #initial position
      q0 = [1.0; 0; 0; 0] #initial attitude
      v0 = [0; 0; 0] #initial velocity
      ω0 = [0; 0; 0] #initial angular rate
      x0 = [r0; q0; v0; ω0; m0]
      prob = ODEProblem(open_loop!, x0, (0.0, 10.0), vehicle_params)
      soln = solve(prob)
      return soln
end

function open_loop!(ẋ,x,params,t)
      #Just a function to test open-loop inputs
      Fmax = params[:Fmax]

      #No inputs
      u0 = zeros(4)

      #Test inputs for thrusters
      u_fx = [Fmax; Fmax; Fmax; Fmax] #Max forward acceleration
      u_τx = [Fmax; 0; Fmax; 0] #Max x-axis (roll) torque
      u_τy = [0; 0; Fmax; Fmax] #Max y-axis (pitch) torque
      u_τz = [Fmax; 0; 0; Fmax] #Max z-axis (yaw) torque

      vehicle_dynamics!(ẋ, x, u_τy, params,t)
end

function vehicle_dynamics!(ẋ,x,u,params,t)

      # Parameters:
      #TODO: Add thruster geometry to params struct
      m_dry = params[:m_dry]
      J_dry = params[:J_dry]
      Isp = params[:Isp]
      g0 = params[:g0]
      r_tank = params[:r_tank]

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      ṙ = x[8:10] #vehicle velocity (inertial frame)
      ω = x[11:13] #vehicle angular velocity (body frame)
      m_fuel = x[14] #current fuel mass

      #Clip inputs (no negative thrust or fuel allowed)
      m_fuel = max(m_fuel, 0.0)
      if m_fuel == 0.0
            u = zeros(4)
      else
            u = max.(u, 0.0)
      end

      m = m_dry + m_fuel #total vehicle mass

      #Total vehicle inertia assuming spherical fuel mass at COM
      J = J_dry + (2*m_fuel*r_tank*r_tank/5)*I

      # Some rotation stuff we're going to use a lot
      R = qtoR(q) #convert quaternion to rotation matrix
      ω̂ = hat(ω) #skew-symmetric cross product matrix

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

      #Torques (body frame)
      τ = Bτ*u - cross(ω,J*ω) #Torques from all thrusters + gyroscopic term

      #Forces (inertial frame)
      F = R*Bf*u #+ m*gravity(r,t) #Forces from all thrusters

      # Output:
      ẋ[1:3] = ṙ #Vehicle velocity
      ẋ[4] = -0.5*q[2:4]'*ω #Quaternion kinematics (scalar part)
      ẋ[5:7] = 0.5*(q[1]*ω + cross(q[2:4], ω)) #Quaternion kinematics (vector part)
      ẋ[8:10] = F/m #Vehicle acceleration
      ẋ[11:13] = J\τ #Vehicle angular acceleration
      ẋ[14] = -sum(u)/(g0*Isp) #Fuel use
end

function tesseract_test(params)
      p = params[:p_slosh] # slosh pivot point
      l = params[:l_slosh]

      r0 = [0; 0; 0]
      q0 = [1.0; 0; 0; 0]
      s0 = r0 + p - [0; 0; l]
      v0 = [0; 0; 0]
      ω0 = [0; 0; 0]
      ṡ0 = v0
      x0 = [r0; q0; s0; v0; ω0; ṡ0]
      prob = ODEProblem(tesseract_open_loop!, x0, (0.0, 20.0), vehicle_params)
      soln = solve(prob)
      plot(soln)
end

function tesseract_open_loop!(ẋ,x,params,t)
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

      tesseract_vehicle_dynamics_slosh!(ẋ, x, u_fx, params,t)
end

function tesseract_vehicle_dynamics_slosh!(ẋ,x,u,params,t)

      # Parameters:
      #TODO: Add thruster geometry to params struct
      m = params[:m] # vehicle mass
      J = params[:J] # vehicle inertia matrix
      p_s = params[:p_slosh] # slosh pivot point
      m_s = params[:m_slosh] # slosh mass
      l_s = params[:l_slosh] # slosh pendulum length

      # Mass matrix for entire model
      M = [m*I zeros(3,6); zeros(3,3) J zeros(3,3); zeros(3,6) m_s*I]
      Minv = inv(M)

      # State:
      r = x[1:3] #vehicle position (inertial frame)
      q = x[4:7]/norm(x[4:7]) #quaternion (body to inertial, scalar first)
      s = x[8:10] #slosh mass position vector (inertial frame)
      ṙ = x[11:13] #vehicle velocity (inertial frame)
      ω = x[14:16] #vehicle angular velocity (body frame)
      ṡ = x[17:19] #slosh mass velocity vector (inertial frame)

      # Inputs:
      θ = u[1] #main engine gimbal angle about body x-axis
      ϕ = u[2] #main engine gimbal angle about body y-axis
      t = u[3:19] #thruster forces (main engine first)

      # Some rotation stuff we're going to use a lot
      R = qtoR(q) #convert quaternion to rotation matrix
      ω̂ = hat(ω) #skew-symmetric cross product matrix

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

      #Torques (body frame)
      τ = Bτ*t - cross(ω,J*ω) #Torques from all thrusters + gyroscopic term

      #Forces (inertial frame)
      F = R*Bf*t #+ m*gravity(r,t) #Forces from all thrusters + gravity

      #Fuel slosh pendulum stuff
      k = 10 # Constraint stabilization gain

      y = r + R*p_s - s #mass to pivot point vector (inertial frame)
      ẏ = ṙ + R*ω̂*p_s - ṡ

      ϕ = y'*y - l_s*l_s; #constraint
      ϕ̇ = 2*y'*ẏ; #derivative of constraint

      c = ẏ'*ẏ + y'*R*ω̂*ω̂*p_s
      G = [y' -y'*R*hat(p_s) -y']

      Fs = [0.0; 0.0; 0.0] #TODO: apply gravity to pendulum

      λ = -(G*Minv*G')\(G*Minv*[F; τ; Fs] + c + k*k*ϕ + 2*k*ϕ̇);

      # Output:
      ẋ[1:3] = ṙ #Vehicle velocity
      ẋ[4] = -0.5*q[2:4]'*ω #Quaternion kinematics (scalar part)
      ẋ[5:7] = 0.5*(q[1]*ω + cross(q[2:4], ω)) #Quaternion kinematics (vector part)
      ẋ[8:10] = ṡ #Slosh mass velocity
      ẋ[11:19] = Minv*(G'*λ + [F; τ; Fs]) #Accelerations
      # ẋ[11:13] = F/m #Vehicle acceleration
      # ẋ[14:16] = J\τ #Vehicle angular acceleration
      # ẋ[17:19] = Fs/m_slosh #Slosh mass acceleration
end

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
