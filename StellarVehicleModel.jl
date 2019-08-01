using LinearAlgebra

vehicle_params = (m_dry=60.0, #dry mass (kg)
                  m_fuel=48.0, #fuel mass (kg)
                  J_dry=[10.0 0 0; 0 30.0 0; 0 0 30.0], #dry inertia (kg*m^2)
                  Isp=285.0, #Isp (s)
                  g0=9.80665, #Standard gravity used to define Isp (m/s^2)
                  Fmax=10.0, #Max thruster force (N)
                  r_tank=0.5, #fuel tank radius (m)
                  τ_thrust=.0625, #minimum valve switching time (s)
                  q_thrust=2 #number of quantization bits for thrusters
                  )

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
