using LinearAlgebra
using DifferentialEquations

function closed_loop(params)
      m_fuel = params[:m_fuel]
      m0 = m_fuel
      r0 = [0; 0; 0] #initial position
      q0 = [1.0; (10*pi/360)*randn(3)] #initial attitude
      q0 = q0/norm(q0)
      v0 = [0; 0; 0] #initial velocity
      ω0 = [0; 0; 0] #initial angular rate
      x0 = [r0; q0; v0; ω0; m0]

      A, qpprob = attitude_tracking_setup(x0, params)
      Fmax = params[:Fmax]
      Fmin = params[:Fmin]
      τ_thrust = params[:τ_thrust]
      q_thrust = params[:q_thrust]
      h = τ_thrust*(2^q_thrust)
      tfinal = 860
      Nt = Int(tfinal/h)
      Nx = 14
      Nu = 4
      thist = collect(h*(0:(Nt-1)))
      xhist = zeros(Nx,Nt)
      uhist = zeros(Nu,Nt)
      xhist[:,1] = x0
      qref = [1.0;0;0;0]
      for k = 1:(Nt-1)
            δu = attitude_tracking(qref, xhist[:,k], A, qpprob, params)
            u = [Fmax;Fmax;Fmax;Fmax] - δu
            if all(u .>= Fmin)
                  uhist[:,k] = u
                  xhist[:,k+1] = rk4_step(xhist[:,k],u,h,params)
            else
                  xhist[:,k+1], uhist[:,k] = quantized_rk_step(xhist[:,k],u,h,q_thrust,params)
            end
      end

      return xhist, uhist, thist
end

function quantized_rk_step(x,u,h,q_thrust,params)
      Fmax = params[:Fmax]

      t_on = round.(u/Fmax,digits=q_thrust,base=2)
      Nsteps = 2^q_thrust
      step = (1/2)^q_thrust
      xn = deepcopy(x)
      for k = 1:Nsteps
            xn = rk4_step(xn,Fmax*(t_on.>=(k*step)),h/Nsteps,params)
      end

      return xn, Fmax*t_on
end

function rk4_step(x,u,h,params)
      ẋ1 = zeros(14)
      ẋ2 = zeros(14)
      ẋ3 = zeros(14)
      ẋ4 = zeros(14)
      vehicle_dynamics!(ẋ1,x,u,params,0)
      vehicle_dynamics!(ẋ2,x+(h/2)*ẋ1,u,params,0)
      vehicle_dynamics!(ẋ3,x+(h/2)*ẋ2,u,params,0)
      vehicle_dynamics!(ẋ4,x+h*ẋ3,u,params,0)
      xn = x + (h/6)*(ẋ1 + 2*ẋ2 + 2*ẋ3 + ẋ4)
      xn[4:7] = xn[4:7]/norm(xn[4:7])
      return xn
end

function simple_test(params)
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
