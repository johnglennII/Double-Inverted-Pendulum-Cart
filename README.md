# Double-Inverted-Pendulum-Cart
*Swing-up and stabilization control of double inverted pendulum on a cart*

## Project Overview
In this project, I derived the system dynamics for an underactuated, highly nonlinear double inverted pendulum on a cart (DIPC) using Lagrangian mechanics and utilized a two-phase control strategy to drive the system from a stable "down-down" initial state to an unstable "up-up" equilibrium. I implemented and compared two different swing-up methods: energy shaping and offline trajectory optimization then used a Linear Model Predictive Controller (MPC) to catch and stabilize the chaotic system once it entered the domain of attraction.

## Key Features
* **Energy Shaping (Passivity-Based Control):** Utilized input-output feedback linearization and passivity-based control to drive the pendulum's energy to the desired equilibrium level, achieving swing-up in 19.1 seconds.
* **Trajectory Optimization (Direct Collocation):** Solved a discrete-time nonlinear programming (NLP) optimization problem using trapezoidal direct collocation to generate an optimal pre-computed trajectory, significantly reducing the swing-up time to just 4.26 seconds.
* **Linear MPC Regulation:** Designed a Linear Model Predictive Controller (LMPC) to stabilize the chaotic dynamics and maintain the desired equilibrium point once the switching conditions were met.
* **Dynamics Modeling:** Derived the full Equations of Motion (EOMs) using Euler-Lagrange equations and converted them into a compact state-space model for advanced control design.

## Visuals
### 1. Control System Comparison

<p align="center">
  <img src="./media/trajectory_optimization/trajectory_optimization_vid.webp" width="49%" alt="Trajectory Optimization Demonstration"/>
  <img src="./media/energy_shaping/energy_shaping_vid.webp" width="49%" alt="Energy Shaping Demonstration"/>
</p>

### 2. Trajectory Optimization: Actual vs Reference (left) Positions, Input (right)

<p align="center">
  <img src="./media/trajectory_optimization/state_traj.png" width="49%" alt="Actual vs Reference Trajectory"/>
  <img src="./media/trajectory_optimization/pos_u.png" width="49%" alt="Positions and Input"/>
</p>

### 3. Energy Shaping: Figures

<p align="center">
  <img src="./media/energy_shaping/theta_u_E.png" width="49%" alt="Theta, u, E"/>
  <img src="./media/energy_shaping/cart_v_u.png" width="49%" alt="Cart, v, u"/>
</p>

## Skills & Software Used
* **Software:** MATLAB, Simulink,
* **Hardware:** N/A
* **Concepts:** Optimal Control, Linear Model Predictive Control (LMPC), Trajectory Optimization (Direct Collocation), Passivity-Based Control, Feedback Linearization, Lagrangian Mechanics
