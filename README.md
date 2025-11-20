# 🌞⚡ Microgrid Optimization Project

This project focuses on the optimization of a hybrid **microgrid** that includes **solar PV**, **battery storage**, an optional **wind turbine**, and **grid interaction**.
The goal is to determine cost-optimal or emission-optimal system design over different planning horizons while considering uncertainty and multiple operating scenarios.

The project was implemented using **Julia**, **JuMP**, and **Gurobi**, following the tasks from the *Introduction to Numerical Optimization* course.

---

## 📂 Project Scope

The optimization problem is developed progressively across several tasks:

1. **Microgrid cost minimization mathematical model (PV + Battery)**
2. **Simulation for different planning horizons**
3. **Solver performance comparison (Simplex vs Barrier)**
4. **Reduced LP formulation (no explicit SOC variable)**
5. **Dual variables interpretation**
6. **Sensitivity analysis of battery price ( \pi_B )**
7. **CO₂ emissions minimization model and sensitivity analysis of amount of CO₂ produced by a powerplant ( \theta_B )**
8. **Wind turbine integration with 3-scenario uncertainty (expected + worst-case cost)**

All results include numerical outputs and plots of power balance, production, consumption, and battery state of charge (SOC).

---

## ⚙️ System Description

The microgrid includes:

* **PV panel array** producing power proportional to irradiance.
* **Battery storage system** with charging/discharging efficiency and state-of-charge constraints.
* **Wind turbine** (in task 8), producing power proportional to wind availability.
* **Grid connection** enabling buying and selling of electricity.
* **24-hour typical day** repeated for any number of years.

All flows are modeled in **kW** and **kWh**, using a 1-hour time step.

---

## 🎯 Optimization Problems

### **1. Cost minimization**

The objective includes:

* PV investment cost
* Battery investment cost
* Wind turbine investment cost (in task 8)
* Cost of electricity bought from the grid
* Revenue from electricity sold to the grid

### **2. CO₂ emissions minimization**

The objective includes:

* Emissions embedded in manufacturing PV and battery
* Emissions from grid electricity
* (Optional) offset from electricity exported to the grid

### **3. Scenario-based cost minimization with wind**

Three wind scenarios, each with probability **1/3**, are evaluated:

* The objective is **sum of scenario costs + worst-case scenario cost**
* Encourages a robust capacity choice under uncertainty


---

## 📊 Key Results

### **Task 2 — Cost minimization with PV & Battery**

* **1-year horizon:**
  No PV or battery is installed. Grid-only solution is optimal due to insufficient time to recover investment costs.
* **5-year horizon:**
  A small PV installation (~0.5 kWp) becomes optimal; still no battery.

### **Task 3 — Solver comparison**

* **Simplex is significantly faster** than the Barrier method.
* Barrier requires solving dense KKT systems → high computational cost.
* Simplex uses sparsity more efficiently.

### **Task 4 — Reduced model**

* Replacing SOC dynamics with a cumulative energy constraint reduces the number of variables.
* **Both Simplex and Barrier become faster** due to:

  * fewer variables,
  * smaller basis pivots (Simplex),
  * smaller KKT systems (Barrier).

### **Task 5 — Dual variable interpretation**

* **λₜ**: marginal cost of supplying 1 kWh of demand in hour *t*.
* **µₜ**: value of shifting energy from hour *t* to *t+1*.
* **ξ**: shadow price enforcing periodic SOC.

### **Task 6 — Sensitivity of battery price**

* For **1 year**, optimal solution never includes a battery → stability interval ([0, ∞)).
* For **5 years**, buying a battery becomes attractive only when price falls below **20.43 €/kWh**.

### **Task 7 — CO₂ minimization**

* Similar conclusions:
  CO₂-optimal solution includes PV only when grid CO₂ factor exceeds a threshold.
* For a 5-year horizon, PV becomes optimal when
  ( θ_G > 0.13\ \text{kgCO₂/kWh} ).

### **Task 8 — Microgrid with Wind Turbine (3-scenario model)**

#### **1-year horizon**

* No investment in PV, battery, or wind turbine.
* All energy purchased from the grid.
* Identical costs across scenarios → no worst-case differentiation.

#### **5-year horizon**

* Optimal solution includes a **0.642 kWp wind turbine**.
* PV still not selected due to low selling price and poor irradiance.
* Battery not selected because:

  * energy surpluses are small,
  * buying/selling price difference is too small for profitable arbitrage.

Wind power helps significantly reduce operating costs, and the robust objective justifies the turbine investment.

---

## 🚀 How to Run

1. Install Julia and required packages:

   ```julia
   using Pkg
   Pkg.add(["JuMP", "Gurobi", "Plots", "Statistics", "Printf", "MathOptInterface", "BenchmarkTools"])
   ```
2. Ensure your Gurobi license is active.
3. Run a simulation example:

   ```julia
   include("microgrid_opt.jl")
   ```
5. Plots and tables will be saved automatically in the same folder.

---

## 📝 Conclusion

This project demonstrates:

* How linear programming can be used to design renewable microgrids.
* How model structure affects solver performance.
* How uncertainty can be handled via scenario-based optimization.
* How economic and environmental objectives lead to different optimal designs.

The extended model with wind scenarios reveals strong benefits of multi-scenario optimization and highlights when renewable investments become economically and environmentally justified.
