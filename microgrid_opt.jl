using JuMP, Gurobi, Plots, Statistics, MathOptInterface, Printf, BenchmarkTools



# ================
# ==== TASK 2 ====
# ================
function simulate_microgrid(n_years::Int; plot_overview::Bool=true, method::Int=0, pi_B::Float64=500.0)
    include("data.jl")


    # PARAMETERS
    P_C = consumption ./ 1000.0  # W → kW  # consumption of the load
    i = irradiance
    T = 24 # 24h in a day
    N_days = 365 * n_years

    eta_PV = 0.86 # PV panels efficiency [%]
    eta_B  = 0.95 # battery efficiency [%]
    pi_PV = 800.0 # cost per installed capacity PV [$/kWp]
    pi_G_plus  = 0.1 # cost per unit of electricity bought [$/kWh]
    pi_G_minus = 0.02 # revenue per unit of electricity sold [$/kWh]
    # Also: pi_B = 500 -> cost per installed capacity of battery [€/kWhp] 
    # sent to the function because of the next tasks


    # MODEL 
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", method) # using simplex or barrier method


    # VARIABLES
    @variable(model, CPV >= 0) # installed PV capacity [kWp]
    @variable(model, E_B >= 0, base_name="E_B") # (energy) capacity of the battery [kWhp]
    @variable(model, P_PV[1:T] >= 0) # the quantity of power the panels produce [kw]
    @variable(model, P_B_ch[1:T] >= 0) # battery charging power [kW]
    @variable(model, P_B_dis[1:T] >= 0) # battery discharging power [kW]
    @variable(model, SOC[1:T+1] >= 0) # battery state of charge level [kWh]
    @variable(model, P_G_plus[1:T] >= 0) # energy bought from the network [kW]
    @variable(model, P_G_minus[1:T] >= 0) # energy sold to the network [kW]


    # CONSTRAINTS
    # constraint on the energy produced by the panels (taking into account possible curtailment)
    for t in 1:T
        @constraint(model, P_PV[t] <= eta_PV * i[t] * CPV)
    end

    # energy balance
    for t in 1:T
        @constraint(model, P_PV[t] + P_G_plus[t] + P_B_dis[t] == P_C[t] + P_G_minus[t] + P_B_ch[t])
    end

    # SOC hour by hour
    for t in 1:T
        @constraint(model, SOC[t+1] == SOC[t] + eta_B * P_B_ch[t] - (1/eta_B) * P_B_dis[t])
    end
    # SOC bounds
    for t in 1:(T+1)
        @constraint(model, 0 <= SOC[t])
        @constraint(model, SOC[t] <= E_B)  
    end
    # SOC periodicity
    @constraint(model, SOC[1] == SOC[T+1])


    # OBJECTIVE
    invest_cost = pi_PV*CPV + pi_B*E_B
    operating_cost = N_days * sum(pi_G_plus*P_G_plus[t] - pi_G_minus*P_G_minus[t] for t in 1:T)
    @objective(model, Min, invest_cost + operating_cost)


    # SENSITIVITY ANALYSIS
    # Save the model to lp file (for sensitivity analysis)
    JuMP.write_to_file(model, "microgrid_pi_B_$(n_years).lp")

    # OPTIMIZATION
    optimize!(model)
    status = termination_status(model) # we check whether the model has returned optimal results to avoid errors


    # DOWNLOADING MODEL RESULTS
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        CPV_val = value(CPV)
        E_B_val = value(E_B)
        P_PV_val = value.(P_PV)
        P_B_ch_val = value.(P_B_ch)
        P_B_dis_val = value.(P_B_dis)
        P_G_plus_val = value.(P_G_plus)
        P_G_minus_val = value.(P_G_minus)
        SOC_val = value.(SOC)
        objective_val = objective_value(model)
    else
        @warn "!!! Solver did not finish optimally, results unavailable !!!"
        CPV_val = NaN
        E_B_val = NaN
        P_PV_val = fill(NaN, length(P_PV))
        P_B_ch_val = fill(NaN, length(P_B_ch))
        P_B_dis_val = fill(NaN, length(P_B_dis))
        P_G_plus_val = fill(NaN, length(P_G_plus))
        P_G_minus_val = fill(NaN, length(P_G_minus))
        SOC_val = fill(NaN, length(SOC))
        objective_val = NaN
    end

    # OUTPUT OF RESULTS AND SUMMARY PLOT
    if plot_overview
        # display results
        println("\n" * "="^55)
        println(" MICROGRID – COST MINIMIZATION RESULTS")
        println("="^55)
        println(" Simulation period : ", n_years, " years")
        println(" Status            : ", termination_status(model))
        println("-"^55)
        @printf(" %-25s : %10.3f kWp\n", "PV capacity (CPV)", CPV_val)
        @printf(" %-25s : %10.3f kWh\n", "Battery capacity (E_B)", E_B_val)
        @printf(" %-25s : %10.2f\n", "Objective value", objective_val)
        println("="^55 * "\n")

        # plot
        hours = 1:(7*T)
        P_PV_week = repeat(P_PV_val, 7)
        P_G_plus_week = repeat(P_G_plus_val, 7)
        P_G_minus_week = repeat(P_G_minus_val, 7)
        P_C_week = repeat(P_C, 7)
        SOC_week = repeat(SOC_val[1:T], 7)

        p = plot(hours, P_PV_week,
            lw=7, color=:orange, label="PV production",
            xlabel="Hour of the week", ylabel="Power [kW]",
            title="Microgrid Power Flow and Battery SOC ($(n_years)_years)",
            legend=:topright, grid=true)

        plot!(p, hours, P_G_plus_week, lw=5, color=:blue, label="Power bought")
        plot!(p, hours, P_G_minus_week, lw=5, color=:red, label="Power sold")
        plot!(p, hours, P_C_week, lw=3, color=:green, label="Load consumption")
        plot!(p, hours, SOC_week, lw=1, color=:purple, label="Battery SOC", ylabel="Power [kW]", right_ylabel="Energy [kWh]")

        savefig(p, "microgrid_overview_$(n_years)_years.png")
        display(p)
        println("Plot saved as: microgrid_overview_$(n_years)_years.png")
    end

    return (CPV_val=CPV_val, E_B_val=E_B_val, objective=objective_val)
end

# RUN SIMULATIONS
simulate_microgrid(1)
simulate_microgrid(5)
#simulate_microgrid(1, plot_overview=true, method=0, pi_B=0.0)
#simulate_microgrid(5, plot_overview=true, method=0, pi_B=20.0)
simulate_microgrid(20)



# ================
# ==== TASK 4 ====
# ================
function simulate_microgrid_reduced(n_years::Int; print_summary::Bool=true, method::Int=0, pi_B::Float64=500.0)
    include("data.jl") 


    # PARAMETERS
    P_C = consumption ./ 1000.0 # W → kW  # consumption of the load
    i = irradiance
    T = 24
    N_days = 365 * n_years

    eta_PV = 0.86 # PV panels efficiency [%]
    eta_B  = 0.95 # battery efficiency [%]
    pi_PV = 800.0 # cost per installed capacity PV [$/kWp]
    pi_G_plus  = 0.1 # cost per unit of electricity bought [$/kWh]
    pi_G_minus = 0.02 # revenue per unit of electricity sold [$/kWh]
    # Also: pi_B = 500 [$/kWp] sent to the function because of the next tasks


    # MODEL WITHOUT EXPLICIT SOC
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", method)  # using simplex or barrier method


    # VARIABLES
    @variable(model, CPV >= 0) # installed PV capacity [kWp]
    @variable(model, E_B >= 0) # (energy) capacity of the battery [kWhp]
    @variable(model, P_PV[1:T] >= 0) # the quantity of power the panels produce [kW]
    @variable(model, P_B_ch[1:T] >= 0) # battery charging power [kW]
    @variable(model, P_B_dis[1:T] >= 0) # battery discharging power [kW]
    @variable(model, P_G_plus[1:T] >= 0) # energy bought from the network [kW]
    @variable(model, P_G_minus[1:T] >= 0) # energy sold to the network [kW]
    @variable(model, SOC0 >= 0) # initial state of charge

    @constraint(model, SOC0 <= E_B)

    # for each hour t: SOC0 + cumulated sum in [0, E_B]
    for t in 1:T
        @constraint(model, SOC0 + sum(eta_B * P_B_ch[k] - (1/eta_B) * P_B_dis[k] for k in 1:t) >= 0)
        @constraint(model, SOC0 + sum(eta_B * P_B_ch[k] - (1/eta_B) * P_B_dis[k] for k in 1:t) <= E_B)
    end

    # periodicity 
    @constraint(model, sum(eta_B * P_B_ch[k] - (1/eta_B) * P_B_dis[k] for k in 1:T) == 0)


    # CONSTRAINTS
    # constraint on the energy produced by the panels (taking into account possible curtailment)
    for t in 1:T
        @constraint(model, P_PV[t] <= eta_PV * i[t] * CPV)
    end

    # energy balance
    for t in 1:T
        @constraint(model, P_PV[t] + P_G_plus[t] + P_B_dis[t] == P_C[t] + P_B_ch[t] + P_G_minus[t])
    end

    # OBJECTIVE
    invest_cost = pi_PV*CPV + pi_B*E_B
    operating_cost = N_days * sum(pi_G_plus*P_G_plus[t] - pi_G_minus*P_G_minus[t] for t in 1:T)
    @objective(model, Min, invest_cost + operating_cost)


    # PTIMIZATION
    optimize!(model)
    status = termination_status(model) # we check whether the model has returned optimal results to avoid errors
    

    # DOWNLOADING MODEL RESULTS
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        CPV_val = value(CPV)
        E_B_val = value(E_B)
        objective_val = objective_value(model)

        # OUTPUT OF RESULTS
        if print_summary
            println("\n" * "="^55)
            println(" MICROGRID – COST MINIMIZATION (REDUCED SOC) ")
            println("="^55)
            println(" Simulation period : ", n_years, " years")
            println(" Status            : ", termination_status(model))
            println("-"^55)
            @printf(" %-25s : %10.3f kWp\n", "PV capacity (CPV)", CPV_val)
            @printf(" %-25s : %10.3f kWh\n", "Battery capacity (E_B)", E_B_val)
            @printf(" %-25s : %10.2f\n", "Objective value", objective_val)
            println("="^55 * "\n")
        end
    else
        @warn "!!! Solver did not finish optimally, results unavailable !!!"
        CPV_val = NaN
        E_B_val = NaN
        objective_val = NaN
    end

    return (CPV_val=CPV_val, E_B_val=E_B_val, objective=objective_val)
end

# RUN SIMULATIONS FOR REDUCED MODEL
simulate_microgrid_reduced(1)
simulate_microgrid_reduced(5)
simulate_microgrid_reduced(20)




# ================
# ==== TASK 3 ====
# ================
function solver_times(years_list::Vector{Int}, method_t::Int; use_reduced::Bool=false)
    times = Float64[]
    
    for y in years_list
        times_trials = Float64[]
        if use_reduced
            elapsed_s = @belapsed simulate_microgrid_reduced($y; print_summary=false, method=$method_t) samples=8 evals=4
            push!(times_trials, elapsed_s)
        else
            elapsed_s = @belapsed simulate_microgrid($y; plot_overview=false, method=$method_t) samples=8 evals=4
            push!(times_trials, elapsed_s)
        end
        avg_time = mean(times_trials)
        push!(times, avg_time)
    end
    
    # CHOSE METHOD
    method_name = method_t == 0 ? "Simplex" : method_t == 2 ? "Barrier" : "Unknown"
    model_type = use_reduced ? "Reduced model" : "Full model"
    
    println("\n" * "="^50)
    println(" SOLVER TIME RESULTS: $model_type, Method = $method_name ($method_t)")
    println("="^50)
    @printf("%-10s | %-10s\n", "Years", "Solve time [ms]")
    println("-"^25)
    for (i, y) in enumerate(years_list)
        @printf("%-10d | %-15.4f\n", y, times[i] * 1000)
    end
    println("="^50 * "\n")

    return times
end

# TEST
# years_to_test = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
# # full model
# times_full_simplex = solver_times(years_to_test, 0; use_reduced=false)
# times_full_barrier = solver_times(years_to_test, 2; use_reduced=false)
# # reduced model
# times_reduced_simplex = solver_times(years_to_test, 0; use_reduced=true)
# times_reduced_barrier= solver_times(years_to_test, 2; use_reduced=true)
# # plot all results together
# p = plot(years_to_test, times_full_simplex, 
#     lw=2, marker=:o, color=:blue, label="Full model - Simplex")

# plot!(p, years_to_test, times_full_barrier, 
#     lw=2, marker=:square, color=:red, label="Full model - Barrier")

# plot!(p, years_to_test, times_reduced_simplex, 
#     lw=2, marker=:diamond, color=:green, label="Reduced model - Simplex")

# plot!(p, years_to_test, times_reduced_barrier, 
#     lw=2, marker=:utriangle, color=:purple, label="Reduced model - Barrier")

# xlabel!("Years of operation")
# ylabel!("Solve time [s]")
# title!("Microgrid LP: Solver Time Comparison (Full vs Reduced)")


# savefig("solver_time_comparison_all.png")
# display(p)
# println("Plot saved as: solver_time_comparison_all.png")



# ================
# ==== TASK 7 ====
# ================
function simulate_microgrid_co2(n_years::Int; theta_G::Float64 = 0.1, plot_overview::Bool=true, find_threshold::Bool=false)

    include("data.jl") 
    P_C = consumption ./ 1000.0 # W → kW  # consumption of the load
    i = irradiance
    T = 24
    N_days = 365 * n_years

    eta_PV = 0.86 # PV panels efficiency [%]
    eta_B  = 0.95 # battery efficiency [%]
    theta_PV = 1000.0 # CO2 emissions of PV panels production [kgCO2/kWp]
    theta_B = 150.0  # CO2 emissions of battery production [kgCO2/kWhp]
    # Also: theta_G = 0.1  -> CO2 emissions of grid electricity production [kgCO2/kWh]
    # sent to the function because of the next tasks

    # MODEL MINIMIZING CO2 EMISSIONS
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", 0) 



    # VARIABLES
    @variable(model, CPV >= 0) # installed PV capacity [kWp]
    @variable(model, E_B >= 0) # (energy) capacity of the battery [kWhp]          
    @variable(model, P_PV[1:T] >= 0) # the quantity of power the panels produce [kW]    
    @variable(model, P_B_ch[1:T] >= 0) # battery charging power [kW]
    @variable(model, P_B_dis[1:T] >= 0) # battery discharging power [kW]
    @variable(model, SOC[1:T+1] >= 0) # battery state of charge level [kWh]
    @variable(model, P_G_plus[1:T] >= 0) # energy bought from the network [kW]
    @variable(model, P_G_minus[1:T] >= 0) # energy sold to the network [kW]

    # CONSTRAINTS
    # constraint on the energy produced by the panels (taking into account possible curtailment)
    for t in 1:T
        @constraint(model, P_PV[t] <= eta_PV * i[t] * CPV)
    end

    # energy balance
    for t in 1:T
        @constraint(model, P_PV[t] + P_G_plus[t] + P_B_dis[t] == P_C[t] + P_G_minus[t] + P_B_ch[t])
    end

    # SOC hour by hour
    for t in 1:T
        @constraint(model, SOC[t+1] == SOC[t] + eta_B * P_B_ch[t] - (1.0/eta_B) * P_B_dis[t])
    end
    # SOC bounds
    for t in 1:(T+1)
        @constraint(model, 0 <= SOC[t])
        @constraint(model, SOC[t] <= E_B)
    end
    
    # SOC periodicity
    @constraint(model, SOC[1] == SOC[T+1])


    # OBJECTIVE
    pv_emissions = theta_PV * CPV 
    battery_emissions = theta_B * E_B

    grid_emissions = N_days * sum(theta_G * P_G_plus[t] for t in 1:T)
    @objective(model, Min, pv_emissions + battery_emissions + grid_emissions)

    # grid_emission_decrease = N_days * sum(theta_G * P_G_minus[t] for t in 1:T) # sold to the network
    # @constraint(model, pv_emissions + battery_emissions + grid_emissions - grid_emission_decrease >= 0)
    #@objective(model, Min, pv_emissions + battery_emissions + grid_emissions - grid_emission_decrease)
    
    JuMP.write_to_file(model, "microgrid_theta_$(n_years).lp")
    # OPTIMIZATION
    optimize!(model)
    status = termination_status(model)

    # DOWNLOADING MODEL RESULTS
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        CPV_val = value(CPV)
        E_B_val  = value(E_B)
        P_PV_val = value.(P_PV)
        P_B_ch_val = value.(P_B_ch)
        P_B_dis_val = value.(P_B_dis)
        P_G_plus_val = value.(P_G_plus)
        P_G_minus_val = value.(P_G_minus)
        SOC_val = value.(SOC)
        objective_val = objective_value(model)  # total CO2 emission

        pv_emissions_val = theta_PV * CPV_val
        battery_emissions_val = theta_B * E_B_val
        grid_emissions_val = N_days * sum(theta_G * P_G_plus_val[t] for t in 1:T)

        #grid_emission_decrease_val = N_days * sum(theta_G * P_G_minus_val[t] for t in 1:T)
        # total_emissions_val = pv_emissions_val + battery_emissions_val + grid_emissions_val - grid_emission_decrease_val

        # OUTPUT OF RESULTS AND SUMMARY PLOT
        if plot_overview
            # display results
            println("\n" * "="^60)
            println(" MICROGRID – CO₂ MINIMIZATION RESULTS")
            println("="^60)
            println(" Simulation period : ", n_years, " years (", N_days, " days)")
            println(" Solver status     : ", status)
            println("-"^60)
            @printf(" %-35s : %10.3f kWp\n", "PV installed capacity", CPV_val)
            @printf(" %-35s : %10.3f kWh\n", "Battery capacity", E_B_val)
            println("-"^60)
            @printf(" %-35s : %10.3f kg CO₂\n", "Total CO₂ emissions (objective)", objective_val)
            @printf("   %-33s : %10.3f kg CO₂\n", "PV manufacturing", pv_emissions_val)
            @printf("   %-33s : %10.3f kg CO₂\n", "Battery manufacturing", battery_emissions_val)
            @printf("   %-33s : %10.3f kg CO₂\n", "Grid electricity (bought)", grid_emissions_val)
            # @printf("   %-33s : %10.3f kg CO₂\n", "Grid electricity (sold)", grid_emission_decrease_val)
            println("="^60 * "\n")

            # plot
            hours = 1:(7*T)
            P_PV_week = repeat(P_PV_val, 7)
            P_G_plus_week = repeat(P_G_plus_val, 7)
            P_G_minus_week = repeat(P_G_minus_val, 7)
            P_C_week = repeat(P_C, 7)
            SOC_week = repeat(SOC_val[1:T], 7)

            p = plot(hours, P_PV_week,
                lw=7, label="PV production",
                xlabel="Hour of the week", ylabel="Power [kW]",
                title="Microgrid Power Flow and Battery SOC (CO2-min) ($(n_years) years)",
                legend=:topright, grid=true)

            plot!(p, hours, P_G_plus_week, lw=5, label="Power bought (grid → microgrid)")
            plot!(p, hours, P_G_minus_week, lw=5, label="Power sold (microgrid → grid)")
            plot!(p, hours, P_C_week, lw=2, label="Load consumption")
            plot!(p, hours, SOC_week, lw=1, label="Battery SOC", ylabel="Power [kW]", right_ylabel="Energy [kWh]")

            savefig(p, "microgrid_co2_overview_$(n_years)_years.png")
            display(p)
            println("Plot saved as: microgrid_co2_overview_$(n_years)_years.png")
        end
    else
        @warn "!!! Solver did not finish optimally, results unavailable !!!"
        CPV_val = NaN
        E_B_val = NaN
        P_PV_val = fill(NaN, length(P_PV))
        P_B_ch_val = fill(NaN, length(P_B_ch))
        P_B_dis_val = fill(NaN, length(P_B_dis))
        P_G_plus_val = fill(NaN, length(P_G_plus))
        P_G_minus_val = fill(NaN, length(P_G_minus))
        SOC_val = fill(NaN, length(SOC))
        objective_val = NaN
    end

    if find_threshold
        upper_threshold = nothing
        for theta in theta_G:0.01:(theta_G*20)
            pv_emissions = theta_PV * CPV 
            battery_emissions = theta_B * E_B
            grid_emissions = N_days * sum(theta * P_G_plus[t] for t in 1:T)
            #grid_emission_decrease = N_days * sum(theta * P_G_minus[t] for t in 1:T)
            # @objective(model, Min, pv_emissions + battery_emissions + grid_emissions - grid_emission_decrease)

            @objective(model, Min, pv_emissions + battery_emissions + grid_emissions)
            optimize!(model)

            CPV_this_theta = value(CPV)
            E_B_this_theta  = value(E_B)

            if CPV_this_theta != CPV_val || E_B_this_theta != E_B_val
                upper_threshold = theta
                break
            end
        end

        lower_threshold = nothing
        for theta in theta_G:-0.01:0.0
            pv_emissions = theta_PV * CPV 
            battery_emissions = theta_B * E_B
            grid_emissions = N_days * sum(theta * P_G_plus[t] for t in 1:T)
            #grid_emission_decrease = N_days * sum(theta * P_G_minus[t] for t in 1:T)
            # @objective(model, Min, pv_emissions + battery_emissions + grid_emissions - grid_emission_decrease)
            
            @objective(model, Min, pv_emissions + battery_emissions + grid_emissions)
            optimize!(model)

            CPV_this_theta = value(CPV)
            E_B_this_theta  = value(E_B)

            if CPV_this_theta != CPV_val || E_B_this_theta != E_B_val
                lower_threshold = theta
                break
            end
        end
        println("\n" * "="^55)
        println("\nSensitivity interval for theta_G where PV and battery decisions do not change:")
        println("Lower bound: ", lower_threshold)
        println("Upper bound: ", upper_threshold)
        println("Reference decisions: CPV =", CPV_val, ", E_B =", E_B_val)
        println("\n" * "="^55)
    end



    return Dict(
        :CPV_kWp => CPV_val,
        :E_B_kWh => E_B_val,
        # :total_CO2_kg => total_emissions_val,
        :pv_CO2_kg => pv_emissions_val,
        :battery_CO2_kg => battery_emissions_val,
        :grid_CO2_kg => grid_emissions_val,
        :P_PV => P_PV_val,
        :P_G_plus => P_G_plus_val,
        :P_G_minus => P_G_minus_val,
        :P_B_ch => P_B_ch_val,
        :P_B_dis => P_B_dis_val,
        :SOC => SOC_val,
        :termination_status => status
    )
end

# RUN SIMULATIONS (CO2)
result1 = simulate_microgrid_co2(1, theta_G=0.1, plot_overview=true, find_threshold=true)
result5 = simulate_microgrid_co2(5, theta_G=0.1, plot_overview=true, find_threshold=true)

#result1m = simulate_microgrid_co2(1, theta_G=0.60, plot_overview=true, find_threshold=false)
#result1m = simulate_microgrid_co2(1, theta_G=0.62, plot_overview=true, find_threshold=false)
# result5m = simulate_microgrid_co2(5, theta_G=0.12, plot_overview=true, find_threshold=false)
#result5m = simulate_microgrid_co2(5, theta_G=0.14, plot_overview=true, find_threshold=false)



# ================
# ==== TASK 8 ====
# ================
function simulate_microgrid_with_wind(n_years::Int; plot_overview::Bool=true, method::Int=0, pi_W::Float64=1500.0)
    include("data.jl")

    # PARAMETERS
    P_C = consumption ./ 1000.0 # W → kW  # consumption of the load
    i = irradiance
    T = 24
    S = 3 # number of scenarios
    wind = vcat(wind_1', wind_2', wind_3') # 3×24 matrix: 3 scenarios, 24 hours
    size(wind) # should be (3,24)
    N_days = 365 * n_years

    eta_PV = 0.86 # PV panels efficiency [%]
    eta_B  = 0.95 # battery efficiency [%]
    pi_PV = 800.0 # cost per installed capacity PV [$/kWp]
    pi_B  = 500.0 # cost per installed capacity of battery [$/kWhp] 
    pi_G_plus  = 0.1 # cost per unit of electricity bought [$/kWh]
    pi_G_minus = 0.02 # revenue per unit of electricity sold [$/kWh]


    # MODEL WITH WIND TURBINE
    model = Model(Gurobi.Optimizer)
    set_silent(model)
    set_optimizer_attribute(model, "Method", method) # using simplex or barrier method


    # VARIABLES
    @variable(model, CPV >= 0) # installed PV capacity [kWp]
    @variable(model, E_B >= 0) # (energy) capacity of the battery [kWhp]
    @variable(model, CW >= 0) # installed wind turbine capacity [kWp]

    # scenario-specific variables
    @variable(model, P_PV[s=1:S, t=1:T] >= 0) # the quantity of power the panels produce [kW]
    @variable(model, P_B_ch[s=1:S, t=1:T] >= 0) # battery charging power [kW]
    @variable(model, P_B_dis[s=1:S, t=1:T] >= 0) # battery discharging power [kW]
    @variable(model, SOC[s=1:S, t=1:T+1] >= 0) # battery state of charge level [kWh]
    @variable(model, P_G_plus[s=1:S, t=1:T] >= 0) # energy bought from the network [kW]
    @variable(model, P_G_minus[s=1:S, t=1:T] >= 0) # energy sold to the network [kW]

    #wind production (cannot be curtailed) as expression
    @expression(model, P_W[s=1:S, t=1:T], wind[s,t] * CW)

    # CONSTRAINTS
    # constraint on the energy produced by the panels (taking into account possible curtailment), same CPV across scenarios
    for s in 1:S, t in 1:T
        @constraint(model, P_PV[s,t] <= eta_PV * i[t] * CPV)
    end

    # energy balance per scenario
    for s in 1:S, t in 1:T
        @constraint(model, P_PV[s,t] + P_W[s,t] + P_G_plus[s,t] + P_B_dis[s,t] == P_C[t] + P_G_minus[s,t] + P_B_ch[s,t])
    end

    # SOC hour by hour per scenario
    for s in 1:S, t in 1:T
        @constraint(model, SOC[s,t+1] == SOC[s,t] + eta_B * P_B_ch[s,t] - (1/eta_B) * P_B_dis[s,t])
    end

    # SOC bounds (same E_B across scenarios)
    for s in 1:S, t in 1:T+1
        @constraint(model, 0 <= SOC[s,t])
        @constraint(model, SOC[s,t] <= E_B)
    end
    # SOC periodicity
    for s in 1:S
        @constraint(model, SOC[s,1] == SOC[s,T+1])
    end

    
    # OBJECTIVE
    invest_cost = pi_PV * CPV + pi_B * E_B + pi_W * CW # shared investment cost
    @expression(model, operating_cost[s=1:S], N_days * sum(pi_G_plus * P_G_plus[s,t] - pi_G_minus * P_G_minus[s,t] for t in 1:T))
    @expression(model, total_cost[s=1:S], invest_cost + operating_cost[s])

    # z models the worst-case cost (max_{s} total_cost[s])
    @variable(model, z >= 0)
    for s in 1:S
        @constraint(model, z >= total_cost[s])
    end

    # objective: sum of scenario total costs + worst-case (z)
    @objective(model, Min, sum(total_cost[s] for s in 1:S) + z)


    # OPTIMIZATION
    optimize!(model)
    status = termination_status(model)


    # DOWNLOADING MODEL RESULTS
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        CPV_val = value(CPV)
        E_B_val  = value(E_B)
        CW_val   = value(CW)
        P_PV_val = value.(P_PV)
        P_B_ch_val = value.(P_B_ch)
        P_B_dis_val = value.(P_B_dis)
        P_G_plus_val = value.(P_G_plus)
        P_G_minus_val = value.(P_G_minus)
        SOC_val = value.(SOC)
        total_cost_vals = [value(total_cost[s]) for s in 1:S]
        operating_cost_vals = [value(operating_cost[s]) for s in 1:S]
        z_val = value(z)
        objective_val = objective_value(model)
        if plot_overview
            # display results
            println("\n" * "="^70)
            println(" MICROGRID WITH WIND — RESULTS")
            println("="^70)
            println(" Simulation period : ", n_years, " years (", N_days, " days)")
            println(" Solver status     : ", status)
            println("-"^70)
            @printf(" %-30s : %10.3f kWp\n", "PV capacity (CPV)", CPV_val)
            @printf(" %-30s : %10.3f kWh\n", "Battery capacity (E_B)", E_B_val)
            @printf(" %-30s : %10.3f kWp\n", "Wind capacity (CW)", CW_val)
            println("-"^70)
            for s in 1:S
                @printf(" Scenario %d: total cost = %12.8f  , operating = %12.8f\n", s, total_cost_vals[s], operating_cost_vals[s])
            end
            @printf(" Worst-case (z)     : %12.3f\n", z_val)
            @printf(" Objective (sum+worst): %12.3f\n", objective_val)
            println("="^70 * "\n")

            # plot
            hours = 1:(7*24)
            P_PV_week = repeat(P_PV_val[1, :], 7)
            P_W_week  = repeat([wind[1,t]*CW_val for t in 1:24], 7)
            P_G_plus_week = repeat(P_G_plus_val[1, :], 7)
            P_G_minus_week = repeat(P_G_minus_val[1, :], 7)
            P_C_week = repeat(P_C, 7)
            SOC_week = repeat(SOC_val[1,1:24], 7)
            
            p = plot(hours, P_W_week, lw=6, color=:cyan, label="Wind production")
            plot!(p, hours, P_PV_week, lw=5, color=:orange, label="PV production")
            plot!(p, hours, P_G_plus_week, lw=4, color=:blue, label="Power bought (grid → microgrid)")
            plot!(p, hours, P_G_minus_week, lw=4, color=:red, label="Power sold (microgrid → grid)")
            plot!(p, hours, P_C_week, lw=3, color=:green, label="Load consumption")
            plot!(p, hours, SOC_week, lw=2, color=:purple, label="Battery SOC", ylabel="Power [kW]", right_ylabel="Energy [kWh]")
            
            xlabel!("Hour of the week")
            ylabel!("Power [kW]")
            title!("Microgrid Power Flow with Wind (Scenario 1) ($(n_years) years)")
            
            savefig(p, "microgrid_wind_overview_$(n_years)_years.png")
            display(p)
            println("Plot saved as: microgrid_wind_overview_$(n_years)_years.png")
        end
    else
        @warn "!!! Solver did not finish optimally, results unavailable !!!"
        CPV_val = NaN
        E_B_val = NaN
        CW_val = NaN
        P_PV_val = fill(NaN, length(P_PV))
        P_B_ch_val = fill(NaN, length(P_B_ch))
        P_B_dis_val = fill(NaN, length(P_B_dis))
        P_G_plus_val = fill(NaN, length(P_G_plus))
        P_G_minus_val = fill(NaN, length(P_G_minus))
        SOC_val = fill(NaN, length(SOC))
        total_cost_vals = fill(NaN, S)
        operating_cost_vals = fill(NaN, S)
        z_val = NaN
        objective_val = NaN
    end


    return Dict(
        :CPV_kWp => CPV_val,
        :E_B_kWh => E_B_val,
        :CW_kWp  => CW_val,
        :total_costs => total_cost_vals,
        :operating_costs => operating_cost_vals,
        :z => z_val,
        :objective => objective_val,
        :P_PV => P_PV_val,
        :P_W_expr => [(wind[s,t]*CW_val) for s in 1:S, t in 1:T],
        :P_G_plus => P_G_plus_val,
        :P_G_minus => P_G_minus_val,
        :P_B_ch => P_B_ch_val,
        :P_B_dis => P_B_dis_val,
        :SOC => SOC_val,
        :termination_status => status
    )
end

# Example call:
simulate_microgrid_with_wind(1)
simulate_microgrid_with_wind(5)
simulate_microgrid_with_wind(8)


# THE END





#--------------------------------------------------------------------------------------------------------

# GETTING SENSITIVITIES FROM SAVED LP FILES IN GUROBI INTERACTIVE SHELL

#--------------------------------------------------------------------------------------------------------

# # SENSITIVITY pi_B FOR 1 YEAR
# lp_path = r"microgrid_pi_B_1.lp" # paste filepath here
# m = read(lp_path)
# print("LP file loaded:", lp_path)
# m.optimize()
# print("Optimization finished. Status:", m.status)
# var_E_B = m.getVarByName("_E_B")
# low_EB  = var_E_B.SAObjLow
# up_EB   = var_E_B.SAObjUp
# z_opt = m.objVal
# S = var_E_B.x
# print("\nAbsolute sensitivity interval for delta pi_B (kgCO2/kWh battery capacity):") 
# print("Lower bound:", low_EB) 
# print("Upper bound:", up_EB)
# print("\nOptimal objective z* =", z_opt)
# print("Optimal E_B* =", S)
# print("Derivative dz/d(pi_B) =", S)
# print("\nSensitivity formula:")
# print(f"z(Δpi_B) = {z_opt} + ({S}) * Δpi_B")



# # SENSITIVITY pi_B FOR 1 YEAR
# lp_path = r"microgrid_pi_B_5.lp" # paste filepath here
# m = read(lp_path)
# print("LP file loaded:", lp_path)
# m.optimize()
# print("Optimization finished. Status:", m.status)
# var_E_B = m.getVarByName("_E_B")
# low_EB  = var_E_B.SAObjLow
# up_EB   = var_E_B.SAObjUp
# z_opt = m.objVal
# S = var_E_B.x
# print("\nAbsolute sensitivity interval for delta pi_B (kgCO2/kWh battery capacity):") 
# print("Lower bound:", low_EB) 
# print("Upper bound:", up_EB)
# print("\nOptimal objective z* =", z_opt)
# print("Optimal E_B* =", S)
# print("Derivative dz/d(pi_B) =", S)
# print("\nSensitivity formula:")
# print(f"z(Δpi_B) = {z_opt} + ({S}) * Δpi_B")