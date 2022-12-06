### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ e5a18d4c-14cd-11ed-36d5-69de0fd02830
# ╠═╡ show_logs = false
begin
	local_dir = joinpath(splitpath(@__FILE__)[1:end-1])
	import Pkg
	Pkg.activate(mktempdir())
	Pkg.add([
		Pkg.PackageSpec(name="CairoMakie"),
		Pkg.PackageSpec(name="MakiePublication"),
		Pkg.PackageSpec(name="ForwardDiff"),
		Pkg.PackageSpec(name="AbstractDifferentiation"),
		Pkg.PackageSpec(name="CSV"),
		Pkg.PackageSpec(name="ComponentArrays"),
		Pkg.PackageSpec(name="DataFrames"),
		Pkg.PackageSpec(name="Optimization"),
		Pkg.PackageSpec(name="OptimizationOptimJL"),
		Pkg.PackageSpec(name="LabelledArrays"),
		Pkg.PackageSpec(name="HypertextLiteral"),
		Pkg.PackageSpec(name="PlutoUI"),
	])
	Pkg.develop([
		Pkg.PackageSpec(path=joinpath(local_dir, "InverseLangevinApproximations")),
		Pkg.PackageSpec(path=joinpath(local_dir, "NonlinearContinua")),
		Pkg.PackageSpec(path=joinpath(local_dir, "Hyperelastics")),
		# Pkg.PackageSpec(path=joinpath(local_dir, "PlutoUI")),
	])
	using PlutoUI, AbstractDifferentiation, ForwardDiff, CSV, ComponentArrays, DataFrames, Optimization, OptimizationOptimJL, InverseLangevinApproximations, LabelledArrays, CairoMakie, MakiePublication
end

# ╔═╡ 2d189645-189f-4886-a6d5-5718a613798f
using Hyperelastics

# ╔═╡ 73ab5774-dc3c-4759-92c4-7f7917c18cbf
HTML("""<center><h1>Vagus <br> Hyperelastic Model Fitting Toolbox</h1></center>
		<center><h2>Upload Uniaxial Test Data</h2></center>
		""")

# ╔═╡ 692b1d0d-2353-4931-b289-490f74988811
md"""
Test Type: $(@bind test_type Select([:Uniaxial, :Biaxial]))
"""

# ╔═╡ 6434a755-90be-45f8-8e1e-cdeba4be244b
@bind data FilePicker()

# ╔═╡ f12538a9-f595-4fae-b76c-078179bc5109
if !isnothing(data) HTML("""<center><h3>Verification Plot</h3></center>""") end

# ╔═╡ d0319d95-f335-48fa-b789-59daf9a0f1a4
if !isnothing(data) HTML("""<center><h2>Select Hyperelastic Model</h2></center>""") end

# ╔═╡ 9343a51e-5002-4489-a55f-12c49f5b8cf3
if !isnothing(data)
md"""
!!! note "Note"
	- When selecting a phenomenological model, be aware that using higher order models may result in overfitting of the data.
	- All moduli in models are in the defined stress units above
"""
end

# ╔═╡ da3634ea-48d7-4d4f-a853-c631a6fa7bf4
if !isnothing(data) html"""<center><h3> Model Information</h3></center>""" end

# ╔═╡ c6e726ab-ea78-4129-a662-338976633cd5
if !isnothing(data) html"""<center><h2> Set initial parameter guess</h2></center>""" end

# ╔═╡ 08d775f2-94fc-4ca8-bcdd-e9535cfd129a
if !isnothing(data)
md"""
Optimizer: $(@bind optimizer Select([:LBFGS, :BFGS, :NelderMead])) - *If parameters are not converging, try using a different optimizer or changing your initial guess*
"""
end

# ╔═╡ 7196aa51-e86d-4f0e-ae40-cc6aa74aa237
md"---"

# ╔═╡ 0dd8b7de-570d-41a7-b83d-d1bbe39c017e
TableOfContents()

# ╔═╡ d495c5e5-bf33-475c-a49a-5c9f8dc13789
set_theme!(MakiePublication.theme_web(width = 1000))

# ╔═╡ 6f061996-be32-493d-80e2-daedec8bb103
exclude = [:HorganMurphy, :KhiemItskov, :GeneralCompressible, :LogarithmicCompressible, :GeneralMooneyRivlin];

# ╔═╡ e0e7407d-fe60-4583-8060-3ba38c22c409
begin
	hyperelastic_models = filter(x -> typeof(getfield(Hyperelastics, x)) <: Union{DataType, UnionAll},names(Hyperelastics))
	hyperelastic_models = filter(x -> !(getfield(Hyperelastics, x) <: Hyperelastics.AbstractDataDrivenHyperelasticModel) && (getfield(Hyperelastics, x) <: Hyperelastics.AbstractHyperelasticModel), hyperelastic_models)
	hyperelastic_models = filter(x -> !(x in exclude), hyperelastic_models)
	map(model->parameters(model()), Base.Fix1(getfield, Hyperelastics).( hyperelastic_models))
end;

# ╔═╡ 2f1fde4b-6bd8-42b4-bf5c-d61006d55f10
if !isnothing(data) @bind model Select(hyperelastic_models) end

# ╔═╡ a75d209e-93cb-4b21-899e-4c567f0dfb09
if !isnothing(data) eval(:(@doc $(getfield(Hyperelastics, model)()))) end

# ╔═╡ 7998136a-de3d-42f9-9028-1172415c8b75
if !isnothing(data)
	df = CSV.read(data["data"], DataFrame);
end;

# ╔═╡ 69068002-ca3a-4e19-9562-6736d3b15dea
if !isnothing(data)
	if test_type == :Uniaxial
md"""
Stress Column: $(@bind stress_column Select(names(df)))
Stretch Column: $(@bind stretch_column Select(names(df)))
Stress Units: $(@bind stress_units TextField())
Test Name: $(@bind test_name TextField())
"""
	elseif test_type == :Biaxial
md"""
Stress-1 Column: $(@bind stress1_column Select(names(df)))
Stress-2 Column: $(@bind stress2_column Select(names(df)))
Stretch-1 Column: $(@bind stretch1_column Select(names(df)))
Stretch-2 Column: $(@bind stretch2_column Select(names(df)))
Stress Units: $(@bind stress_units TextField())
Test Name: $(@bind test_name TextField())
"""
	end
end

# ╔═╡ 2607b1b6-9c9c-482f-b38b-35e83a57f5d3
if !isnothing(data)
	if test_type == :Uniaxial
	f, ax, p = scatter(
		df[!, stretch_column],
		df[!, stress_column],
		axis = (xlabel = "Stretch", ylabel = "Stress [$stress_units]"),
		label = test_name*" - Experimental"
	)
	elseif test_type == :Biaxial
	f, ax, p = scatter(
		df[!, stretch1_column],
		df[!, stress1_column],
		axis = (xlabel = "Stretch", ylabel = "Stress [$stress_units]"),
		label = test_name*" - Experimental - 1"
	)
	scatter!(
		df[!, stretch2_column],
		df[!, stress2_column],
		label = test_name * " - Experimental - 2"
	)
	end
	axislegend(position = :lt)
	f
end

# ╔═╡ 12256359-1dca-4a71-a225-66994e2dfd66
if !isnothing(data)
	if test_type == :Uniaxial
		he_data = HyperelasticUniaxialTest(df[!, stretch_column],df[!, stress_column],  name = test_name);
	elseif test_type == :Biaxial
		he_data = HyperelasticBiaxialTest(
			df[!, stretch1_column],
			df[!, stretch2_column],
			df[!, stress1_column],
			df[!, stress2_column],
			name = test_name);
	end
	map(model->Hyperelastics.parameter_bounds(model(), he_data), Base.Fix1(getfield, Hyperelastics).( hyperelastic_models))
end;

# ╔═╡ 4d6f03c0-203a-4536-8ca2-c3dd77182ce6
function set_parameters(model,data)
	ψ = getfield(Hyperelastics, model)()
	ps = Hyperelastics.parameters(ψ)
	bounds = Hyperelastics.parameter_bounds(ψ, data)
	if isnothing(bounds.lb)
		lb = Dict(ps.=>-Inf)
	else
		lb = Dict{Symbol, Float64}(pairs(bounds.lb))
	end
	if isnothing(bounds.ub)
		ub = Dict(ps.=>Inf)
	else
		ub = Dict(pairs(bounds.ub))
	end

	return PlutoUI.combine() do Child
		inputs = [
			md"""$(string(p)) [ $(round(lb[p], digits = 3)) to $(round(ub[p], digits = 3)) ] $(Child(string(p), TextField()))"""
			for p in ps
		]
		md"""
		$(inputs)
		Fit Model: $(@bind fit_model CheckBox(default = false))
		"""
	end
end;

# ╔═╡ 703091d0-bf33-4baf-b75e-43e01b42ec0b
if !isnothing(data)
	@bind ps set_parameters(model, he_data)
end

# ╔═╡ d0713eb0-fe75-4ea4-bf20-2d4e9b722da5
if !isnothing(data)
			parsed, p₀ = try
				pair_ps = map(x->x.first => parse.(Float64, replace.(replace.(split(x.second, ","), " "=>""), ","=>"")), collect(pairs(ps)))
				ps = []
				for i in eachindex(pair_ps)
					if length(pair_ps[i].second) == 1
						push!(ps, pair_ps[i].first => pair_ps[i].second[1])
					else
						push!(ps, pair_ps[i])
					end
				end
				true, ComponentVector(NamedTuple(ps))
			catch
				false, nothing
		end
end;

# ╔═╡ 1018d35f-42e9-4970-8a5f-f5cc6e951cbc
begin
	if !isnothing(data)
		if parsed && fit_model
			ψ = getfield(Hyperelastics, Symbol(model))()
			heprob = HyperelasticProblem(ψ, he_data, p₀)
			opt = getfield(OptimizationOptimJL, optimizer)()
			solution = solve(heprob, opt)
			sol = NamedTuple(solution.u)
		end
	end
end;

# ╔═╡ 0fa152b1-462a-4f34-9753-13ef6ef63071
begin
	if !isnothing(data)
		if !(isnothing(ps))
			if @isdefined sol
			str_table = let
				_str = "<span>"
				columns = string.(parameters(getfield(Hyperelastics, model)()))
				####################             HTML        #######################
				_str *= """<center><h2> Final Parameters</h2></center>"""
				_str *= """<table>"""
				for column ∈ columns
					_str *=	"""<th>$(replace(column, "_" => " "))</th>"""
				end
				_str *= """<tr>"""
				for column ∈ columns
					_str *= """<td>$(round.(getfield(sol, Symbol(column)), sigdigits = 6))</input></td>"""
				end
				_str *= """</tr>"""
				_str *= """<tfoot><tr><td></td>"""
				_str *= """<td></td></tr></tfoot>"""
				_str *= """</table>"""
				_str *= """</span>"""
			end
			HTML(str_table)
			end
		end
	end
end

# ╔═╡ 1345476c-ee08-4233-8506-0ebc94a2bec5
let
	if !isnothing(data)
		if parsed && fit_model
			if @isdefined sol
		ψ = getfield(Hyperelastics, Symbol(model))()
		ŷ = predict(ψ, he_data, sol)
		if test_type == :Uniaxial
			Δs₁₃ = getindex.(ŷ.data.s, 1)
			λ₁ = getindex.(ŷ.data.λ, 1)
			s₁ = getindex.(he_data.data.s, 1)
			f = Figure()
			ax = CairoMakie.Axis(f,xlabel = "Stretch", xticks = 1:maximum(df[!, stretch_column]), ylabel = "Stress [$stress_units]")
			s1 = scatter!(
				ax,
				λ₁,s₁,
			)
			l1 = lines!(
				ax,
				λ₁,
				Δs₁₃,
				color = MakiePublication.seaborn_muted()[2],
			)
			axislegend(ax, [[s1], [l1]], [test_name*" - Experimental", split(string(typeof(ψ)), ".")[2]], position = :lt)
			f[1,1] = ax
			f
		elseif test_type == :Biaxial
			Δs₁₃ = getindex.(ŷ.data.s, 1)
			Δs₂₃ = getindex.(ŷ.data.s, 2)
			λ₁ = getindex.(ŷ.data.λ, 1)
			λ₂ = getindex.(ŷ.data.λ, 2)
			s₁ = getindex.(he_data.data.s, 1)
			s₂ = getindex.(he_data.data.s, 2)
			fig = Figure()
			ax = CairoMakie.Axis(fig,xlabel = "Stretch", ylabel = "Stress [$stress_units]")
			s1 = scatter!(
				ax,
				λ₁,s₁,
			)
			s2 = scatter!(
				ax,
				λ₂,s₂,
			)
			l1 = lines!(
				ax,
				λ₁,
				Δs₁₃,
				# color = MakiePublication.seaborn_muted()[2],
			)
			l2 = lines!(
				ax,
				λ₂,
				Δs₂₃,
				# color = MakiePublication.seaborn_muted()[2],
			)
			axislegend(ax, [[s1], [s2], [l1], [l2]], [
				test_name*" - Experimental - 1",
				test_name*" - Experimental - 2",
				split(string(typeof(ψ)), ".")[2]*" - 1",
				split(string(typeof(ψ)), ".")[2]*" - 2",
			], position = :lt)
			fig[1,1] = ax
			fig
		end

		end
		end
	end
end

# ╔═╡ 86f7e74f-c0a9-4561-85b9-f3ed6559facc
function ShearModulus(ψ, ps; adb = AD.ForwardDiffBackend())
	s(x) = AD.gradient(adb,x->StrainEnergyDensity(ψ, x, ps), x)[1][1]-AD.gradient(adb,x->StrainEnergyDensity(ψ, x, ps), x)[1][3]*x[3]
	AD.gradient(adb,y->s(y)[1], [1.0, 1.0, 1.0])[1][1]
end;

# ╔═╡ 8ea07dab-06dc-456d-9769-5e9c3980a777
ElasticModulus(ψ, ps) = ShearModulus(ψ, ps)*3;

# ╔═╡ 9441279c-49d9-4640-aca5-4576e6ee29ed
if !isnothing(data) && fit_model && @isdefined sol
	if parsed && !isnothing(data)
		HTML("""
		<center><h2> Other Values </h2></center>
		Small Strain Shear Modulus: $(round(ShearModulus(ψ, sol), digits = 3)) $(stress_units)
		<br>
		Small Strain Elastic Modulus: $(round(ElasticModulus(ψ, sol), digits = 3)) $(stress_units)
		""")
	end
end

# ╔═╡ bcf0c08c-cc7a-4785-a87b-2be47633eb85
function model_note(ψ::Gent)
	return (
	μ = "Small strain shear modulus",
	Jₘ = "Limiting Stretch Invariant"
	)
end;

# ╔═╡ Cell order:
# ╟─73ab5774-dc3c-4759-92c4-7f7917c18cbf
# ╟─692b1d0d-2353-4931-b289-490f74988811
# ╟─6434a755-90be-45f8-8e1e-cdeba4be244b
# ╟─69068002-ca3a-4e19-9562-6736d3b15dea
# ╟─f12538a9-f595-4fae-b76c-078179bc5109
# ╟─2607b1b6-9c9c-482f-b38b-35e83a57f5d3
# ╟─d0319d95-f335-48fa-b789-59daf9a0f1a4
# ╟─9343a51e-5002-4489-a55f-12c49f5b8cf3
# ╟─2f1fde4b-6bd8-42b4-bf5c-d61006d55f10
# ╟─da3634ea-48d7-4d4f-a853-c631a6fa7bf4
# ╟─a75d209e-93cb-4b21-899e-4c567f0dfb09
# ╟─c6e726ab-ea78-4129-a662-338976633cd5
# ╟─703091d0-bf33-4baf-b75e-43e01b42ec0b
# ╟─08d775f2-94fc-4ca8-bcdd-e9535cfd129a
# ╟─1018d35f-42e9-4970-8a5f-f5cc6e951cbc
# ╟─0fa152b1-462a-4f34-9753-13ef6ef63071
# ╟─1345476c-ee08-4233-8506-0ebc94a2bec5
# ╟─9441279c-49d9-4640-aca5-4576e6ee29ed
# ╟─7196aa51-e86d-4f0e-ae40-cc6aa74aa237
# ╟─e5a18d4c-14cd-11ed-36d5-69de0fd02830
# ╟─0dd8b7de-570d-41a7-b83d-d1bbe39c017e
# ╟─2d189645-189f-4886-a6d5-5718a613798f
# ╟─d495c5e5-bf33-475c-a49a-5c9f8dc13789
# ╟─6f061996-be32-493d-80e2-daedec8bb103
# ╟─e0e7407d-fe60-4583-8060-3ba38c22c409
# ╟─7998136a-de3d-42f9-9028-1172415c8b75
# ╟─12256359-1dca-4a71-a225-66994e2dfd66
# ╟─4d6f03c0-203a-4536-8ca2-c3dd77182ce6
# ╟─d0713eb0-fe75-4ea4-bf20-2d4e9b722da5
# ╟─86f7e74f-c0a9-4561-85b9-f3ed6559facc
# ╟─8ea07dab-06dc-456d-9769-5e9c3980a777
# ╟─bcf0c08c-cc7a-4785-a87b-2be47633eb85
