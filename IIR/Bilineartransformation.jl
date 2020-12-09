### A Pluto.jl notebook ###
# v0.12.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 35f938a4-3a51-11eb-0de7-294bd4dc1657
begin
	using Pkg
	Pkg.add("DSP")
	Pkg.add("LaTeXStrings")
	Pkg.add("PlutoUI")
	
	using Plots, DSP, LaTeXStrings, PlutoUI
end

# ╔═╡ d3d460ca-3a50-11eb-12b9-7dbe361e0d8c
md"# Julia-Notebook zur Veranschaulichung des IIR-Filterentwurfs mit Bilineartransformation

Guido Dietl, 9. Dezember 2020

## Initialisierung"

# ╔═╡ bcda031a-3a51-11eb-1526-e395952c371b
md"""
## IIR-Filterentwurf mit Bilineartransformation

Designmethoden:
* **Butterworth**-Filter: `Butterworth`
* **Tschebyscheff** mit Welligkeit im Durchlassbereich: `Chebyshev1`
* **Tschebyscheff** mit Welligkeit im Sperrbereich: `Chebyshev2`
* **Elliptisches** oder Cauer-Filter: `Elliptic`

Einstellung der Parameter:
* Filterordnung $N$ des analogen Musterfilters: $(@bind N Slider(2:20))
* normierte Grenzfrequenz $f_{1}/f_s$: $(@bind f1_norm Slider(.05:.01:.2))
* normierte Grenzfrequenz $f_{2}/f_s$:$(@bind f2_norm Slider(0.25:.01:.45))
* Welligkeit $\Delta H_\text{pass}$ im Durchlassbereich: $(@bind ΔHpass Slider(0.01:.01:0.2))
* Welligkeit $\Delta H_\text{stop}$ im Sperrbereich: $(@bind ΔHstop Slider(0.01:.01:0.2))
* Designmethode: $(@bind design_type Select(["Butterworth", "Chebyshev1", "Chebyshev2", "Elliptic"]))
* Filtertyp: $(@bind filter_type Select(["TP", "HP", "BP", "BS"]))
"""

# ╔═╡ 46a789d8-3a5b-11eb-2206-edbb9eb63eca
md"Werte der Parameter:
* Filterordnung $N=$ $(N)
* normierte Grenzfrequenz $f_1/f_s=$ $(f1_norm)
* normierte Grenzfrequenz $f_2/f_s=$ $(f2_norm)
* Welligkeit im Durchlassbereich $\Delta H_\mathrm{pass}=$ $(ΔHpass)
* Welligkeit im Sperrbereich $\Delta H_\mathrm{stop}=$ $(ΔHstop)

Amplitudengang:"

# ╔═╡ a16e9ce4-3a5b-11eb-1c47-03dee9040484
begin
	ΔHpassdB = -20*log10(1-ΔHpass)
    ΔHstopdB = -20*log10(ΔHstop)
    if filter_type == "TP"
        responsetype = Lowpass(f1_norm*2)
    elseif filter_type == "HP"
        responsetype = Highpass(f1_norm*2)
    elseif filter_type == "BP"
        responsetype = Bandpass(f1_norm*2, f2_norm*2)
    elseif filter_type == "BS"
        responsetype = Bandstop(f1_norm*2, f2_norm*2)
    end
    if design_type == "Butterworth"
        designmethod = Butterworth(N)
    elseif design_type == "Chebyshev1"
        designmethod = Chebyshev1(N, ΔHpassdB)
    elseif design_type == "Chebyshev2"
        designmethod = Chebyshev2(N, ΔHstopdB)
    elseif design_type == "Elliptic"
        designmethod = Elliptic(N, ΔHpassdB, ΔHstopdB)
    end
    filter = digitalfilter(responsetype, designmethod)
    f_norm = 0:.001:.5
    H = freqz(filter, f_norm*2*π)
    plot(f_norm, abs.(H), label=false, lw=2, c = :blue)
    vline!([f1_norm], label = L"f = f_1", c = :green)
    if filter_type == "BP" || filter_type == "BS"
        vline!([f2_norm], label = L"f = f_2", c = :green)
    end
    xlims!((0,.5))
    ylims!((0,1.2))
    xlabel!(L"f/f_s")
    ylabel!(L"H^\prime(f)")
end

# ╔═╡ 2516391c-3a5c-11eb-3fc6-29f504263840
md"Filterkoeffizienten:"

# ╔═╡ 1b703d5e-3a5c-11eb-3578-0712c4939fc2
filtercoefficients = [coefb(filter)'; coefa(filter)']'

# ╔═╡ Cell order:
# ╟─d3d460ca-3a50-11eb-12b9-7dbe361e0d8c
# ╠═35f938a4-3a51-11eb-0de7-294bd4dc1657
# ╟─bcda031a-3a51-11eb-1526-e395952c371b
# ╟─46a789d8-3a5b-11eb-2206-edbb9eb63eca
# ╠═a16e9ce4-3a5b-11eb-1c47-03dee9040484
# ╟─2516391c-3a5c-11eb-3fc6-29f504263840
# ╠═1b703d5e-3a5c-11eb-3578-0712c4939fc2
