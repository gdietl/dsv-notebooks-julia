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

# ╔═╡ acb44b30-3a01-11eb-00ad-c72d740adc29
begin
	using Pkg
	Pkg.add("DSP")
	Pkg.add("LaTeXStrings")
	Pkg.add("PlutoUI")

	using Plots, DSP, LaTeXStrings, PlutoUI
end

# ╔═╡ 79cc6450-3a01-11eb-085d-2dee3a616165
md"# Julia-Notebook zur Veranschaulichung der Optimalmethode

Guido Dietl, 9. Dezember 2020

## Initialisierung"

# ╔═╡ ba70a08c-3a01-11eb-362f-9388281db67b
md"""## Parks-McClellan-Methode

Die Parks-McClellan-Methode findet das optimale FIR-Filter unter Verwendung einer Tschebyscheff-Approximation. Es entspricht dem Remez-Verfahren, das auch für IIR-Filter angewandt werden kann.

Einstellung der Parameter:
* Filterordnung $N$: $(@bind N Slider(10:100))
* normierte Grenzfrequenz $f_{1}/f_s$: $(@bind f1_norm Slider(.05:.01:.2))
* normierter Übergangsbereich $\Delta f_1/f_s$: $(@bind Δf1_norm Slider(0.01:.01:.1))
* normierte Grenzfrequenz $f_{2}/f_s$: $(@bind f2_norm Slider(0.25:.01:.45))
* normierter Übergangsbereich $\Delta f_{2}/f_s$: $(@bind Δf2_norm Slider(0.01:.01:.1))
* Filtertyp: $(@bind filter_type Select(["TP", "HP", "BP", "BS"]))


Des weiteren wird angenommen, dass $H^\prime(f)=1$ im Durchlassbereich und $H^\prime(f)=0$ im Sperrbereich gilt."""

# ╔═╡ 30a4a298-3a05-11eb-24ba-73efa31073a2
md"Werte der Paramter:
* Filterordnung $N=$ $(N)
* normierte Grenzfrequenz $f_{1}/f_s=$ $(f1_norm)
* normierter Übergangsbereich $\Delta f_1/f_s=$ $(Δf1_norm)
* normierte Grenzfrequenz $f_{2}/f_s=$ $(f2_norm)
* normierter Übergangsbereich $\Delta f_2/f_s=$ $(Δf2_norm)

Frequenzgang:"

# ╔═╡ 91ef5e0a-3a08-11eb-2fdc-2b3fd3af83ec
begin
	n = 0:N
    if filter_type == "TP"
        filter = remez(N, [(0, (f1_norm-Δf1_norm/2)) => 1, ((f1_norm+Δf1_norm/2), .5) => 0])
    elseif filter_type == "HP"
        filter = remez(N, [(0, (f1_norm-Δf1_norm/2)) => 0, ((f1_norm+Δf1_norm/2), .5) => 1])
    elseif filter_type == "BP"
        filter = remez(N, [(0, (f1_norm-Δf1_norm/2)) => 0, ((f1_norm+Δf1_norm/2), (f2_norm-Δf2_norm/2)) => 1, ((f2_norm+Δf2_norm/2), .5) => 0])
    elseif filter_type == "BS"
        filter = remez(N, [(0, (f1_norm-Δf1_norm/2)) => 1, ((f1_norm+Δf1_norm/2), (f2_norm-Δf2_norm/2)) => 0, ((f2_norm+Δf2_norm/2), .5) => 1])
    end
    P = PolynomialRatio(filter, [1])
    f_norm = 0:.001:.5
    H = freqz(P, f_norm*2*π)
    plot(f_norm, 20*log10.(abs.(H)), label=false, lw=2, linecolor = :blue, legend = :bottomright)
    vline!([f1_norm-Δf1_norm/2], label = L"f = f_1-\Delta f_1/2", linecolor = :green)
    vline!([f1_norm+Δf1_norm/2], label = L"f = f_1+\Delta f_1/2", linecolor = :lightgreen)
    if filter_type == "BP" || filter_type == "BS"
        vline!([f2_norm-Δf2_norm/2], label = L"f = f_2-\Delta f_2/2", linecolor = :red)
        vline!([f2_norm+Δf2_norm/2], label = L"f = f_2+\Delta f_2/2", linecolor = :orange)
    end
    xlims!((0,.5))
    ylims!((-150,5))
    xlabel!(L"f/f_s")
    ylabel!(L"20 \log_{10}H^\prime(f)")
end

# ╔═╡ Cell order:
# ╟─79cc6450-3a01-11eb-085d-2dee3a616165
# ╠═acb44b30-3a01-11eb-00ad-c72d740adc29
# ╟─ba70a08c-3a01-11eb-362f-9388281db67b
# ╟─30a4a298-3a05-11eb-24ba-73efa31073a2
# ╠═91ef5e0a-3a08-11eb-2fdc-2b3fd3af83ec
