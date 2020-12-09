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
md"# Julia-Notebook zur Veranschaulichung der Pol-Nullstellen-Methode

Guido Dietl, 9. Dezember 2020

## Initialisierung"

# ╔═╡ 5fca3ed8-3a51-11eb-317c-3defe6947409
md"Definition nützliche Funktionen:"

# ╔═╡ 7ac2f1f8-3a51-11eb-1738-cb7b9ed4c0dc
function circle(r::Float64)
   	θ = LinRange(0, 2*π, 500)
   	sin.(θ), cos.(θ)
end

# ╔═╡ 09486782-3a52-11eb-2a43-c76c672bdb1c
function zplane(z::Array, p::Array)
   	plot(circle(1.0), seriestype = [:shape,], c = :black, fillalpha = .2, aspect_ratio = 1, label=false, legend=:bottomright)
   	scatter!(real(z), imag(z), c = :blue, label="Nullstellen")
   	xlabel!(L"\mathrm{Re}{z}")
   	ylabel!(L"\mathrm{Im}{z}")
   	scatter!(real(p), imag(p), c = :red, label="Polstellen")
end

# ╔═╡ bcda031a-3a51-11eb-1526-e395952c371b
md"## Bandpass

Parameter:
* Abtastfrequenz $f_s$
* Bandbreite $B$
* Mittenfrequenz $f_0$
* ein gegebener Wert von $H^\prime(f)=H\left(e^{j2\pi f/ f_s}\right)$, z.B., $H^\prime(f_0)=1$

Nullstellen: $z_{0,1}=1$ und $z_{0,2}=-1$

Polstellen: $z_{\infty,1/2}=r e^{\pm j 2 \pi \frac{f_0}{f_s}}$ mit $r=1-\frac{B}{f_s}\pi$

Funktion zur Ermittlung der Filterkoeffzienten (Annahme: $H^\prime(f_0)=1$):"

# ╔═╡ eaf8f466-3a51-11eb-3f00-0b3b7d530063
function pole_zero_method_BP(f0_norm::Float64, B_norm::Float64)
    # Nullstellen
    z = [1, -1]

    # Pole
    r = 1 - B_norm*π
    z0 = exp(im*2*π*f0_norm)
    p = [r*z0, r*conj(z0)]

    # Berechnung des Koeffizienten b_0 (gain k)
    k = (z0-p[1])*(z0-p[2])/(z0-z[1])/(z0-z[2])

    # Filterkoeffizienten
    filter = ZeroPoleGain(z, p, k)
    
    return z, p, filter
    
end

# ╔═╡ 3a25bf44-3a52-11eb-2e33-0bf8ebad32a1
md"Einstellung der Parameter:

* normierte Mittenfrequenz $f_0/f_s$: $(@bind f0_norm Slider(0:.01:.5))
* normierte Bandbreite $B/f_s$: $(@bind B_norm Slider(0.01:.01:.5))"

# ╔═╡ f04b8a2e-3a52-11eb-136b-456bbe94f7b6
md"Werte der Parameter:
* normierte Mittenfrequenz $f_0/f_s=$ $(f0_norm)
* normierte Bandbreite $B/f_s=$ $(B_norm)

Amplitudengang und Pol-Nullstellen-Diagramm:"

# ╔═╡ eae71bcc-3a53-11eb-044c-9fb89d7e346a
md"## Bandsperre

Parameter:
* Abtastfrequenz $f_s$
* Bandbreite $B$
* Mittenfrequenz $f_0$
* ein gegebener Wert von $H^\prime(f)=H\left(e^{j2\pi f/ f_s}\right)$, z.B., $H^\prime(0)=1$

Nullstellen: $z_{0,1/2}=e^{\pm j 2 \pi \frac{f_0}{f_s}}$

Polstellen: $z_{\infty,1/2}=r e^{\pm j 2 \pi \frac{f_0}{f_s}}$ mit $r=1-\frac{B}{f_s}\pi$

Funktion zur Ermittlung der Filterkoeffzienten (Annahme: $H^\prime(0)=1$):"

# ╔═╡ 11023d02-3a54-11eb-0388-851a2e950a2c
function pole_zero_method_BS(f0_norm::Float64, B_norm::Float64)
    # Nullstellen
    z0 = exp(im*2*π*f0_norm)
    z = [z0, conj(z0)]

    # Pole
    r = 1 - B_norm*π
    p = [r*z0, r*conj(z0)]

    # Berechnung des Koeffizienten b_0 (gain k)
    k = (1-p[1])*(1-p[2])/(1-z[1])/(1-z[2])

    # Filterkoeffizienten
    filter = ZeroPoleGain(z, p, k)
    return z, p, filter
    
end

# ╔═╡ 0d73e6c6-3a55-11eb-1620-a3d4cc12c9fb
function plot_spectrum_zplane(f0_norm::Float64, B_norm::Float64, filter_type)
	f_norm = 0:.001:.5
	if filter_type == "BP"
    	z, p, filter = pole_zero_method_BP(f0_norm, B_norm)
	elseif filter_type =="BS"
		z, p, filter = pole_zero_method_BS(f0_norm, B_norm)
	end
    H = freqz(filter, f_norm*2*π)
    a = plot(f_norm, 20*log10.(abs.(H)), label=false, lw=2, c = :blue, legend=:bottomright)
    vline!([f0_norm-B_norm/2], label = L"f = f_0-B/2", c = :green)
    vline!([f0_norm+B_norm/2], label = L"f = f_0+B/2", c = :green)
    xlims!((0,.5))
    ylims!((-50,5))
    xlabel!(L"f/f_s")
    ylabel!(L"20 \log_{10}H^\prime(f)")
    b = zplane(z, p)
    plot(a, b)
end

# ╔═╡ 003785d2-3a53-11eb-0198-55120fd7b69d
plot_spectrum_zplane(f0_norm, B_norm, "BP")

# ╔═╡ 2596c882-3a54-11eb-3fcd-5bdfdefaf573
md"Werte der Parameter:
* normierte Mittenfrequenz $f_0/f_s=$ $(f0_norm)
* normierte Bandbreite $B/f_s=$ $(B_norm)

Amplitudengang und Pol-Nullstellen-Diagramm:"

# ╔═╡ 48559c68-3a54-11eb-2785-e582a5512075
plot_spectrum_zplane(f0_norm, B_norm, "BS")

# ╔═╡ Cell order:
# ╟─d3d460ca-3a50-11eb-12b9-7dbe361e0d8c
# ╠═35f938a4-3a51-11eb-0de7-294bd4dc1657
# ╟─5fca3ed8-3a51-11eb-317c-3defe6947409
# ╠═7ac2f1f8-3a51-11eb-1738-cb7b9ed4c0dc
# ╠═09486782-3a52-11eb-2a43-c76c672bdb1c
# ╠═0d73e6c6-3a55-11eb-1620-a3d4cc12c9fb
# ╟─bcda031a-3a51-11eb-1526-e395952c371b
# ╠═eaf8f466-3a51-11eb-3f00-0b3b7d530063
# ╟─3a25bf44-3a52-11eb-2e33-0bf8ebad32a1
# ╟─f04b8a2e-3a52-11eb-136b-456bbe94f7b6
# ╠═003785d2-3a53-11eb-0198-55120fd7b69d
# ╟─eae71bcc-3a53-11eb-044c-9fb89d7e346a
# ╠═11023d02-3a54-11eb-0388-851a2e950a2c
# ╟─2596c882-3a54-11eb-3fcd-5bdfdefaf573
# ╠═48559c68-3a54-11eb-2785-e582a5512075
