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

# ╔═╡ 93a36792-3ac5-11eb-2653-85d2ba3306f2
begin
	using Pkg
	Pkg.add("DSP")
	Pkg.add("LaTeXStrings")
	
	using Plots, DSP, LaTeXStrings, PlutoUI
end

# ╔═╡ 6f854664-3ac5-11eb-2675-9301f51b27b9
md"# Julia-Notebook zur Veranschaulichung der Fenstermethode

Guido Dietl, 10. Dezember 2020

## Initialisierung"

# ╔═╡ d0193936-3ac5-11eb-1aa2-892211210245
md"## Abgetastete, gefensterte und verschobene Impulsantwort eines idealen analogen Tiefpassfilters

Impulsantwort (mit $\mathrm{sinc}(x)=\sin(\pi x)/(\pi x)$):

$$h[n]=\frac{\Omega_g}{\pi}\mathrm{sinc}\frac{\Omega_g\left(n-\frac{N}{2}\right)}{\pi}\cdot w[n]$$

mit normierter Grenzkreisfrequenz $\Omega_g=2 \pi f_g / f_s$ (Abtastfrequenz $f_s$), Filterordnung $N$ und Fensterfunktion $w[n]$.

Fensterfunktionen:
* **Rechteck**:
  
  $$w[n]=\begin{cases}1,& 0\leq n \leq N\\0,& \text{otherwise}\end{cases}$$

* **Hann** (auch fälschlicherweise als **Hanning** bezeichnet, mit *hanning* wurde, abgeleitet von *to hann*, im Ursprungspaper die Anwendung des Hann-Fensters bezeichnet): Raised-Cosine-Fenster mit Nullen an den Endpunkten

  $$w[n]=\begin{cases}\cos^2\left(\pi\left(\frac{n}{N}-0.5\right)\right),& 0\leq n \leq N\\0,& \text{otherwise}\end{cases}$$

* **Hamming**: keine Nullen an den Eckpunkten und kleinere Flankensteilheit gegenüber Hann-Fenster, aber entfernt erste Nebenkeule

  $$w[n]=\begin{cases}0.54+0.46\cos\left(2\pi\left(\frac{n}{N}-0.5\right)\right),& 0\leq n \leq N\\0,& \text{otherwise}\end{cases}$$

* **Tukey**: flacher Verlauf um die Mitte, gleich Null an den Endpunkten und sinusoidaler Übergang mit Parameter $a\in[0,1]$ dazwischen (Grenzfall $a=0$ entspricht Rechteckfenster, Grenzfall $a=1$ Hann-Fenster)

  $$w[n]=\begin{cases}\frac{1}{2}\left(1+\cos\left(\pi\frac{N\frac{1-a}{2}-n}{N\frac{1-a}{2}}\right)\right),& 0\leq n < N\frac{1-a}{2}\\
    1,& N\frac{1-a}{2} \leq n \leq N\frac{1+a}{2}\\
    \frac{1}{2}\left(1+\cos\left(\pi\frac{n-N\frac{1+a}{2}}{N\frac{1-a}{2}}\right)\right),& N\frac{1+a}{2} < n \leq N\\0,& \text{otherwise}\end{cases}$$

* **Bartlett**: Dreiecksverlauf mit Nullen an den Endpunkten

  $$w[n]=\begin{cases}1-\left|\frac{2n}{N}-1\right|,& 0\leq n \leq N\\0,& \text{otherwise}\end{cases}$$

* **Kaiser**: basiert auf einer Besselfunktion mit Parameter $b$, größere Werte von $a$ erhöhen die Sperrdämpfung, verbreitern allerdings auch die Hauptkeule, typischer Wert: $b=3$

* weitere, hier nicht betrachtete Fensterfunktionen: **Blackman**, **Gauss**, **Dreieck**, **Lanczos**, **Cosinus**"

# ╔═╡ fda49044-3ac5-11eb-0c24-93d13c82705d
h(n, Ωg, N) = Ωg/π * sinc(Ωg*(n-N/2)/π)

# ╔═╡ 06afbe16-3ac6-11eb-23aa-49aeee7aeefb
md"""
Einstellung der Parameter:
* Filterordnung $N$: $(@bind N Slider(10:100))
* normierte Grenzkreisfrequenz $\Omega_g$: $(@bind Ωg Slider(0:.1:π))
* Tukey-Parameter $a$: $(@bind a_tukey Slider(0:.1:1))
* Kaiser-Parameter $b$: $(@bind a_kaiser Slider(.5:.1:6))
* Fensterfunktion: $(@bind win_type Select(["Rechteck", "Hann", "Hamming", "Bartlett", "Tukey", "Kaiser"]))
"""

# ╔═╡ e23de5b6-3ac6-11eb-245f-a37826b0c7be
md"Werte der Parameter:
* Filterordnung $N=$ $(N)
* normierte Grenzkreisfrequenz $\Omega_g=$ $(Ωg)
* Tukey-Parameter $a=$ $(a_tukey)
* Kaiser-Parameter $b=$ $(a_kaiser)
"

# ╔═╡ 13b88218-3ac7-11eb-3aa2-877e5f56ced9
begin
	if win_type == "Rechteck"
        window = rect(N+1)
    elseif win_type == "Hann"
        window = hanning(N+1)
    elseif win_type == "Hamming"
        window = hamming(N+1)
    elseif win_type == "Bartlett"
        window = bartlett(N+1)
    elseif win_type == "Tukey"
        window = tukey(N+1,a_tukey)
    elseif win_type == "Kaiser"
        window = kaiser(N+1,a_kaiser)
    end
    n = 0:N
    plot(n, window, label=false, lw=2, c = :lightblue)
    plot!(n, window, label=L"w[n]", lw=3, seriestype = :scatter, markercolor = :blue)
    plot!(n,h.(n, Ωg, N).*window, label = false, lw=2, c = :lightgreen)
    plot!(n,h.(n, Ωg, N).*window, label = L"h[n]", lw=3, seriestype = :scatter, markercolor = :green)
    ylims!((-.25,1.1))
    xlabel!(L"n")
end

# ╔═╡ 456426d6-3ac8-11eb-1ad0-01c5472e0a38
md"## Resultierender Frequenzgang

Der Frequenzgang $H^\prime(f)$ ermittelt sich aus der Übertragungsfunktion $H(z)$ wie folgt:

$$H^\prime(f) = H\left(e^{j 2 \pi f / f_s}\right).$$"

# ╔═╡ a24eaf38-3ac8-11eb-3252-e3cb8cd356ae
md"Normierte Grenzfrequenz $f_g/f_s=\Omega_g/(2\pi)$:"

# ╔═╡ ad82d87c-3ac8-11eb-0d58-37e409631f85
fg_norm = Ωg/2/π

# ╔═╡ 5a249128-3ac8-11eb-05d5-93e0a00cb89a
begin
	P = PolynomialRatio(h.(n,fg_norm*2*π,N).*window, [1])
    f_norm = 0:.001:.5
    H = freqz(P, f_norm*2*π)
    plot(f_norm, 20*log10.(abs.(H)), label=false, lw=2, linecolor = :blue)
    vline!([fg_norm], label = L"f = f_g", linecolor = :green)
    xlims!((0,.5))
    ylims!((-150,5))
    xlabel!(L"f/f_s")
    ylabel!(L"20 \log_{10}H^\prime(f)")
end

# ╔═╡ f03e11b4-3ac8-11eb-2574-fbef1b2c52c6
md"**Nachteil der Fenstermethode:**

Man kann mit der Ordnung $N$ und der Wahl der Fensterfunktion $w[n]$ den Frequenzgang dahingehend beeinflussen, dass er die Entwurfskriterien (maximale Welligkeit im Durchlassbereich, minimale Sperrdämpfung im Sperrbereich, obere Grenzfrequenz des Durchlassbereichs, untere Grenzfrequenz des Sperrbereichs) erfüllt. Damit kann die optimale Wahl von $N$ und $w[n]$ empirisch ermittelt werden.

**Optimalmethode:**

Im Gegensatz dazu ermittelt die Optimalmethode (siehe gesondertes Julia-Notebook) die Filterkoeffizienten so, dass bei gegebenen $N$, die Entwurfskriterien bestmöglich erfüllt werden. Des weiteren gibt es Algorithmen, die das kleinste $N$ ermitteln, für das die Entwurfskriterien bei Anwendung der Optimalmethode gerade noch erfüllt sind.

## Transformation in Hochpass, Bandpass und Bandsperre

Es folgen die Transformationsgleichungen zur Umwandlung des Tiefpassfilters in ein Hochpass-, Bandpass- bzw. Bandsperren-Filter.

* **Tiefpass** mit Grenzfrequenz $f_1:=f_g$: $h_{TP}[n]=h[n]$
* **Hochpass** mit Grenzfrequenz $f_1:=f_g$ ($N$ gerade): $h_{HP}[n]=\begin{cases}1-h[0],& n=\frac{N}{2}\\-h[n], & n\neq \frac{N}{2}\end{cases}$
* **Bandpass** mit unterer Grenzfrequenz $f_1$ und oberer Grenzfrequenz $f_2$: $h_{BP}[n]=\left.h[n]\right|_ {f_g=f_2} - \left.h[n]\right|_ {f_g=f_1}$
* **Bandsperre** mit unterer Grenzfrequenz $f_1$ und oberer Grenzfrequenz $f_2$: $h_{BS}[n]=\left.h[n]\right|_ {f_g=f_1}+\left.h_{HP}[n]\right|_ {f_g=f_2}$"

# ╔═╡ 245a50cc-3ac9-11eb-1a99-8f9aefe5c530
function filter_HP(n, Ωg, N)
    if n==N/2
        1 - Ωg/π
    else
        - h(n, Ωg, N)
    end
end

# ╔═╡ 2b3a5ed4-3ac9-11eb-289f-7912c7febe41
begin
	filter_TP(n, Ωg, N) = h(n, Ωg, N)
	filter_BP(n, Ω1, Ω2, N) = h(n, Ω2, N)-h(n, Ω1, N)
	filter_BS(n, Ω1, Ω2, N) = h(n, Ω1, N)+filter_HP(n, Ω2, N)
end

# ╔═╡ 93eec210-3ac9-11eb-2a8a-55ef961a8420
md"""
Einstellung der Parameter:
* Filterordnung $N$: $(@bind N_ Slider(10:2:100))
* erste normierte Grenzfrequenz $f_1/f_g$: $(@bind f1_norm Slider(0.1:.01:.25))
* zweite normierte Grenzfrequenz $f_2/f_g$: $(@bind f2_norm Slider(0.25:.01:.5))
* Filtertyp: $(@bind filter_type Select(["TP", "HP", "BP", "BS"]))
"""

# ╔═╡ 9c39d176-3ac9-11eb-3749-1b48657314f7
md"Werte der Parameter:
* Filterordnung $N=$ $(N_)
* erste normierte Grenzfrequenz $f_1/f_g=$ $(f1_norm)
* zweite normierte Grenzfrequenz $f_2/f_g=$ $(f2_norm)

Frequenzgang:"

# ╔═╡ bac3a23e-3ac9-11eb-0217-0b1d6540bb00
begin
	n_ = 0:N_
    if filter_type == "TP"
        P_ = PolynomialRatio(filter_TP.(n_, f1_norm*2*π, N), [1])
    elseif filter_type == "HP"
        P_ = PolynomialRatio(filter_HP.(n_, f1_norm*2*π, N), [1])
    elseif filter_type == "BP"
        P_ = PolynomialRatio(filter_BP.(n_, f1_norm*2*π, f2_norm*2*π, N), [1])
    elseif filter_type == "BS"
        P_ = PolynomialRatio(filter_BS.(n_, f1_norm*2*π, f2_norm*2*π, N), [1])
    end
    H_ = freqz(P_, f_norm*2*π)
    plot(f_norm, 20*log10.(abs.(H_)), label=false, lw=2, linecolor = :blue)
    vline!([f1_norm], label = L"f = f_1", linecolor = :green)
    if filter_type == "BP" || filter_type == "BS"
        vline!([f2_norm], label = L"f = f_2", linecolor = :green)
    end
    xlims!((0,.5))
    ylims!((-150,5))
    xlabel!(L"f/f_s")
    ylabel!(L"20 \log_{10}H^\prime(f)")
end

# ╔═╡ 4a97e1f0-3ac9-11eb-39bd-c5af7384d883
md"## Hinweis

Das Paket `DSP.jl` besitzt bereits fertige Funktionen zur Berechnung der Filterkoeffizienten nach der Fenstermethode, die hier, bis auf die Fensterfunktionen $w[n]$, nicht verwendet wurden. Den Funktionen `Lowpass`, `Highpass`, `Bandpass`, `Bandstop`, `FIRWindow`, der zahlreichen, bereits verwendeten Fensterfunktionen und `digitalfilter` berechnen die Filterkoeffizienten direkt.

Folgendes Beispiel berechnet die Filterkoeffizienten einer Bandsperre mit $f_1=100$, $f_2=200$, $f_s=1000$, $N+1=201$ und einem Hann-Fensters mit diesen Funktionen und stellt das resultierende Spektrum graphisch dar."

# ╔═╡ 63afdef6-3ac9-11eb-2dd7-2995958ed61e
begin
	filter_length=201
	f1_norm_=.1
	f2_norm_=.2
	responsetype = Bandstop(f1_norm_*2, f2_norm_*2)
	designmethod = FIRWindow(hanning(filter_length))
	b = digitalfilter(responsetype, designmethod)
end

# ╔═╡ 8416f39e-3ac9-11eb-2506-33d6015bfe5b
begin
	P__ = PolynomialRatio(b, [1])
	H__ = freqz(P__, f_norm*2*π)
	plot(f_norm, 20*log10.(abs.(H__)), label=false, lw=2, linecolor = :blue)
	vline!([f1_norm_], label = L"f = f_1", linecolor = :green)
	vline!([f2_norm_], label = L"f = f_2", linecolor = :green)
	xlims!((0,.5))
	ylims!((-100,5))
	xlabel!(L"f/f_s")
	ylabel!(L"20 \log_{10}H^\prime(f)")
end

# ╔═╡ Cell order:
# ╟─6f854664-3ac5-11eb-2675-9301f51b27b9
# ╠═93a36792-3ac5-11eb-2653-85d2ba3306f2
# ╟─d0193936-3ac5-11eb-1aa2-892211210245
# ╠═fda49044-3ac5-11eb-0c24-93d13c82705d
# ╟─06afbe16-3ac6-11eb-23aa-49aeee7aeefb
# ╟─e23de5b6-3ac6-11eb-245f-a37826b0c7be
# ╟─13b88218-3ac7-11eb-3aa2-877e5f56ced9
# ╟─456426d6-3ac8-11eb-1ad0-01c5472e0a38
# ╟─a24eaf38-3ac8-11eb-3252-e3cb8cd356ae
# ╠═ad82d87c-3ac8-11eb-0d58-37e409631f85
# ╠═5a249128-3ac8-11eb-05d5-93e0a00cb89a
# ╟─f03e11b4-3ac8-11eb-2574-fbef1b2c52c6
# ╠═245a50cc-3ac9-11eb-1a99-8f9aefe5c530
# ╠═2b3a5ed4-3ac9-11eb-289f-7912c7febe41
# ╟─93eec210-3ac9-11eb-2a8a-55ef961a8420
# ╟─9c39d176-3ac9-11eb-3749-1b48657314f7
# ╠═bac3a23e-3ac9-11eb-0217-0b1d6540bb00
# ╟─4a97e1f0-3ac9-11eb-39bd-c5af7384d883
# ╠═63afdef6-3ac9-11eb-2dd7-2995958ed61e
# ╠═8416f39e-3ac9-11eb-2506-33d6015bfe5b
