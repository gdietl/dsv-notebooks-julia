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

# ╔═╡ 42fdc084-3a5d-11eb-3764-67389120ce6f
begin
	using Pkg
	Pkg.add("FFTW")
	Pkg.add("DSP")
	Pkg.add("LaTeXStrings")
	Pkg.add("LinearAlgebra")
	Pkg.add("PlutoUI")

	using Plots, FFTW, DSP, LaTeXStrings, LinearAlgebra, PlutoUI
end

# ╔═╡ 0cf506b6-3a5d-11eb-1113-195b58449428
md"# Spektralanalyse

Guido Dietl, 9. Dezember 2020

In diesem Julia-Notebook werden populäre Vefahren zur Schätzung des Leistungsdichte- bzw. Amplitudenspektrums untersucht. Der Schwerpunkt liegt hier in der Anwendung und weniger in der Herleitung bzw. Beschreibung der Verfahren.

## Initialisierung

Zunächst müssen zur Ausführung der Julia-Skripte die benötigten Bibliotheken geladen werden. Das Hinzügen der Pakete ist für die Ausführung auf Binder notwendig, kann aber entfallen, wenn die Pakete bereits installiert wurden."

# ╔═╡ 8fd94040-3a5d-11eb-291b-b12462e55d97
md"Die Abtastfrequenz sei im Folgenden als $f_\mathrm{s}=100\,\mathrm{kHz}$ angenommen. Die folgenden Skripte können auch für jede andere Abtastfrequenz ausgeführt werden. Gerne kann zum Ausprobieren die Abtastfrequenz `fs` in der folgenden Zeile geändert werden."

# ╔═╡ 9b0e1c06-3a5d-11eb-151e-b79dfcbbb363
fs=100e3;

# ╔═╡ cf539192-3a5d-11eb-2cc7-6552b1a4bea3
md"""
## Theoretische Grundlagen

### Zeitdiskrete Fouriertransformation eines abgetasteten Signals

Methoden zur Schätzungen des Leistungsdichtespektrums $\Phi_x(f)$ eines Zeitsignals $x(t)\in\mathbb{R}$ sind numerisch nur mit Hilfe der Abtastung des Signals möglich. Wir tasten somit $x(t)$ mit der Abtastfrequenz $f_\mathrm{s}=1/T_\mathrm{s}$ ab und erhalten das zeitdiskreten Signals

$$x[k]=x(kT_\mathrm{s}), k\in\mathbb{Z}.$$ 

Zur Berechnung der diskreten Fouriertransformation $X[\ell]\in\mathbb{C}$ muss $x[k]$ $N$-periodisch sein, d.h. $x[k]=x[k+N]$. Es reicht daher lediglich die $N$ Werte von $x[k]$, $k\in\{0,2,\dots,N-1\}$, einer Periode zu betrachten. Ist $x[k]$ nicht periodisch, dann wird lediglich ein Ausschnitt von $x[k]$ mit $N$ Werten betrachtet. Man sagt, dass Signal $x[k]$ wird mit einer Rechteckfunktion der Länge $N$ gefenstert. Zur Berechnung von $X[\ell]$ denkt man sich dieses Fenster dann periodisch fortgesetzt. Beachte, dass $X[\ell]$ ebenfalls $N$-periodisch ist, d.h. $X[\ell]=X[\ell+N]$. Auch hier reicht es aus lediglich $N$ Werte zu betrachten. Die diskrete Fouriertransformierte berechnet sich schließlich zu

$$X[\ell]=\sum_{k=0}^{N-1} x[k]e^{-j k\ell \frac{2 \pi}{N}},\ \ell\in\{0,1,\dots,N-1\}.$$

Entspricht $N$ einer Potenz von $2$, d.h. $N=2^i$, $i\in\mathbb{N}$, dann kann die diskreten Fouriertransformation $X[\ell]$ effizient mit Hilfe der *Fast-Fourier-Transform* (FFT) implementiert werden. Die Ordnung der Komplexität wird mit diesem Verfahren von ${\cal O}(N^2)$ auf ${\cal O}(N\log_2 N)$ reduziert.

### Zusammenhang zwischen der zeitdiskreten Fouriertransformation und der Fourierreihe

Der Zusammenhang zwischen der diskreten Fouriertransformation $X[\ell]$ und den Fourierkoeffizienten, d.h. den Amplituden $a_\ell\in\mathbb{R}_0^+$ und Phasen $\alpha_\ell\in[0,2\pi[$, der Fourierreihe eines $T$-periodischen Signals

$$x(t)=a_0+\sum_{\ell=1}^\infty a_\ell \cos\left(2 \pi \ell \frac{t}{T} + \alpha_\ell\right),$$

kann wie folgt zusammengefasst werden (siehe Vorlesungsskript):

* Wert bei $\ell=0$: $X[0]\approx Na_0$,
* Werte für $\ell\in\left\{0,1,\dots,\left\lfloor\frac{N}{2}\right\rfloor\right\}$: $X[\ell]\approx \frac{N}{2}a_\ell e^{j\alpha_\ell}$ und
* Werte für $\ell\in\left\{\left\lfloor\frac{N}{2}\right\rfloor+1, \left\lfloor\frac{N}{2}\right\rfloor+2,\dots,N-1\right\}$: $X[\ell]\approx \left(X\left[N-\ell\right]\right)^*$.

Der letzte Zusammenhang gilt aufgrund der Achsensymmetrie des Realteils und der Punktsymmetrie des Imaginärteils bzw. der Achsensymmetrie des Betrags und der Punktsymmetrie der Phase von $X[\ell]$. Die Ursache dieser Symmetrien liegt in der Tatsache, dass $x[k]\in\mathbb{R}$ reellwertig ist.

Als Amplitudenspektrum $A_x[\ell]$ definieren wir im Folgenden über die Amplituden der Fourierreihe, d.h.

$$a_\ell\approx A_x[\ell]:=\begin{cases}
\frac{X[0]}{N}, & \ell=0, \\
\frac{2|X[\ell]|}{N}, & \mathrm{sonst}.
\end{cases}$$

Eine Periode $T$ enthält $N$ Abtastwerte im Abstand von $1/f_\mathrm{s}$. Es gilt damit $T=N/f_\mathrm{s}$. Der Abstand zweier Werte des Spektrums entspricht damit einer Frequenzdifferenz von $f_\mathrm{s}$ und es folgt:

$$f=\ell f_\mathrm{s}.$$

## Beispiel eines stochastischen Signals zur Untersuchung

Um die Verfahren zur Spektralanalyse untersuchen zu können, definieren wir zunächst ein stochastisches Signal, das aus zwei unterschiedlich gewichteten Sinussignalen, überlagert mit additiv weissem Gauß'schen Rauschen, besteht, d.h.
$$x(t)=a_1\sin\left(2 \pi f_1 t + \varphi_1\right) + a_2\sin\left(2 \pi f_2 t + \varphi_2\right) + n(t).$$
Dabei seien die Amplituden exemplarisch $a_1=1$ bzw. $a_2=0{,}2$, die Frequenzen $f_1=1{,}5\,\mathrm{kHz}$ bzw. $f_2=1\,\mathrm{kHz}$ und die Phasen $\varphi_1=0$ bzw. $\varphi_2=\frac{\pi}{2}$. Die beiden Sinussignale besitzen also Frequenzen, die sehr nahe beieinander liegen und sehr unterschiedliche Amplituden aufweisen. Die Werte des Rauschsignals $n(t)$ folgen einer Normalverteilung mit Mittelwert $m_x=0$ und Varianz $\sigma_x^2=1$. Benachbarte Werte von $n(t)$ sind nicht korreliert.

Folgende Funktion `signal_gen` erzeugt für einen gegebenen Vektor `t` von Abtastzeitpunkten $t$ den dazugehörigen Vektor `x` von Signalwerte, die den Abtastwerten von $x(t)$ entsprechen. Um später auch auf das unverrauschte Signal zugreifen zu können, wird dieses als Vektor `x0` mit ausgegebenen.
"""

# ╔═╡ 98c8a738-3a5e-11eb-0dd0-47174081d391
function signal_gen(t)
    a1 = 1.
    f1 = 1.5e3
    phi1 = 0.
    a2 = .2
    f2 = 1e3
    phi2 = pi/2
    x0 = a1*sin.(2*pi*f1*t.+phi1) .+ a2*sin.(2*pi*f2*t.+phi2)
    # Additiver weisser Gauß'scher Rauschvektor mit Mittelwert Null und Varianz 1
    # (Länge entspricht der Länge von t)
    n = randn(length(t))
    x = x0 .+ n
    return x, x0
end

# ╔═╡ a638d6ba-3a5e-11eb-2d05-7919d0c608b2
md"Um sich das Signal $x(t)$ besser vorstellen zu können, erzeugt folgender Code eine möglich Realisierung. Durch wiederholte Ausführung des Codes können weitere Realisierungen generiert werden. Das Signal $x(t)$ wird für $t/\mathrm{s}\in[0,1]$ dargestellt. Das zweite Diagramm zeigt einen kleineren Zeitausschnitt."

# ╔═╡ b43ba904-3a5e-11eb-3441-b5ed01c305e2
begin
	# Erzeuge Zeit- und Signalvektor
	t = LinRange(0,1, Int(fs))
	x, x0 = signal_gen(t)
	# Ausgabe des Ergbnisses
	plot(plot(t,x, label=L"x(t)"), plot!(t,x0, label=L"x(t)-n(t)"),
    plot(t[100:500],x[100:500], label=L"x(t)"), plot!(t[100:500],x0[100:500], 			label=L"x(t)-n(t)"))
	xlabel!(L"t/s")
end

# ╔═╡ e70c2c1e-3a5e-11eb-1ea9-c9341955442d
md"## Schätzung des Amplitudenspektrums

Populäre Methoden zur Schätzung des Leistungsdichte- bzw. Leistungsspektrums sind das Periodogramm, das Bartlett- und das Welch-Verfahren. Das Periodogramm basiert auf den Werten der FFT, wobei die Werte derart gewichtet werden, dass sie einer Schätzung der Leistungsdichte entsprechen. Der Zusammenhang zwischen dem Leistungsdichtespektrum $\text{PSD}_x[\ell]$ und dem oben beschriebenen Amplitudenspektrum $A_x[\ell]$ ergibt sich zu
$$\text{PSD}_x[\ell]=\frac{\left(A_x[\ell]\right)^2}{2}\cdot\frac{N}{f_s}\ \mathrm{bzw.}\ A_x[\ell]=\sqrt{\frac{2f_s}{N}\text{PSD}_x[\ell]},$$
unter der Annahme von sinusförmigen Basisfunktionen. 

Die Bartlett- und die Welch-Methode betrachten eine Mittelung von Periodogrammen, um den Einfluß des Rauschens im Falle von stochastischen Signalen zu minimieren. Dabei findet die Mittelung über kleine Zeitausschnitte, so genannte Zeitfenster, statt. Die Welch-Methode unterscheidet sich von der Bartlett-Methode in der Tatsache, dass sie eine Überlappung der Zeitfenster zulässt.

Die Funktionen `periodogram`, `welch_pgram`, `freq` und `power` des Julia-Pakets `DSP.jl` berechnen ein Leistungsdichtespektrum zu einem gegebenen Signalvektor (die Dokumentation von `power` ist irreführend, da hier von Leistung anstelle von Leistungsdichte gesprochen wird, obwohl eindeutig die Leistungsdichte berechnet wird). Die FFT-Funktion ist als `fft` im Paket `FFTW.jl` enthalten. Sie sollen im folgenden Programmcode verwendet werden, um das Amplitudenspektrum des obigen stochastischen Signals zu bestimmen. Dabei spielen folgende Parameter eine entscheidende Rolle:

* die zur Untersuchung des Signals verwendete Abtastfrequenz $f_\mathrm{sd}=f_\mathrm{s}/\delta$ (Parameter `fsd`), $\delta\in\mathbb{N}$  (Parameter `δ`), die durch *Downsampling* des gegebenen Signals mit der ursprünglichen Abtastfrequenz $f_\mathrm{s}$ entsteht ($\delta$ nennt man den *Downsampling*-Faktor),
* Anzahl $N\in\mathbb{N}$ (Parameter `N`) der Abtastwerte pro Zeitfenster,
* Größe $N_\mathrm{FFT}=2^{n_\mathrm{FFT}}$ (Parameter `Nfft`), $n_\mathrm{FFT}\in\mathbb{N}$ (Parameter `nfft`), der FFT (Falls $N_\mathrm{FFT}>N$, dann werden nicht vorhandene Werte mit Nullen aufgefüllt. Dieses Vorgehen bezeichnet man als Zero-Padding und dient dazu, die Anzahl der Werte im Spektrum zu erhöhen, ohne dabei allerdings die Auflösung zu verbessern. Es handelt sich also ausschließlich um eine Interpolation.) und
* der Parameter `wtype`: Fensterfunktion zur Gewichtung der Werte eines Zeitfensters (siehe <a href=https://de.wikipedia.org/wiki/Fensterfunktion>hier</a> für Details zu den verschiedenen Fensterfunktionen).

Beachten Sie, dass im folgenden Programmcode die Fensterfunktion lediglich im Periodogramm und im Welch-Verfahren angewendet wird. Das Ergebnis der FFT basiert stets auf dem Rechteckfenster."

# ╔═╡ 28714c48-3a5f-11eb-0c02-5b07adb0fa50
md"""
Einstellung der Parameter:
* Downsampling-Faktor $\delta$: $(@bind δ Slider(1:1:15))
* Anzahl $N$ der Abtastwerte pro Zeitfenster: $(@bind N Slider(100:10:10000))
* Größe $N_\mathrm{FFT}$ der FFT: $(@bind nfft Slider(1:1:16))
* additive Störgröße: $(@bind noise Select(["ohne Rauschen", "mit Rauschen"]))
* Fensterfunktion: $(@bind wtype Select(["Rechteck", "Hann", "Hamming", "Bartlett"]))
"""

# ╔═╡ e35ce606-3a60-11eb-03ae-3dde234369f4
md"Werte der Parameter:
* Downsampling-Faktor $\delta=$ $(δ)
* Anzahl der Abtastwerte: $N=$ $(N)
* Größe der FFT: $N_\mathrm{FFT}=$ $(nfft)"

# ╔═╡ 1deba4d0-3a61-11eb-10a9-ade8a64ccca8
begin
	if wtype == "Rechteck"
        window = rect(N)
    elseif wtype == "Hann"
        window = hanning(N)
    elseif wtype == "Hamming"
        window = hamming(N)
    elseif wtype == "Bartlett"
        window = bartlett(N)
    end
    # Normalisiere alle Fenster auf den Wert sqrt(N)
    window = window/norm(window)*sqrt(N)
    Nfft=2^nfft
    # Downsampling
    if noise=="mit Rauschen"
        xsd = x[1:δ:end]
    else
        xsd = x0[1:δ:end]
    end
    fsd = fs/δ
    
    # Berechnung des Amplitudenspektrums FFT-basiert
    Xfft = abs.(fft([xsd[1:N].*window; zeros(max(0,Nfft-N), 1)]))/N
    Xfft[2:end] = Xfft[2:end]*2 # Umrechnung in Koeffizienten der Fourierreihe (Faktor 2 für k>0)
    
    # Berechnung des Leistungsspektrums mit Periodogramm
    p = periodogram(xsd[1:N]; onesided=true, nfft=max(N,Nfft), fs=fsd, window=window)
    f_P = freq(p)
    Pxx_P = power(p)
    
    # Berechnung des Leistungsspektrums mit Welch-Methode
    w = welch_pgram(xsd, N, 0; onesided=true, nfft=max(N,Nfft), fs=fsd, window=window)
    f_W = freq(w)
    Pxx_W = power(w)
    
    f_norm = LinRange(0, fsd, length(Xfft))
    plot(f_norm, Xfft, yscale=:log10, label="FFT-basiert", lw=2, linecolor = :blue)
    #= Darstellung des Amplitudenspektrums (Wurzel des Leistungsdichtespektrums mal Wurzel der spektralen Breite fsd/N eines Frequenzbins
    mal Wurzel von 2 falls sinusförmige Größen angenommen werden)=#
    plot!(f_P, sqrt.(2*Pxx_P*fsd/N), yscale=:log10, label="Periodogramm", lw=2, linecolor = :red)
    plot!(f_W, sqrt.(2*Pxx_W*fsd/N), yscale=:log10, label="Welch", lw=2, linecolor = :green, grid = (:all, :grey, :dash, 0.5, 0.9), minorgrid = (:all, :grey, :dash, 0.7, 0.9))
    xlims!((.5e3,2e3))
    ylims!((.001,1.2))
    xlabel!(L"f/\mathrm{Hz}")
    ylabel!("Amplitudenspektrum")
end

# ╔═╡ 34f7d5ae-3a5f-11eb-04da-2b689e926175
md"## Aufgaben

Spielen Sie mit den Parametern und versuchen Sie folgende Fragen zu beantworten:
* Wie ändert sich das Spektrum mit der Abtastfrequenz? Begründen Sie das Verhalten.
* Wie ändert sich das Spektrum mit der Länge des Zeitfensters? Begründen Sie das Verhalten.
* Wie wirkt die Gewichtung des Zeitfensters auf das Spektrum aus? Wann ist welches Fenster zu bevorzugen?
* Welchen Effekt hat Zero-Padding?"

# ╔═╡ Cell order:
# ╟─0cf506b6-3a5d-11eb-1113-195b58449428
# ╠═42fdc084-3a5d-11eb-3764-67389120ce6f
# ╟─8fd94040-3a5d-11eb-291b-b12462e55d97
# ╠═9b0e1c06-3a5d-11eb-151e-b79dfcbbb363
# ╟─cf539192-3a5d-11eb-2cc7-6552b1a4bea3
# ╠═98c8a738-3a5e-11eb-0dd0-47174081d391
# ╟─a638d6ba-3a5e-11eb-2d05-7919d0c608b2
# ╠═b43ba904-3a5e-11eb-3441-b5ed01c305e2
# ╟─e70c2c1e-3a5e-11eb-1ea9-c9341955442d
# ╟─28714c48-3a5f-11eb-0c02-5b07adb0fa50
# ╟─e35ce606-3a60-11eb-03ae-3dde234369f4
# ╠═1deba4d0-3a61-11eb-10a9-ade8a64ccca8
# ╟─34f7d5ae-3a5f-11eb-04da-2b689e926175
