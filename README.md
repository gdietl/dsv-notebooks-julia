#dsv-notebooks-julia

Julia-Notebooks zum besseren Verständnis der digitalen Signalverarbeitung

Öffne Jupyterlab im Browser mit:
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/gdietl/dsv-notebooks-julia/main?urlpath=lab)

Repository ist derzeit nicht öffentlich, da es bei der Ausführungder Notebooks mit mybinder.org Probleme mit WebIO gibt. Das längerfristige Ziel ist es allerdings, die Notebooks über mybinder.org zugänglich zu machen.

Notebooks müssen derzeit mit `jupyter-lab *.ipynb` lokal ausgeführt werden.

Falls bei der Ausführung auf mybinder.org zusätzliche Packete benötigt werden, so können Sie zu `Project.toml` folgendermaßen hinzugefügt werden:
```
]
activate .
add PACKAGE_NAME
```
