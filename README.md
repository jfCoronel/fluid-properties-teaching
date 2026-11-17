# fluid-properties-teaching

Scripts and examples accompanying the article

> Coronel, J.F., Pérez-Lombard, L., Velázquez Alonso, D. *From steam tables to the browser: fluid-property evaluation in engineering education.* Submitted to *Education for Chemical Engineers*.

## What is here

| Path | Content |
|---|---|
| `analysis.py` | Emulates the manual methods with [CoolProp](http://www.coolprop.org) 6.4.1: double linear interpolation in superheated-vapour tables of water and R-134a built on textbook-like grids (4000 random states per fluid), and a standard R-134a vapour-compression cycle solved with the equation of state, by table interpolation and with simulated chart-reading errors (Monte Carlo, 20 000 readings). |
| `figs.py` | Draws Figures 1–3 of the article (PNG, TIFF and PDF). |
| `results/analysis.json` | The results reported in the article, for comparison. |
| `examples/` | The three [fProperties](https://fproperties.jfcoronel.org) problems shown in the article, in the fProperties problem format. |

## Reproducing the results

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python analysis.py
.venv/bin/python figs.py figures
```

`analysis.py` uses a fixed random seed, so the `analysis.json` it writes should match `results/analysis.json`. `figs.py` uses Arial; matplotlib falls back to another font if it is not installed.

## The fProperties examples

Each example opens directly in the browser:

| Example | Link |
|---|---|
| Vapour-compression cycle with R-134a | https://fproperties.jfcoronel.org/#ex=r134a-cycle |
| Student submission with two inconsistent processes | https://fproperties.jfcoronel.org/#ex=r134a-student-errors |
| Mixing and cooling coil on the psychrometric chart | https://fproperties.jfcoronel.org/#ex=ahu-mixing-coil |

The JSON files in `examples/` can also be loaded with the import button of fProperties.

## Licence

MIT. See [LICENSE](LICENSE).
