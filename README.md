# covid19model-reanalysis
This code is derived from https://github.com/ImperialCollegeLondon/covid19model
(only the `nature` directory).
The original code is the one that was used in Flaxman, Mishra, Gandy et al.
"Estimating the effects of non-pharmaceutical interventions on COVID-19 in
Europe," Nature, 2020. https://www.nature.com/articles/s41586-020-2405-7.

The following commands were used to produce our results:

```
Rscript run-model.r --full --fixseed flaxman
Rscript run-model.r --full --fixseed gomes

Rscript run-model.r --full --fixseed --forecast 14 flaxman
Rscript run-model.r --full --fixseed --forecast 14 gomes

Rscript run-model.r --full --fixseed --scaleifr 0.275 flaxman
Rscript run-model.r --full --fixseed --scaleifr 0.275 gomes

Rscript run-model.r --full --fixseed --scaleifr 0.275 --forecast 14 flaxman
Rscript run-model.r --full --fixseed --scaleifr 0.275 --forecast 14 gomes

Rscript run-model.r --full --fixseed --gradeffect 0.7 flaxman
Rscript run-model.r --full --fixseed --gradeffect 0.7 gomes

Rscript run-model.r --full --fixseed flaxman-halfnormal
Rscript run-model.r --full --fixseed gomes-halfnormal

Rscript run-model.r --full --fixseed flaxman-nolastintervention
Rscript run-model.r --full --fixseed gomes-nolastintervention
```
