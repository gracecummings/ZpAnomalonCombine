# The Combine Incantations for the ZpAnomalon Analysis

These are some (hopefully all) of the scripts used to do combine tests for the ZpAnomalon Analysis. Theses are executed with a Combine environment set-up as in the [Combine Docs](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#setting-up-the-environment-and-installation) (visited 2022-09-16).

These scripts have built in name-parsing for the ZpAnomalon datacards, so if you are not using a datacard naming convention like

```
Run2_161718_ZllHbbMET_datacard_mumu_Zp5000ND800NS200_Zptcut100.0_Hptcut300.0_metcut75.0_btagwp0.8.txt
```

the automated output naming will not make sense. The script uses `mumu_` and `_Zpt` in the string parsing to get the signal sample name.

## Signal Injection Tests

The B2G Combine Review has very specific requests for the signal injections tests. It requires that four different signal strengths are injected:

1. 0 injected signal (background only)
2. Median expected limit
3. 1 sigma up expected limit
4. 1 sigmal down expected limit

If you already know the values you want to inject, you can use the `generateSignalInjectionFiles.py` script. This calls Combine to generate the toys, and run the FitDiagnostics incantation.```

```
python generateSignalInjectionFiles.py -d YOURDATACARD.txt -i EXTRA_NAMING_STRING -exps SIG_To_Inject_fb -t NTOYS
```

This outputs the standard fitDiagnostc output of combined, with the specific naming.

If you just want all four points done at one time, you can use `runSignalInjectionTestAtDiffR.py`. This runs the AsymptoticLimits method in Combine to find the expected limits, and then calls `generateSignalInjectionFiles.py` for you with these points.

```
python runSignalInjectionTestAtDiffR.py -d YOURDATACARD.txt -i EXTRA_NAMING_STRING -t NTOYS
```