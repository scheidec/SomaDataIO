
# Submission 6.6.1

- This is a submission to an existing CRAN package.

- It contains a bug fix to an existing function.


## Reverse Dependencies

- There are no reverse dependencies.

- One package, `SomaScan.db`, has a reverse *Suggests*.

- It is an associated `BioConductor` package which we also maintain.


## R CMD check results
```
0 errors | 0 warnings | 0 notes
```

### Notes

0 `notes` were displayed during
`devtools::check(remote = TRUE, manual = TRUE)`

## Package Size

- Package tarball from `R CMD build` is 4.2Mb


## Example Timings from Win-builder

```
name	                  user	system	elapsed
Col.Meta	              0.29	0.00	0.29	
SeqId	                  0.06	0.00	0.06	
SomaDataIO-package	    0.15	0.00	0.16	
SomaScanObjects	        0.40	0.01	0.43	
adat-helpers	          0.06	0.00	0.06	
adat2eSet	              0.56	0.03	0.60	
addClass	              0.00 	0.00	0.00	
calcOutlierMap	        4.78	0.31	5.10	
calc_eLOD	              0.03	0.00	0.03	
cleanNames	            0.00	0.00	0.00	
diffAdats	              0.42	0.03	0.46	
getAnalyteInfo	        0.15	0.00	0.16	
getAnalytes	            0.02	0.00	0.01	
getOutlierIds	          2.03	0.07	2.09	
groupGenerics	          0.51	0.02	0.53	
is_intact_attr	        0.02	0.00	0.01	
is_seqFormat	          0.01	0.00	0.02	
lift_adat	              0.42	0.03	0.47	
loadAdatsAsList	        1.72	0.05	1.77	
medianNormalize	        0.00	0.00	0.00	
merge_clin	            0.14	0.00	0.14	
parseHeader	            0.12	0.00	0.12	
pivotExpressionSet	    0.13	0.00	0.12	
plot.Map              	2.07	0.11	2.19	
preProcessAdat	        2.19	0.05	2.25	
read_adat	              1.71	0.01	1.73	
read_annotations	      0.00	0.00	0.00	
reverseMedianNormalize	0.00	0.00	0.00	
rownames	              0.01	0.00	0.02	
soma_adat	              0.36	0.03	0.39	
transform	              0.00	0.00	0.00	
updateColMeta	          0.00	0.00	0.00	
write_adat	            0.37	0.06	0.44	
```
