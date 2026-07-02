### Capabilities

A web-server application using TurboGears to carry out an 
in silico enzymatic digest of a user-provided protein sequence


Users should be able to specify:
```
- min and max length
- min and max molecular weight
- number of missed cleavages
- specific enzyme
```

Output is a table of peptides, with their:
```- length
- molecular weight
- number of missed cleavages
- amino acids to the left and right of eaach peptide in the protein sequence
```
