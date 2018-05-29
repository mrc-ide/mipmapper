
### TODO Functionality list:

---
  
1. ~~Dummy dataset creation function, for ease of examples/tests etc.~~
2. Dataset checking function. Something that looks for needed variables, and 
either identifies them intelligently or the user provides their locations.
3. Bit more option and control over plot_pca options
4. Thre is a poorly written generate_dummy_data function that hopes to generate
dummy data at large scale, where you can control the parameters about number of
samples, countries etc. At the moment the data is nonsense, and has no structure
or PCA utility, so might be nice to come up with imporvements to make faking
large clustered genetic data.
5. Shiny interface for loading a dataset, sliders for filtering, and then 
displayed plots. 
6. Fst permutation test