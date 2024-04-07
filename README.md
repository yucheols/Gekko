# Assessing the establishment potentials for the two gecko species introduced to South Korea (Gekko japonicus and G. swinhonis)
MaxEnt ecological niche modeling workflow based on the SDMtune package (Vignali et al. 2020) and niche analyses using the packages ENMTools (Warren et al. 2021) and ecospat (Di Cola et al. 2017).



## Study background
![Figure_1](https://github.com/yucheols/Gekko/assets/85914125/7687a0a5-a940-4038-b85c-b759949b2af6)

- Currently, there are two gecko species known from South Korea: Gekko japonicus (blue dots) and G. swinhonis (pink dots). Both species are not native to South Korea.
- G. japonicus was known to be present in South Korea from at least 1885. The presence of G. swinhonis is first reported in 2021.
- In South Korea, G. japonicus is found mostly in the southern coastal regions. On the other hand, G. swinhonis is only found in several locations in northwestern coastal region of South Korea.
- Both gecko species were found in sympatry in this northwestern coastal site, implying that G. japonicus may have additional suitable habitats in South Korea.
- Therefore, we aimed to map out potential suitable habitats of the two gecko species across South Korea using ecological niche modeling. Furthermore, we aimed to gain insights into the patterns of partial range overlaps of these two species given the environmental conditions of South Korea.  

## Ecological niche modeling
![image](https://github.com/yucheols/Gekko/assets/85914125/3ccb9d09-bdfa-484d-a5d6-dbfa1ced0c70)

- We modeled the current and future suitable habitats of two gecko species in South Korea using the MaxEnt algorithm.
- The model predictions under current environment conditions suggest that there are broad areas of suitable, but currently unoccupied, habitats in South Korea for both species.

## Niche quantification
![Figure_6](https://github.com/yucheols/Gekko/assets/85914125/81b89675-745e-4280-9395-4ffe22cde992)

- We quantified and compared the niches of two gecko species using the R packaged ENMTools and ecospat.
- To do so, we conducted both niche identity test and symmetric background test.
- The results suggest non-identical niches attributable to underlying environmental differences in the ranges of the two gecko species.

## Dataset
- The "Dataset" folder contains all the data necessary to conduct the analyses done in this paper.
- The "bg" subfolder contains the background points for each species that were used for MaxEnt modeling.
- The "envs" subfolder contains the raster layers used for model calibration.
- The "occs" subfolder contains the spatially thinned occurrence points for each species that were used for MaxEnt modeling.
- The "proj_envs" subfolder contains the future climate layers used for model transfer.

## Citation
A research article associated with this project is currently in press in the journal NeoBiota.

```
I-K Park, Y Shin, H-J Baek, J Kim, D-I Kim, M Seok, Y Oh, D Park. 2024. Establishment potential for the two gecko
species adapted to different climates, Gekko japonicus and G. swinhonis, introduced to South Korea. NeoBiota.
```
