# Recipes of Monte Carlo simulations

These are the four exemplified Monte Carlo simulations in Figure 4 of [Fieremans and Lee, NeuroImage 2018](https://doi.org/10.1016/j.neuroimage.2018.06.046) with more details in Supplmentary Information

**Sept 21st, 2018: Update Densely packed cylinders**

* **Example 1 (Figure 4, point 1):** Free diffusion in 2d
* **Example 2 (Figure 4, point 2):** Check short-time limit of diffusion in a geometry composed of randomly packed impermeable cylinders (2d)
* **Example 3 (Figure 4, point 3):** Check against analytical formulas for diffusion within an impermeable non-absorbing cylinder (2d)
* **Example 4-5 (Figure 4, point 5):** Calculate membrane's permeability by starting diffusing particles from the center of a permeable cylinder (2d)
* **Analytical solution** of time-dependent diffisivity and kurtosis between parallel planes, inside cylinders, and inside spheres.
* **Densely packed cylinders:** Generation of randomly packed cylinders with the freedom of tuning axonal water fraction, inner diameter distribution, and g-ratio.
* **Densely packed spheres:** Generation of randomly packed spheres with the freedom of tuning cellular water fraction, inner diameter distribution, and g-ratio.

These are good exercises if you just start your own MC simulation codes.
Some results can suprise you, even if you are well experienced!!

## References
1. If you use this code in your study, please cite:
Els Fieremans, Hong-Hsi Lee, [Physical and numerical phantoms for the validation of brain microstructural MRI: A cookbook](https://doi.org/10.1016/j.neuroimage.2018.06.046), NeuroImage 2018

2. If you use this code to generate your own packing geometry, please **also** cite:
Aleksandar Donev, Salvatore Torquato, Frank H. Stillinger, [Neighbor list collision-driven molecular dynamics simulation for nonspherical hard particles. I. Algorithmic details](https://doi.org/10.1016/j.jcp.2004.08.014), Journal of Computational Physics 2005

3. If you use the axonal diameter histogram of corpus callosum in this code, pleas **also** cite:
Francisco Aboitiz, Arnold B. Scheibel, Robin S. Fisher, Eran Zaidel, [Fiber composition of the human corpus callosum](https://doi.org/10.1016/0006-8993(92)90178-C), Brain Research 1992

## Authors
Hong-Hsi Lee - [our page](http://www.diffusion-mri.com/people/hong-hsi-lee)

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/leehhtw/monte-carlo-simulation-recipes/blob/example1/LICENSE) file for details
