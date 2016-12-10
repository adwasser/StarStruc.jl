# StarStruc.jl
Stellar structure code in Julia

It's a glass cannon, but the 4 Msun case runs on the order of a minute.

To install from Julia, use
```julia
Pkg.clone("https://github.com/adwasser/StarStruc.jl.git")
```

StarStruc exports the following names:
* shootf,
* shootf!
* profiles 
* Msun
* Rsun
* Star
* starplot

To run a model and plot the resulting structure profiles:
```julia
using StarStruc
star = Star(4 * Msun)
shootf!(star)
starplot(star)
```
