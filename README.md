# StarStruc.jl
Stellar structure code in Julia

It's a glass cannon, but the 4 Msun case runs in less than a minute.

```julia
using StarStruc
star = Star(4 * Msun)
shootf!(star)
starplot(star)
```
