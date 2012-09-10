library(sensR)
#################################
## profBinom uses:
## pSeq <- seq(from = 1e-8, to = 1-1e-8, length = nProf)
## when cumputing the likelihood function. 


#################################
## Tests of triangle family object:
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
print(cbind(delta.lim,
            tri.linkinv=triangle()$linkinv(delta.lim),
            muEta=triangle()$mu.eta(delta.lim)), digits=7)

pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-9, 1/3 + 1e-7,
        4:9/10, 1-1e-7, 1-1e-9, 1, 2, Inf)
print(cbind(pc, triangle()$linkfun(pc)), digits=7)

1 - (pcc <- triangle()$linkinv(14))
triangle()$linkfun(pcc)
1 - (pcc <- triangle()$linkinv(14 + .1))
triangle()$linkfun(pcc)
1 - (pcc <- triangle()$linkinv(15))
triangle()$linkfun(pcc)

#################################
## Tests of 3-AFC family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values: 
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-4, 1/3 + 1e-3,
        4:9/10, 1-1e-7, 1-1e-8, 1, 2, Inf)

print(cbind(delta.lim,
            linkinv= threeAFC()$linkinv(delta.lim),
            muEta= threeAFC()$mu.eta(delta.lim)
            ), digits=7)
print(cbind(pc, threeAFC()$linkfun(pc)), digits=7)

## Going smooth toward zero:
threeAFC()$mu.eta(8 - 10^(0:-7))
diff(threeAFC()$mu.eta(8 - 10^(0:-7)))

## Upper limit for linkinv:
1 - threeAFC()$linkinv(9-.1) # 3.102338e-10
1 - threeAFC()$linkinv(9) # 0
1 - threeAFC()$linkinv(9+.1) # 0
## Lower limit for linkinv:
threeAFC()$linkinv(1e-9) - 1/3 # 1.791433e-10
threeAFC()$linkinv(1e-10) - 1/3 # 0

## The upper limit for linkfun:
threeAFC()$linkfun(1 - 1e-7) ## > Inf ~ 7.53
threeAFC()$linkfun(1 - 1.1e-8) ## < Inf ~ 8.0819
threeAFC()$linkfun(1 - 1e-8) ## Inf
## Lower limit for linkfun:
threeAFC()$linkfun(1/3 + 1e-5) # 6.103516e-05
threeAFC()$linkfun(1/3 + 1e-6) # 0

#################################
## Tests of tetrad family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values: 
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:3/10, 1/3, 1/3 + 1e-4, 1/3 + 1e-3,
        4:9/10, 1-1e-7, 1-1e-8, 1, 2, Inf)

print(cbind(delta.lim,
            linkinv= tetrad()$linkinv(delta.lim),
            muEta= tetrad()$mu.eta(delta.lim)
            ), digits=7)
print(cbind(pc, tetrad()$linkfun(pc)), digits=7)

## Going smooth toward zero:
tetrad()$mu.eta(8 - 10^(0:-7))
diff(tetrad()$mu.eta(8 - 10^(0:-7)))

## Upper limit for linkinv:
1 - tetrad()$linkinv(9-.1) # 5.737801e-10
1 - tetrad()$linkinv(9) # 0
1 - tetrad()$linkinv(9+.1) # 0
## Lower limit for linkinv:
tetrad()$linkinv(1e-4) - 1/3 # 2.043666e-09
tetrad()$linkinv(1e-5) - 1/3 # 2.242808e-10
tetrad()$linkinv(1e-6) - 1/3 # 2.060869e-10
tetrad()$linkinv(1e-7) - 1/3 # 2.059049e-10
tetrad()$linkinv(1e-8) - 1/3 # 0

## The upper limit for linkfun:
tetrad()$linkfun(1 - 1e-7) ## > Inf ~ 7.706576
tetrad()$linkfun(1 - 1.1e-8) ## < Inf ~ 8.235212
tetrad()$linkfun(1 - 1e-8) ## Inf
## Lower limit for linkfun:
tetrad()$linkfun(1/3 + 1e-8) # 0.0002441406
tetrad()$linkfun(1/3 + 1e-9) # 6.103516e-05
tetrad()$linkfun(1/3 + 1e-10) # 0
tetrad()$linkfun(1/3) # 0

#################################
## Tests of duo-trio family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values: 
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
pc <- c(-Inf, -1, 0:5/10, 0.5 + 1e-9, 0.5 + 1e-7,
        6:9/10, 1-1e-9, 1-1e-10, 1, 2, Inf)
print(cbind(delta.lim,
            DT.linkinv=duotrio()$linkinv(delta.lim),
            DT.mu.eta=duotrio()$mu.eta(delta.lim)), digits=7)
print(cbind(pc, duotrio()$linkfun(pc)), digits=7)

## linkinv - upper limit:
1 - duotrio()$linkinv(19)
1 - duotrio()$linkinv(20)
## eta > 20 -> mu = 1

## linkinv - lower limit:
.5 - duotrio()$linkinv(1e-8)
.5 - duotrio()$linkinv(1e-7)

## mu.eta - upper limit:
duotrio()$mu.eta(94)
duotrio()$mu.eta(95)

## mu.eta - lower limit
duotrio()$mu.eta(0)
duotrio()$mu.eta(1e-15)
duotrio()$mu.eta(1e-9)

## linkfun - upper limit
duotrio()$linkfun(1 - 1e-8)
duotrio()$linkfun(1 - 1e-9)
duotrio()$linkfun(1 - 1.000001e-10)
duotrio()$linkfun(1 - 1e-10)

## linkfun - lower limit:
duotrio()$linkfun(.5)
duotrio()$linkfun(.5 + 1e-10)
duotrio()$linkfun(.5 + 1.7e-10)
duotrio()$linkfun(.5 + 1.8e-10)
duotrio()$linkfun(.5 + 5e-10)
duotrio()$linkfun(.5 + 1e-9)
duotrio()$linkfun(.5 + 1e-8)
duotrio()$linkfun(.6)

#################################
## Tests of the 2-AFC family object:

## Evaluating linkinv, mu.eta and linkfun for a range of values: 
delta.lim <- c(-Inf, -1, 0, 1e-3, 1e-2, .1, 1:20, 30, Inf)
print(cbind(delta.lim,
            DT.linkinv=twoAFC()$linkinv(delta.lim),
            DT.mu.eta=twoAFC()$mu.eta(delta.lim)), digits=7)

pc <- c(-Inf, -1, 0:5/10, 0.5 + 1e-9, 0.5 + 1e-7,
        6:9/10, 1-1e-9, 1-1e-10, 1, 2, Inf)
print(cbind(pc, twoAFC()$linkfun(pc)), digits=7)

## linkinv - upper limit:
1 - twoAFC()$linkinv(11)
1 - twoAFC()$linkinv(12)
## eta > 11 -> mu = 1

## linkinv - lower limit:
.5 - twoAFC()$linkinv(1e-16)
.5 - twoAFC()$linkinv(1e-15)

## mu.eta - upper limit:
twoAFC()$mu.eta(54)
twoAFC()$mu.eta(55)

## mu.eta - lower limit
twoAFC()$mu.eta(0)

## linkfun - upper limit
twoAFC()$linkfun(1 - 1e-16)
twoAFC()$linkfun(1 - 1e-17)

## linkfun - lower limit:
twoAFC()$linkfun(.5)
twoAFC()$linkfun(.5 + 1e-17)
twoAFC()$linkfun(.5 + 1e-16)
