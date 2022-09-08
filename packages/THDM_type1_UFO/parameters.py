# This file was automatically created by FeynRules 2.3.24
# Mathematica version: 10.4.1 for Linux x86 (64-bit) (April 11, 2016)
# Date: Fri 7 Jul 2017 10:14:54



from object_library import all_parameters, Parameter


from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot

# This is a default parameter object representing 0.
ZERO = Parameter(name = 'ZERO',
                 nature = 'internal',
                 type = 'real',
                 value = '0.0',
                 texname = '0')

# User-defined parameters.
l2 = Parameter(name = 'l2',
               nature = 'external',
               type = 'real',
               value = 0.5,
               texname = '\\lambda _2',
               lhablock = 'Higgs',
               lhacode = [ 1 ])

l3 = Parameter(name = 'l3',
               nature = 'external',
               type = 'real',
               value = 1,
               texname = '\\lambda _3',
               lhablock = 'Higgs',
               lhacode = [ 2 ])

lR7 = Parameter(name = 'lR7',
                nature = 'external',
                type = 'real',
                value = 0.1,
                texname = '\\text{lR7}',
                lhablock = 'Higgs',
                lhacode = [ 3 ])

aEWM1 = Parameter(name = 'aEWM1',
                  nature = 'external',
                  type = 'real',
                  value = 127.9,
                  texname = '\\text{aEWM1}',
                  lhablock = 'SMINPUTS',
                  lhacode = [ 1 ])

Gf = Parameter(name = 'Gf',
               nature = 'external',
               type = 'real',
               value = 0.000011663900000000002,
               texname = '\\text{Gf}',
               lhablock = 'SMINPUTS',
               lhacode = [ 2 ])

aS = Parameter(name = 'aS',
               nature = 'external',
               type = 'real',
               value = 0.118,
               texname = '\\text{aS}',
               lhablock = 'SMINPUTS',
               lhacode = [ 3 ])

ymdo = Parameter(name = 'ymdo',
                 nature = 'external',
                 type = 'real',
                 value = 0.00504,
                 texname = '\\text{ymdo}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 1 ])

ymup = Parameter(name = 'ymup',
                 nature = 'external',
                 type = 'real',
                 value = 0.00255,
                 texname = '\\text{ymup}',
                 lhablock = 'YUKAWA',
                 lhacode = [ 2 ])

yms = Parameter(name = 'yms',
                nature = 'external',
                type = 'real',
                value = 0.101,
                texname = '\\text{yms}',
                lhablock = 'YUKAWA',
                lhacode = [ 3 ])

ymc = Parameter(name = 'ymc',
                nature = 'external',
                type = 'real',
                value = 1.27,
                texname = '\\text{ymc}',
                lhablock = 'YUKAWA',
                lhacode = [ 4 ])

ymb = Parameter(name = 'ymb',
                nature = 'external',
                type = 'real',
                value = 4.7,
                texname = '\\text{ymb}',
                lhablock = 'YUKAWA',
                lhacode = [ 5 ])

ymt = Parameter(name = 'ymt',
                nature = 'external',
                type = 'real',
                value = 172,
                texname = '\\text{ymt}',
                lhablock = 'YUKAWA',
                lhacode = [ 6 ])

yme = Parameter(name = 'yme',
                nature = 'external',
                type = 'real',
                value = 0.000511,
                texname = '\\text{yme}',
                lhablock = 'YUKAWA',
                lhacode = [ 11 ])

ymm = Parameter(name = 'ymm',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{ymm}',
                lhablock = 'YUKAWA',
                lhacode = [ 13 ])

ymtau = Parameter(name = 'ymtau',
                  nature = 'external',
                  type = 'real',
                  value = 1.777,
                  texname = '\\text{ymtau}',
                  lhablock = 'YUKAWA',
                  lhacode = [ 15 ])

TB = Parameter(name = 'TB',
               nature = 'external',
               type = 'real',
               value = 1.53748,
               texname = '\\text{TB}',
               lhablock = 'FRBlock',
               lhacode = [ 1 ])

sinbma = Parameter(name = 'sinbma',
                   nature = 'external',
                   type = 'real',
                   value = 0.1,
                   texname = '\\text{sinbma}',
                   lhablock = 'FRBlock',
                   lhacode = [ 2 ])

MZ = Parameter(name = 'MZ',
               nature = 'external',
               type = 'real',
               value = 91.1876,
               texname = '\\text{MZ}',
               lhablock = 'MASS',
               lhacode = [ 23 ])

Me = Parameter(name = 'Me',
               nature = 'external',
               type = 'real',
               value = 0.000511,
               texname = '\\text{Me}',
               lhablock = 'MASS',
               lhacode = [ 11 ])

MMU = Parameter(name = 'MMU',
                nature = 'external',
                type = 'real',
                value = 0.10566,
                texname = '\\text{MMU}',
                lhablock = 'MASS',
                lhacode = [ 13 ])

MTA = Parameter(name = 'MTA',
                nature = 'external',
                type = 'real',
                value = 1.777,
                texname = '\\text{MTA}',
                lhablock = 'MASS',
                lhacode = [ 15 ])

MU = Parameter(name = 'MU',
               nature = 'external',
               type = 'real',
               value = 0.00255,
               texname = 'M',
               lhablock = 'MASS',
               lhacode = [ 2 ])

MC = Parameter(name = 'MC',
               nature = 'external',
               type = 'real',
               value = 1.27,
               texname = '\\text{MC}',
               lhablock = 'MASS',
               lhacode = [ 4 ])

MT = Parameter(name = 'MT',
               nature = 'external',
               type = 'real',
               value = 172,
               texname = '\\text{MT}',
               lhablock = 'MASS',
               lhacode = [ 6 ])

MD = Parameter(name = 'MD',
               nature = 'external',
               type = 'real',
               value = 0.00504,
               texname = '\\text{MD}',
               lhablock = 'MASS',
               lhacode = [ 1 ])

MS = Parameter(name = 'MS',
               nature = 'external',
               type = 'real',
               value = 0.101,
               texname = '\\text{MS}',
               lhablock = 'MASS',
               lhacode = [ 3 ])

MB = Parameter(name = 'MB',
               nature = 'external',
               type = 'real',
               value = 4.7,
               texname = '\\text{MB}',
               lhablock = 'MASS',
               lhacode = [ 5 ])

mhc = Parameter(name = 'mhc',
                nature = 'external',
                type = 'real',
                value = 150,
                texname = '\\text{mhc}',
                lhablock = 'MASS',
                lhacode = [ 37 ])

mh1 = Parameter(name = 'mh1',
                nature = 'external',
                type = 'real',
                value = 120,
                texname = '\\text{mh1}',
                lhablock = 'MASS',
                lhacode = [ 25 ])

mh2 = Parameter(name = 'mh2',
                nature = 'external',
                type = 'real',
                value = 130,
                texname = '\\text{mh2}',
                lhablock = 'MASS',
                lhacode = [ 35 ])

mh3 = Parameter(name = 'mh3',
                nature = 'external',
                type = 'real',
                value = 140,
                texname = '\\text{mh3}',
                lhablock = 'MASS',
                lhacode = [ 36 ])

WZ = Parameter(name = 'WZ',
               nature = 'external',
               type = 'real',
               value = 2.4952,
               texname = '\\text{WZ}',
               lhablock = 'DECAY',
               lhacode = [ 23 ])

WW = Parameter(name = 'WW',
               nature = 'external',
               type = 'real',
               value = 2.085,
               texname = '\\text{WW}',
               lhablock = 'DECAY',
               lhacode = [ 24 ])

WT = Parameter(name = 'WT',
               nature = 'external',
               type = 'real',
               value = 1.50833649,
               texname = '\\text{WT}',
               lhablock = 'DECAY',
               lhacode = [ 6 ])

whc = Parameter(name = 'whc',
                nature = 'external',
                type = 'real',
                value = 1,
                texname = '\\text{whc}',
                lhablock = 'DECAY',
                lhacode = [ 37 ])

Wh1 = Parameter(name = 'Wh1',
                nature = 'external',
                type = 'real',
                value = 1,
                texname = '\\text{Wh1}',
                lhablock = 'DECAY',
                lhacode = [ 25 ])

Wh2 = Parameter(name = 'Wh2',
                nature = 'external',
                type = 'real',
                value = 1,
                texname = '\\text{Wh2}',
                lhablock = 'DECAY',
                lhacode = [ 35 ])

Wh3 = Parameter(name = 'Wh3',
                nature = 'external',
                type = 'real',
                value = 1,
                texname = '\\text{Wh3}',
                lhablock = 'DECAY',
                lhacode = [ 36 ])

CB = Parameter(name = 'CB',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1/(1 + TB**2))',
               texname = '\\text{CB}')

cosbma = Parameter(name = 'cosbma',
                   nature = 'internal',
                   type = 'real',
                   value = 'cmath.sqrt(1 - sinbma**2)',
                   texname = '\\text{cosbma}')

aEW = Parameter(name = 'aEW',
                nature = 'internal',
                type = 'real',
                value = '1/aEWM1',
                texname = '\\text{aEW}')

G = Parameter(name = 'G',
              nature = 'internal',
              type = 'real',
              value = '2*cmath.sqrt(aS)*cmath.sqrt(cmath.pi)',
              texname = 'G')

l7 = Parameter(name = 'l7',
               nature = 'internal',
               type = 'complex',
               value = 'lR7',
               texname = '\\lambda _7')

lI5 = Parameter(name = 'lI5',
                nature = 'internal',
                type = 'real',
                value = '0',
                texname = '\\text{lI5}')

lI6 = Parameter(name = 'lI6',
                nature = 'internal',
                type = 'real',
                value = '0',
                texname = '\\text{lI6}')

MW = Parameter(name = 'MW',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(MZ**2/2. + cmath.sqrt(MZ**4/4. - (aEW*cmath.pi*MZ**2)/(Gf*cmath.sqrt(2))))',
               texname = '\\text{MW}')

SB = Parameter(name = 'SB',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - CB**2)',
               texname = '\\text{SB}')

ee = Parameter(name = 'ee',
               nature = 'internal',
               type = 'real',
               value = '2*cmath.sqrt(aEW)*cmath.sqrt(cmath.pi)',
               texname = 'e')

CA = Parameter(name = 'CA',
               nature = 'internal',
               type = 'real',
               value = 'CB*cosbma + SB*sinbma',
               texname = '\\text{CA}')

sw2 = Parameter(name = 'sw2',
                nature = 'internal',
                type = 'real',
                value = '1 - MW**2/MZ**2',
                texname = '\\text{sw2}')

cw = Parameter(name = 'cw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - sw2)',
               texname = 'c_w')

SA = Parameter(name = 'SA',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(1 - CA**2)',
               texname = '\\text{SA}')

sw = Parameter(name = 'sw',
               nature = 'internal',
               type = 'real',
               value = 'cmath.sqrt(sw2)',
               texname = 's_w')

mu2 = Parameter(name = 'mu2',
                nature = 'internal',
                type = 'real',
                value = 'mhc**2 - (2*l3*MW**2*sw**2)/ee**2',
                texname = '\\text{mu2}')

l1 = Parameter(name = 'l1',
               nature = 'internal',
               type = 'real',
               value = '-(ee**2*(-(mh1**2*sinbma**2) - mh2**2*(1 - sinbma**2)))/(8.*MW**2*sw**2)',
               texname = '\\lambda _1')

l4 = Parameter(name = 'l4',
               nature = 'internal',
               type = 'real',
               value = '(ee**2*(2*mh1**2 + 2*mh2**2 + 4*mh3**2 - 8*mhc**2 + 2*(-mh1**2 + mh2**2)*cmath.cos(2*(cmath.pi/2. - cmath.asin(sinbma)))))/(16.*MW**2*sw**2)',
               texname = '\\lambda _4')

lR5 = Parameter(name = 'lR5',
                nature = 'internal',
                type = 'real',
                value = '(ee**2*(2*(mh1**2 + mh2**2 - 2*mh3**2) - 2*(mh1 - mh2)*(mh1 + mh2)*cmath.cos(2*(cmath.pi/2. - cmath.asin(sinbma)))))/(32.*MW**2*sw**2)',
                texname = '\\text{lR5}')

lR6 = Parameter(name = 'lR6',
                nature = 'internal',
                type = 'real',
                value = '(ee**2*(-mh1**2 + mh2**2)*sinbma*cmath.sqrt(1 - sinbma**2))/(4.*MW**2*sw**2)',
                texname = '\\text{lR6}')

yb = Parameter(name = 'yb',
               nature = 'internal',
               type = 'real',
               value = '(ee*ymb)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{yb}')

yc = Parameter(name = 'yc',
               nature = 'internal',
               type = 'real',
               value = '(ee*ymc)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{yc}')

ydo = Parameter(name = 'ydo',
                nature = 'internal',
                type = 'real',
                value = '(ee*ymdo)/(MW*sw*cmath.sqrt(2))',
                texname = '\\text{ydo}')

ye = Parameter(name = 'ye',
               nature = 'internal',
               type = 'real',
               value = '(ee*yme)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{ye}')

ym = Parameter(name = 'ym',
               nature = 'internal',
               type = 'real',
               value = '(ee*ymm)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{ym}')

ys = Parameter(name = 'ys',
               nature = 'internal',
               type = 'real',
               value = '(ee*yms)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{ys}')

yt = Parameter(name = 'yt',
               nature = 'internal',
               type = 'real',
               value = '(ee*ymt)/(MW*sw*cmath.sqrt(2))',
               texname = '\\text{yt}')

ytau = Parameter(name = 'ytau',
                 nature = 'internal',
                 type = 'real',
                 value = '(ee*ymtau)/(MW*sw*cmath.sqrt(2))',
                 texname = '\\text{ytau}')

yup = Parameter(name = 'yup',
                nature = 'internal',
                type = 'real',
                value = '(ee*ymup)/(MW*sw*cmath.sqrt(2))',
                texname = '\\text{yup}')

mu1 = Parameter(name = 'mu1',
                nature = 'internal',
                type = 'real',
                value = '(-4*l1*MW**2*sw**2)/ee**2',
                texname = '\\text{mu1}')

l5 = Parameter(name = 'l5',
               nature = 'internal',
               type = 'complex',
               value = 'complex(0,1)*lI5 + lR5',
               texname = '\\lambda _5')

l6 = Parameter(name = 'l6',
               nature = 'internal',
               type = 'complex',
               value = 'complex(0,1)*lI6 + lR6',
               texname = '\\lambda _6')

mu3 = Parameter(name = 'mu3',
                nature = 'internal',
                type = 'complex',
                value = '(-2*l6*MW**2*sw**2)/ee**2',
                texname = '\\text{mu3}')

