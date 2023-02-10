# This file was automatically created by FeynRules 2.3.24
# Mathematica version: 10.4.1 for Linux x86 (64-bit) (April 11, 2016)
# Date: Fri 7 Jul 2017 10:14:54


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-(ee*complex(0,1))/3.',
                order = {'QED':1})

GC_2 = Coupling(name = 'GC_2',
                value = '(2*ee*complex(0,1))/3.',
                order = {'QED':1})

GC_3 = Coupling(name = 'GC_3',
                value = '-(ee*complex(0,1))',
                order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = 'ee*complex(0,1)',
                order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = 'ee**2*complex(0,1)',
                order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '2*ee**2*complex(0,1)',
                order = {'QED':2})

GC_7 = Coupling(name = 'GC_7',
                value = '-ee**2/(2.*cw)',
                order = {'QED':2})

GC_8 = Coupling(name = 'GC_8',
                value = 'ee**2/(2.*cw)',
                order = {'QED':2})

GC_9 = Coupling(name = 'GC_9',
                value = '-(cosbma*ee**2*complex(0,1))/(2.*cw)',
                order = {'QED':2})

GC_10 = Coupling(name = 'GC_10',
                 value = '-G',
                 order = {'QCD':1})

GC_11 = Coupling(name = 'GC_11',
                 value = 'complex(0,1)*G',
                 order = {'QCD':1})

GC_12 = Coupling(name = 'GC_12',
                 value = 'complex(0,1)*G**2',
                 order = {'QCD':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '-2*complex(0,1)*l1',
                 order = {'QED':2})

GC_14 = Coupling(name = 'GC_14',
                 value = '-4*complex(0,1)*l1',
                 order = {'QED':2})

GC_15 = Coupling(name = 'GC_15',
                 value = '-6*complex(0,1)*l1',
                 order = {'QED':2})

GC_16 = Coupling(name = 'GC_16',
                 value = '-2*complex(0,1)*l2',
                 order = {'QED':2})

GC_17 = Coupling(name = 'GC_17',
                 value = '-4*complex(0,1)*l2',
                 order = {'QED':2})

GC_18 = Coupling(name = 'GC_18',
                 value = '-6*complex(0,1)*l2',
                 order = {'QED':2})

GC_19 = Coupling(name = 'GC_19',
                 value = '-(complex(0,1)*l3)',
                 order = {'QED':2})

GC_20 = Coupling(name = 'GC_20',
                 value = '-(complex(0,1)*l3) - complex(0,1)*l4',
                 order = {'QED':2})

GC_21 = Coupling(name = 'GC_21',
                 value = '-(complex(0,1)*l4)/2. - complex(0,1)*l5',
                 order = {'QED':2})

GC_22 = Coupling(name = 'GC_22',
                 value = '-(complex(0,1)*l3) - complex(0,1)*l4 - 2*complex(0,1)*l5',
                 order = {'QED':2})

GC_23 = Coupling(name = 'GC_23',
                 value = '-4*complex(0,1)*l5',
                 order = {'QED':2})

GC_24 = Coupling(name = 'GC_24',
                 value = '(cosbma*l4)/2. - cosbma*l5',
                 order = {'QED':2})

GC_25 = Coupling(name = 'GC_25',
                 value = '-(cosbma*l4)/2. + cosbma*l5',
                 order = {'QED':2})

GC_26 = Coupling(name = 'GC_26',
                 value = '-(complex(0,1)*l6)',
                 order = {'QED':2})

GC_27 = Coupling(name = 'GC_27',
                 value = '-2*complex(0,1)*l6',
                 order = {'QED':2})

GC_28 = Coupling(name = 'GC_28',
                 value = '-3*complex(0,1)*l6',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '-(complex(0,1)*l7)',
                 order = {'QED':2})

GC_30 = Coupling(name = 'GC_30',
                 value = '-2*complex(0,1)*l7',
                 order = {'QED':2})

GC_31 = Coupling(name = 'GC_31',
                 value = '-3*complex(0,1)*l7',
                 order = {'QED':2})

GC_32 = Coupling(name = 'GC_32',
                 value = '-(ee*MW)',
                 order = {'QED':1})

GC_33 = Coupling(name = 'GC_33',
                 value = 'ee*complex(0,1)*MW',
                 order = {'QED':1})

GC_34 = Coupling(name = 'GC_34',
                 value = 'ee*MW',
                 order = {'QED':1})

GC_35 = Coupling(name = 'GC_35',
                 value = '-(ee**2*complex(0,1)*sinbma)/(2.*cw)',
                 order = {'QED':2})

GC_36 = Coupling(name = 'GC_36',
                 value = '(ee**2*complex(0,1)*sinbma)/(2.*cw)',
                 order = {'QED':2})

GC_37 = Coupling(name = 'GC_37',
                 value = '(l4*sinbma)/2. - l5*sinbma',
                 order = {'QED':2})

GC_38 = Coupling(name = 'GC_38',
                 value = '-(l4*sinbma)/2. + l5*sinbma',
                 order = {'QED':2})

GC_39 = Coupling(name = 'GC_39',
                 value = '-(cosbma**2*complex(0,1)*l3) - 2*cosbma*complex(0,1)*l6*sinbma - 2*complex(0,1)*l1*sinbma**2',
                 order = {'QED':2})

GC_40 = Coupling(name = 'GC_40',
                 value = '-(cosbma**2*complex(0,1)*l3) - cosbma**2*complex(0,1)*l4 + 2*cosbma**2*complex(0,1)*l5 - 2*cosbma*complex(0,1)*l6*sinbma - 2*complex(0,1)*l1*sinbma**2',
                 order = {'QED':2})

GC_41 = Coupling(name = 'GC_41',
                 value = '-(cosbma**2*complex(0,1)*l3) + 2*cosbma*complex(0,1)*l7*sinbma - 2*complex(0,1)*l2*sinbma**2',
                 order = {'QED':2})

GC_42 = Coupling(name = 'GC_42',
                 value = '-(cosbma**2*complex(0,1)*l3) - cosbma**2*complex(0,1)*l4 + 2*cosbma**2*complex(0,1)*l5 + 2*cosbma*complex(0,1)*l7*sinbma - 2*complex(0,1)*l2*sinbma**2',
                 order = {'QED':2})

GC_43 = Coupling(name = 'GC_43',
                 value = '-2*cosbma**2*complex(0,1)*l1 + 2*cosbma*complex(0,1)*l6*sinbma - complex(0,1)*l3*sinbma**2',
                 order = {'QED':2})

GC_44 = Coupling(name = 'GC_44',
                 value = '-2*cosbma**2*complex(0,1)*l2 - 2*cosbma*complex(0,1)*l7*sinbma - complex(0,1)*l3*sinbma**2',
                 order = {'QED':2})

GC_45 = Coupling(name = 'GC_45',
                 value = '-(cosbma**2*complex(0,1)*l4)/2. - cosbma**2*complex(0,1)*l5 - cosbma*complex(0,1)*l6*sinbma + cosbma*complex(0,1)*l7*sinbma + (complex(0,1)*l4*sinbma**2)/2. + complex(0,1)*l5*sinbma**2',
                 order = {'QED':2})

GC_46 = Coupling(name = 'GC_46',
                 value = '-2*cosbma**2*complex(0,1)*l5 - cosbma*complex(0,1)*l6*sinbma + cosbma*complex(0,1)*l7*sinbma + 2*complex(0,1)*l5*sinbma**2',
                 order = {'QED':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '-2*cosbma**2*complex(0,1)*l1 + 2*cosbma*complex(0,1)*l6*sinbma - complex(0,1)*l3*sinbma**2 - complex(0,1)*l4*sinbma**2 + 2*complex(0,1)*l5*sinbma**2',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '-2*cosbma**2*complex(0,1)*l2 - 2*cosbma*complex(0,1)*l7*sinbma - complex(0,1)*l3*sinbma**2 - complex(0,1)*l4*sinbma**2 + 2*complex(0,1)*l5*sinbma**2',
                 order = {'QED':2})

GC_49 = Coupling(name = 'GC_49',
                 value = '-(cosbma**2*complex(0,1)*l7) - cosbma*complex(0,1)*l4*sinbma - 2*cosbma*complex(0,1)*l5*sinbma - complex(0,1)*l6*sinbma**2',
                 order = {'QED':2})

GC_50 = Coupling(name = 'GC_50',
                 value = '-(cosbma**2*complex(0,1)*l7) - 4*cosbma*complex(0,1)*l5*sinbma - complex(0,1)*l6*sinbma**2',
                 order = {'QED':2})

GC_51 = Coupling(name = 'GC_51',
                 value = '-(cosbma**2*complex(0,1)*l6) - 2*cosbma*complex(0,1)*l1*sinbma + cosbma*complex(0,1)*l3*sinbma + complex(0,1)*l6*sinbma**2',
                 order = {'QED':2})

GC_52 = Coupling(name = 'GC_52',
                 value = '-(cosbma**2*complex(0,1)*l6) - 2*cosbma*complex(0,1)*l1*sinbma + cosbma*complex(0,1)*l3*sinbma + cosbma*complex(0,1)*l4*sinbma - 2*cosbma*complex(0,1)*l5*sinbma + complex(0,1)*l6*sinbma**2',
                 order = {'QED':2})

GC_53 = Coupling(name = 'GC_53',
                 value = '-(cosbma**2*complex(0,1)*l6) + cosbma*complex(0,1)*l4*sinbma + 2*cosbma*complex(0,1)*l5*sinbma - complex(0,1)*l7*sinbma**2',
                 order = {'QED':2})

GC_54 = Coupling(name = 'GC_54',
                 value = '-(cosbma**2*complex(0,1)*l6) + 4*cosbma*complex(0,1)*l5*sinbma - complex(0,1)*l7*sinbma**2',
                 order = {'QED':2})

GC_55 = Coupling(name = 'GC_55',
                 value = '-(cosbma**2*complex(0,1)*l7) + 2*cosbma*complex(0,1)*l2*sinbma - cosbma*complex(0,1)*l3*sinbma + complex(0,1)*l7*sinbma**2',
                 order = {'QED':2})

GC_56 = Coupling(name = 'GC_56',
                 value = '-(cosbma**2*complex(0,1)*l7) + 2*cosbma*complex(0,1)*l2*sinbma - cosbma*complex(0,1)*l3*sinbma - cosbma*complex(0,1)*l4*sinbma + 2*cosbma*complex(0,1)*l5*sinbma + complex(0,1)*l7*sinbma**2',
                 order = {'QED':2})

GC_57 = Coupling(name = 'GC_57',
                 value = '-6*cosbma**4*complex(0,1)*l2 - 12*cosbma**3*complex(0,1)*l7*sinbma - 6*cosbma**2*complex(0,1)*l3*sinbma**2 - 6*cosbma**2*complex(0,1)*l4*sinbma**2 - 12*cosbma**2*complex(0,1)*l5*sinbma**2 - 12*cosbma*complex(0,1)*l6*sinbma**3 - 6*complex(0,1)*l1*sinbma**4',
                 order = {'QED':2})

GC_58 = Coupling(name = 'GC_58',
                 value = '-6*cosbma**4*complex(0,1)*l1 + 12*cosbma**3*complex(0,1)*l6*sinbma - 6*cosbma**2*complex(0,1)*l3*sinbma**2 - 6*cosbma**2*complex(0,1)*l4*sinbma**2 - 12*cosbma**2*complex(0,1)*l5*sinbma**2 + 12*cosbma*complex(0,1)*l7*sinbma**3 - 6*complex(0,1)*l2*sinbma**4',
                 order = {'QED':2})

GC_59 = Coupling(name = 'GC_59',
                 value = '-(cosbma**4*complex(0,1)*l3) - cosbma**4*complex(0,1)*l4 - 2*cosbma**4*complex(0,1)*l5 - 6*cosbma**3*complex(0,1)*l6*sinbma + 6*cosbma**3*complex(0,1)*l7*sinbma - 6*cosbma**2*complex(0,1)*l1*sinbma**2 - 6*cosbma**2*complex(0,1)*l2*sinbma**2 + 4*cosbma**2*complex(0,1)*l3*sinbma**2 + 4*cosbma**2*complex(0,1)*l4*sinbma**2 + 8*cosbma**2*complex(0,1)*l5*sinbma**2 + 6*cosbma*complex(0,1)*l6*sinbma**3 - 6*cosbma*complex(0,1)*l7*sinbma**3 - complex(0,1)*l3*sinbma**4 - complex(0,1)*l4*sinbma**4 - 2*complex(0,1)*l5*sinbma**4',
                 order = {'QED':2})

GC_60 = Coupling(name = 'GC_60',
                 value = '-3*cosbma**4*complex(0,1)*l7 + 6*cosbma**3*complex(0,1)*l2*sinbma - 3*cosbma**3*complex(0,1)*l3*sinbma - 3*cosbma**3*complex(0,1)*l4*sinbma - 6*cosbma**3*complex(0,1)*l5*sinbma - 9*cosbma**2*complex(0,1)*l6*sinbma**2 + 9*cosbma**2*complex(0,1)*l7*sinbma**2 - 6*cosbma*complex(0,1)*l1*sinbma**3 + 3*cosbma*complex(0,1)*l3*sinbma**3 + 3*cosbma*complex(0,1)*l4*sinbma**3 + 6*cosbma*complex(0,1)*l5*sinbma**3 + 3*complex(0,1)*l6*sinbma**4',
                 order = {'QED':2})

GC_61 = Coupling(name = 'GC_61',
                 value = '-3*cosbma**4*complex(0,1)*l6 - 6*cosbma**3*complex(0,1)*l1*sinbma + 3*cosbma**3*complex(0,1)*l3*sinbma + 3*cosbma**3*complex(0,1)*l4*sinbma + 6*cosbma**3*complex(0,1)*l5*sinbma + 9*cosbma**2*complex(0,1)*l6*sinbma**2 - 9*cosbma**2*complex(0,1)*l7*sinbma**2 + 6*cosbma*complex(0,1)*l2*sinbma**3 - 3*cosbma*complex(0,1)*l3*sinbma**3 - 3*cosbma*complex(0,1)*l4*sinbma**3 - 6*cosbma*complex(0,1)*l5*sinbma**3 + 3*complex(0,1)*l7*sinbma**4',
                 order = {'QED':2})

GC_62 = Coupling(name = 'GC_62',
                 value = '(cosbma**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sinbma**2)/(2.*sw**2)',
                 order = {'QED':2})

GC_63 = Coupling(name = 'GC_63',
                 value = '(ee**2*complex(0,1))/(2.*sw**2)',
                 order = {'QED':2})

GC_64 = Coupling(name = 'GC_64',
                 value = '-((ee**2*complex(0,1))/sw**2)',
                 order = {'QED':2})

GC_65 = Coupling(name = 'GC_65',
                 value = '(cw**2*ee**2*complex(0,1))/sw**2',
                 order = {'QED':2})

GC_66 = Coupling(name = 'GC_66',
                 value = 'ee/(2.*sw)',
                 order = {'QED':1})

GC_67 = Coupling(name = 'GC_67',
                 value = '(ee*complex(0,1))/(sw*cmath.sqrt(2))',
                 order = {'QED':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '-(cosbma*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '(cosbma*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '-(cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_71 = Coupling(name = 'GC_71',
                 value = '(cw*ee*complex(0,1))/(2.*sw)',
                 order = {'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = '-((cw*ee*complex(0,1))/sw)',
                 order = {'QED':1})

GC_73 = Coupling(name = 'GC_73',
                 value = '(cw*ee*complex(0,1))/sw',
                 order = {'QED':1})

GC_74 = Coupling(name = 'GC_74',
                 value = '-ee**2/(2.*sw)',
                 order = {'QED':2})

GC_75 = Coupling(name = 'GC_75',
                 value = 'ee**2/(2.*sw)',
                 order = {'QED':2})

GC_76 = Coupling(name = 'GC_76',
                 value = '(cosbma*ee**2*complex(0,1))/(2.*sw)',
                 order = {'QED':2})

GC_77 = Coupling(name = 'GC_77',
                 value = '(-2*cw*ee**2*complex(0,1))/sw',
                 order = {'QED':2})

GC_78 = Coupling(name = 'GC_78',
                 value = '-(ee*MW)/(2.*sw)',
                 order = {'QED':1})

GC_79 = Coupling(name = 'GC_79',
                 value = '(ee*MW)/(2.*sw)',
                 order = {'QED':1})

GC_80 = Coupling(name = 'GC_80',
                 value = '-(cosbma*ee*complex(0,1)*MW)/(2.*sw)',
                 order = {'QED':1})

GC_81 = Coupling(name = 'GC_81',
                 value = '(cosbma*ee*complex(0,1)*MW)/sw',
                 order = {'QED':1})

GC_82 = Coupling(name = 'GC_82',
                 value = '-(ee*complex(0,1)*sinbma)/(2.*sw)',
                 order = {'QED':1})

GC_83 = Coupling(name = 'GC_83',
                 value = '(ee*complex(0,1)*sinbma)/(2.*sw)',
                 order = {'QED':1})

GC_84 = Coupling(name = 'GC_84',
                 value = '-(ee**2*complex(0,1)*sinbma)/(2.*sw)',
                 order = {'QED':2})

GC_85 = Coupling(name = 'GC_85',
                 value = '(ee**2*complex(0,1)*sinbma)/(2.*sw)',
                 order = {'QED':2})

GC_86 = Coupling(name = 'GC_86',
                 value = '-(ee*complex(0,1)*MW*sinbma)/(2.*sw)',
                 order = {'QED':1})

GC_87 = Coupling(name = 'GC_87',
                 value = '(ee*complex(0,1)*MW*sinbma)/sw',
                 order = {'QED':1})

GC_88 = Coupling(name = 'GC_88',
                 value = '-(ee*complex(0,1)*sw)/(6.*cw)',
                 order = {'QED':1})

GC_89 = Coupling(name = 'GC_89',
                 value = '(ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_90 = Coupling(name = 'GC_90',
                 value = '-((ee*complex(0,1)*MW*sw)/cw)',
                 order = {'QED':1})

GC_91 = Coupling(name = 'GC_91',
                 value = '-(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_92 = Coupling(name = 'GC_92',
                 value = '(cw*ee*complex(0,1))/(2.*sw) + (ee*complex(0,1)*sw)/(2.*cw)',
                 order = {'QED':1})

GC_93 = Coupling(name = 'GC_93',
                 value = '-(cosbma*cw*ee)/(2.*sw) - (cosbma*ee*sw)/(2.*cw)',
                 order = {'QED':1})

GC_94 = Coupling(name = 'GC_94',
                 value = '(cosbma*cw*ee)/(2.*sw) + (cosbma*ee*sw)/(2.*cw)',
                 order = {'QED':1})

GC_95 = Coupling(name = 'GC_95',
                 value = '(cw*ee**2*complex(0,1))/sw - (ee**2*complex(0,1)*sw)/cw',
                 order = {'QED':2})

GC_96 = Coupling(name = 'GC_96',
                 value = '-(cw*ee*MW)/(2.*sw) - (ee*MW*sw)/(2.*cw)',
                 order = {'QED':1})

GC_97 = Coupling(name = 'GC_97',
                 value = '(cw*ee*MW)/(2.*sw) - (ee*MW*sw)/(2.*cw)',
                 order = {'QED':1})

GC_98 = Coupling(name = 'GC_98',
                 value = '-(cw*ee*MW)/(2.*sw) + (ee*MW*sw)/(2.*cw)',
                 order = {'QED':1})

GC_99 = Coupling(name = 'GC_99',
                 value = '(cw*ee*MW)/(2.*sw) + (ee*MW*sw)/(2.*cw)',
                 order = {'QED':1})

GC_100 = Coupling(name = 'GC_100',
                  value = '(l4*MW*sw)/ee - (2*l5*MW*sw)/ee',
                  order = {'QED':1})

GC_101 = Coupling(name = 'GC_101',
                  value = '-((l4*MW*sw)/ee) + (2*l5*MW*sw)/ee',
                  order = {'QED':1})

GC_102 = Coupling(name = 'GC_102',
                  value = '-(cw*ee*sinbma)/(2.*sw) - (ee*sinbma*sw)/(2.*cw)',
                  order = {'QED':1})

GC_103 = Coupling(name = 'GC_103',
                  value = '(-2*cosbma*complex(0,1)*l6*MW*sw)/ee - (4*complex(0,1)*l1*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_104 = Coupling(name = 'GC_104',
                  value = '(-2*cosbma*complex(0,1)*l7*MW*sw)/ee - (2*complex(0,1)*l3*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_105 = Coupling(name = 'GC_105',
                  value = '(-2*cosbma*complex(0,1)*l6*MW*sw)/ee + (complex(0,1)*l4*MW*sinbma*sw)/ee + (2*complex(0,1)*l5*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_106 = Coupling(name = 'GC_106',
                  value = '(-2*cosbma*complex(0,1)*l6*MW*sw)/ee + (4*complex(0,1)*l5*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_107 = Coupling(name = 'GC_107',
                  value = '(-2*cosbma*complex(0,1)*l7*MW*sw)/ee - (2*complex(0,1)*l3*MW*sinbma*sw)/ee - (2*complex(0,1)*l4*MW*sinbma*sw)/ee + (4*complex(0,1)*l5*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_108 = Coupling(name = 'GC_108',
                  value = '-((cosbma*complex(0,1)*l4*MW*sw)/ee) - (2*cosbma*complex(0,1)*l5*MW*sw)/ee - (2*complex(0,1)*l6*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_109 = Coupling(name = 'GC_109',
                  value = '(-4*cosbma*complex(0,1)*l5*MW*sw)/ee - (2*complex(0,1)*l6*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_110 = Coupling(name = 'GC_110',
                  value = '(-4*cosbma*complex(0,1)*l1*MW*sw)/ee + (2*complex(0,1)*l6*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_111 = Coupling(name = 'GC_111',
                  value = '(-2*cosbma*complex(0,1)*l3*MW*sw)/ee + (2*complex(0,1)*l7*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_112 = Coupling(name = 'GC_112',
                  value = '(-2*cosbma*complex(0,1)*l3*MW*sw)/ee - (2*cosbma*complex(0,1)*l4*MW*sw)/ee + (4*cosbma*complex(0,1)*l5*MW*sw)/ee + (2*complex(0,1)*l7*MW*sinbma*sw)/ee',
                  order = {'QED':1})

GC_113 = Coupling(name = 'GC_113',
                  value = '(-6*cosbma**3*complex(0,1)*l7*MW*sw)/ee - (6*cosbma**2*complex(0,1)*l3*MW*sinbma*sw)/ee - (6*cosbma**2*complex(0,1)*l4*MW*sinbma*sw)/ee - (12*cosbma**2*complex(0,1)*l5*MW*sinbma*sw)/ee - (18*cosbma*complex(0,1)*l6*MW*sinbma**2*sw)/ee - (12*complex(0,1)*l1*MW*sinbma**3*sw)/ee',
                  order = {'QED':1})

GC_114 = Coupling(name = 'GC_114',
                  value = '(-6*cosbma**3*complex(0,1)*l6*MW*sw)/ee - (12*cosbma**2*complex(0,1)*l1*MW*sinbma*sw)/ee + (4*cosbma**2*complex(0,1)*l3*MW*sinbma*sw)/ee + (4*cosbma**2*complex(0,1)*l4*MW*sinbma*sw)/ee + (8*cosbma**2*complex(0,1)*l5*MW*sinbma*sw)/ee + (12*cosbma*complex(0,1)*l6*MW*sinbma**2*sw)/ee - (6*cosbma*complex(0,1)*l7*MW*sinbma**2*sw)/ee - (2*complex(0,1)*l3*MW*sinbma**3*sw)/ee - (2*complex(0,1)*l4*MW*sinbma**3*sw)/ee - (4*complex(0,1)*l5*MW*sinbma**3*sw)/ee',
                  order = {'QED':1})

GC_115 = Coupling(name = 'GC_115',
                  value = '(-2*cosbma**3*complex(0,1)*l3*MW*sw)/ee - (2*cosbma**3*complex(0,1)*l4*MW*sw)/ee - (4*cosbma**3*complex(0,1)*l5*MW*sw)/ee - (12*cosbma**2*complex(0,1)*l6*MW*sinbma*sw)/ee + (6*cosbma**2*complex(0,1)*l7*MW*sinbma*sw)/ee - (12*cosbma*complex(0,1)*l1*MW*sinbma**2*sw)/ee + (4*cosbma*complex(0,1)*l3*MW*sinbma**2*sw)/ee + (4*cosbma*complex(0,1)*l4*MW*sinbma**2*sw)/ee + (8*cosbma*complex(0,1)*l5*MW*sinbma**2*sw)/ee + (6*complex(0,1)*l6*MW*sinbma**3*sw)/ee',
                  order = {'QED':1})

GC_116 = Coupling(name = 'GC_116',
                  value = '(-12*cosbma**3*complex(0,1)*l1*MW*sw)/ee + (18*cosbma**2*complex(0,1)*l6*MW*sinbma*sw)/ee - (6*cosbma*complex(0,1)*l3*MW*sinbma**2*sw)/ee - (6*cosbma*complex(0,1)*l4*MW*sinbma**2*sw)/ee - (12*cosbma*complex(0,1)*l5*MW*sinbma**2*sw)/ee + (6*complex(0,1)*l7*MW*sinbma**3*sw)/ee',
                  order = {'QED':1})

GC_117 = Coupling(name = 'GC_117',
                  value = '-(ee**2*complex(0,1)) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                  order = {'QED':2})

GC_118 = Coupling(name = 'GC_118',
                  value = 'ee**2*complex(0,1) + (cw**2*ee**2*complex(0,1))/(2.*sw**2) + (ee**2*complex(0,1)*sw**2)/(2.*cw**2)',
                  order = {'QED':2})

GC_119 = Coupling(name = 'GC_119',
                  value = 'cosbma**2*ee**2*complex(0,1) + ee**2*complex(0,1)*sinbma**2 + (cosbma**2*cw**2*ee**2*complex(0,1))/(2.*sw**2) + (cw**2*ee**2*complex(0,1)*sinbma**2)/(2.*sw**2) + (cosbma**2*ee**2*complex(0,1)*sw**2)/(2.*cw**2) + (ee**2*complex(0,1)*sinbma**2*sw**2)/(2.*cw**2)',
                  order = {'QED':2})

GC_120 = Coupling(name = 'GC_120',
                  value = '-(cosbma*cw**2*ee*complex(0,1)*MW)/(2.*sw) - cosbma*ee*complex(0,1)*MW*sw - (cosbma*ee*complex(0,1)*MW*sw**3)/(2.*cw**2)',
                  order = {'QED':1})

GC_121 = Coupling(name = 'GC_121',
                  value = '(cosbma*cw**2*ee*complex(0,1)*MW)/sw + 2*cosbma*ee*complex(0,1)*MW*sw + (cosbma*ee*complex(0,1)*MW*sw**3)/cw**2',
                  order = {'QED':1})

GC_122 = Coupling(name = 'GC_122',
                  value = '-(cw**2*ee*complex(0,1)*MW*sinbma)/(2.*sw) - ee*complex(0,1)*MW*sinbma*sw - (ee*complex(0,1)*MW*sinbma*sw**3)/(2.*cw**2)',
                  order = {'QED':1})

GC_123 = Coupling(name = 'GC_123',
                  value = '(cw**2*ee*complex(0,1)*MW*sinbma)/sw + 2*ee*complex(0,1)*MW*sinbma*sw + (ee*complex(0,1)*MW*sinbma*sw**3)/cw**2',
                  order = {'QED':1})

GC_124 = Coupling(name = 'GC_124',
                  value = '-(complex(0,1)*yb)',
                  order = {'QED':1})

GC_125 = Coupling(name = 'GC_125',
                  value = '-(yb/cmath.sqrt(2))',
                  order = {'QED':1})

GC_126 = Coupling(name = 'GC_126',
                  value = '-((cosbma*complex(0,1)*yb)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_127 = Coupling(name = 'GC_127',
                  value = '-((complex(0,1)*sinbma*yb)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_128 = Coupling(name = 'GC_128',
                  value = 'complex(0,1)*yc',
                  order = {'QED':1})

GC_129 = Coupling(name = 'GC_129',
                  value = 'yc/cmath.sqrt(2)',
                  order = {'QED':1})

GC_130 = Coupling(name = 'GC_130',
                  value = '-(complex(0,1)*ydo)',
                  order = {'QED':1})

GC_131 = Coupling(name = 'GC_131',
                  value = '-(ydo/cmath.sqrt(2))',
                  order = {'QED':1})

GC_132 = Coupling(name = 'GC_132',
                  value = '-(complex(0,1)*ye)',
                  order = {'QED':1})

GC_133 = Coupling(name = 'GC_133',
                  value = '-(ye/cmath.sqrt(2))',
                  order = {'QED':1})

GC_134 = Coupling(name = 'GC_134',
                  value = '-(complex(0,1)*ym)',
                  order = {'QED':1})

GC_135 = Coupling(name = 'GC_135',
                  value = '-(ym/cmath.sqrt(2))',
                  order = {'QED':1})

GC_136 = Coupling(name = 'GC_136',
                  value = '(ee*ymb)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YB':1})

GC_137 = Coupling(name = 'GC_137',
                  value = '-((ee*complex(0,1)*ymb)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1,'YB':1})

GC_138 = Coupling(name = 'GC_138',
                  value = '-(cosbma*ee*complex(0,1)*ymb)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YB':1})

GC_139 = Coupling(name = 'GC_139',
                  value = '(ee*complex(0,1)*sinbma*ymb)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YB':1})

GC_140 = Coupling(name = 'GC_140',
                  value = '-(ee*ymc)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_141 = Coupling(name = 'GC_141',
                  value = '(ee*complex(0,1)*ymc)/(MW*sw*TB*cmath.sqrt(2))',
                  order = {'QED':1})

GC_142 = Coupling(name = 'GC_142',
                  value = '-(cosbma*ee*complex(0,1)*ymc)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*yc)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_143 = Coupling(name = 'GC_143',
                  value = '(ee*complex(0,1)*sinbma*ymc)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*yc)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_144 = Coupling(name = 'GC_144',
                  value = '(ee*ymdo)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_145 = Coupling(name = 'GC_145',
                  value = '-((ee*complex(0,1)*ymdo)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_146 = Coupling(name = 'GC_146',
                  value = '-(cosbma*ee*complex(0,1)*ymdo)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*ydo)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_147 = Coupling(name = 'GC_147',
                  value = '(ee*complex(0,1)*sinbma*ymdo)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*ydo)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_148 = Coupling(name = 'GC_148',
                  value = '(ee*yme)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_149 = Coupling(name = 'GC_149',
                  value = '-((ee*complex(0,1)*yme)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_150 = Coupling(name = 'GC_150',
                  value = '-(cosbma*ee*complex(0,1)*yme)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*ye)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_151 = Coupling(name = 'GC_151',
                  value = '(ee*complex(0,1)*sinbma*yme)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*ye)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_152 = Coupling(name = 'GC_152',
                  value = '(ee*ymm)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_153 = Coupling(name = 'GC_153',
                  value = '-((ee*complex(0,1)*ymm)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_154 = Coupling(name = 'GC_154',
                  value = '-(cosbma*ee*complex(0,1)*ymm)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*ym)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_155 = Coupling(name = 'GC_155',
                  value = '(ee*complex(0,1)*sinbma*ymm)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*ym)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_156 = Coupling(name = 'GC_156',
                  value = '(ee*yms)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_157 = Coupling(name = 'GC_157',
                  value = '-((ee*complex(0,1)*yms)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_158 = Coupling(name = 'GC_158',
                  value = '-(ee*ymt)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YT':1})

GC_159 = Coupling(name = 'GC_159',
                  value = '(ee*complex(0,1)*ymt)/(MW*sw*TB*cmath.sqrt(2))',
                  order = {'QED':1,'YT':1})

GC_160 = Coupling(name = 'GC_160',
                  value = '-(cosbma*ee*complex(0,1)*ymt)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YT':1})

GC_161 = Coupling(name = 'GC_161',
                  value = '(ee*complex(0,1)*sinbma*ymt)/(2.*MW*sw*TB)',
                  order = {'QED':1,'YT':1})

GC_162 = Coupling(name = 'GC_162',
                  value = '(ee*ymtau)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_163 = Coupling(name = 'GC_163',
                  value = '-((ee*complex(0,1)*ymtau)/(MW*sw*TB*cmath.sqrt(2)))',
                  order = {'QED':1})

GC_164 = Coupling(name = 'GC_164',
                  value = '-(ee*ymup)/(2.*MW*sw*TB)',
                  order = {'QED':1})

GC_165 = Coupling(name = 'GC_165',
                  value = '(ee*complex(0,1)*ymup)/(MW*sw*TB*cmath.sqrt(2))',
                  order = {'QED':1})

GC_166 = Coupling(name = 'GC_166',
                  value = '-(complex(0,1)*ys)',
                  order = {'QED':1})

GC_167 = Coupling(name = 'GC_167',
                  value = '-(ys/cmath.sqrt(2))',
                  order = {'QED':1})

GC_168 = Coupling(name = 'GC_168',
                  value = '(ee*complex(0,1)*sinbma*yms)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*ys)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_169 = Coupling(name = 'GC_169',
                  value = '-(cosbma*ee*complex(0,1)*yms)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*ys)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_170 = Coupling(name = 'GC_170',
                  value = 'complex(0,1)*yt',
                  order = {'QED':1})

GC_171 = Coupling(name = 'GC_171',
                  value = 'yt/cmath.sqrt(2)',
                  order = {'QED':1})

GC_172 = Coupling(name = 'GC_172',
                  value = '-((cosbma*complex(0,1)*yt)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_173 = Coupling(name = 'GC_173',
                  value = '-((complex(0,1)*sinbma*yt)/cmath.sqrt(2))',
                  order = {'QED':1})

GC_174 = Coupling(name = 'GC_174',
                  value = '-(complex(0,1)*ytau)',
                  order = {'QED':1})

GC_175 = Coupling(name = 'GC_175',
                  value = '-(ytau/cmath.sqrt(2))',
                  order = {'QED':1})

GC_176 = Coupling(name = 'GC_176',
                  value = '(ee*complex(0,1)*sinbma*ymtau)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*ytau)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_177 = Coupling(name = 'GC_177',
                  value = '-(cosbma*ee*complex(0,1)*ymtau)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*ytau)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_178 = Coupling(name = 'GC_178',
                  value = 'complex(0,1)*yup',
                  order = {'QED':1})

GC_179 = Coupling(name = 'GC_179',
                  value = 'yup/cmath.sqrt(2)',
                  order = {'QED':1})

GC_180 = Coupling(name = 'GC_180',
                  value = '(ee*complex(0,1)*sinbma*ymup)/(2.*MW*sw*TB) - (cosbma*complex(0,1)*yup)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_181 = Coupling(name = 'GC_181',
                  value = '-(cosbma*ee*complex(0,1)*ymup)/(2.*MW*sw*TB) - (complex(0,1)*sinbma*yup)/cmath.sqrt(2)',
                  order = {'QED':1})

