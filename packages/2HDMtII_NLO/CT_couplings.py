# This file was automatically created by FeynRules 2.3.13
# Mathematica version: 10.1.0  for Mac OS X x86 (64-bit) (March 24, 2015)
# Date: Tue 8 Dec 2015 11:44:25


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



R2GC_212_1 = Coupling(name = 'R2GC_212_1',
                      value = '-(cw*complex(0,1)*G**2*g1)/(9.*cmath.pi**2)',
                      order = {'QCD':2,'QED':1})

R2GC_213_2 = Coupling(name = 'R2GC_213_2',
                      value = '(complex(0,1)*G**2*g1*sw)/(9.*cmath.pi**2)',
                      order = {'QCD':2,'QED':1})

R2GC_214_3 = Coupling(name = 'R2GC_214_3',
                      value = '(cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2)',
                      order = {'QCD':2,'QED':1})

R2GC_215_4 = Coupling(name = 'R2GC_215_4',
                      value = '-(complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2)',
                      order = {'QCD':2,'QED':1})

R2GC_220_5 = Coupling(name = 'R2GC_220_5',
                      value = '-(complex(0,1)*G**2*MB**2)/(8.*cmath.pi**2)',
                      order = {'QCD':2})

R2GC_220_6 = Coupling(name = 'R2GC_220_6',
                      value = '-(complex(0,1)*G**2*MT**2)/(8.*cmath.pi**2)',
                      order = {'QCD':2})

R2GC_221_7 = Coupling(name = 'R2GC_221_7',
                      value = '-(complex(0,1)*G**2*tanbeta**2*TH3x3**2*ymb**2)/(8.*cmath.pi**2*vev**2)',
                      order = {'QCD':2,'QED':2})

R2GC_221_8 = Coupling(name = 'R2GC_221_8',
                      value = '-(complex(0,1)*G**2*TH3x3**2*ymt**2)/(8.*cmath.pi**2*tanbeta**2*vev**2)',
                      order = {'QCD':2,'QED':2})

R2GC_222_9 = Coupling(name = 'R2GC_222_9',
                      value = '-(complex(0,1)*G**2*MB*tanbeta*TH2x1*ymb)/(8.*cmath.pi**2*vev) - (complex(0,1)*G**2*MB*TH1x1*yb)/(8.*cmath.pi**2*cmath.sqrt(2))',
                      order = {'QCD':2,'QED':1})

R2GC_222_10 = Coupling(name = 'R2GC_222_10',
                       value = '(complex(0,1)*G**2*MT*TH2x1*ymt)/(8.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*MT*TH1x1*yt)/(8.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_223_11 = Coupling(name = 'R2GC_223_11',
                       value = '-(complex(0,1)*G**2*MB*tanbeta*TH2x2*ymb)/(8.*cmath.pi**2*vev) - (complex(0,1)*G**2*MB*TH1x2*yb)/(8.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_223_12 = Coupling(name = 'R2GC_223_12',
                       value = '(complex(0,1)*G**2*MT*TH2x2*ymt)/(8.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*MT*TH1x2*yt)/(8.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_224_13 = Coupling(name = 'R2GC_224_13',
                       value = '-(complex(0,1)*G**2*tanbeta*TH3x3*yb*ymb)/(8.*cmath.pi**2*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_224_14 = Coupling(name = 'R2GC_224_14',
                       value = '(complex(0,1)*G**2*TH3x3*ymt*yt)/(8.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_225_15 = Coupling(name = 'R2GC_225_15',
                       value = '-(complex(0,1)*G**2*yb**2)/(16.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_225_16 = Coupling(name = 'R2GC_225_16',
                       value = '-(complex(0,1)*G**2*yt**2)/(16.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_226_17 = Coupling(name = 'R2GC_226_17',
                       value = '-(complex(0,1)*G**2*TH1x1**2*yb**2)/(16.*cmath.pi**2) - (complex(0,1)*G**2*tanbeta**2*TH2x1**2*ymb**2)/(8.*cmath.pi**2*vev**2) - (complex(0,1)*G**2*tanbeta*TH1x1*TH2x1*yb*ymb)/(4.*cmath.pi**2*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_226_18 = Coupling(name = 'R2GC_226_18',
                       value = '-(complex(0,1)*G**2*TH2x1**2*ymt**2)/(8.*cmath.pi**2*tanbeta**2*vev**2) - (complex(0,1)*G**2*TH1x1**2*yt**2)/(16.*cmath.pi**2) + (complex(0,1)*G**2*TH1x1*TH2x1*ymt*yt)/(4.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_227_19 = Coupling(name = 'R2GC_227_19',
                       value = '-(complex(0,1)*G**2*TH1x1*TH1x2*yb**2)/(16.*cmath.pi**2) - (complex(0,1)*G**2*tanbeta**2*TH2x1*TH2x2*ymb**2)/(8.*cmath.pi**2*vev**2) - (complex(0,1)*G**2*tanbeta*TH1x2*TH2x1*yb*ymb)/(8.*cmath.pi**2*vev*cmath.sqrt(2)) - (complex(0,1)*G**2*tanbeta*TH1x1*TH2x2*yb*ymb)/(8.*cmath.pi**2*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_227_20 = Coupling(name = 'R2GC_227_20',
                       value = '-(complex(0,1)*G**2*TH2x1*TH2x2*ymt**2)/(8.*cmath.pi**2*tanbeta**2*vev**2) - (complex(0,1)*G**2*TH1x1*TH1x2*yt**2)/(16.*cmath.pi**2) + (complex(0,1)*G**2*TH1x2*TH2x1*ymt*yt)/(8.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) + (complex(0,1)*G**2*TH1x1*TH2x2*ymt*yt)/(8.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_228_21 = Coupling(name = 'R2GC_228_21',
                       value = '-(complex(0,1)*G**2*TH1x2**2*yb**2)/(16.*cmath.pi**2) - (complex(0,1)*G**2*tanbeta**2*TH2x2**2*ymb**2)/(8.*cmath.pi**2*vev**2) - (complex(0,1)*G**2*tanbeta*TH1x2*TH2x2*yb*ymb)/(4.*cmath.pi**2*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_228_22 = Coupling(name = 'R2GC_228_22',
                       value = '-(complex(0,1)*G**2*TH2x2**2*ymt**2)/(8.*cmath.pi**2*tanbeta**2*vev**2) - (complex(0,1)*G**2*TH1x2**2*yt**2)/(16.*cmath.pi**2) + (complex(0,1)*G**2*TH1x2*TH2x2*ymt*yt)/(4.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_229_23 = Coupling(name = 'R2GC_229_23',
                       value = '-(cw*complex(0,1)*G**3*gw)/(192.*cmath.pi**2) + (complex(0,1)*G**3*g1*sw)/(576.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_229_24 = Coupling(name = 'R2GC_229_24',
                       value = '(cw*complex(0,1)*G**3*gw)/(192.*cmath.pi**2) - (5*complex(0,1)*G**3*g1*sw)/(576.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_230_25 = Coupling(name = 'R2GC_230_25',
                       value = '(-3*cw*complex(0,1)*G**3*gw)/(64.*cmath.pi**2) - (3*complex(0,1)*G**3*g1*sw)/(64.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_230_26 = Coupling(name = 'R2GC_230_26',
                       value = '(3*cw*complex(0,1)*G**3*gw)/(64.*cmath.pi**2) + (3*complex(0,1)*G**3*g1*sw)/(64.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_231_27 = Coupling(name = 'R2GC_231_27',
                       value = '-(cw*complex(0,1)*G**3*g1)/(576.*cmath.pi**2) - (complex(0,1)*G**3*gw*sw)/(192.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_231_28 = Coupling(name = 'R2GC_231_28',
                       value = '(5*cw*complex(0,1)*G**3*g1)/(576.*cmath.pi**2) + (complex(0,1)*G**3*gw*sw)/(192.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_232_29 = Coupling(name = 'R2GC_232_29',
                       value = '(-3*cw*complex(0,1)*G**3*g1)/(64.*cmath.pi**2) + (3*complex(0,1)*G**3*gw*sw)/(64.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_232_30 = Coupling(name = 'R2GC_232_30',
                       value = '(3*cw*complex(0,1)*G**3*g1)/(64.*cmath.pi**2) - (3*complex(0,1)*G**3*gw*sw)/(64.*cmath.pi**2)',
                       order = {'QCD':3,'QED':1})

R2GC_233_31 = Coupling(name = 'R2GC_233_31',
                       value = '(cw**2*complex(0,1)*G**2*gw**2)/(192.*cmath.pi**2) + (cw*complex(0,1)*G**2*g1*gw*sw)/(288.*cmath.pi**2) + (5*complex(0,1)*G**2*g1**2*sw**2)/(1728.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_233_32 = Coupling(name = 'R2GC_233_32',
                       value = '(cw**2*complex(0,1)*G**2*gw**2)/(192.*cmath.pi**2) - (cw*complex(0,1)*G**2*g1*gw*sw)/(288.*cmath.pi**2) + (17*complex(0,1)*G**2*g1**2*sw**2)/(1728.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_234_33 = Coupling(name = 'R2GC_234_33',
                       value = '-(cw**2*complex(0,1)*G**2*g1*gw)/(576.*cmath.pi**2) - (5*cw*complex(0,1)*G**2*g1**2*sw)/(1728.*cmath.pi**2) + (cw*complex(0,1)*G**2*gw**2*sw)/(192.*cmath.pi**2) + (complex(0,1)*G**2*g1*gw*sw**2)/(576.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_234_34 = Coupling(name = 'R2GC_234_34',
                       value = '(cw**2*complex(0,1)*G**2*g1*gw)/(576.*cmath.pi**2) - (17*cw*complex(0,1)*G**2*g1**2*sw)/(1728.*cmath.pi**2) + (cw*complex(0,1)*G**2*gw**2*sw)/(192.*cmath.pi**2) - (complex(0,1)*G**2*g1*gw*sw**2)/(576.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_235_35 = Coupling(name = 'R2GC_235_35',
                       value = '(5*cw**2*complex(0,1)*G**2*g1**2)/(1728.*cmath.pi**2) - (cw*complex(0,1)*G**2*g1*gw*sw)/(288.*cmath.pi**2) + (complex(0,1)*G**2*gw**2*sw**2)/(192.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_235_36 = Coupling(name = 'R2GC_235_36',
                       value = '(17*cw**2*complex(0,1)*G**2*g1**2)/(1728.*cmath.pi**2) + (cw*complex(0,1)*G**2*g1*gw*sw)/(288.*cmath.pi**2) + (complex(0,1)*G**2*gw**2*sw**2)/(192.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_236_37 = Coupling(name = 'R2GC_236_37',
                       value = '-(cw*complex(0,1)*G**2*gw)/(12.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_237_38 = Coupling(name = 'R2GC_237_38',
                       value = '-(cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(12.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_238_39 = Coupling(name = 'R2GC_238_39',
                       value = '(cw*complex(0,1)*G**2*gw)/(12.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_239_40 = Coupling(name = 'R2GC_239_40',
                       value = '-(cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(12.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_242_41 = Coupling(name = 'R2GC_242_41',
                       value = '-(complex(0,1)*G**2*tanbeta**2*ymb**2)/(8.*cmath.pi**2*vev**2) - (complex(0,1)*G**2*ymt**2)/(8.*cmath.pi**2*tanbeta**2*vev**2)',
                       order = {'QCD':2,'QED':2})

R2GC_243_42 = Coupling(name = 'R2GC_243_42',
                       value = '-(complex(0,1)*G**2*tanbeta*yb*ymb)/(8.*cmath.pi**2*vev*cmath.sqrt(2)) + (complex(0,1)*G**2*ymt*yt)/(8.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':2})

R2GC_244_43 = Coupling(name = 'R2GC_244_43',
                       value = '-(complex(0,1)*G**2*yb**2)/(16.*cmath.pi**2) - (complex(0,1)*G**2*yt**2)/(16.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_245_44 = Coupling(name = 'R2GC_245_44',
                       value = '(complex(0,1)*G**2*gw**2)/(96.*cmath.pi**2)',
                       order = {'QCD':2,'QED':2})

R2GC_254_45 = Coupling(name = 'R2GC_254_45',
                       value = '-G**4/(192.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_254_46 = Coupling(name = 'R2GC_254_46',
                       value = 'G**4/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_255_47 = Coupling(name = 'R2GC_255_47',
                       value = '-(complex(0,1)*G**4)/(192.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_255_48 = Coupling(name = 'R2GC_255_48',
                       value = '(complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_256_49 = Coupling(name = 'R2GC_256_49',
                       value = '(complex(0,1)*G**4)/(192.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_256_50 = Coupling(name = 'R2GC_256_50',
                       value = '-(complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_257_51 = Coupling(name = 'R2GC_257_51',
                       value = '-(complex(0,1)*G**4)/(48.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_258_52 = Coupling(name = 'R2GC_258_52',
                       value = '(complex(0,1)*G**4)/(288.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_258_53 = Coupling(name = 'R2GC_258_53',
                       value = '-(complex(0,1)*G**4)/(32.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_259_54 = Coupling(name = 'R2GC_259_54',
                       value = '-(complex(0,1)*G**4)/(16.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_259_55 = Coupling(name = 'R2GC_259_55',
                       value = '(complex(0,1)*G**4)/(4.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_260_56 = Coupling(name = 'R2GC_260_56',
                       value = '(-3*complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_260_57 = Coupling(name = 'R2GC_260_57',
                       value = '(-23*complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_261_58 = Coupling(name = 'R2GC_261_58',
                       value = '-(complex(0,1)*G**3)/(6.*cmath.pi**2)',
                       order = {'QCD':3})

R2GC_262_59 = Coupling(name = 'R2GC_262_59',
                       value = '(complex(0,1)*G**2)/(12.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_275_60 = Coupling(name = 'R2GC_275_60',
                       value = '-(complex(0,1)*G**2*gw)/(6.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_280_61 = Coupling(name = 'R2GC_280_61',
                       value = '(complex(0,1)*G**2*MB)/(6.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_284_62 = Coupling(name = 'R2GC_284_62',
                       value = '-(G**2*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_285_63 = Coupling(name = 'R2GC_285_63',
                       value = '(complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(3.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x1*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_286_64 = Coupling(name = 'R2GC_286_64',
                       value = '(complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(3.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x2*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_287_65 = Coupling(name = 'R2GC_287_65',
                       value = '-(G**2*tanbeta*TH3x3*ymb)/(3.*cmath.pi**2*vev)',
                       order = {'QCD':2,'QED':1})

R2GC_288_66 = Coupling(name = 'R2GC_288_66',
                       value = '(complex(0,1)*G**2)/(48.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_288_67 = Coupling(name = 'R2GC_288_67',
                       value = '(3*complex(0,1)*G**2)/(32.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_289_68 = Coupling(name = 'R2GC_289_68',
                       value = '-(complex(0,1)*G**2)/(16.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_290_69 = Coupling(name = 'R2GC_290_69',
                       value = 'G**3/(24.*cmath.pi**2)',
                       order = {'QCD':3})

R2GC_290_70 = Coupling(name = 'R2GC_290_70',
                       value = '(11*G**3)/(64.*cmath.pi**2)',
                       order = {'QCD':3})

R2GC_291_71 = Coupling(name = 'R2GC_291_71',
                       value = '(5*complex(0,1)*G**4)/(48.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_291_72 = Coupling(name = 'R2GC_291_72',
                       value = '(19*complex(0,1)*G**4)/(32.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_292_73 = Coupling(name = 'R2GC_292_73',
                       value = '(23*complex(0,1)*G**4)/(192.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_293_74 = Coupling(name = 'R2GC_293_74',
                       value = '(31*complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_294_75 = Coupling(name = 'R2GC_294_75',
                       value = '(-17*complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_295_76 = Coupling(name = 'R2GC_295_76',
                       value = '(-7*complex(0,1)*G**4)/(32.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_296_77 = Coupling(name = 'R2GC_296_77',
                       value = '(7*complex(0,1)*G**4)/(64.*cmath.pi**2)',
                       order = {'QCD':4})

R2GC_301_78 = Coupling(name = 'R2GC_301_78',
                       value = '(complex(0,1)*G**2*MT)/(6.*cmath.pi**2)',
                       order = {'QCD':2})

R2GC_305_79 = Coupling(name = 'R2GC_305_79',
                       value = '(complex(0,1)*G**2*yb)/(3.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_306_80 = Coupling(name = 'R2GC_306_80',
                       value = '(complex(0,1)*G**2*tanbeta*ymb*cmath.sqrt(2))/(3.*cmath.pi**2*vev)',
                       order = {'QCD':2,'QED':1})

R2GC_307_81 = Coupling(name = 'R2GC_307_81',
                       value = '(complex(0,1)*G**2*ymt*cmath.sqrt(2))/(3.*cmath.pi**2*tanbeta*vev)',
                       order = {'QCD':2,'QED':1})

R2GC_308_82 = Coupling(name = 'R2GC_308_82',
                       value = '-(G**2*TH3x3*ymt)/(3.*cmath.pi**2*tanbeta*vev)',
                       order = {'QCD':2,'QED':1})

R2GC_309_83 = Coupling(name = 'R2GC_309_83',
                       value = '-(complex(0,1)*G**2*yt)/(3.*cmath.pi**2)',
                       order = {'QCD':2,'QED':1})

R2GC_310_84 = Coupling(name = 'R2GC_310_84',
                       value = '(G**2*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_311_85 = Coupling(name = 'R2GC_311_85',
                       value = '-(complex(0,1)*G**2*TH2x1*ymt)/(3.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x1*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

R2GC_312_86 = Coupling(name = 'R2GC_312_86',
                       value = '-(complex(0,1)*G**2*TH2x2*ymt)/(3.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x2*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',
                       order = {'QCD':2,'QED':1})

UVGC_248_1 = Coupling(name = 'UVGC_248_1',
                      value = {-1:'(51*G**3)/(128.*cmath.pi**2)'},
                      order = {'QCD':3})

UVGC_249_2 = Coupling(name = 'UVGC_249_2',
                      value = {-1:'G**3/(128.*cmath.pi**2)'},
                      order = {'QCD':3})

UVGC_250_3 = Coupling(name = 'UVGC_250_3',
                      value = {-1:'-(complex(0,1)*G**2)/(12.*cmath.pi**2)'},
                      order = {'QCD':2})

UVGC_254_4 = Coupling(name = 'UVGC_254_4',
                      value = {-1:'(3*G**4)/(512.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_254_5 = Coupling(name = 'UVGC_254_5',
                      value = {-1:'(-3*G**4)/(512.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_255_6 = Coupling(name = 'UVGC_255_6',
                      value = {-1:'(3*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_255_7 = Coupling(name = 'UVGC_255_7',
                      value = {-1:'(-3*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_257_8 = Coupling(name = 'UVGC_257_8',
                      value = {-1:'-(complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_257_9 = Coupling(name = 'UVGC_257_9',
                      value = {-1:'(complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                      order = {'QCD':4})

UVGC_258_10 = Coupling(name = 'UVGC_258_10',
                       value = {-1:'(-3*complex(0,1)*G**4)/(256.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_258_11 = Coupling(name = 'UVGC_258_11',
                       value = {-1:'(3*complex(0,1)*G**4)/(256.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_259_12 = Coupling(name = 'UVGC_259_12',
                       value = {-1:'-(complex(0,1)*G**4)/(24.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_259_13 = Coupling(name = 'UVGC_259_13',
                       value = {-1:'(47*complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_260_14 = Coupling(name = 'UVGC_260_14',
                       value = {-1:'(-253*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_260_15 = Coupling(name = 'UVGC_260_15',
                       value = {-1:'(5*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_261_16 = Coupling(name = 'UVGC_261_16',
                       value = {-1:'(-13*complex(0,1)*G**3)/(48.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_262_17 = Coupling(name = 'UVGC_262_17',
                       value = {-1:'(complex(0,1)*G**2)/(12.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_263_18 = Coupling(name = 'UVGC_263_18',
                       value = {-1:'( 0 if MB else (complex(0,1)*G**3)/(48.*cmath.pi**2) )'},
                       order = {'QCD':3})

UVGC_263_19 = Coupling(name = 'UVGC_263_19',
                       value = {-1:'(complex(0,1)*G**3)/(48.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_263_20 = Coupling(name = 'UVGC_263_20',
                       value = {-1:'(-19*complex(0,1)*G**3)/(128.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_263_21 = Coupling(name = 'UVGC_263_21',
                       value = {-1:'-(complex(0,1)*G**3)/(128.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_263_22 = Coupling(name = 'UVGC_263_22',
                       value = {-1:'( 0 if MT else (complex(0,1)*G**3)/(48.*cmath.pi**2) )'},
                       order = {'QCD':3})

UVGC_263_23 = Coupling(name = 'UVGC_263_23',
                       value = {-1:'(complex(0,1)*G**3)/(12.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_275_24 = Coupling(name = 'UVGC_275_24',
                       value = {-1:'(complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_275_25 = Coupling(name = 'UVGC_275_25',
                       value = {-1:'-(complex(0,1)*G**2*gw)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_277_26 = Coupling(name = 'UVGC_277_26',
                       value = {-1:'( (complex(0,1)*G**2)/(6.*cmath.pi**2) if MB else -(complex(0,1)*G**2)/(12.*cmath.pi**2) ) + (complex(0,1)*G**2)/(12.*cmath.pi**2)',0:'( (5*complex(0,1)*G**2)/(12.*cmath.pi**2) - (complex(0,1)*G**2*reglog(MB/MU_R))/(2.*cmath.pi**2) if MB else (complex(0,1)*G**2)/(12.*cmath.pi**2) ) - (complex(0,1)*G**2)/(12.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_278_27 = Coupling(name = 'UVGC_278_27',
                       value = {-1:'( -(complex(0,1)*G**3)/(6.*cmath.pi**2) if MB else (complex(0,1)*G**3)/(12.*cmath.pi**2) )',0:'( (-5*complex(0,1)*G**3)/(12.*cmath.pi**2) + (complex(0,1)*G**3*reglog(MB/MU_R))/(2.*cmath.pi**2) if MB else -(complex(0,1)*G**3)/(12.*cmath.pi**2) ) + (complex(0,1)*G**3)/(12.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_279_28 = Coupling(name = 'UVGC_279_28',
                       value = {-1:'( (cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2) if MB else -(cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2)',0:'( (5*cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) - (cw*complex(0,1)*G**2*g1*reglog(MB/MU_R))/(6.*cmath.pi**2) if MB else (cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_280_29 = Coupling(name = 'UVGC_280_29',
                       value = {-1:'( (complex(0,1)*G**2*MB)/(6.*cmath.pi**2) if MB else -(complex(0,1)*G**2*MB)/(12.*cmath.pi**2) ) + (complex(0,1)*G**2*MB)/(3.*cmath.pi**2)',0:'( (3*complex(0,1)*G**2*MB)/(4.*cmath.pi**2) - (complex(0,1)*G**2*MB*reglog(MB/MU_R))/cmath.pi**2 if MB else (complex(0,1)*G**2*MB)/(12.*cmath.pi**2) ) - (complex(0,1)*G**2*MB)/(12.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_281_30 = Coupling(name = 'UVGC_281_30',
                       value = {-1:'( (cw*complex(0,1)*G**2*gw)/(12.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2) if MB else -(cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2)',0:'( (5*cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (5*complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) - (cw*complex(0,1)*G**2*gw*reglog(MB/MU_R))/(4.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_282_31 = Coupling(name = 'UVGC_282_31',
                       value = {-1:'( -(complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2) if MB else (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2) ) - (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2)',0:'( (-5*complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw*reglog(MB/MU_R))/(6.*cmath.pi**2) if MB else -(complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2) ) + (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_283_32 = Coupling(name = 'UVGC_283_32',
                       value = {-1:'( -(cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(12.*cmath.pi**2) if MB else (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2)',0:'( (-5*cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) + (5*complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) + (cw*complex(0,1)*G**2*g1*reglog(MB/MU_R))/(12.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw*reglog(MB/MU_R))/(4.*cmath.pi**2) if MB else -(cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_284_33 = Coupling(name = 'UVGC_284_33',
                       value = {-1:'( -(G**2*yb)/(6.*cmath.pi**2*cmath.sqrt(2)) if MB else (G**2*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (G**2*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (-3*G**2*yb)/(4.*cmath.pi**2*cmath.sqrt(2)) + (G**2*yb*reglog(MB/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MB else -(G**2*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (G**2*yb)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_285_34 = Coupling(name = 'UVGC_285_34',
                       value = {-1:'( (complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(6.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x1*yb)/(6.*cmath.pi**2*cmath.sqrt(2)) if MB else -(complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(12.*cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x1*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(3.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x1*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (3*complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(4.*cmath.pi**2*vev) + (3*complex(0,1)*G**2*TH1x1*yb)/(4.*cmath.pi**2*cmath.sqrt(2)) - (complex(0,1)*G**2*tanbeta*TH2x1*ymb*reglog(MB/MU_R))/(cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x1*yb*reglog(MB/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MB else (complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(12.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x1*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (complex(0,1)*G**2*tanbeta*TH2x1*ymb)/(12.*cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x1*yb)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_286_35 = Coupling(name = 'UVGC_286_35',
                       value = {-1:'( (complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(6.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x2*yb)/(6.*cmath.pi**2*cmath.sqrt(2)) if MB else -(complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(12.*cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x2*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(3.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x2*yb)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (3*complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(4.*cmath.pi**2*vev) + (3*complex(0,1)*G**2*TH1x2*yb)/(4.*cmath.pi**2*cmath.sqrt(2)) - (complex(0,1)*G**2*tanbeta*TH2x2*ymb*reglog(MB/MU_R))/(cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x2*yb*reglog(MB/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MB else (complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(12.*cmath.pi**2*vev) + (complex(0,1)*G**2*TH1x2*yb)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (complex(0,1)*G**2*tanbeta*TH2x2*ymb)/(12.*cmath.pi**2*vev) - (complex(0,1)*G**2*TH1x2*yb)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_287_36 = Coupling(name = 'UVGC_287_36',
                       value = {-1:'( -(G**2*tanbeta*TH3x3*ymb)/(6.*cmath.pi**2*vev) if MB else (G**2*tanbeta*TH3x3*ymb)/(12.*cmath.pi**2*vev) ) - (G**2*tanbeta*TH3x3*ymb)/(3.*cmath.pi**2*vev)',0:'( (-3*G**2*tanbeta*TH3x3*ymb)/(4.*cmath.pi**2*vev) + (G**2*tanbeta*TH3x3*ymb*reglog(MB/MU_R))/(cmath.pi**2*vev) if MB else -(G**2*tanbeta*TH3x3*ymb)/(12.*cmath.pi**2*vev) ) + (G**2*tanbeta*TH3x3*ymb)/(12.*cmath.pi**2*vev)'},
                       order = {'QCD':2,'QED':1})

UVGC_288_37 = Coupling(name = 'UVGC_288_37',
                       value = {-1:'( 0 if MB else -(complex(0,1)*G**2)/(24.*cmath.pi**2) ) + (complex(0,1)*G**2)/(24.*cmath.pi**2)',0:'( -(complex(0,1)*G**2*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':2})

UVGC_288_38 = Coupling(name = 'UVGC_288_38',
                       value = {-1:'( 0 if MT else -(complex(0,1)*G**2)/(24.*cmath.pi**2) ) + (complex(0,1)*G**2)/(24.*cmath.pi**2)',0:'( -(complex(0,1)*G**2*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':2})

UVGC_289_39 = Coupling(name = 'UVGC_289_39',
                       value = {-1:'( 0 if MB else (complex(0,1)*G**2)/(24.*cmath.pi**2) ) - (complex(0,1)*G**2)/(24.*cmath.pi**2)',0:'( (complex(0,1)*G**2*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':2})

UVGC_289_40 = Coupling(name = 'UVGC_289_40',
                       value = {-1:'(3*complex(0,1)*G**2)/(64.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_289_41 = Coupling(name = 'UVGC_289_41',
                       value = {-1:'(-3*complex(0,1)*G**2)/(64.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_289_42 = Coupling(name = 'UVGC_289_42',
                       value = {-1:'( 0 if MT else (complex(0,1)*G**2)/(24.*cmath.pi**2) ) - (complex(0,1)*G**2)/(24.*cmath.pi**2)',0:'( (complex(0,1)*G**2*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':2})

UVGC_290_43 = Coupling(name = 'UVGC_290_43',
                       value = {-1:'( 0 if MB else -G**3/(16.*cmath.pi**2) ) + G**3/(24.*cmath.pi**2)',0:'( -(G**3*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':3})

UVGC_290_44 = Coupling(name = 'UVGC_290_44',
                       value = {-1:'-G**3/(48.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_290_45 = Coupling(name = 'UVGC_290_45',
                       value = {-1:'( 0 if MT else -G**3/(16.*cmath.pi**2) ) + G**3/(24.*cmath.pi**2)',0:'( -(G**3*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':3})

UVGC_291_46 = Coupling(name = 'UVGC_291_46',
                       value = {-1:'( 0 if MB else -(complex(0,1)*G**4)/(12.*cmath.pi**2) ) + (complex(0,1)*G**4)/(12.*cmath.pi**2)',0:'( -(complex(0,1)*G**4*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':4})

UVGC_291_47 = Coupling(name = 'UVGC_291_47',
                       value = {-1:'(147*complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_291_48 = Coupling(name = 'UVGC_291_48',
                       value = {-1:'(3*complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_291_49 = Coupling(name = 'UVGC_291_49',
                       value = {-1:'( 0 if MT else -(complex(0,1)*G**4)/(12.*cmath.pi**2) ) + (complex(0,1)*G**4)/(12.*cmath.pi**2)',0:'( -(complex(0,1)*G**4*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':4})

UVGC_292_50 = Coupling(name = 'UVGC_292_50',
                       value = {-1:'(147*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_292_51 = Coupling(name = 'UVGC_292_51',
                       value = {-1:'(21*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_293_52 = Coupling(name = 'UVGC_293_52',
                       value = {-1:'( 0 if MB else -(complex(0,1)*G**4)/(12.*cmath.pi**2) )',0:'( -(complex(0,1)*G**4*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':4})

UVGC_293_53 = Coupling(name = 'UVGC_293_53',
                       value = {-1:'-(complex(0,1)*G**4)/(12.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_293_54 = Coupling(name = 'UVGC_293_54',
                       value = {-1:'(523*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_293_55 = Coupling(name = 'UVGC_293_55',
                       value = {-1:'(13*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_293_56 = Coupling(name = 'UVGC_293_56',
                       value = {-1:'( 0 if MT else -(complex(0,1)*G**4)/(12.*cmath.pi**2) )',0:'( -(complex(0,1)*G**4*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':4})

UVGC_294_57 = Coupling(name = 'UVGC_294_57',
                       value = {-1:'( 0 if MB else (complex(0,1)*G**4)/(12.*cmath.pi**2) ) - (complex(0,1)*G**4)/(24.*cmath.pi**2)',0:'( (complex(0,1)*G**4*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':4})

UVGC_294_58 = Coupling(name = 'UVGC_294_58',
                       value = {-1:'(complex(0,1)*G**4)/(24.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_294_59 = Coupling(name = 'UVGC_294_59',
                       value = {-1:'(-341*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_294_60 = Coupling(name = 'UVGC_294_60',
                       value = {-1:'(-11*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_294_61 = Coupling(name = 'UVGC_294_61',
                       value = {-1:'( 0 if MT else (complex(0,1)*G**4)/(12.*cmath.pi**2) ) - (complex(0,1)*G**4)/(24.*cmath.pi**2)',0:'( (complex(0,1)*G**4*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':4})

UVGC_295_62 = Coupling(name = 'UVGC_295_62',
                       value = {-1:'(-83*complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_295_63 = Coupling(name = 'UVGC_295_63',
                       value = {-1:'(-5*complex(0,1)*G**4)/(128.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_296_64 = Coupling(name = 'UVGC_296_64',
                       value = {-1:'( 0 if MB else (complex(0,1)*G**4)/(12.*cmath.pi**2) )',0:'( (complex(0,1)*G**4*reglog(MB/MU_R))/(12.*cmath.pi**2) if MB else 0 )'},
                       order = {'QCD':4})

UVGC_296_65 = Coupling(name = 'UVGC_296_65',
                       value = {-1:'(complex(0,1)*G**4)/(12.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_296_66 = Coupling(name = 'UVGC_296_66',
                       value = {-1:'(-85*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_296_67 = Coupling(name = 'UVGC_296_67',
                       value = {-1:'(-19*complex(0,1)*G**4)/(512.*cmath.pi**2)'},
                       order = {'QCD':4})

UVGC_296_68 = Coupling(name = 'UVGC_296_68',
                       value = {-1:'( 0 if MT else (complex(0,1)*G**4)/(12.*cmath.pi**2) )',0:'( (complex(0,1)*G**4*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else 0 )'},
                       order = {'QCD':4})

UVGC_297_69 = Coupling(name = 'UVGC_297_69',
                       value = {-1:'( (complex(0,1)*G**2)/(6.*cmath.pi**2) if MT else -(complex(0,1)*G**2)/(12.*cmath.pi**2) ) + (complex(0,1)*G**2)/(12.*cmath.pi**2)',0:'( (5*complex(0,1)*G**2)/(12.*cmath.pi**2) - (complex(0,1)*G**2*reglog(MT/MU_R))/(2.*cmath.pi**2) if MT else (complex(0,1)*G**2)/(12.*cmath.pi**2) ) - (complex(0,1)*G**2)/(12.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_298_70 = Coupling(name = 'UVGC_298_70',
                       value = {-1:'( -(complex(0,1)*G**3)/(6.*cmath.pi**2) if MT else (complex(0,1)*G**3)/(12.*cmath.pi**2) )',0:'( (-5*complex(0,1)*G**3)/(12.*cmath.pi**2) + (complex(0,1)*G**3*reglog(MT/MU_R))/(2.*cmath.pi**2) if MT else -(complex(0,1)*G**3)/(12.*cmath.pi**2) ) + (complex(0,1)*G**3)/(12.*cmath.pi**2)'},
                       order = {'QCD':3})

UVGC_299_71 = Coupling(name = 'UVGC_299_71',
                       value = {-1:'( -(cw*complex(0,1)*G**2*g1)/(9.*cmath.pi**2) if MT else (cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2)',0:'( (-5*cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2) + (cw*complex(0,1)*G**2*g1*reglog(MT/MU_R))/(3.*cmath.pi**2) if MT else -(cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*g1)/(18.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_300_72 = Coupling(name = 'UVGC_300_72',
                       value = {-1:'( -(complex(0,1)*G**2*gw)/(12.*cmath.pi**2*cmath.sqrt(2)) if MB else (complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) )',0:'( (-5*complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) + (complex(0,1)*G**2*gw*reglog(MB/MU_R))/(4.*cmath.pi**2*cmath.sqrt(2)) if MB else -(complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_300_73 = Coupling(name = 'UVGC_300_73',
                       value = {-1:'( -(complex(0,1)*G**2*gw)/(12.*cmath.pi**2*cmath.sqrt(2)) if MT else (complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) )',0:'( (-5*complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) + (complex(0,1)*G**2*gw*reglog(MT/MU_R))/(4.*cmath.pi**2*cmath.sqrt(2)) if MT else -(complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*gw)/(24.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_301_74 = Coupling(name = 'UVGC_301_74',
                       value = {-1:'( (complex(0,1)*G**2*MT)/(6.*cmath.pi**2) if MT else -(complex(0,1)*G**2*MT)/(12.*cmath.pi**2) ) + (complex(0,1)*G**2*MT)/(3.*cmath.pi**2)',0:'( (3*complex(0,1)*G**2*MT)/(4.*cmath.pi**2) - (complex(0,1)*G**2*MT*reglog(MT/MU_R))/cmath.pi**2 if MT else (complex(0,1)*G**2*MT)/(12.*cmath.pi**2) ) - (complex(0,1)*G**2*MT)/(12.*cmath.pi**2)'},
                       order = {'QCD':2})

UVGC_302_75 = Coupling(name = 'UVGC_302_75',
                       value = {-1:'( -(cw*complex(0,1)*G**2*gw)/(12.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(36.*cmath.pi**2) if MT else (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2)',0:'( (-5*cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (5*complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) + (cw*complex(0,1)*G**2*gw*reglog(MT/MU_R))/(4.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw*reglog(MT/MU_R))/(12.*cmath.pi**2) if MT else -(cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) + (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*gw)/(24.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw)/(72.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_303_76 = Coupling(name = 'UVGC_303_76',
                       value = {-1:'( (complex(0,1)*G**2*g1*sw)/(9.*cmath.pi**2) if MT else -(complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2) ) + (complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2)',0:'( (5*complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2) - (complex(0,1)*G**2*g1*sw*reglog(MT/MU_R))/(3.*cmath.pi**2) if MT else (complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2) ) - (complex(0,1)*G**2*g1*sw)/(18.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_304_77 = Coupling(name = 'UVGC_304_77',
                       value = {-1:'( -(cw*complex(0,1)*G**2*g1)/(36.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(12.*cmath.pi**2) if MT else (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) ) - (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2)',0:'( (-5*cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) - (5*complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) + (cw*complex(0,1)*G**2*g1*reglog(MT/MU_R))/(12.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw*reglog(MT/MU_R))/(4.*cmath.pi**2) if MT else -(cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) - (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2) ) + (cw*complex(0,1)*G**2*g1)/(72.*cmath.pi**2) + (complex(0,1)*G**2*gw*sw)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_305_78 = Coupling(name = 'UVGC_305_78',
                       value = {-1:'( (complex(0,1)*G**2*yb)/(12.*cmath.pi**2) if MB else -(complex(0,1)*G**2*yb)/(24.*cmath.pi**2) )',0:'( (13*complex(0,1)*G**2*yb)/(24.*cmath.pi**2) - (3*complex(0,1)*G**2*yb*reglog(MB/MU_R))/(4.*cmath.pi**2) if MB else (complex(0,1)*G**2*yb)/(24.*cmath.pi**2) ) - (complex(0,1)*G**2*yb)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_305_79 = Coupling(name = 'UVGC_305_79',
                       value = {-1:'( (complex(0,1)*G**2*yb)/(12.*cmath.pi**2) if MT else -(complex(0,1)*G**2*yb)/(24.*cmath.pi**2) )',0:'( (5*complex(0,1)*G**2*yb)/(24.*cmath.pi**2) - (complex(0,1)*G**2*yb*reglog(MT/MU_R))/(4.*cmath.pi**2) if MT else (complex(0,1)*G**2*yb)/(24.*cmath.pi**2) ) - (complex(0,1)*G**2*yb)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_305_80 = Coupling(name = 'UVGC_305_80',
                       value = {-1:'(complex(0,1)*G**2*yb)/(3.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_306_81 = Coupling(name = 'UVGC_306_81',
                       value = {-1:'( (complex(0,1)*G**2*tanbeta*ymb)/(6.*cmath.pi**2*vev*cmath.sqrt(2)) if MB else -(complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) )',0:'( (13*complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) - (3*complex(0,1)*G**2*tanbeta*ymb*reglog(MB/MU_R))/(2.*cmath.pi**2*vev*cmath.sqrt(2)) if MB else (complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) ) - (complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_306_82 = Coupling(name = 'UVGC_306_82',
                       value = {-1:'( (complex(0,1)*G**2*tanbeta*ymb)/(6.*cmath.pi**2*vev*cmath.sqrt(2)) if MT else -(complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) )',0:'( (5*complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) - (complex(0,1)*G**2*tanbeta*ymb*reglog(MT/MU_R))/(2.*cmath.pi**2*vev*cmath.sqrt(2)) if MT else (complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2)) ) - (complex(0,1)*G**2*tanbeta*ymb)/(12.*cmath.pi**2*vev*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_306_83 = Coupling(name = 'UVGC_306_83',
                       value = {-1:'(complex(0,1)*G**2*tanbeta*ymb*cmath.sqrt(2))/(3.*cmath.pi**2*vev)'},
                       order = {'QCD':2,'QED':1})

UVGC_307_84 = Coupling(name = 'UVGC_307_84',
                       value = {-1:'( (complex(0,1)*G**2*ymt)/(6.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) if MB else -(complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) )',0:'( (5*complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) - (complex(0,1)*G**2*ymt*reglog(MB/MU_R))/(2.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) if MB else (complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) ) - (complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_307_85 = Coupling(name = 'UVGC_307_85',
                       value = {-1:'( (complex(0,1)*G**2*ymt)/(6.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) if MT else -(complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) )',0:'( (13*complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) - (3*complex(0,1)*G**2*ymt*reglog(MT/MU_R))/(2.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) if MT else (complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2)) ) - (complex(0,1)*G**2*ymt)/(12.*cmath.pi**2*tanbeta*vev*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_307_86 = Coupling(name = 'UVGC_307_86',
                       value = {-1:'(complex(0,1)*G**2*ymt*cmath.sqrt(2))/(3.*cmath.pi**2*tanbeta*vev)'},
                       order = {'QCD':2,'QED':1})

UVGC_308_87 = Coupling(name = 'UVGC_308_87',
                       value = {-1:'( -(G**2*TH3x3*ymt)/(6.*cmath.pi**2*tanbeta*vev) if MT else (G**2*TH3x3*ymt)/(12.*cmath.pi**2*tanbeta*vev) ) - (G**2*TH3x3*ymt)/(3.*cmath.pi**2*tanbeta*vev)',0:'( (-3*G**2*TH3x3*ymt)/(4.*cmath.pi**2*tanbeta*vev) + (G**2*TH3x3*ymt*reglog(MT/MU_R))/(cmath.pi**2*tanbeta*vev) if MT else -(G**2*TH3x3*ymt)/(12.*cmath.pi**2*tanbeta*vev) ) + (G**2*TH3x3*ymt)/(12.*cmath.pi**2*tanbeta*vev)'},
                       order = {'QCD':2,'QED':1})

UVGC_309_88 = Coupling(name = 'UVGC_309_88',
                       value = {-1:'( -(complex(0,1)*G**2*yt)/(12.*cmath.pi**2) if MB else (complex(0,1)*G**2*yt)/(24.*cmath.pi**2) )',0:'( (-5*complex(0,1)*G**2*yt)/(24.*cmath.pi**2) + (complex(0,1)*G**2*yt*reglog(MB/MU_R))/(4.*cmath.pi**2) if MB else -(complex(0,1)*G**2*yt)/(24.*cmath.pi**2) ) + (complex(0,1)*G**2*yt)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_309_89 = Coupling(name = 'UVGC_309_89',
                       value = {-1:'( -(complex(0,1)*G**2*yt)/(12.*cmath.pi**2) if MT else (complex(0,1)*G**2*yt)/(24.*cmath.pi**2) )',0:'( (-13*complex(0,1)*G**2*yt)/(24.*cmath.pi**2) + (3*complex(0,1)*G**2*yt*reglog(MT/MU_R))/(4.*cmath.pi**2) if MT else -(complex(0,1)*G**2*yt)/(24.*cmath.pi**2) ) + (complex(0,1)*G**2*yt)/(24.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_309_90 = Coupling(name = 'UVGC_309_90',
                       value = {-1:'-(complex(0,1)*G**2*yt)/(3.*cmath.pi**2)'},
                       order = {'QCD':2,'QED':1})

UVGC_310_91 = Coupling(name = 'UVGC_310_91',
                       value = {-1:'( (G**2*yt)/(6.*cmath.pi**2*cmath.sqrt(2)) if MT else -(G**2*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (G**2*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (3*G**2*yt)/(4.*cmath.pi**2*cmath.sqrt(2)) - (G**2*yt*reglog(MT/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MT else (G**2*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (G**2*yt)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_311_92 = Coupling(name = 'UVGC_311_92',
                       value = {-1:'( -(complex(0,1)*G**2*TH2x1*ymt)/(6.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x1*yt)/(6.*cmath.pi**2*cmath.sqrt(2)) if MT else (complex(0,1)*G**2*TH2x1*ymt)/(12.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x1*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (complex(0,1)*G**2*TH2x1*ymt)/(3.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x1*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (-3*complex(0,1)*G**2*TH2x1*ymt)/(4.*cmath.pi**2*tanbeta*vev) + (3*complex(0,1)*G**2*TH1x1*yt)/(4.*cmath.pi**2*cmath.sqrt(2)) + (complex(0,1)*G**2*TH2x1*ymt*reglog(MT/MU_R))/(cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x1*yt*reglog(MT/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MT else -(complex(0,1)*G**2*TH2x1*ymt)/(12.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x1*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*TH2x1*ymt)/(12.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x1*yt)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

UVGC_312_93 = Coupling(name = 'UVGC_312_93',
                       value = {-1:'( -(complex(0,1)*G**2*TH2x2*ymt)/(6.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x2*yt)/(6.*cmath.pi**2*cmath.sqrt(2)) if MT else (complex(0,1)*G**2*TH2x2*ymt)/(12.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x2*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) - (complex(0,1)*G**2*TH2x2*ymt)/(3.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x2*yt)/(3.*cmath.pi**2*cmath.sqrt(2))',0:'( (-3*complex(0,1)*G**2*TH2x2*ymt)/(4.*cmath.pi**2*tanbeta*vev) + (3*complex(0,1)*G**2*TH1x2*yt)/(4.*cmath.pi**2*cmath.sqrt(2)) + (complex(0,1)*G**2*TH2x2*ymt*reglog(MT/MU_R))/(cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x2*yt*reglog(MT/MU_R))/(cmath.pi**2*cmath.sqrt(2)) if MT else -(complex(0,1)*G**2*TH2x2*ymt)/(12.*cmath.pi**2*tanbeta*vev) + (complex(0,1)*G**2*TH1x2*yt)/(12.*cmath.pi**2*cmath.sqrt(2)) ) + (complex(0,1)*G**2*TH2x2*ymt)/(12.*cmath.pi**2*tanbeta*vev) - (complex(0,1)*G**2*TH1x2*yt)/(12.*cmath.pi**2*cmath.sqrt(2))'},
                       order = {'QCD':2,'QED':1})

