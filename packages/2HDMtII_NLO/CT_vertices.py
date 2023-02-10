# This file was automatically created by FeynRules 2.3.13
# Mathematica version: 10.1.0  for Mac OS X x86 (64-bit) (March 24, 2015)
# Date: Tue 8 Dec 2015 11:44:25


from object_library import all_vertices, all_CTvertices, Vertex, CTVertex
import particles as P
import CT_couplings as C
import lorentz as L


V_1 = CTVertex(name = 'V_1',
               type = 'R2',
               particles = [ P.g, P.g, P.g ],
               color = [ 'f(1,2,3)' ],
               lorentz = [ L.VVV8 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ] ],
               couplings = {(0,0,0):C.R2GC_290_69,(0,0,1):C.R2GC_290_70})

V_2 = CTVertex(name = 'V_2',
               type = 'R2',
               particles = [ P.g, P.g, P.g, P.g ],
               color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
               lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5, L.VVVV8 ],
               loop_particles = [ [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ] ],
               couplings = {(2,0,0):C.R2GC_256_49,(2,0,1):C.R2GC_256_50,(0,0,0):C.R2GC_256_49,(0,0,1):C.R2GC_256_50,(4,0,0):C.R2GC_254_45,(4,0,1):C.R2GC_254_46,(3,0,0):C.R2GC_254_45,(3,0,1):C.R2GC_254_46,(8,0,0):C.R2GC_255_47,(8,0,1):C.R2GC_255_48,(7,0,0):C.R2GC_260_56,(7,0,1):C.R2GC_294_75,(6,0,0):C.R2GC_259_54,(6,0,1):C.R2GC_295_76,(5,0,0):C.R2GC_254_45,(5,0,1):C.R2GC_254_46,(1,0,0):C.R2GC_254_45,(1,0,1):C.R2GC_254_46,(11,3,0):C.R2GC_258_52,(11,3,1):C.R2GC_258_53,(10,3,0):C.R2GC_258_52,(10,3,1):C.R2GC_258_53,(9,3,1):C.R2GC_257_51,(2,1,0):C.R2GC_256_49,(2,1,1):C.R2GC_256_50,(0,1,0):C.R2GC_256_49,(0,1,1):C.R2GC_256_50,(4,1,0):C.R2GC_254_45,(4,1,1):C.R2GC_254_46,(3,1,0):C.R2GC_254_45,(3,1,1):C.R2GC_254_46,(8,1,0):C.R2GC_255_47,(8,1,1):C.R2GC_296_77,(6,1,0):C.R2GC_291_71,(6,1,1):C.R2GC_291_72,(7,1,0):C.R2GC_260_56,(7,1,1):C.R2GC_260_57,(5,1,0):C.R2GC_254_45,(5,1,1):C.R2GC_254_46,(1,1,0):C.R2GC_254_45,(1,1,1):C.R2GC_254_46,(2,2,0):C.R2GC_256_49,(2,2,1):C.R2GC_256_50,(0,2,0):C.R2GC_256_49,(0,2,1):C.R2GC_256_50,(4,2,0):C.R2GC_254_45,(4,2,1):C.R2GC_254_46,(3,2,0):C.R2GC_254_45,(3,2,1):C.R2GC_254_46,(8,2,0):C.R2GC_255_47,(8,2,1):C.R2GC_293_74,(6,2,0):C.R2GC_259_54,(6,2,1):C.R2GC_259_55,(7,2,0):C.R2GC_292_73,(7,2,1):C.R2GC_256_50,(5,2,0):C.R2GC_254_45,(5,2,1):C.R2GC_254_46,(1,2,0):C.R2GC_254_45,(1,2,1):C.R2GC_254_46})

V_3 = CTVertex(name = 'V_3',
               type = 'R2',
               particles = [ P.t__tilde__, P.b, P.h__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS5 ],
               loop_particles = [ [ [P.b, P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_306_80,(0,1,0):C.R2GC_307_81})

V_4 = CTVertex(name = 'V_4',
               type = 'R2',
               particles = [ P.b__tilde__, P.t, P.h__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS5 ],
               loop_particles = [ [ [P.b, P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_307_81,(0,1,0):C.R2GC_306_80})

V_5 = CTVertex(name = 'V_5',
               type = 'R2',
               particles = [ P.t__tilde__, P.b, P.G__plus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS5 ],
               loop_particles = [ [ [P.b, P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_305_79,(0,1,0):C.R2GC_309_83})

V_6 = CTVertex(name = 'V_6',
               type = 'R2',
               particles = [ P.b__tilde__, P.b, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1 ],
               loop_particles = [ [ [P.b, P.g] ] ],
               couplings = {(0,0,0):C.R2GC_284_62})

V_7 = CTVertex(name = 'V_7',
               type = 'R2',
               particles = [ P.t__tilde__, P.t, P.G0 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS1 ],
               loop_particles = [ [ [P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_310_84})

V_8 = CTVertex(name = 'V_8',
               type = 'R2',
               particles = [ P.b__tilde__, P.t, P.G__minus__ ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS3, L.FFS5 ],
               loop_particles = [ [ [P.b, P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_309_83,(0,1,0):C.R2GC_305_79})

V_9 = CTVertex(name = 'V_9',
               type = 'R2',
               particles = [ P.t__tilde__, P.t, P.h1 ],
               color = [ 'Identity(1,2)' ],
               lorentz = [ L.FFS2 ],
               loop_particles = [ [ [P.g, P.t] ] ],
               couplings = {(0,0,0):C.R2GC_311_85})

V_10 = CTVertex(name = 'V_10',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_312_86})

V_11 = CTVertex(name = 'V_11',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_308_82})

V_12 = CTVertex(name = 'V_12',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_285_63})

V_13 = CTVertex(name = 'V_13',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_286_64})

V_14 = CTVertex(name = 'V_14',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_287_65})

V_15 = CTVertex(name = 'V_15',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_237_38,(0,1,0):C.R2GC_212_1})

V_16 = CTVertex(name = 'V_16',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_237_38,(0,1,0):C.R2GC_212_1})

V_17 = CTVertex(name = 'V_17',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_237_38,(0,1,0):C.R2GC_212_1})

V_18 = CTVertex(name = 'V_18',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_19 = CTVertex(name = 'V_19',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_20 = CTVertex(name = 'V_20',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_21 = CTVertex(name = 'V_21',
                type = 'R2',
                particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_22 = CTVertex(name = 'V_22',
                type = 'R2',
                particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_23 = CTVertex(name = 'V_23',
                type = 'R2',
                particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_24 = CTVertex(name = 'V_24',
                type = 'R2',
                particles = [ P.u__tilde__, P.u, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_236_37,(0,1,0):C.R2GC_213_2})

V_25 = CTVertex(name = 'V_25',
                type = 'R2',
                particles = [ P.c__tilde__, P.c, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_236_37,(0,1,0):C.R2GC_213_2})

V_26 = CTVertex(name = 'V_26',
                type = 'R2',
                particles = [ P.t__tilde__, P.t, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_236_37,(0,1,0):C.R2GC_213_2})

V_27 = CTVertex(name = 'V_27',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_239_40,(0,1,0):C.R2GC_214_3})

V_28 = CTVertex(name = 'V_28',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_239_40,(0,1,0):C.R2GC_214_3})

V_29 = CTVertex(name = 'V_29',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_239_40,(0,1,0):C.R2GC_214_3})

V_30 = CTVertex(name = 'V_30',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_31 = CTVertex(name = 'V_31',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_32 = CTVertex(name = 'V_32',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_261_58})

V_33 = CTVertex(name = 'V_33',
                type = 'R2',
                particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_34 = CTVertex(name = 'V_34',
                type = 'R2',
                particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_35 = CTVertex(name = 'V_35',
                type = 'R2',
                particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.b, P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_275_60})

V_36 = CTVertex(name = 'V_36',
                type = 'R2',
                particles = [ P.d__tilde__, P.d, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_238_39,(0,1,0):C.R2GC_215_4})

V_37 = CTVertex(name = 'V_37',
                type = 'R2',
                particles = [ P.s__tilde__, P.s, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_238_39,(0,1,0):C.R2GC_215_4})

V_38 = CTVertex(name = 'V_38',
                type = 'R2',
                particles = [ P.b__tilde__, P.b, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_238_39,(0,1,0):C.R2GC_215_4})

V_39 = CTVertex(name = 'V_39',
                type = 'R2',
                particles = [ P.u__tilde__, P.u ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_262_59})

V_40 = CTVertex(name = 'V_40',
                type = 'R2',
                particles = [ P.c__tilde__, P.c ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_262_59})

V_41 = CTVertex(name = 'V_41',
                type = 'R2',
                particles = [ P.t__tilde__, P.t ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_301_78,(0,2,0):C.R2GC_301_78,(0,1,0):C.R2GC_262_59,(0,3,0):C.R2GC_262_59})

V_42 = CTVertex(name = 'V_42',
                type = 'R2',
                particles = [ P.d__tilde__, P.d ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_262_59})

V_43 = CTVertex(name = 'V_43',
                type = 'R2',
                particles = [ P.s__tilde__, P.s ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.R2GC_262_59})

V_44 = CTVertex(name = 'V_44',
                type = 'R2',
                particles = [ P.b__tilde__, P.b ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.R2GC_280_61,(0,2,0):C.R2GC_280_61,(0,1,0):C.R2GC_262_59,(0,3,0):C.R2GC_262_59})

V_45 = CTVertex(name = 'V_45',
                type = 'R2',
                particles = [ P.g, P.g ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VV1, L.VV2, L.VV3 ],
                loop_particles = [ [ [P.b] ], [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.g] ], [ [P.t] ] ],
                couplings = {(0,0,2):C.R2GC_289_68,(0,1,0):C.R2GC_220_5,(0,1,3):C.R2GC_220_6,(0,2,1):C.R2GC_288_66,(0,2,2):C.R2GC_288_67})

V_46 = CTVertex(name = 'V_46',
                type = 'R2',
                particles = [ P.g, P.g, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_222_9,(0,0,1):C.R2GC_222_10})

V_47 = CTVertex(name = 'V_47',
                type = 'R2',
                particles = [ P.g, P.g, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_223_11,(0,0,1):C.R2GC_223_12})

V_48 = CTVertex(name = 'V_48',
                type = 'R2',
                particles = [ P.g, P.g, P.W__minus__, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b, P.t], [P.c, P.s], [P.d, P.u] ] ],
                couplings = {(0,0,0):C.R2GC_245_44,(0,1,0):C.R2GC_245_44,(0,2,0):C.R2GC_245_44})

V_49 = CTVertex(name = 'V_49',
                type = 'R2',
                particles = [ P.g, P.g, P.Z, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                couplings = {(0,0,0):C.R2GC_233_31,(0,0,1):C.R2GC_233_32,(0,1,0):C.R2GC_233_31,(0,1,1):C.R2GC_233_32,(0,2,0):C.R2GC_233_31,(0,2,1):C.R2GC_233_32})

V_50 = CTVertex(name = 'V_50',
                type = 'R2',
                particles = [ P.a, P.g, P.g, P.Z ],
                color = [ 'Identity(2,3)' ],
                lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                couplings = {(0,0,0):C.R2GC_234_33,(0,0,1):C.R2GC_234_34,(0,1,0):C.R2GC_234_33,(0,1,1):C.R2GC_234_34,(0,2,0):C.R2GC_234_33,(0,2,1):C.R2GC_234_34})

V_51 = CTVertex(name = 'V_51',
                type = 'R2',
                particles = [ P.a, P.a, P.g, P.g ],
                color = [ 'Identity(3,4)' ],
                lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                couplings = {(0,0,0):C.R2GC_235_35,(0,0,1):C.R2GC_235_36,(0,1,0):C.R2GC_235_35,(0,1,1):C.R2GC_235_36,(0,2,0):C.R2GC_235_35,(0,2,1):C.R2GC_235_36})

V_52 = CTVertex(name = 'V_52',
                type = 'R2',
                particles = [ P.g, P.g, P.g, P.Z ],
                color = [ 'd(1,2,3)', 'f(1,2,3)' ],
                lorentz = [ L.VVVV1, L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                couplings = {(1,0,0):C.R2GC_230_25,(1,0,1):C.R2GC_230_26,(0,1,0):C.R2GC_229_23,(0,1,1):C.R2GC_229_24,(0,2,0):C.R2GC_229_23,(0,2,1):C.R2GC_229_24,(0,3,0):C.R2GC_229_23,(0,3,1):C.R2GC_229_24})

V_53 = CTVertex(name = 'V_53',
                type = 'R2',
                particles = [ P.a, P.g, P.g, P.g ],
                color = [ 'd(2,3,4)', 'f(2,3,4)' ],
                lorentz = [ L.VVVV1, L.VVVV2, L.VVVV3, L.VVVV5 ],
                loop_particles = [ [ [P.b], [P.d], [P.s] ], [ [P.c], [P.t], [P.u] ] ],
                couplings = {(1,0,0):C.R2GC_232_29,(1,0,1):C.R2GC_232_30,(0,1,0):C.R2GC_231_27,(0,1,1):C.R2GC_231_28,(0,2,0):C.R2GC_231_27,(0,2,1):C.R2GC_231_28,(0,3,0):C.R2GC_231_27,(0,3,1):C.R2GC_231_28})

V_54 = CTVertex(name = 'V_54',
                type = 'R2',
                particles = [ P.g, P.g, P.G0, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_225_15,(0,0,1):C.R2GC_225_16})

V_55 = CTVertex(name = 'V_55',
                type = 'R2',
                particles = [ P.g, P.g, P.G__minus__, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_244_43})

V_56 = CTVertex(name = 'V_56',
                type = 'R2',
                particles = [ P.g, P.g, P.G__minus__, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_243_42})

V_57 = CTVertex(name = 'V_57',
                type = 'R2',
                particles = [ P.g, P.g, P.G__plus__, P.h__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_243_42})

V_58 = CTVertex(name = 'V_58',
                type = 'R2',
                particles = [ P.g, P.g, P.h__minus__, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b, P.t] ] ],
                couplings = {(0,0,0):C.R2GC_242_41})

V_59 = CTVertex(name = 'V_59',
                type = 'R2',
                particles = [ P.g, P.g, P.h1, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_226_17,(0,0,1):C.R2GC_226_18})

V_60 = CTVertex(name = 'V_60',
                type = 'R2',
                particles = [ P.g, P.g, P.h1, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_227_19,(0,0,1):C.R2GC_227_20})

V_61 = CTVertex(name = 'V_61',
                type = 'R2',
                particles = [ P.g, P.g, P.h2, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_228_21,(0,0,1):C.R2GC_228_22})

V_62 = CTVertex(name = 'V_62',
                type = 'R2',
                particles = [ P.g, P.g, P.G0, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_224_13,(0,0,1):C.R2GC_224_14})

V_63 = CTVertex(name = 'V_63',
                type = 'R2',
                particles = [ P.g, P.g, P.h3, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.VVSS1 ],
                loop_particles = [ [ [P.b] ], [ [P.t] ] ],
                couplings = {(0,0,0):C.R2GC_221_7,(0,0,1):C.R2GC_221_8})

V_64 = CTVertex(name = 'V_64',
                type = 'UV',
                particles = [ P.g, P.g, P.g ],
                color = [ 'f(1,2,3)' ],
                lorentz = [ L.VVV7, L.VVV8, L.VVV9 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                couplings = {(0,1,0):C.UVGC_290_43,(0,1,1):C.UVGC_290_44,(0,1,4):C.UVGC_290_45,(0,2,2):C.UVGC_248_1,(0,0,3):C.UVGC_249_2})

V_65 = CTVertex(name = 'V_65',
                type = 'UV',
                particles = [ P.g, P.g, P.g, P.g ],
                color = [ 'd(-1,1,3)*d(-1,2,4)', 'd(-1,1,3)*f(-1,2,4)', 'd(-1,1,4)*d(-1,2,3)', 'd(-1,1,4)*f(-1,2,3)', 'd(-1,2,3)*f(-1,1,4)', 'd(-1,2,4)*f(-1,1,3)', 'f(-1,1,2)*f(-1,3,4)', 'f(-1,1,3)*f(-1,2,4)', 'f(-1,1,4)*f(-1,2,3)', 'Identity(1,2)*Identity(3,4)', 'Identity(1,3)*Identity(2,4)', 'Identity(1,4)*Identity(2,3)' ],
                lorentz = [ L.VVVV2, L.VVVV3, L.VVVV5, L.VVVV8 ],
                loop_particles = [ [ [P.b] ], [ [P.b], [P.c], [P.d], [P.s], [P.t], [P.u] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                couplings = {(2,0,3):C.UVGC_255_7,(2,0,4):C.UVGC_255_6,(0,0,3):C.UVGC_255_7,(0,0,4):C.UVGC_255_6,(4,0,3):C.UVGC_254_4,(4,0,4):C.UVGC_254_5,(3,0,3):C.UVGC_254_4,(3,0,4):C.UVGC_254_5,(8,0,3):C.UVGC_255_6,(8,0,4):C.UVGC_255_7,(7,0,0):C.UVGC_294_57,(7,0,2):C.UVGC_294_58,(7,0,3):C.UVGC_294_59,(7,0,4):C.UVGC_294_60,(7,0,5):C.UVGC_294_61,(6,0,0):C.UVGC_294_57,(6,0,2):C.UVGC_294_58,(6,0,3):C.UVGC_295_62,(6,0,4):C.UVGC_295_63,(6,0,5):C.UVGC_294_61,(5,0,3):C.UVGC_254_4,(5,0,4):C.UVGC_254_5,(1,0,3):C.UVGC_254_4,(1,0,4):C.UVGC_254_5,(11,3,3):C.UVGC_258_10,(11,3,4):C.UVGC_258_11,(10,3,3):C.UVGC_258_10,(10,3,4):C.UVGC_258_11,(9,3,3):C.UVGC_257_8,(9,3,4):C.UVGC_257_9,(2,1,3):C.UVGC_255_7,(2,1,4):C.UVGC_255_6,(0,1,3):C.UVGC_255_7,(0,1,4):C.UVGC_255_6,(4,1,3):C.UVGC_254_4,(4,1,4):C.UVGC_254_5,(3,1,3):C.UVGC_254_4,(3,1,4):C.UVGC_254_5,(8,1,0):C.UVGC_296_64,(8,1,2):C.UVGC_296_65,(8,1,3):C.UVGC_296_66,(8,1,4):C.UVGC_296_67,(8,1,5):C.UVGC_296_68,(6,1,0):C.UVGC_291_46,(6,1,3):C.UVGC_291_47,(6,1,4):C.UVGC_291_48,(6,1,5):C.UVGC_291_49,(7,1,1):C.UVGC_259_12,(7,1,3):C.UVGC_260_14,(7,1,4):C.UVGC_260_15,(5,1,3):C.UVGC_254_4,(5,1,4):C.UVGC_254_5,(1,1,3):C.UVGC_254_4,(1,1,4):C.UVGC_254_5,(2,2,3):C.UVGC_255_7,(2,2,4):C.UVGC_255_6,(0,2,3):C.UVGC_255_7,(0,2,4):C.UVGC_255_6,(4,2,3):C.UVGC_254_4,(4,2,4):C.UVGC_254_5,(3,2,3):C.UVGC_254_4,(3,2,4):C.UVGC_254_5,(8,2,0):C.UVGC_293_52,(8,2,2):C.UVGC_293_53,(8,2,3):C.UVGC_293_54,(8,2,4):C.UVGC_293_55,(8,2,5):C.UVGC_293_56,(6,2,1):C.UVGC_259_12,(6,2,3):C.UVGC_259_13,(6,2,4):C.UVGC_257_8,(7,2,0):C.UVGC_291_46,(7,2,3):C.UVGC_292_50,(7,2,4):C.UVGC_292_51,(7,2,5):C.UVGC_291_49,(5,2,3):C.UVGC_254_4,(5,2,4):C.UVGC_254_5,(1,2,3):C.UVGC_254_4,(1,2,4):C.UVGC_254_5})

V_66 = CTVertex(name = 'V_66',
                type = 'UV',
                particles = [ P.t__tilde__, P.b, P.h__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS5 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_306_81,(0,0,2):C.UVGC_306_82,(0,0,1):C.UVGC_306_83,(0,1,0):C.UVGC_307_84,(0,1,2):C.UVGC_307_85,(0,1,1):C.UVGC_307_86})

V_67 = CTVertex(name = 'V_67',
                type = 'UV',
                particles = [ P.b__tilde__, P.t, P.h__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS5 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_307_84,(0,0,2):C.UVGC_307_85,(0,0,1):C.UVGC_307_86,(0,1,0):C.UVGC_306_81,(0,1,2):C.UVGC_306_82,(0,1,1):C.UVGC_306_83})

V_68 = CTVertex(name = 'V_68',
                type = 'UV',
                particles = [ P.t__tilde__, P.b, P.G__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS5 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_305_78,(0,0,2):C.UVGC_305_79,(0,0,1):C.UVGC_305_80,(0,1,0):C.UVGC_309_88,(0,1,2):C.UVGC_309_89,(0,1,1):C.UVGC_309_90})

V_69 = CTVertex(name = 'V_69',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_284_33})

V_70 = CTVertex(name = 'V_70',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.G0 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_310_91})

V_71 = CTVertex(name = 'V_71',
                type = 'UV',
                particles = [ P.b__tilde__, P.t, P.G__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS3, L.FFS5 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_309_88,(0,0,2):C.UVGC_309_89,(0,0,1):C.UVGC_309_90,(0,1,0):C.UVGC_305_78,(0,1,2):C.UVGC_305_79,(0,1,1):C.UVGC_305_80})

V_72 = CTVertex(name = 'V_72',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_311_92})

V_73 = CTVertex(name = 'V_73',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_312_93})

V_74 = CTVertex(name = 'V_74',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_308_87})

V_75 = CTVertex(name = 'V_75',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.h1 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_285_34})

V_76 = CTVertex(name = 'V_76',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.h2 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS2 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_286_35})

V_77 = CTVertex(name = 'V_77',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.h3 ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFS1 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_287_36})

V_78 = CTVertex(name = 'V_78',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_304_77,(0,1,0):C.UVGC_299_71})

V_79 = CTVertex(name = 'V_79',
                type = 'UV',
                particles = [ P.u__tilde__, P.u, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.u] ], [ [P.t] ] ],
                couplings = {(0,0,4):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,1):C.UVGC_263_19,(0,1,2):C.UVGC_263_20,(0,1,3):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,4):C.UVGC_263_23,(0,2,0):C.UVGC_263_18,(0,2,1):C.UVGC_263_19,(0,2,2):C.UVGC_263_20,(0,2,3):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,4):C.UVGC_263_23})

V_80 = CTVertex(name = 'V_80',
                type = 'UV',
                particles = [ P.c__tilde__, P.c, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.c, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                couplings = {(0,0,2):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,1):C.UVGC_263_19,(0,1,3):C.UVGC_263_20,(0,1,4):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,2):C.UVGC_263_23,(0,2,0):C.UVGC_263_18,(0,2,1):C.UVGC_263_19,(0,2,3):C.UVGC_263_20,(0,2,4):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,2):C.UVGC_263_23})

V_81 = CTVertex(name = 'V_81',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.t] ], [ [P.t] ] ],
                couplings = {(0,0,4):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,1):C.UVGC_263_19,(0,1,2):C.UVGC_263_20,(0,1,3):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,4):C.UVGC_298_70,(0,2,0):C.UVGC_263_18,(0,2,1):C.UVGC_263_19,(0,2,2):C.UVGC_263_20,(0,2,3):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,4):C.UVGC_298_70})

V_82 = CTVertex(name = 'V_82',
                type = 'UV',
                particles = [ P.d__tilde__, P.u, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.UVGC_275_24,(0,0,1):C.UVGC_275_25})

V_83 = CTVertex(name = 'V_83',
                type = 'UV',
                particles = [ P.s__tilde__, P.c, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.UVGC_275_24,(0,0,1):C.UVGC_275_25})

V_84 = CTVertex(name = 'V_84',
                type = 'UV',
                particles = [ P.b__tilde__, P.t, P.W__minus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_300_72,(0,0,2):C.UVGC_300_73,(0,0,1):C.UVGC_275_25})

V_85 = CTVertex(name = 'V_85',
                type = 'UV',
                particles = [ P.t__tilde__, P.t, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_302_75,(0,1,0):C.UVGC_303_76})

V_86 = CTVertex(name = 'V_86',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.a ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_283_32,(0,1,0):C.UVGC_279_28})

V_87 = CTVertex(name = 'V_87',
                type = 'UV',
                particles = [ P.d__tilde__, P.d, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.d, P.g] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                couplings = {(0,0,2):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,1):C.UVGC_263_19,(0,1,3):C.UVGC_263_20,(0,1,4):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,2):C.UVGC_263_23,(0,2,0):C.UVGC_263_18,(0,2,1):C.UVGC_263_19,(0,2,3):C.UVGC_263_20,(0,2,4):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,2):C.UVGC_263_23})

V_88 = CTVertex(name = 'V_88',
                type = 'UV',
                particles = [ P.s__tilde__, P.s, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.g, P.s] ], [ [P.t] ] ],
                couplings = {(0,0,4):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,1):C.UVGC_263_19,(0,1,2):C.UVGC_263_20,(0,1,3):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,4):C.UVGC_263_23,(0,2,0):C.UVGC_263_18,(0,2,1):C.UVGC_263_19,(0,2,2):C.UVGC_263_20,(0,2,3):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,4):C.UVGC_263_23})

V_89 = CTVertex(name = 'V_89',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.g ],
                color = [ 'T(3,2,1)' ],
                lorentz = [ L.FFV1, L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b] ], [ [P.b, P.g] ], [ [P.c], [P.d], [P.s], [P.u] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                couplings = {(0,0,1):C.UVGC_261_16,(0,1,0):C.UVGC_263_18,(0,1,2):C.UVGC_263_19,(0,1,3):C.UVGC_263_20,(0,1,4):C.UVGC_263_21,(0,1,5):C.UVGC_263_22,(0,1,1):C.UVGC_278_27,(0,2,0):C.UVGC_263_18,(0,2,2):C.UVGC_263_19,(0,2,3):C.UVGC_263_20,(0,2,4):C.UVGC_263_21,(0,2,5):C.UVGC_263_22,(0,2,1):C.UVGC_278_27})

V_90 = CTVertex(name = 'V_90',
                type = 'UV',
                particles = [ P.u__tilde__, P.d, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.d, P.g], [P.g, P.u] ], [ [P.d, P.g, P.u] ] ],
                couplings = {(0,0,0):C.UVGC_275_24,(0,0,1):C.UVGC_275_25})

V_91 = CTVertex(name = 'V_91',
                type = 'UV',
                particles = [ P.c__tilde__, P.s, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.c, P.g], [P.g, P.s] ], [ [P.c, P.g, P.s] ] ],
                couplings = {(0,0,0):C.UVGC_275_24,(0,0,1):C.UVGC_275_25})

V_92 = CTVertex(name = 'V_92',
                type = 'UV',
                particles = [ P.t__tilde__, P.b, P.W__plus__ ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2 ],
                loop_particles = [ [ [P.b, P.g] ], [ [P.b, P.g, P.t] ], [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_300_72,(0,0,2):C.UVGC_300_73,(0,0,1):C.UVGC_275_25})

V_93 = CTVertex(name = 'V_93',
                type = 'UV',
                particles = [ P.b__tilde__, P.b, P.Z ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FFV2, L.FFV3 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_281_30,(0,1,0):C.UVGC_282_31})

V_94 = CTVertex(name = 'V_94',
                type = 'UV',
                particles = [ P.u__tilde__, P.u ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1, L.FF3, L.FF5 ],
                loop_particles = [ [ [P.g, P.u] ] ],
                couplings = {(0,0,0):C.UVGC_262_17,(0,1,0):C.UVGC_250_3,(0,2,0):C.UVGC_250_3})

V_95 = CTVertex(name = 'V_95',
                type = 'UV',
                particles = [ P.c__tilde__, P.c ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1, L.FF3, L.FF5 ],
                loop_particles = [ [ [P.c, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_262_17,(0,1,0):C.UVGC_250_3,(0,2,0):C.UVGC_250_3})

V_96 = CTVertex(name = 'V_96',
                type = 'UV',
                particles = [ P.t__tilde__, P.t ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                loop_particles = [ [ [P.g, P.t] ] ],
                couplings = {(0,0,0):C.UVGC_301_74,(0,2,0):C.UVGC_301_74,(0,1,0):C.UVGC_297_69,(0,3,0):C.UVGC_297_69})

V_97 = CTVertex(name = 'V_97',
                type = 'UV',
                particles = [ P.d__tilde__, P.d ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1, L.FF3, L.FF5 ],
                loop_particles = [ [ [P.d, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_262_17,(0,1,0):C.UVGC_250_3,(0,2,0):C.UVGC_250_3})

V_98 = CTVertex(name = 'V_98',
                type = 'UV',
                particles = [ P.s__tilde__, P.s ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF1, L.FF3, L.FF5 ],
                loop_particles = [ [ [P.g, P.s] ] ],
                couplings = {(0,0,0):C.UVGC_262_17,(0,1,0):C.UVGC_250_3,(0,2,0):C.UVGC_250_3})

V_99 = CTVertex(name = 'V_99',
                type = 'UV',
                particles = [ P.b__tilde__, P.b ],
                color = [ 'Identity(1,2)' ],
                lorentz = [ L.FF2, L.FF3, L.FF4, L.FF5 ],
                loop_particles = [ [ [P.b, P.g] ] ],
                couplings = {(0,0,0):C.UVGC_280_29,(0,2,0):C.UVGC_280_29,(0,1,0):C.UVGC_277_26,(0,3,0):C.UVGC_277_26})

V_100 = CTVertex(name = 'V_100',
                 type = 'UV',
                 particles = [ P.g, P.g ],
                 color = [ 'Identity(1,2)' ],
                 lorentz = [ L.VV1, L.VV3 ],
                 loop_particles = [ [ [P.b] ], [ [P.g] ], [ [P.ghG] ], [ [P.t] ] ],
                 couplings = {(0,0,0):C.UVGC_289_39,(0,0,1):C.UVGC_289_40,(0,0,2):C.UVGC_289_41,(0,0,3):C.UVGC_289_42,(0,1,0):C.UVGC_288_37,(0,1,3):C.UVGC_288_38})

