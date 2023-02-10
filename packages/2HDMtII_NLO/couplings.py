# This file was automatically created by FeynRules 2.3.13
# Mathematica version: 10.1.0  for Mac OS X x86 (64-bit) (March 24, 2015)
# Date: Tue 8 Dec 2015 11:44:24


from object_library import all_couplings, Coupling

from function_library import complexconjugate, re, im, csc, sec, acsc, asec, cot



GC_1 = Coupling(name = 'GC_1',
                value = '-G',
                order = {'QCD':1})

GC_10 = Coupling(name = 'GC_10',
                 value = '-(cw*complex(0,1)*gw)',
                 order = {'QED':1})

GC_100 = Coupling(name = 'GC_100',
                  value = '-2*complex(0,1)*l1*TH1x1*TH1x2 - complex(0,1)*l6*TH1x2*TH2x1 - complex(0,1)*l6*TH1x1*TH2x2 - complex(0,1)*l3*TH2x1*TH2x2',
                  order = {'QED':2})

GC_101 = Coupling(name = 'GC_101',
                  value = '-2*complex(0,1)*l1*TH1x1*TH1x2 - complex(0,1)*l6*TH1x2*TH2x1 - complex(0,1)*l6*TH1x1*TH2x2 - complex(0,1)*l3*TH2x1*TH2x2 - complex(0,1)*l4*TH2x1*TH2x2 + 2*complex(0,1)*l5*TH2x1*TH2x2',
                  order = {'QED':2})

GC_102 = Coupling(name = 'GC_102',
                  value = '-(complex(0,1)*l6*TH1x1*TH1x2) - (complex(0,1)*l4*TH1x2*TH2x1)/2. - complex(0,1)*l5*TH1x2*TH2x1 - (complex(0,1)*l4*TH1x1*TH2x2)/2. - complex(0,1)*l5*TH1x1*TH2x2 - complex(0,1)*l7*TH2x1*TH2x2',
                  order = {'QED':2})

GC_103 = Coupling(name = 'GC_103',
                  value = '(cw**2*complex(0,1)*gw**2*TH1x1*TH1x2)/2. + cw*complex(0,1)*g1*gw*sw*TH1x1*TH1x2 + (complex(0,1)*g1**2*sw**2*TH1x1*TH1x2)/2. + (cw**2*complex(0,1)*gw**2*TH2x1*TH2x2)/2. + cw*complex(0,1)*g1*gw*sw*TH2x1*TH2x2 + (complex(0,1)*g1**2*sw**2*TH2x1*TH2x2)/2.',
                  order = {'QED':2})

GC_104 = Coupling(name = 'GC_104',
                  value = '-(cw**2*complex(0,1)*g1*gw*TH1x1*TH1x2)/2. - (cw*complex(0,1)*g1**2*sw*TH1x1*TH1x2)/2. + (cw*complex(0,1)*gw**2*sw*TH1x1*TH1x2)/2. + (complex(0,1)*g1*gw*sw**2*TH1x1*TH1x2)/2. - (cw**2*complex(0,1)*g1*gw*TH2x1*TH2x2)/2. - (cw*complex(0,1)*g1**2*sw*TH2x1*TH2x2)/2. + (cw*complex(0,1)*gw**2*sw*TH2x1*TH2x2)/2. + (complex(0,1)*g1*gw*sw**2*TH2x1*TH2x2)/2.',
                  order = {'QED':2})

GC_105 = Coupling(name = 'GC_105',
                  value = '(cw**2*complex(0,1)*g1**2*TH1x1*TH1x2)/2. - cw*complex(0,1)*g1*gw*sw*TH1x1*TH1x2 + (complex(0,1)*gw**2*sw**2*TH1x1*TH1x2)/2. + (cw**2*complex(0,1)*g1**2*TH2x1*TH2x2)/2. - cw*complex(0,1)*g1*gw*sw*TH2x1*TH2x2 + (complex(0,1)*gw**2*sw**2*TH2x1*TH2x2)/2.',
                  order = {'QED':2})

GC_106 = Coupling(name = 'GC_106',
                  value = '-6*complex(0,1)*l1*TH1x1**3*TH1x2 - 9*complex(0,1)*l6*TH1x1**2*TH1x2*TH2x1 - 3*complex(0,1)*l3*TH1x1*TH1x2*TH2x1**2 - 3*complex(0,1)*l4*TH1x1*TH1x2*TH2x1**2 - 6*complex(0,1)*l5*TH1x1*TH1x2*TH2x1**2 - 3*complex(0,1)*l7*TH1x2*TH2x1**3 - 3*complex(0,1)*l6*TH1x1**3*TH2x2 - 3*complex(0,1)*l3*TH1x1**2*TH2x1*TH2x2 - 3*complex(0,1)*l4*TH1x1**2*TH2x1*TH2x2 - 6*complex(0,1)*l5*TH1x1**2*TH2x1*TH2x2 - 9*complex(0,1)*l7*TH1x1*TH2x1**2*TH2x2 - 6*complex(0,1)*l2*TH2x1**3*TH2x2',
                  order = {'QED':2})

GC_107 = Coupling(name = 'GC_107',
                  value = '(complex(0,1)*gw**2*TH1x2**2)/2. + (complex(0,1)*gw**2*TH2x2**2)/2.',
                  order = {'QED':2})

GC_108 = Coupling(name = 'GC_108',
                  value = '-(complex(0,1)*l3*TH1x2**2) - 2*complex(0,1)*l7*TH1x2*TH2x2 - 2*complex(0,1)*l2*TH2x2**2',
                  order = {'QED':2})

GC_109 = Coupling(name = 'GC_109',
                  value = '-2*complex(0,1)*l1*TH1x2**2 - 2*complex(0,1)*l6*TH1x2*TH2x2 - complex(0,1)*l3*TH2x2**2',
                  order = {'QED':2})

GC_11 = Coupling(name = 'GC_11',
                 value = 'cw*complex(0,1)*gw',
                 order = {'QED':1})

GC_110 = Coupling(name = 'GC_110',
                  value = '-2*complex(0,1)*l1*TH1x2**2 - 2*complex(0,1)*l6*TH1x2*TH2x2 - complex(0,1)*l3*TH2x2**2 - complex(0,1)*l4*TH2x2**2 + 2*complex(0,1)*l5*TH2x2**2',
                  order = {'QED':2})

GC_111 = Coupling(name = 'GC_111',
                  value = '-(complex(0,1)*l6*TH1x2**2) - complex(0,1)*l4*TH1x2*TH2x2 - 2*complex(0,1)*l5*TH1x2*TH2x2 - complex(0,1)*l7*TH2x2**2',
                  order = {'QED':2})

GC_112 = Coupling(name = 'GC_112',
                  value = '(cw**2*complex(0,1)*gw**2*TH1x2**2)/2. + cw*complex(0,1)*g1*gw*sw*TH1x2**2 + (complex(0,1)*g1**2*sw**2*TH1x2**2)/2. + (cw**2*complex(0,1)*gw**2*TH2x2**2)/2. + cw*complex(0,1)*g1*gw*sw*TH2x2**2 + (complex(0,1)*g1**2*sw**2*TH2x2**2)/2.',
                  order = {'QED':2})

GC_113 = Coupling(name = 'GC_113',
                  value = '-(cw**2*complex(0,1)*g1*gw*TH1x2**2)/2. - (cw*complex(0,1)*g1**2*sw*TH1x2**2)/2. + (cw*complex(0,1)*gw**2*sw*TH1x2**2)/2. + (complex(0,1)*g1*gw*sw**2*TH1x2**2)/2. - (cw**2*complex(0,1)*g1*gw*TH2x2**2)/2. - (cw*complex(0,1)*g1**2*sw*TH2x2**2)/2. + (cw*complex(0,1)*gw**2*sw*TH2x2**2)/2. + (complex(0,1)*g1*gw*sw**2*TH2x2**2)/2.',
                  order = {'QED':2})

GC_114 = Coupling(name = 'GC_114',
                  value = '(cw**2*complex(0,1)*g1**2*TH1x2**2)/2. - cw*complex(0,1)*g1*gw*sw*TH1x2**2 + (complex(0,1)*gw**2*sw**2*TH1x2**2)/2. + (cw**2*complex(0,1)*g1**2*TH2x2**2)/2. - cw*complex(0,1)*g1*gw*sw*TH2x2**2 + (complex(0,1)*gw**2*sw**2*TH2x2**2)/2.',
                  order = {'QED':2})

GC_115 = Coupling(name = 'GC_115',
                  value = '-6*complex(0,1)*l1*TH1x1**2*TH1x2**2 - 6*complex(0,1)*l6*TH1x1*TH1x2**2*TH2x1 - complex(0,1)*l3*TH1x2**2*TH2x1**2 - complex(0,1)*l4*TH1x2**2*TH2x1**2 - 2*complex(0,1)*l5*TH1x2**2*TH2x1**2 - 6*complex(0,1)*l6*TH1x1**2*TH1x2*TH2x2 - 4*complex(0,1)*l3*TH1x1*TH1x2*TH2x1*TH2x2 - 4*complex(0,1)*l4*TH1x1*TH1x2*TH2x1*TH2x2 - 8*complex(0,1)*l5*TH1x1*TH1x2*TH2x1*TH2x2 - 6*complex(0,1)*l7*TH1x2*TH2x1**2*TH2x2 - complex(0,1)*l3*TH1x1**2*TH2x2**2 - complex(0,1)*l4*TH1x1**2*TH2x2**2 - 2*complex(0,1)*l5*TH1x1**2*TH2x2**2 - 6*complex(0,1)*l7*TH1x1*TH2x1*TH2x2**2 - 6*complex(0,1)*l2*TH2x1**2*TH2x2**2',
                  order = {'QED':2})

GC_116 = Coupling(name = 'GC_116',
                  value = '-6*complex(0,1)*l1*TH1x1*TH1x2**3 - 3*complex(0,1)*l6*TH1x2**3*TH2x1 - 9*complex(0,1)*l6*TH1x1*TH1x2**2*TH2x2 - 3*complex(0,1)*l3*TH1x2**2*TH2x1*TH2x2 - 3*complex(0,1)*l4*TH1x2**2*TH2x1*TH2x2 - 6*complex(0,1)*l5*TH1x2**2*TH2x1*TH2x2 - 3*complex(0,1)*l3*TH1x1*TH1x2*TH2x2**2 - 3*complex(0,1)*l4*TH1x1*TH1x2*TH2x2**2 - 6*complex(0,1)*l5*TH1x1*TH1x2*TH2x2**2 - 9*complex(0,1)*l7*TH1x2*TH2x1*TH2x2**2 - 3*complex(0,1)*l7*TH1x1*TH2x2**3 - 6*complex(0,1)*l2*TH2x1*TH2x2**3',
                  order = {'QED':2})

GC_117 = Coupling(name = 'GC_117',
                  value = '-6*complex(0,1)*l1*TH1x2**4 - 12*complex(0,1)*l6*TH1x2**3*TH2x2 - 6*complex(0,1)*l3*TH1x2**2*TH2x2**2 - 6*complex(0,1)*l4*TH1x2**2*TH2x2**2 - 12*complex(0,1)*l5*TH1x2**2*TH2x2**2 - 12*complex(0,1)*l7*TH1x2*TH2x2**3 - 6*complex(0,1)*l2*TH2x2**4',
                  order = {'QED':2})

GC_118 = Coupling(name = 'GC_118',
                  value = '-(gw*TH3x3)/2.',
                  order = {'QED':1})

GC_119 = Coupling(name = 'GC_119',
                  value = '(gw*TH3x3)/2.',
                  order = {'QED':1})

GC_12 = Coupling(name = 'GC_12',
                 value = '-(cw*g1*gw)/2.',
                 order = {'QED':2})

GC_120 = Coupling(name = 'GC_120',
                  value = '-(cw*g1*gw*TH3x3)/2.',
                  order = {'QED':2})

GC_121 = Coupling(name = 'GC_121',
                  value = '(cw*g1*gw*TH3x3)/2.',
                  order = {'QED':2})

GC_122 = Coupling(name = 'GC_122',
                  value = '-(complex(0,1)*l6*TH3x3)',
                  order = {'QED':2})

GC_123 = Coupling(name = 'GC_123',
                  value = '-3*complex(0,1)*l6*TH3x3',
                  order = {'QED':2})

GC_124 = Coupling(name = 'GC_124',
                  value = '-(complex(0,1)*l7*TH3x3)',
                  order = {'QED':2})

GC_125 = Coupling(name = 'GC_125',
                  value = '-(g1*gw*sw*TH3x3)/2.',
                  order = {'QED':2})

GC_126 = Coupling(name = 'GC_126',
                  value = '(g1*gw*sw*TH3x3)/2.',
                  order = {'QED':2})

GC_127 = Coupling(name = 'GC_127',
                  value = '(complex(0,1)*gw**2*TH3x3**2)/2.',
                  order = {'QED':2})

GC_128 = Coupling(name = 'GC_128',
                  value = '-2*complex(0,1)*l2*TH3x3**2',
                  order = {'QED':2})

GC_129 = Coupling(name = 'GC_129',
                  value = '-(complex(0,1)*l3*TH3x3**2)',
                  order = {'QED':2})

GC_13 = Coupling(name = 'GC_13',
                 value = '(cw*g1*gw)/2.',
                 order = {'QED':2})

GC_130 = Coupling(name = 'GC_130',
                  value = '-(complex(0,1)*l7*TH3x3**2)',
                  order = {'QED':2})

GC_131 = Coupling(name = 'GC_131',
                  value = '-3*complex(0,1)*l7*TH3x3**3',
                  order = {'QED':2})

GC_132 = Coupling(name = 'GC_132',
                  value = '-6*complex(0,1)*l2*TH3x3**4',
                  order = {'QED':2})

GC_133 = Coupling(name = 'GC_133',
                  value = '-(complex(0,1)*l4*TH3x3)/2. - complex(0,1)*l5*TH3x3',
                  order = {'QED':2})

GC_134 = Coupling(name = 'GC_134',
                  value = '(l4*TH1x1*TH3x3)/2. - l5*TH1x1*TH3x3',
                  order = {'QED':2})

GC_135 = Coupling(name = 'GC_135',
                  value = '-(l4*TH1x1*TH3x3)/2. + l5*TH1x1*TH3x3',
                  order = {'QED':2})

GC_136 = Coupling(name = 'GC_136',
                  value = '(l4*TH1x2*TH3x3)/2. - l5*TH1x2*TH3x3',
                  order = {'QED':2})

GC_137 = Coupling(name = 'GC_137',
                  value = '-(l4*TH1x2*TH3x3)/2. + l5*TH1x2*TH3x3',
                  order = {'QED':2})

GC_138 = Coupling(name = 'GC_138',
                  value = '-(cw*gw*TH2x1*TH3x3)/2. - (g1*sw*TH2x1*TH3x3)/2.',
                  order = {'QED':1})

GC_139 = Coupling(name = 'GC_139',
                  value = '(cw*gw*TH2x1*TH3x3)/2. + (g1*sw*TH2x1*TH3x3)/2.',
                  order = {'QED':1})

GC_14 = Coupling(name = 'GC_14',
                 value = '(complex(0,1)*gw**2)/2.',
                 order = {'QED':2})

GC_140 = Coupling(name = 'GC_140',
                  value = '-(cw*g1*TH2x1*TH3x3)/2. + (gw*sw*TH2x1*TH3x3)/2.',
                  order = {'QED':1})

GC_141 = Coupling(name = 'GC_141',
                  value = '-(complex(0,1)*l6*TH1x1**2*TH3x3) - 4*complex(0,1)*l5*TH1x1*TH2x1*TH3x3 - complex(0,1)*l7*TH2x1**2*TH3x3',
                  order = {'QED':2})

GC_142 = Coupling(name = 'GC_142',
                  value = '-(cw*gw*TH2x2*TH3x3)/2. - (g1*sw*TH2x2*TH3x3)/2.',
                  order = {'QED':1})

GC_143 = Coupling(name = 'GC_143',
                  value = '(cw*gw*TH2x2*TH3x3)/2. + (g1*sw*TH2x2*TH3x3)/2.',
                  order = {'QED':1})

GC_144 = Coupling(name = 'GC_144',
                  value = '-(cw*g1*TH2x2*TH3x3)/2. + (gw*sw*TH2x2*TH3x3)/2.',
                  order = {'QED':1})

GC_145 = Coupling(name = 'GC_145',
                  value = '-(complex(0,1)*l6*TH1x1*TH1x2*TH3x3) - 2*complex(0,1)*l5*TH1x2*TH2x1*TH3x3 - 2*complex(0,1)*l5*TH1x1*TH2x2*TH3x3 - complex(0,1)*l7*TH2x1*TH2x2*TH3x3',
                  order = {'QED':2})

GC_146 = Coupling(name = 'GC_146',
                  value = '-(complex(0,1)*l6*TH1x2**2*TH3x3) - 4*complex(0,1)*l5*TH1x2*TH2x2*TH3x3 - complex(0,1)*l7*TH2x2**2*TH3x3',
                  order = {'QED':2})

GC_147 = Coupling(name = 'GC_147',
                  value = '-(complex(0,1)*l3*TH3x3**2) - complex(0,1)*l4*TH3x3**2 - 2*complex(0,1)*l5*TH3x3**2',
                  order = {'QED':2})

GC_148 = Coupling(name = 'GC_148',
                  value = '(cw**2*complex(0,1)*gw**2*TH3x3**2)/2. + cw*complex(0,1)*g1*gw*sw*TH3x3**2 + (complex(0,1)*g1**2*sw**2*TH3x3**2)/2.',
                  order = {'QED':2})

GC_149 = Coupling(name = 'GC_149',
                  value = '-(cw**2*complex(0,1)*g1*gw*TH3x3**2)/2. - (cw*complex(0,1)*g1**2*sw*TH3x3**2)/2. + (cw*complex(0,1)*gw**2*sw*TH3x3**2)/2. + (complex(0,1)*g1*gw*sw**2*TH3x3**2)/2.',
                  order = {'QED':2})

GC_15 = Coupling(name = 'GC_15',
                 value = '-(complex(0,1)*gw**2)',
                 order = {'QED':2})

GC_150 = Coupling(name = 'GC_150',
                  value = '(cw**2*complex(0,1)*g1**2*TH3x3**2)/2. - cw*complex(0,1)*g1*gw*sw*TH3x3**2 + (complex(0,1)*gw**2*sw**2*TH3x3**2)/2.',
                  order = {'QED':2})

GC_151 = Coupling(name = 'GC_151',
                  value = '-(complex(0,1)*l3*TH1x1**2*TH3x3**2) - complex(0,1)*l4*TH1x1**2*TH3x3**2 + 2*complex(0,1)*l5*TH1x1**2*TH3x3**2 - 2*complex(0,1)*l7*TH1x1*TH2x1*TH3x3**2 - 2*complex(0,1)*l2*TH2x1**2*TH3x3**2',
                  order = {'QED':2})

GC_152 = Coupling(name = 'GC_152',
                  value = '-(complex(0,1)*l3*TH1x1*TH1x2*TH3x3**2) - complex(0,1)*l4*TH1x1*TH1x2*TH3x3**2 + 2*complex(0,1)*l5*TH1x1*TH1x2*TH3x3**2 - complex(0,1)*l7*TH1x2*TH2x1*TH3x3**2 - complex(0,1)*l7*TH1x1*TH2x2*TH3x3**2 - 2*complex(0,1)*l2*TH2x1*TH2x2*TH3x3**2',
                  order = {'QED':2})

GC_153 = Coupling(name = 'GC_153',
                  value = '-(complex(0,1)*l3*TH1x2**2*TH3x3**2) - complex(0,1)*l4*TH1x2**2*TH3x3**2 + 2*complex(0,1)*l5*TH1x2**2*TH3x3**2 - 2*complex(0,1)*l7*TH1x2*TH2x2*TH3x3**2 - 2*complex(0,1)*l2*TH2x2**2*TH3x3**2',
                  order = {'QED':2})

GC_154 = Coupling(name = 'GC_154',
                  value = '(cw*complex(0,1)*g1*gw*vev)/2.',
                  order = {'QED':1})

GC_155 = Coupling(name = 'GC_155',
                  value = '-(gw**2*vev)/4.',
                  order = {'QED':1})

GC_156 = Coupling(name = 'GC_156',
                  value = '(gw**2*vev)/4.',
                  order = {'QED':1})

GC_157 = Coupling(name = 'GC_157',
                  value = '-(complex(0,1)*g1*gw*sw*vev)/2.',
                  order = {'QED':1})

GC_158 = Coupling(name = 'GC_158',
                  value = '-(complex(0,1)*gw**2*TH1x1*vev)/4.',
                  order = {'QED':1})

##################Coupling for WWh1################################################################
GC_159 = Coupling(name = 'GC_159',
                  value = '(complex(0,1)*gw**2*TH1x1*vev)/2.',
                  order = {'QED':1})
###################################################################################################
GC_16 = Coupling(name = 'GC_16',
                 value = '2*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_160 = Coupling(name = 'GC_160',
                  value = '-(complex(0,1)*gw**2*TH1x2*vev)/4.',
                  order = {'QED':1})

#################Coupling for WWh2#################################################################
GC_161 = Coupling(name = 'GC_161',
                  value = '(complex(0,1)*gw**2*TH1x2*vev)/2.',
                  order = {'QED':1})
##################################################################################################

GC_162 = Coupling(name = 'GC_162',
                  value = '-(cw*gw**2*vev)/4. - (g1*gw*sw*vev)/4.',
                  order = {'QED':1})

GC_163 = Coupling(name = 'GC_163',
                  value = '(cw*gw**2*vev)/4. - (g1*gw*sw*vev)/4.',
                  order = {'QED':1})

GC_164 = Coupling(name = 'GC_164',
                  value = '-(cw*gw**2*vev)/4. + (g1*gw*sw*vev)/4.',
                  order = {'QED':1})

GC_165 = Coupling(name = 'GC_165',
                  value = '(cw*gw**2*vev)/4. + (g1*gw*sw*vev)/4.',
                  order = {'QED':1})

GC_166 = Coupling(name = 'GC_166',
                  value = '-(cw*g1*gw*vev)/4. - (gw**2*sw*vev)/4.',
                  order = {'QED':1})

GC_167 = Coupling(name = 'GC_167',
                  value = '(cw*g1*gw*vev)/4. - (gw**2*sw*vev)/4.',
                  order = {'QED':1})

GC_168 = Coupling(name = 'GC_168',
                  value = '-(cw*g1*gw*vev)/4. + (gw**2*sw*vev)/4.',
                  order = {'QED':1})

GC_169 = Coupling(name = 'GC_169',
                  value = '(cw*g1*gw*vev)/4. + (gw**2*sw*vev)/4.',
                  order = {'QED':1})

GC_17 = Coupling(name = 'GC_17',
                 value = 'cw**2*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_170 = Coupling(name = 'GC_170',
                  value = '-(cw**2*complex(0,1)*gw**2*TH1x1*vev)/4. - (cw*complex(0,1)*g1*gw*sw*TH1x1*vev)/2. - (complex(0,1)*g1**2*sw**2*TH1x1*vev)/4.',
                  order = {'QED':1})

GC_171 = Coupling(name = 'GC_171',
                  value = '(cw**2*complex(0,1)*gw**2*TH1x1*vev)/2. + cw*complex(0,1)*g1*gw*sw*TH1x1*vev + (complex(0,1)*g1**2*sw**2*TH1x1*vev)/2.',
                  order = {'QED':1})

GC_172 = Coupling(name = 'GC_172',
                  value = '(cw**2*complex(0,1)*g1*gw*TH1x1*vev)/4. + (cw*complex(0,1)*g1**2*sw*TH1x1*vev)/4. - (cw*complex(0,1)*gw**2*sw*TH1x1*vev)/4. - (complex(0,1)*g1*gw*sw**2*TH1x1*vev)/4.',
                  order = {'QED':1})

GC_173 = Coupling(name = 'GC_173',
                  value = '-(cw**2*complex(0,1)*g1*gw*TH1x1*vev)/2. - (cw*complex(0,1)*g1**2*sw*TH1x1*vev)/2. + (cw*complex(0,1)*gw**2*sw*TH1x1*vev)/2. + (complex(0,1)*g1*gw*sw**2*TH1x1*vev)/2.',
                  order = {'QED':1})

GC_174 = Coupling(name = 'GC_174',
                  value = '-(cw**2*complex(0,1)*g1**2*TH1x1*vev)/4. + (cw*complex(0,1)*g1*gw*sw*TH1x1*vev)/2. - (complex(0,1)*gw**2*sw**2*TH1x1*vev)/4.',
                  order = {'QED':1})

GC_175 = Coupling(name = 'GC_175',
                  value = '(cw**2*complex(0,1)*g1**2*TH1x1*vev)/2. - cw*complex(0,1)*g1*gw*sw*TH1x1*vev + (complex(0,1)*gw**2*sw**2*TH1x1*vev)/2.',
                  order = {'QED':1})

GC_176 = Coupling(name = 'GC_176',
                  value = '-(cw**2*complex(0,1)*gw**2*TH1x2*vev)/4. - (cw*complex(0,1)*g1*gw*sw*TH1x2*vev)/2. - (complex(0,1)*g1**2*sw**2*TH1x2*vev)/4.',
                  order = {'QED':1})

GC_177 = Coupling(name = 'GC_177',
                  value = '(cw**2*complex(0,1)*gw**2*TH1x2*vev)/2. + cw*complex(0,1)*g1*gw*sw*TH1x2*vev + (complex(0,1)*g1**2*sw**2*TH1x2*vev)/2.',
                  order = {'QED':1})

GC_178 = Coupling(name = 'GC_178',
                  value = '(cw**2*complex(0,1)*g1*gw*TH1x2*vev)/4. + (cw*complex(0,1)*g1**2*sw*TH1x2*vev)/4. - (cw*complex(0,1)*gw**2*sw*TH1x2*vev)/4. - (complex(0,1)*g1*gw*sw**2*TH1x2*vev)/4.',
                  order = {'QED':1})

GC_179 = Coupling(name = 'GC_179',
                  value = '-(cw**2*complex(0,1)*g1*gw*TH1x2*vev)/2. - (cw*complex(0,1)*g1**2*sw*TH1x2*vev)/2. + (cw*complex(0,1)*gw**2*sw*TH1x2*vev)/2. + (complex(0,1)*g1*gw*sw**2*TH1x2*vev)/2.',
                  order = {'QED':1})

GC_18 = Coupling(name = 'GC_18',
                 value = '-2*cw**2*complex(0,1)*gw**2',
                 order = {'QED':2})

GC_180 = Coupling(name = 'GC_180',
                  value = '-(cw**2*complex(0,1)*g1**2*TH1x2*vev)/4. + (cw*complex(0,1)*g1*gw*sw*TH1x2*vev)/2. - (complex(0,1)*gw**2*sw**2*TH1x2*vev)/4.',
                  order = {'QED':1})

GC_181 = Coupling(name = 'GC_181',
                  value = '(cw**2*complex(0,1)*g1**2*TH1x2*vev)/2. - cw*complex(0,1)*g1*gw*sw*TH1x2*vev + (complex(0,1)*gw**2*sw**2*TH1x2*vev)/2.',
                  order = {'QED':1})

GC_182 = Coupling(name = 'GC_182',
                  value = '-(complex(0,1)*l6*TH1x1*vev) - (complex(0,1)*l4*TH2x1*vev)/2. - complex(0,1)*l5*TH2x1*vev',
                  order = {'QED':1})

GC_183 = Coupling(name = 'GC_183',
                  value = '-2*complex(0,1)*l1*TH1x1*vev - complex(0,1)*l6*TH2x1*vev',
                  order = {'QED':1})

GC_184 = Coupling(name = 'GC_184',
                  value = '-(complex(0,1)*l3*TH1x1*vev) - complex(0,1)*l7*TH2x1*vev',
                  order = {'QED':1})

GC_185 = Coupling(name = 'GC_185',
                  value = '-6*complex(0,1)*l1*TH1x1**3*vev - 9*complex(0,1)*l6*TH1x1**2*TH2x1*vev - 3*complex(0,1)*l3*TH1x1*TH2x1**2*vev - 3*complex(0,1)*l4*TH1x1*TH2x1**2*vev - 6*complex(0,1)*l5*TH1x1*TH2x1**2*vev - 3*complex(0,1)*l7*TH2x1**3*vev',
                  order = {'QED':1})

GC_186 = Coupling(name = 'GC_186',
                  value = '-(complex(0,1)*l6*TH1x2*vev) - (complex(0,1)*l4*TH2x2*vev)/2. - complex(0,1)*l5*TH2x2*vev',
                  order = {'QED':1})

GC_187 = Coupling(name = 'GC_187',
                  value = '-2*complex(0,1)*l1*TH1x2*vev - complex(0,1)*l6*TH2x2*vev',
                  order = {'QED':1})

GC_188 = Coupling(name = 'GC_188',
                  value = '-(complex(0,1)*l3*TH1x2*vev) - complex(0,1)*l7*TH2x2*vev',
                  order = {'QED':1})

GC_189 = Coupling(name = 'GC_189',
                  value = '-6*complex(0,1)*l1*TH1x1**2*TH1x2*vev - 6*complex(0,1)*l6*TH1x1*TH1x2*TH2x1*vev - complex(0,1)*l3*TH1x2*TH2x1**2*vev - complex(0,1)*l4*TH1x2*TH2x1**2*vev - 2*complex(0,1)*l5*TH1x2*TH2x1**2*vev - 3*complex(0,1)*l6*TH1x1**2*TH2x2*vev - 2*complex(0,1)*l3*TH1x1*TH2x1*TH2x2*vev - 2*complex(0,1)*l4*TH1x1*TH2x1*TH2x2*vev - 4*complex(0,1)*l5*TH1x1*TH2x1*TH2x2*vev - 3*complex(0,1)*l7*TH2x1**2*TH2x2*vev',
                  order = {'QED':1})

GC_19 = Coupling(name = 'GC_19',
                 value = '-(complex(0,1)*I1a33)',
                 order = {'QED':1})

GC_190 = Coupling(name = 'GC_190',
                  value = '-6*complex(0,1)*l1*TH1x1*TH1x2**2*vev - 3*complex(0,1)*l6*TH1x2**2*TH2x1*vev - 6*complex(0,1)*l6*TH1x1*TH1x2*TH2x2*vev - 2*complex(0,1)*l3*TH1x2*TH2x1*TH2x2*vev - 2*complex(0,1)*l4*TH1x2*TH2x1*TH2x2*vev - 4*complex(0,1)*l5*TH1x2*TH2x1*TH2x2*vev - complex(0,1)*l3*TH1x1*TH2x2**2*vev - complex(0,1)*l4*TH1x1*TH2x2**2*vev - 2*complex(0,1)*l5*TH1x1*TH2x2**2*vev - 3*complex(0,1)*l7*TH2x1*TH2x2**2*vev',
                  order = {'QED':1})

GC_191 = Coupling(name = 'GC_191',
                  value = '-6*complex(0,1)*l1*TH1x2**3*vev - 9*complex(0,1)*l6*TH1x2**2*TH2x2*vev - 3*complex(0,1)*l3*TH1x2*TH2x2**2*vev - 3*complex(0,1)*l4*TH1x2*TH2x2**2*vev - 6*complex(0,1)*l5*TH1x2*TH2x2**2*vev - 3*complex(0,1)*l7*TH2x2**3*vev',
                  order = {'QED':1})

GC_192 = Coupling(name = 'GC_192',
                  value = '(l4*TH3x3*vev)/2. - l5*TH3x3*vev',
                  order = {'QED':1})

GC_193 = Coupling(name = 'GC_193',
                  value = '-(l4*TH3x3*vev)/2. + l5*TH3x3*vev',
                  order = {'QED':1})

GC_194 = Coupling(name = 'GC_194',
                  value = '-(complex(0,1)*l6*TH1x1*TH3x3*vev) - 2*complex(0,1)*l5*TH2x1*TH3x3*vev',
                  order = {'QED':1})

GC_195 = Coupling(name = 'GC_195',
                  value = '-(complex(0,1)*l6*TH1x2*TH3x3*vev) - 2*complex(0,1)*l5*TH2x2*TH3x3*vev',
                  order = {'QED':1})

GC_196 = Coupling(name = 'GC_196',
                  value = '-(complex(0,1)*l3*TH1x1*TH3x3**2*vev) - complex(0,1)*l4*TH1x1*TH3x3**2*vev + 2*complex(0,1)*l5*TH1x1*TH3x3**2*vev - complex(0,1)*l7*TH2x1*TH3x3**2*vev',
                  order = {'QED':1})

GC_197 = Coupling(name = 'GC_197',
                  value = '-(complex(0,1)*l3*TH1x2*TH3x3**2*vev) - complex(0,1)*l4*TH1x2*TH3x3**2*vev + 2*complex(0,1)*l5*TH1x2*TH3x3**2*vev - complex(0,1)*l7*TH2x2*TH3x3**2*vev',
                  order = {'QED':1})

GC_198 = Coupling(name = 'GC_198',
                  value = '-(yb/cmath.sqrt(2))',
                  order = {'QED':1})

##################Coupling for bbh3########################################################################
GC_199 = Coupling(name = 'GC_199',
                  value = '-((tanbeta*TH3x3*ymb)/vev)',
                  order = {'QED':1})
##########################################################################################################

GC_2 = Coupling(name = 'GC_2',
                value = 'complex(0,1)*G',
                order = {'QCD':1})

GC_20 = Coupling(name = 'GC_20',
                 value = 'complex(0,1)*I2a33',
                 order = {'QED':1})

#################Coupling for bbh1#########################################################################
GC_200 = Coupling(name = 'GC_200',
                  value = '-((complex(0,1)*tanbeta*TH2x1*ymb)/vev) - (complex(0,1)*TH1x1*yb)/cmath.sqrt(2)',
                  order = {'QED':1})
###########################################################################################################
#################Coupling for bbh2#########################################################################
GC_201 = Coupling(name = 'GC_201',
                  value = '-((complex(0,1)*tanbeta*TH2x2*ymb)/vev) - (complex(0,1)*TH1x2*yb)/cmath.sqrt(2)',
                  order = {'QED':1})
###########################################################################################################

#################Coupling for tth3#################################################################################
GC_202 = Coupling(name = 'GC_202',
                  value = '-((TH3x3*ymt)/(tanbeta*vev))',
                  order = {'QED':1})
###################################################################################################################

GC_203 = Coupling(name = 'GC_203',
                  value = '-((complex(0,1)*tanbeta*ymtau*cmath.sqrt(2))/vev)',
                  order = {'QED':1})

GC_204 = Coupling(name = 'GC_204',
                  value = '-((tanbeta*TH3x3*ymtau)/vev)',
                  order = {'QED':1})

GC_205 = Coupling(name = 'GC_205',
                  value = 'yt/cmath.sqrt(2)',
                  order = {'QED':1})

#################Coupling for tth1##################################################################################
GC_206 = Coupling(name = 'GC_206',
                  value = '(complex(0,1)*TH2x1*ymt)/(tanbeta*vev) - (complex(0,1)*TH1x1*yt)/cmath.sqrt(2)',
                  order = {'QED':1})
####################################################################################################################
#################Coupling for tth2#################################################################################
GC_207 = Coupling(name = 'GC_207',
                  value = '(complex(0,1)*TH2x2*ymt)/(tanbeta*vev) - (complex(0,1)*TH1x2*yt)/cmath.sqrt(2)',
                  order = {'QED':1})
####################################################################################################################

GC_208 = Coupling(name = 'GC_208',
                  value = '-(complex(0,1)*ytau)',
                  order = {'QED':1})

GC_209 = Coupling(name = 'GC_209',
                  value = '-(ytau/cmath.sqrt(2))',
                  order = {'QED':1})

GC_21 = Coupling(name = 'GC_21',
                 value = 'complex(0,1)*I3a33',
                 order = {'QED':1})

GC_210 = Coupling(name = 'GC_210',
                  value = '-((complex(0,1)*tanbeta*TH2x1*ymtau)/vev) - (complex(0,1)*TH1x1*ytau)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_211 = Coupling(name = 'GC_211',
                  value = '-((complex(0,1)*tanbeta*TH2x2*ymtau)/vev) - (complex(0,1)*TH1x2*ytau)/cmath.sqrt(2)',
                  order = {'QED':1})

GC_22 = Coupling(name = 'GC_22',
                 value = '-(complex(0,1)*I4a33)',
                 order = {'QED':1})

GC_23 = Coupling(name = 'GC_23',
                 value = '-(complex(0,1)*I5a33)',
                 order = {'QED':1})

GC_24 = Coupling(name = 'GC_24',
                 value = 'complex(0,1)*I6a33',
                 order = {'QED':1})

GC_25 = Coupling(name = 'GC_25',
                 value = 'complex(0,1)*I7a33',
                 order = {'QED':1})

GC_26 = Coupling(name = 'GC_26',
                 value = '-(complex(0,1)*I8a33)',
                 order = {'QED':1})

GC_27 = Coupling(name = 'GC_27',
                 value = '-2*complex(0,1)*l1',
                 order = {'QED':2})

GC_28 = Coupling(name = 'GC_28',
                 value = '-4*complex(0,1)*l1',
                 order = {'QED':2})

GC_29 = Coupling(name = 'GC_29',
                 value = '-6*complex(0,1)*l1',
                 order = {'QED':2})

GC_3 = Coupling(name = 'GC_3',
                value = 'complex(0,1)*G**2',
                order = {'QCD':2})

GC_30 = Coupling(name = 'GC_30',
                 value = '-4*complex(0,1)*l2',
                 order = {'QED':2})

GC_31 = Coupling(name = 'GC_31',
                 value = '-(complex(0,1)*l3)',
                 order = {'QED':2})

GC_32 = Coupling(name = 'GC_32',
                 value = '-(complex(0,1)*l3) - complex(0,1)*l4',
                 order = {'QED':2})

GC_33 = Coupling(name = 'GC_33',
                 value = '-4*complex(0,1)*l5',
                 order = {'QED':2})

GC_34 = Coupling(name = 'GC_34',
                 value = '-(complex(0,1)*l6)',
                 order = {'QED':2})

GC_35 = Coupling(name = 'GC_35',
                 value = '-2*complex(0,1)*l6',
                 order = {'QED':2})

GC_36 = Coupling(name = 'GC_36',
                 value = '-2*complex(0,1)*l7',
                 order = {'QED':2})

GC_37 = Coupling(name = 'GC_37',
                 value = '(complex(0,1)*g1*sw)/3.',
                 order = {'QED':1})

GC_38 = Coupling(name = 'GC_38',
                 value = '(-2*complex(0,1)*g1*sw)/3.',
                 order = {'QED':1})

GC_39 = Coupling(name = 'GC_39',
                 value = 'complex(0,1)*g1*sw',
                 order = {'QED':1})

GC_4 = Coupling(name = 'GC_4',
                value = '-(cw*complex(0,1)*g1)/3.',
                order = {'QED':1})

GC_40 = Coupling(name = 'GC_40',
                 value = '-(complex(0,1)*gw*sw)',
                 order = {'QED':1})

GC_41 = Coupling(name = 'GC_41',
                 value = 'complex(0,1)*gw*sw',
                 order = {'QED':1})

GC_42 = Coupling(name = 'GC_42',
                 value = '-(g1*gw*sw)/2.',
                 order = {'QED':2})

GC_43 = Coupling(name = 'GC_43',
                 value = '(g1*gw*sw)/2.',
                 order = {'QED':2})

GC_44 = Coupling(name = 'GC_44',
                 value = 'cw*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_45 = Coupling(name = 'GC_45',
                 value = '-2*cw*complex(0,1)*gw**2*sw',
                 order = {'QED':2})

GC_46 = Coupling(name = 'GC_46',
                 value = 'complex(0,1)*gw**2*sw**2',
                 order = {'QED':2})

GC_47 = Coupling(name = 'GC_47',
                 value = '-2*complex(0,1)*gw**2*sw**2',
                 order = {'QED':2})

GC_48 = Coupling(name = 'GC_48',
                 value = '-(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw)/6.',
                 order = {'QED':1})

GC_49 = Coupling(name = 'GC_49',
                 value = '(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw)/6.',
                 order = {'QED':1})

GC_5 = Coupling(name = 'GC_5',
                value = '(2*cw*complex(0,1)*g1)/3.',
                order = {'QED':1})

GC_50 = Coupling(name = 'GC_50',
                 value = '(cw*complex(0,1)*gw)/2. - (complex(0,1)*g1*sw)/2.',
                 order = {'QED':1})

GC_51 = Coupling(name = 'GC_51',
                 value = '-(cw*complex(0,1)*gw)/2. + (complex(0,1)*g1*sw)/2.',
                 order = {'QED':1})

GC_52 = Coupling(name = 'GC_52',
                 value = '(cw*complex(0,1)*gw)/2. + (complex(0,1)*g1*sw)/2.',
                 order = {'QED':1})

GC_53 = Coupling(name = 'GC_53',
                 value = '(cw*complex(0,1)*g1)/6. - (complex(0,1)*gw*sw)/2.',
                 order = {'QED':1})

GC_54 = Coupling(name = 'GC_54',
                 value = '-(cw*complex(0,1)*g1)/2. - (complex(0,1)*gw*sw)/2.',
                 order = {'QED':1})

GC_55 = Coupling(name = 'GC_55',
                 value = '(cw*complex(0,1)*g1)/6. + (complex(0,1)*gw*sw)/2.',
                 order = {'QED':1})

GC_56 = Coupling(name = 'GC_56',
                 value = '-(cw*complex(0,1)*g1)/2. + (complex(0,1)*gw*sw)/2.',
                 order = {'QED':1})

GC_57 = Coupling(name = 'GC_57',
                 value = '(cw**2*complex(0,1)*gw**2)/2. - cw*complex(0,1)*g1*gw*sw + (complex(0,1)*g1**2*sw**2)/2.',
                 order = {'QED':2})

GC_58 = Coupling(name = 'GC_58',
                 value = '(cw**2*complex(0,1)*gw**2)/2. + cw*complex(0,1)*g1*gw*sw + (complex(0,1)*g1**2*sw**2)/2.',
                 order = {'QED':2})

GC_59 = Coupling(name = 'GC_59',
                 value = '(cw**2*complex(0,1)*g1*gw)/2. - (cw*complex(0,1)*g1**2*sw)/2. + (cw*complex(0,1)*gw**2*sw)/2. - (complex(0,1)*g1*gw*sw**2)/2.',
                 order = {'QED':2})

GC_6 = Coupling(name = 'GC_6',
                value = '-(cw*complex(0,1)*g1)',
                order = {'QED':1})

GC_60 = Coupling(name = 'GC_60',
                 value = '-(cw**2*complex(0,1)*g1*gw)/2. - (cw*complex(0,1)*g1**2*sw)/2. + (cw*complex(0,1)*gw**2*sw)/2. + (complex(0,1)*g1*gw*sw**2)/2.',
                 order = {'QED':2})

GC_61 = Coupling(name = 'GC_61',
                 value = '(cw**2*complex(0,1)*g1**2)/2. - cw*complex(0,1)*g1*gw*sw + (complex(0,1)*gw**2*sw**2)/2.',
                 order = {'QED':2})

GC_62 = Coupling(name = 'GC_62',
                 value = '(cw**2*complex(0,1)*g1**2)/2. + cw*complex(0,1)*g1*gw*sw + (complex(0,1)*gw**2*sw**2)/2.',
                 order = {'QED':2})

GC_63 = Coupling(name = 'GC_63',
                 value = '-(complex(0,1)*gw*TH1x1)/2.',
                 order = {'QED':1})

GC_64 = Coupling(name = 'GC_64',
                 value = '(complex(0,1)*gw*TH1x1)/2.',
                 order = {'QED':1})

GC_65 = Coupling(name = 'GC_65',
                 value = '(cw*complex(0,1)*g1*gw*TH1x1)/2.',
                 order = {'QED':2})

GC_66 = Coupling(name = 'GC_66',
                 value = '-(complex(0,1)*g1*gw*sw*TH1x1)/2.',
                 order = {'QED':2})

GC_67 = Coupling(name = 'GC_67',
                 value = '-(cw*gw*TH1x1)/2. - (g1*sw*TH1x1)/2.',
                 order = {'QED':1})

GC_68 = Coupling(name = 'GC_68',
                 value = '(cw*gw*TH1x1)/2. + (g1*sw*TH1x1)/2.',
                 order = {'QED':1})

GC_69 = Coupling(name = 'GC_69',
                 value = '(cw*g1*TH1x1)/2. - (gw*sw*TH1x1)/2.',
                 order = {'QED':1})

GC_7 = Coupling(name = 'GC_7',
                value = '-gw/2.',
                order = {'QED':1})

GC_70 = Coupling(name = 'GC_70',
                 value = '-(complex(0,1)*gw*TH1x2)/2.',
                 order = {'QED':1})

GC_71 = Coupling(name = 'GC_71',
                 value = '(complex(0,1)*gw*TH1x2)/2.',
                 order = {'QED':1})

GC_72 = Coupling(name = 'GC_72',
                 value = '(cw*complex(0,1)*g1*gw*TH1x2)/2.',
                 order = {'QED':2})

GC_73 = Coupling(name = 'GC_73',
                 value = '-(complex(0,1)*g1*gw*sw*TH1x2)/2.',
                 order = {'QED':2})

GC_74 = Coupling(name = 'GC_74',
                 value = '-(cw*gw*TH1x2)/2. - (g1*sw*TH1x2)/2.',
                 order = {'QED':1})

GC_75 = Coupling(name = 'GC_75',
                 value = '(cw*gw*TH1x2)/2. + (g1*sw*TH1x2)/2.',
                 order = {'QED':1})

GC_76 = Coupling(name = 'GC_76',
                 value = '(cw*g1*TH1x2)/2. - (gw*sw*TH1x2)/2.',
                 order = {'QED':1})

GC_77 = Coupling(name = 'GC_77',
                 value = '-(complex(0,1)*gw*TH2x1)/2.',
                 order = {'QED':1})

GC_78 = Coupling(name = 'GC_78',
                 value = '(complex(0,1)*gw*TH2x1)/2.',
                 order = {'QED':1})

GC_79 = Coupling(name = 'GC_79',
                 value = '(cw*complex(0,1)*g1*gw*TH2x1)/2.',
                 order = {'QED':2})

GC_8 = Coupling(name = 'GC_8',
                value = 'gw/2.',
                order = {'QED':1})

GC_80 = Coupling(name = 'GC_80',
                 value = '-(complex(0,1)*g1*gw*sw*TH2x1)/2.',
                 order = {'QED':2})

GC_81 = Coupling(name = 'GC_81',
                 value = '(l4*TH2x1)/2. - l5*TH2x1',
                 order = {'QED':2})

GC_82 = Coupling(name = 'GC_82',
                 value = '-(l4*TH2x1)/2. + l5*TH2x1',
                 order = {'QED':2})

GC_83 = Coupling(name = 'GC_83',
                 value = '(complex(0,1)*gw**2*TH1x1**2)/2. + (complex(0,1)*gw**2*TH2x1**2)/2.',
                 order = {'QED':2})

GC_84 = Coupling(name = 'GC_84',
                 value = '-(complex(0,1)*l3*TH1x1**2) - 2*complex(0,1)*l7*TH1x1*TH2x1 - 2*complex(0,1)*l2*TH2x1**2',
                 order = {'QED':2})

GC_85 = Coupling(name = 'GC_85',
                 value = '-2*complex(0,1)*l1*TH1x1**2 - 2*complex(0,1)*l6*TH1x1*TH2x1 - complex(0,1)*l3*TH2x1**2',
                 order = {'QED':2})

GC_86 = Coupling(name = 'GC_86',
                 value = '-2*complex(0,1)*l1*TH1x1**2 - 2*complex(0,1)*l6*TH1x1*TH2x1 - complex(0,1)*l3*TH2x1**2 - complex(0,1)*l4*TH2x1**2 + 2*complex(0,1)*l5*TH2x1**2',
                 order = {'QED':2})

GC_87 = Coupling(name = 'GC_87',
                 value = '-(complex(0,1)*l6*TH1x1**2) - complex(0,1)*l4*TH1x1*TH2x1 - 2*complex(0,1)*l5*TH1x1*TH2x1 - complex(0,1)*l7*TH2x1**2',
                 order = {'QED':2})

GC_88 = Coupling(name = 'GC_88',
                 value = '(cw**2*complex(0,1)*gw**2*TH1x1**2)/2. + cw*complex(0,1)*g1*gw*sw*TH1x1**2 + (complex(0,1)*g1**2*sw**2*TH1x1**2)/2. + (cw**2*complex(0,1)*gw**2*TH2x1**2)/2. + cw*complex(0,1)*g1*gw*sw*TH2x1**2 + (complex(0,1)*g1**2*sw**2*TH2x1**2)/2.',
                 order = {'QED':2})

GC_89 = Coupling(name = 'GC_89',
                 value = '-(cw**2*complex(0,1)*g1*gw*TH1x1**2)/2. - (cw*complex(0,1)*g1**2*sw*TH1x1**2)/2. + (cw*complex(0,1)*gw**2*sw*TH1x1**2)/2. + (complex(0,1)*g1*gw*sw**2*TH1x1**2)/2. - (cw**2*complex(0,1)*g1*gw*TH2x1**2)/2. - (cw*complex(0,1)*g1**2*sw*TH2x1**2)/2. + (cw*complex(0,1)*gw**2*sw*TH2x1**2)/2. + (complex(0,1)*g1*gw*sw**2*TH2x1**2)/2.',
                 order = {'QED':2})

GC_9 = Coupling(name = 'GC_9',
                value = '(complex(0,1)*gw)/cmath.sqrt(2)',
                order = {'QED':1})

GC_90 = Coupling(name = 'GC_90',
                 value = '(cw**2*complex(0,1)*g1**2*TH1x1**2)/2. - cw*complex(0,1)*g1*gw*sw*TH1x1**2 + (complex(0,1)*gw**2*sw**2*TH1x1**2)/2. + (cw**2*complex(0,1)*g1**2*TH2x1**2)/2. - cw*complex(0,1)*g1*gw*sw*TH2x1**2 + (complex(0,1)*gw**2*sw**2*TH2x1**2)/2.',
                 order = {'QED':2})

GC_91 = Coupling(name = 'GC_91',
                 value = '-6*complex(0,1)*l1*TH1x1**4 - 12*complex(0,1)*l6*TH1x1**3*TH2x1 - 6*complex(0,1)*l3*TH1x1**2*TH2x1**2 - 6*complex(0,1)*l4*TH1x1**2*TH2x1**2 - 12*complex(0,1)*l5*TH1x1**2*TH2x1**2 - 12*complex(0,1)*l7*TH1x1*TH2x1**3 - 6*complex(0,1)*l2*TH2x1**4',
                 order = {'QED':2})

GC_92 = Coupling(name = 'GC_92',
                 value = '-(complex(0,1)*gw*TH2x2)/2.',
                 order = {'QED':1})

GC_93 = Coupling(name = 'GC_93',
                 value = '(complex(0,1)*gw*TH2x2)/2.',
                 order = {'QED':1})

GC_94 = Coupling(name = 'GC_94',
                 value = '(cw*complex(0,1)*g1*gw*TH2x2)/2.',
                 order = {'QED':2})

GC_95 = Coupling(name = 'GC_95',
                 value = '-(complex(0,1)*g1*gw*sw*TH2x2)/2.',
                 order = {'QED':2})

GC_96 = Coupling(name = 'GC_96',
                 value = '(l4*TH2x2)/2. - l5*TH2x2',
                 order = {'QED':2})

GC_97 = Coupling(name = 'GC_97',
                 value = '-(l4*TH2x2)/2. + l5*TH2x2',
                 order = {'QED':2})

GC_98 = Coupling(name = 'GC_98',
                 value = '(complex(0,1)*gw**2*TH1x1*TH1x2)/2. + (complex(0,1)*gw**2*TH2x1*TH2x2)/2.',
                 order = {'QED':2})

GC_99 = Coupling(name = 'GC_99',
                 value = '-(complex(0,1)*l3*TH1x1*TH1x2) - complex(0,1)*l7*TH1x2*TH2x1 - complex(0,1)*l7*TH1x1*TH2x2 - 2*complex(0,1)*l2*TH2x1*TH2x2',
                 order = {'QED':2})

