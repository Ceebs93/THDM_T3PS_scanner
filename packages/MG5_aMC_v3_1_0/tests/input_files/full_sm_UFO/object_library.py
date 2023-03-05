##
##
## Feynrules Header
##
##
##
##
##

from __future__ import absolute_import
import cmath
import six


class FRBaseClass(object):
    """The class from which all FeynRules classes are derived."""

    require_args = []

    def __init__(self, *args, **options):
        assert(len(self.require_args) == len (args))
    
        for i, name in enumerate(self.require_args):
            setattr(self, name, args[i])
    
        for (option, value) in options.items():
            setattr(self, option, value)

    def get(self, name):
        return getattr(self, name)
    
    def set(self, name, value):
        setattr(self, name, value)
        
    def get_all(self):
        """Return a dictionary containing all the information of the object"""
        return self.__dict__

    def __str__(self):
        return self.name

    def nice_string(self):
        """ return string with the full information """
        return '\n'.join(['%s \t: %s' %(name, value) for name, value in self.__dict__.items()])

    def __repr__(self):
        replacements = [
            ('+','__plus__'),
            ('-','__minus__'),
            ('@','__at__'),
            ('!','__exclam__'),
            ('?','__quest__'),
            ('*','__star__'),
            ('~','__tilde__')
            ]
        text = self.name
        for orig,sub in replacements:
            text = text.replace(orig,sub)
        return text



all_particles = []

    

class Particle(FRBaseClass):
    """A standard Particle"""

    require_args=['pdg_code', 'name', 'antiname', 'spin', 'color', 'mass', 'width', 'texname', 'antitexname', 'line', 'charge']

    require_args_all = ['pdg_code', 'name', 'antiname', 'spin', 'color', 'mass', 'width', 'texname', 'antitexname', 'line', 'charge', 'propagating', 'goldstoneboson']

    def __init__(self, pdg_code, name, antiname, spin, color, mass, width, texname,
                 antitexname, line, charge , propagating=True, goldstoneboson=False, **options):

        args= (pdg_code, name, antiname, spin, color, mass, width, texname,
                 antitexname, line, float(charge))

        FRBaseClass.__init__(self, *args,  **options)

        global all_particles
        all_particles.append(self)

        self.propagating = propagating
        self.goldstoneboson= goldstoneboson

        self.selfconjugate = (name == antiname)

    def anti(self):
        if self.selfconjugate:
            raise Exception('%s has no anti particle.' % self.name) 
        outdic = {}
        for k,v in six.iteritems(self.__dict__):
            if k not in self.require_args_all:                
                outdic[k] = -v
        if self.color in [1,8]:
            newcolor = self.color
        else:
            newcolor = -self.color
                
        return Particle(-self.pdg_code, self.antiname, self.name, self.spin, newcolor, self.mass, self.width,
                        self.antitexname, self.texname, self.line, -self.charge, self.propagating, self.goldstoneboson, **outdic)



all_parameters = []

class Parameter(FRBaseClass):

    require_args=['name', 'nature', 'type', 'value', 'texname']

    def __init__(self, name, nature, type, value, texname, lhablock=None, lhacode=None):

        args = (name,nature,type,value,texname)

        FRBaseClass.__init__(self, *args)

        args=(name,nature,type,value,texname)

        global all_parameters
        all_parameters.append(self)

        if (lhablock is None or lhacode is None)  and nature == 'external':
            raise Exception('Need LHA information for external parameter "%s".' % name)
        self.lhablock = lhablock
        self.lhacode = lhacode

all_vertices = []

class Vertex(FRBaseClass):

    require_args=['name', 'particles', 'color', 'lorentz', 'couplings']

    def __init__(self, name, particles, color, lorentz, couplings, **opt):
 
        args = (name, particles, color, lorentz, couplings)

        FRBaseClass.__init__(self, *args, **opt)

        args=(particles,color,lorentz,couplings)

        global all_vertices
        all_vertices.append(self)

all_couplings = []

class Coupling(FRBaseClass):

    require_args=['name', 'value', 'order']

    def __init__(self, name, value, order, **opt):

        args =(name, value, order)	
        FRBaseClass.__init__(self, *args, **opt)
        global all_couplings
        all_couplings.append(self)
  


all_lorentz = []

class Lorentz(FRBaseClass):

    require_args=['name','spins','structure']
    
    def __init__(self, name, spins, structure='external', **opt):
        args = (name, spins, structure)
        FRBaseClass.__init__(self, *args, **opt)

        global all_lorentz
        all_lorentz.append(self)


all_functions = []

class Function(object):

    def __init__(self, name, arguments, expression):

        global all_functions
        all_functions.append(self)

        self.name = name
        self.arguments = arguments
        self.expr = expression
    
    def __call__(self, *opt):

        for i, arg in enumerate(self.arguments):
            exec('%s = %s' % (arg, opt[i] ))

        return eval(self.expr)


        
