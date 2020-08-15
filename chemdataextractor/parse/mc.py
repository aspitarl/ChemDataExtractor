# -*- coding: utf-8 -*-
"""
chemdataextractor.parse.mc
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Measured Concentration text parser.

"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import re

from chemdataextractor.parse.cem import cem
from chemdataextractor.parse.common import slash, optdelim
from chemdataextractor.utils import first
from ..model import Compound, MeasuredConcentration
from .actions import merge, join
from .base import BaseParser
from .elements import W, R, I, Optional, OneOrMore, SkipTo

log = logging.getLogger(__name__)

# TODO: Get a better understanding of the pipeline here and replace this with a more efficient implementation
# TODO: rewrite this to replace all hyphen-like characters with same hyphen to standardize output. 
def replace_to(tokens, start, result):
    """Replace the word 'to' with a hyphen in to_range"""
    for e in result:
        for child in e.iter():
            if child.text.startswith('to'):
                child.text = '-'
    return result

si_prefixes = '^[YZEPTGMkhcm\u03BCunpfzyad]?'
floats = R('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')('value')#Matches floats + scientific notation
joined_range = R('[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?[\-–−~∼˜][-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?')
# TODO: ranges with scientific notation, or wait I think they might already be covered
spaced_range = (floats + R('[\-–−~∼˜]') + floats).add_action(merge)
to_range = (floats + I('to') + floats).add_action(replace_to).add_action(merge)

values = (joined_range | spaced_range | to_range | floats)('value')

opt_creatinine = Optional(W('creatinine') | W('Cr'))

molar_unit =  R('^[YZEPTGMkhcm\u03BCunpfzyad]?M$')
mass_unit =  R('^[YZEPTGMkhcm\u03BCunpfzyad]?(mol|g)$')#I know that moles are a unit of quantity, don't @ me
volume_unit =  R('^[YZEPTGMkhcm\u03BCunpfzyad]?(mol|l|L)$')
units = ((molar_unit | mass_unit + slash + volume_unit).add_action(merge) + opt_creatinine)('units').add_action(join)

mc = (values + units)('mc')

mediums = (I('saliva') | I('urine') | I('blood') | I('serum') | I('plasma'))('medium')


infix_words = W('was') | W('were') | W('of') | W('concentration') | W('concentrations') | W('normal') | W('in') | W('level') | W('levels') | W('increased') | W('to')
# Maybe have one set of infix words that doesn't depend on context and another that does? In order to not catch treatment chemicals?
# Or maybe that's a job for a separate parser for treatment chems
# Either way, if the list gets too big then replace it with multiple more fleshed out infix phrases
left_bracket_parentheses = R('[[(]')
right_bracket_parentheses = R('[)\]]')
normal_range = (left_bracket_parentheses.hide() + ((W('normal') + W('range')) | W('NR')) + optdelim + values)('normal_range')
mc_phrase_cem_first = (Optional(mediums) + cem + Optional(OneOrMore(infix_words)).hide() + Optional(left_bracket_parentheses).hide() + mc + Optional(normal_range))('mc_phrase')
mc_phrase_cem_last = (mc + Optional(OneOrMore(infix_words)).hide() +  Optional(mediums) + cem + Optional(SkipTo(W('in') + mediums)) + Optional(normal_range))('mc_phrase')

mc_phrase = mc_phrase_cem_first | mc_phrase_cem_last
#TODO: Ideally, instead of an infix word list, I could perhaps set it up to be sensitive to certain kinds of
# tags or parse ahead until finding a match for things like medium, etc. But that's for the iteration phase
#TODO : make records retrieved by McParser merge with other records
class McParser(BaseParser):
    """"""
    root = mc_phrase

    def interpret(self, result, start, end):
        compound = Compound(
            measured_concentrations=[
                MeasuredConcentration(
                    value=first(result.xpath('./mc/value/text()')),
                    units=first(result.xpath('./mc/units/text()')),
                    location=first(result.xpath('./medium/text()')),
                    normal_concentration = first(result.xpath('./normal_range/value/text()'))
                )
            ]
        )
        cem_el = first(result.xpath('./cem'))
        if cem_el is not None:
            compound.names = cem_el.xpath('./name/text()')
            compound.labels = cem_el.xpath('./label/text()')
        yield compound