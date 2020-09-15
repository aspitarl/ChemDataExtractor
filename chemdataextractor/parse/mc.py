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

from chemdataextractor.parse.cem import cem
from chemdataextractor.parse.common import slash, delim, optdelim
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
floats = R('[><~]?\d*\.?\d+([eE][-+]?\d+)?')('value')#Matches floats + scientific notation
joined_range = R('\d*\.?\d+([eE][-+]?\d+)?[\-–−~∼˜][-+]?\d*\.?\d+([eE][-+]?\d+)?')
spaced_range = (floats + R('[\-–−~∼˜]') + floats).add_action(merge)
to_range = (floats + I('to') + floats).add_action(replace_to).add_action(merge)
stddev_range = (floats + R('±') + floats).add_action(merge)

values = (joined_range | spaced_range | to_range | floats | stddev_range)('value')

opt_creatinine = Optional(W('creatinine') | W('Cr'))

molar_unit =  R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]?M$')
mass_unit =  R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]?(mol|g)$')#I know that moles are a unit of quantity, don't @ me
volume_unit =  R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]?(mol|l|L)$')
units = ((molar_unit | mass_unit + slash + volume_unit).add_action(merge) + opt_creatinine)('units').add_action(join)

left_bracket_parentheses = R('[[(]')
right_bracket_parentheses = R('[)\]]')

mc = (Optional(left_bracket_parentheses).hide() + Optional(cem + optdelim) + values + units + Optional(right_bracket_parentheses).hide())('mc')

mediums = (I('saliva') | I('urine') | I('blood') | I('serum') | I('plasma'))('medium')


infix_words = W('was') | W('were') | W('of') | W('concentration') | W('concentrations') | W('normal') | W('in') | W('level') | W('levels') | W('increased') | W('to') | W('are')
infix_words_opt = Optional(OneOrMore(infix_words)).hide()
# Maybe have one set of infix words that doesn't depend on context and another that does? In order to not catch treatment chemicals?
# Or maybe that's a job for a separate parser for treatment chems
# Either way, if the list gets too big then replace it with multiple more fleshed out infix phrases

#def extract_cem_from_emia(tokens, start, result):
    
normal_range = (left_bracket_parentheses.hide() + ((W('normal') + W('range')) | W('NR')) + optdelim + values)('normal_range')
mc_phrase_cem_first = (Optional(mediums) + cem + infix_words_opt + mc + Optional(normal_range))('mc_phrase')
mc_phrase_cem_last = (mc + infix_words_opt +  Optional(mediums) + cem + Optional(SkipTo(W('in') + mediums)) + Optional(normal_range))('mc_phrase')
mc_phrase_list = (mediums + infix_words_opt + OneOrMore(cem + mc + delim))('mc_phrase')
mc_phrase_and = (mc + cem + W('and').hide() + mc + cem + W('in').hide() + mediums)('mc_phrase')

mc_phrase =  mc_phrase_and | mc_phrase_cem_first | mc_phrase_cem_last | mc_phrase_list
#TODO: Ideally, instead of an infix word list, I could perhaps set it up to be sensitive to certain kinds of
# tags or parse ahead until finding a match for things like medium, etc. But that's for the iteration phase
class McParser(BaseParser):
    """"""
    root = mc_phrase

    def interpret(self, result, start, end):
        for cem_el, mc in zip(result.xpath('/cem'), result.xpath('./mc')):
            compound = Compound(
                measured_concentrations=[
                    MeasuredConcentration(
                        value=first(mc.xpath('./value/text()')),
                        units=first(mc.xpath('./units/text()')),
                        location=first(result.xpath('./medium/text()')),
                        normal_concentration = first(result.xpath('./normal_range/value/text()'))
                        )
                    ]
                )
            if cem_el is not None:
                compound.names = cem_el.xpath('./name/text()')
                compound.labels = cem_el.xpath('./label/text()')
            yield compound