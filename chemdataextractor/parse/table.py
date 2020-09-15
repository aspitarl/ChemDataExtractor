# -*- coding: utf-8 -*-
"""
chemdataextractor.parse.table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
import logging
import re
from lxml.builder import E

from .common import delim, optdelim, slash, hyph
from ..utils import first
from ..model import Compound, UvvisSpectrum, UvvisPeak, QuantumYield, FluorescenceLifetime, MeltingPoint, GlassTransition
from ..model import ElectrochemicalPotential, IrSpectrum, IrPeak, MeasuredConcentration
from .actions import join, merge, fix_whitespace
from .base import BaseParser
from .cem import chemical_label, label_before_name, chemical_name, chemical_label_phrase, solvent_name, lenient_chemical_label, lenient_name, prefixed_label
from .elements import R, I, W, Optional, ZeroOrMore, Any, OneOrMore, Start, End, Group, Not, First, SkipTo

log = logging.getLogger(__name__)


delims = ZeroOrMore(delim)
minus = R('^[\-–−‒]$')


name_blacklist = R('^([\d\.]+)$')   
name_token_blacklist = R('^(al|I|II|control|controls|range)$', re.I)#Currently only works if token at end of cell
#TODO
'''I spent a lot of time trying to find a solution that would let me disqualify an entire cell if I found a blacklisted
substring anywhere in the cell, but couldn't. Either I'm missing something or it isn't possible. Must investigate further'''
#: Compound identifier column heading
compound_heading = R('(^|\b)(comp((oun)?d)?|molecule|ligand|oligomer|complex|dye|porphyrin|substance|sample|material|catalyst|acronym|isomer|(co)?polymer|chromophore|species|quinone|ether|diene|adduct|acid|radical|monomer|amine|analyte|product|system|(photo)?sensitiser|phthalocyanine|MPc)(e?s)?($|\b)', re.I)
solvent_heading = R('(^|\b)(solvent)s?($|\b)', re.I)
solvent_in_heading = Group(solvent_name)('cem')
solvent_cell = Group(solvent_name | chemical_name)('cem')
compound_cell = Group(
    (Start() + chemical_label + End())('cem') |
    (Start() + lenient_chemical_label + End())('cem') |
    chemical_label_phrase('cem') |
    (Not((Start() + OneOrMore(name_blacklist) + End()) | SkipTo(name_token_blacklist)) + OneOrMore(Any())('name').add_action(join).add_action(fix_whitespace) + Optional(W('(').hide() + chemical_label + W(')').hide()))('cem') |
    label_before_name
)('cem_phrase')

uvvis_emi_title = (
    I('emission') + R('max(ima)?') |
    W('λ') + Optional(I('max')) + Optional(W(',')) + R('em(i(ssion)?)?', re.I) |
    R('em(i(ssion)?)?', re.I) + W('λ') + Optional(I('max')) + Optional(W(','))
)
uvvis_abs_title = (
    I('absorption') + R('max(ima)?') |
    W('λ') + OneOrMore(R('^(a|sol)?max$', re.I) | R('abs(or[bp]tion)?', re.I) | I('a') | W(',')) |
    R('uv([-/]?vis)?', re.I)
)
extinction_title = Optional(R('^10\d$') | W('10') + minus + R('^\d$')).hide() + W('ε') + Optional(I('max'))
uvvis_units = (W('nm') | R('^eV[\-–−‒]1$') | W('eV') + minus + W('1'))('uvvis_units').add_action(merge)
multiplier = Optional(I('×')) + (R('^10–?[34]$') | (W('10') + minus + R('^[345]$')))

extinction_units = (
    (Optional(multiplier + delims) + (
        I('M') + minus + I('1') + I('cm') + minus + I('1') |
        I('M') + minus + I('1') + I('cm') + minus + I('1') |
        I('dm3') + I('mol') + minus + I('1') + I('cm') + minus + I('1') |
        I('l') + I('mol') + minus + I('1') + I('cm') + minus + I('1') |
        I('l') + I('cm') + minus + I('1') + I('mol') + minus + I('1')
    )) | multiplier
)('extinction_units').add_action(join)


ir_title = (
    R('^(FT-?)?IR$') + Optional(I('absorption'))
)
ir_units = Optional(W('/')).hide() + (
    R('^\[?cm[-–−]1\]?$') |
    W('cm') + R('^[-–−]$') + W('1')
)('ir_units').add_action(merge)
ir_heading = (OneOrMore(ir_title.hide()) + ZeroOrMore(delims.hide() + ir_units))('ir_heading')
ir_value = (R('^\d{3,5}(\.\d{1,2})?$'))('value')
peak_strength = R('^(sh(oulder)?|br(oad)?)$')('strength')
ir_peak = (
    ir_value + Optional(W('(').hide()) + Optional(peak_strength) + Optional(W(')').hide())
)('ir_peak')
ir_cell = (
    ir_peak + ZeroOrMore(W(',').hide() + ir_peak)
)('ir_cell')

# TODO: (photoluminescence|fluorescence) quantum yield
quantum_yield_title = (R('^(Φ|ϕ)(fl?|pl|ze|t|l|lum)?$', re.I) + Optional(R('^(fl?|pl|ze|t|l|lum)$', re.I)))('quantum_yield_type').add_action(merge)  #  + ZeroOrMore(Any())
quantum_yield_units = W('%')('quantum_yield_units')
quantum_yield_heading = Group(Start() + quantum_yield_title + delims.hide() + Optional(quantum_yield_units) + delims.hide() + End())('quantum_yield_heading')
quantum_yield_value = (Optional(R('^[~∼\<\>]$')) + ((W('10') + minus + R('^\d$')) | R('^(100(\.0+)?|\d\d?(\.\d+)?)$')) + Optional(W('±') + R('^\d+(\.\d+)?$')))('quantum_yield_value').add_action(merge)
quantum_yield_cell = (quantum_yield_value + Optional(quantum_yield_units))('quantum_yield_cell')


def split_uvvis_shape(tokens, start, result):
    """"""
    if result[0].text.endswith('sh') or result[0].text.endswith('br'):
        result.append(E('shape', result[0].text[-2:]))
        result[0].text = result[0].text[:-2]



uvvis_emi_heading = (OneOrMore(uvvis_emi_title.hide()))('uvvis_emi_heading')
uvvis_abs_heading = (OneOrMore(uvvis_abs_title.hide()) + ZeroOrMore(delims.hide() + (uvvis_units | extinction_title.hide() | extinction_units)))('uvvis_abs_heading')
uvvis_abs_disallowed = I('emission')
extinction_heading = (extinction_title.hide() + delims.hide() + Optional(extinction_units))('extinction_heading')
uvvis_value = (R('^\d{3,4}(\.\d{1,2})?(sh|br)?$'))('value').add_action(split_uvvis_shape)
peak_shape = R('^(sh(oulder)?|br(oad)?)$')('shape')
extinction_value = (
    R('^\d+\.\d+$') + Optional(W('±') + R('^\d+\.\d+$')) + Optional(W('×') + R('10\d+')) |  # Scientific notation
    R('^\d{1,3}$') + R('^\d\d\d$') |  # RSC often inserts spaces within values instead of commas
    R('^\d{1,2},?\d{3,3}$')

)('extinction').add_action(merge)


uvvis_abs_emi_quantum_yield_heading = (
    OneOrMore(uvvis_abs_title.hide()) +
    Optional(Optional(delims.hide()) + uvvis_units('uvvis_abs_units') + Optional(delims.hide())) +
    OneOrMore(uvvis_emi_title.hide()) +
    Optional(Optional(delims.hide()) + uvvis_units + Optional(delims.hide())) +
    Optional(delims.hide()) + quantum_yield_title.hide() + Optional(delims.hide()) +
    Optional(Optional(delims.hide()) + quantum_yield_units + Optional(delims.hide()))
)('uvvis_emi_quantum_yield_heading')

uvvis_abs_emi_quantum_yield_cell = (
    uvvis_value('uvvis_abs_value') + delims.hide() + uvvis_value + delims.hide() + quantum_yield_value + Optional(quantum_yield_units)
)('uvvis_emi_quantum_yield_cell')


uvvis_emi_quantum_yield_heading = (
    OneOrMore(uvvis_emi_title.hide()) +
    Optional(Optional(delims.hide()) + uvvis_units + Optional(delims.hide())) +
    Optional(delims.hide()) + quantum_yield_title.hide() + Optional(delims.hide()) +
    Optional(Optional(delims.hide()) + quantum_yield_units + Optional(delims.hide()))
)('uvvis_emi_quantum_yield_heading')

uvvis_emi_quantum_yield_cell = (
    uvvis_value + delims.hide() + quantum_yield_value + Optional(quantum_yield_units)
)('uvvis_emi_quantum_yield_cell')

uvvis_abs_peak = (
    uvvis_value + Optional(peak_shape) + Optional(W('(').hide() + extinction_value + W(')').hide())
)('uvvis_abs_peak')

uvvis_abs_cell = (
    uvvis_abs_peak + ZeroOrMore(W(',').hide() + uvvis_abs_peak)
)('uvvis_abs_cell')

extinction_cell = (
    extinction_value + ZeroOrMore(W(',').hide() + extinction_value)
)('uvvis_abs_cell')

uvvis_emi_peak = (
    uvvis_value + Optional(peak_shape)
)('uvvis_emi_peak')

uvvis_emi_cell = (
    uvvis_emi_peak + ZeroOrMore(W(',').hide() + uvvis_emi_peak)
)('uvvis_emi_cell')


fluorescence_lifetime_title = W('τ') + R('^(e|f|ave|avg|0)$', re.I)
fluorescence_lifetime_units = (W('ns') | W('μ') + W('s'))('fluorescence_lifetime_units').add_action(merge)
fluorescence_lifetime_heading = (fluorescence_lifetime_title.hide() + delims.hide() + Optional(fluorescence_lifetime_units))('fluorescence_lifetime_heading')
fluorescence_lifetime_value = (Optional(R('^[~∼\<\>]$')) + R('^\d+(\.\d+)?$'))('fluorescence_lifetime_value').add_action(merge)
fluorescence_lifetime_cell = (
    fluorescence_lifetime_value + ZeroOrMore(W(',').hide() + fluorescence_lifetime_value)
)('fluorescence_lifetime_cell')

electrochemical_potential_title = ((R('^E(ox|red)1?$', re.I) | W('E') + R('^(ox|red)1?$')) + Optional(W('/') + W('2')))('electrochemical_potential_type').add_action(merge)
electrochemical_potential_units = (W('V'))('electrochemical_potential_units').add_action(merge)
electrochemical_potential_heading = (electrochemical_potential_title + delims.hide() + Optional(electrochemical_potential_units))('electrochemical_potential_heading')
electrochemical_potential_value = (Optional(R('^[~∼\<\>]$')) + Optional(minus) + R('^\d+(\.\d+)?$'))('electrochemical_potential_value').add_action(merge)
electrochemical_potential_cell = (
    electrochemical_potential_value + ZeroOrMore(delims.hide() + electrochemical_potential_value)
)('electrochemical_potential_cell')

subject_phrase_preposition = (I('of') | I('for')) + chemical_name
# TODO: Come back to this. Hopefully what's here works for most cases, but it fails for L-2-HG. Maybe see if it wouldn't be irresponsible to match \w-\d-\w+
subject_phrase_concentration = ((chemical_label | chemical_name | prefixed_label) + R('(concentration)s?', re.I).hide())
subject_phrase = (subject_phrase_preposition | subject_phrase_concentration)('subject_phrase')
solvent_phrase = (I('in') + (solvent_name | chemical_name))('solvent_phrase')

temp_range = (Optional(R('^[\-–−]$')) + (R('^[\+\-–−]?\d+(\.\d+)?[\-–−]\d+(\.\d+)?$') | (R('^[\+\-–−]?\d+(\.\d+)?$') + R('^[\-–−]$') + R('^[\+\-–−]?\d+(\.\d+)?$'))))('temperature').add_action(merge)
temp_value = (Optional(R('^[\-–−]$')) + R('^[\+\-–−]?\d+(\.\d+)?$') + Optional(W('±') + R('^\d+(\.\d+)?$')))('temperature').add_action(merge)
temp_word = (I('room') + R('^temp(erature)?$') | R('^r\.?t\.?$', re.I))('temperature').add_action(merge)
temp = (temp_range | temp_value | temp_word)('value')
temp_units = (W('°') + R('[CFK]') | W('K'))('units').add_action(merge)
temp_with_units = (temp + temp_units)('temp')
temp_with_optional_units = (temp + Optional(temp_units))('temp')

temp_phrase = (I('at') + temp_with_units)('temp_phrase')

melting_point_title = R('^T(melt|m\.p|m)$', re.I) | W('T') + R('^(melt|m\.p|m)?$')
melting_point_heading = (melting_point_title.hide() + delims.hide() + Optional(temp_units))('melting_point_heading')
melting_point_cell = (
    temp_with_optional_units + ZeroOrMore(delims.hide() + temp_with_optional_units)
)('melting_point_cell')

glass_transition_title = R('^T(g\.)$', re.I) | W('T') + R('^(g\.)?$')
glass_transition_heading = (glass_transition_title.hide() + delims.hide() + Optional(temp_units))('glass_transition_heading')
glass_transition_cell = (
    temp_with_optional_units + ZeroOrMore(delims.hide() + temp_with_optional_units)
)('glass_transition_cell')

def strip_superscript(tokens, start, result):
    result[0].text = re.sub('[a-z]$', '', result[0].text)
    
    
mc_opt_creatinine = Optional(R('(creatinine|cr|creat)', re.I))
mc_float_non_scientific = R('^[><~]?\d+\.?\d*$')('value')#Removed ability to match scientific notation because it was causing issues, may reimplement later
mc_float_scientific = (R('^[><~]?\d+\.?\d*E$', re.I) + R('^[-–−+]$') + R('^\d+[a-z]?$').add_action(strip_superscript))('value').add_action(merge)
mc_float = mc_float_scientific | mc_float_non_scientific
mc_range = (mc_float + R('[\-–−~∼˜]') + mc_float)('value').add_action(merge)
mc_value_stddev = (mc_float + W('±') + mc_float)('value').add_action(join)
mc_values = (mc_range | mc_value_stddev | mc_float)

mc_molar_unit =  Optional(R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]$')) + R('^[YZEPTGkhcm\u03BC\u00B5unpfzyad]?M$')
mc_mass_unit =  Optional(R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]$')) + R('^[YZEPTGkhcm\u03BC\u00B5unpfzyad]?(mol|g|M)$', re.I)
mc_volume_unit =  Optional(R('^[YZEPTGMkhcm\u03BC\u00B5unpfzyad]$')) + R('^[YZEPTGkhcm\u03BC\u00B5unpfzyad]?(mol|l|L|M)$', re.I)
mc_units = ((mc_molar_unit | (mc_mass_unit + slash + mc_volume_unit)).add_action(merge) + mc_opt_creatinine)('units').add_action(join)

mc_compound_heading = R('(^|\b)(comp((oun)?d)?|substance|analyte|metabolite|metabolic)(e?s)?($|\b)', re.I)
mc_compound = (lenient_chemical_label | lenient_name)('cem')
mc_compound_cell = (mc_compound + OneOrMore(optdelim) + mc_units)('mc_cem_phrase') 

mc_mediums = (R('^(saliva|urine|blood|serum|csf|plasma)$', re.I))('medium')
mc_value_words = R('(admission|concentration)', re.I)
mc_value_heading = Group(mc_units | mc_mediums | mc_value_words)('mc_value_heading')

mc_value = (mc_values + optdelim + Optional(mc_units))('mc_value_cell')
mc_value_n_range = (mc_float + delim + R('\d+').hide() + OneOrMore(delim) + mc_range)('mc_value_cell')
mc_value_range = (mc_float + delim + mc_range)('mc_value_cell')
mc_value_cell = mc_value_n_range | mc_value_range | mc_value
#TODO figure out how to handle units contained in chemical name cells


caption_context = Group(subject_phrase | solvent_phrase | temp_phrase | mc_units)('caption_context')

class CompoundHeadingParser(BaseParser):
    """"""
    root = compound_heading
    
    def interpret(self, result, start, end):
        """"""
        yield Compound()


class SolventHeadingParser(BaseParser):
    """"""
    root = solvent_heading

    def interpret(self, result, start, end):
        """"""
        yield Compound()


class UvvisAbsDisallowedHeadingParser(BaseParser):
    """"""
    root = uvvis_abs_disallowed

    def interpret(self, result, start, end):
        """"""
        yield Compound()


class SolventInHeadingParser(BaseParser):
    """"""
    root = solvent_in_heading

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        solvent = first(result.xpath('./name/text()'))
        if solvent is not None:
            context = {'solvent': solvent}
            c.melting_points = [MeltingPoint(**context)]
            c.glass_transitions = [GlassTransition(**context)]
            c.quantum_yields = [QuantumYield(**context)]
            c.fluorescence_lifetimes = [FluorescenceLifetime(**context)]
            c.electrochemical_potentials = [ElectrochemicalPotential(**context)]
            c.uvvis_spectra = [UvvisSpectrum(**context)]
        if c.serialize():
            yield c


class TempInHeadingParser(BaseParser):
    """"""
    root = temp_with_units

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        context = {
            'temperature': first(result.xpath('./value/text()')),
            'temperature_units': first(result.xpath('./units/text()'))
        }
        c.quantum_yields = [QuantumYield(**context)]
        c.fluorescence_lifetimes = [FluorescenceLifetime(**context)]
        c.electrochemical_potentials = [ElectrochemicalPotential(**context)]
        c.uvvis_spectra = [UvvisSpectrum(**context)]
        yield c


class SolventCellParser(BaseParser):
    """"""
    root = solvent_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        solvent = first(result.xpath('./name/text()'))
        if solvent is not None:
            context = {'solvent': solvent}
            c.melting_points = [MeltingPoint(**context)]
            c.glass_transitions = [GlassTransition(**context)]
            c.quantum_yields = [QuantumYield(**context)]
            c.fluorescence_lifetimes = [FluorescenceLifetime(**context)]
            c.electrochemical_potentials = [ElectrochemicalPotential(**context)]
            c.uvvis_spectra = [UvvisSpectrum(**context)]
        if c.serialize():
            yield c


class CompoundCellParser(BaseParser):
    """"""
    root = compound_cell

    def interpret(self, result, start, end):
        for cem_el in result.xpath('./cem'):
            c = Compound(
                names=cem_el.xpath('./name/text()'),
                labels=cem_el.xpath('./label/text()')
            )
            yield c


class UvvisEmiHeadingParser(BaseParser):
    """"""
    root = uvvis_emi_heading

    def interpret(self, result, start, end):
        """"""
        uvvis_units = first(result.xpath('./uvvis_units/text()'))
        c = Compound()
        # TODO: Emission peaks
        yield c


class UvvisAbsHeadingParser(BaseParser):
    """"""
    root = uvvis_abs_heading

    def interpret(self, result, start, end):
        """"""
        uvvis_units = first(result.xpath('./uvvis_units/text()'))
        extinction_units = first(result.xpath('./extinction_units/text()'))
        c = Compound()
        if uvvis_units or extinction_units:
            c.uvvis_spectra.append(
                UvvisSpectrum(peaks=[UvvisPeak(units=uvvis_units, extinction_units=extinction_units)])
            )
        yield c


class ExtinctionHeadingParser(BaseParser):
    """"""
    root = extinction_heading

    def interpret(self, result, start, end):
        """"""
        extinction_units = first(result.xpath('./extinction_units/text()'))
        c = Compound()
        if extinction_units:
            c.uvvis_spectra.append(
                UvvisSpectrum(peaks=[UvvisPeak(extinction_units=extinction_units)])
            )
        yield c


class IrHeadingParser(BaseParser):
    """"""
    root = ir_heading

    def interpret(self, result, start, end):
        """"""
        ir_units = first(result.xpath('./ir_units/text()'))
        c = Compound()
        if ir_units:
            c.ir_spectra.append(
                IrSpectrum(peaks=[IrPeak(units=ir_units)])
            )
        yield c


class IrCellParser(BaseParser):
    """"""
    root = ir_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        ir = IrSpectrum()
        for peak in result.xpath('./ir_peak'):
            ir.peaks.append(
                IrPeak(
                    value=first(peak.xpath('./value/text()')),
                    strength=first(peak.xpath('./strength/text()'))
                )
            )
        if ir.peaks:
            c.ir_spectra.append(ir)
            yield c


class QuantumYieldHeadingParser(BaseParser):
    """"""
    root = quantum_yield_heading

    def interpret(self, result, start, end):
        """"""
        c = Compound(
            quantum_yields=[
                QuantumYield(
                    type=first(result.xpath('./quantum_yield_type/text()')),
                    units=first(result.xpath('./quantum_yield_units/text()'))
                )
            ]
        )
        yield c


class QuantumYieldCellParser(BaseParser):
    """"""
    root = quantum_yield_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        qy = QuantumYield(
            value=first(result.xpath('./quantum_yield_value/text()')),
            units=first(result.xpath('./quantum_yield_units/text()'))
        )
        if qy.value:
            c.quantum_yields.append(qy)
            yield c


class UvvisEmiCellParser(BaseParser):
    """"""
    root = uvvis_emi_cell

    def interpret(self, result, start, end):
        """"""
        # TODO: Emission peaks
        return
        yield


class UvvisAbsCellParser(BaseParser):
    """"""
    root = uvvis_abs_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        uvvis = UvvisSpectrum()
        for peak in result.xpath('./uvvis_abs_peak'):
            uvvis.peaks.append(
                UvvisPeak(
                    value=first(peak.xpath('./value/text()')),
                    extinction=first(peak.xpath('./extinction/text()')),
                    shape=first(peak.xpath('./shape/text()'))
                )
            )
        if uvvis.peaks:
            c.uvvis_spectra.append(uvvis)
            yield c


class ExtinctionCellParser(BaseParser):
    """"""
    root = extinction_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        uvvis = UvvisSpectrum()
        for value in result.xpath('./extinction/text()'):
            uvvis.peaks.append(
                UvvisPeak(
                    extinction=value,
                )
            )
        if uvvis.peaks:
            c.uvvis_spectra.append(uvvis)
            yield c


class UvvisAbsEmiQuantumYieldHeadingParser(BaseParser):
    """"""
    root = uvvis_abs_emi_quantum_yield_heading

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        abs_units = first(result.xpath('./uvvis_abs_units/text()'))
        if abs_units:
            c.uvvis_spectra.append(
                UvvisSpectrum(peaks=[UvvisPeak(units=abs_units)])
            )
        qy_units = first(result.xpath('./quantum_yield_units/text()'))
        if qy_units:
            c.quantum_yields.append(
                QuantumYield(units=qy_units)
            )

        yield c


class UvvisAbsEmiQuantumYieldCellParser(BaseParser):
    """"""
    root = uvvis_abs_emi_quantum_yield_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        uvvis = UvvisSpectrum()
        for value in result.xpath('./uvvis_abs_value/text()'):
            uvvis.peaks.append(
                UvvisPeak(
                    value=value,
                )
            )
        if uvvis.peaks:
            c.uvvis_spectra.append(uvvis)
        qy = QuantumYield(
            value=first(result.xpath('./quantum_yield_value/text()'))
        )
        if qy.value:
            c.quantum_yields.append(qy)

        if c.quantum_yields or c.uvvis_spectra:
            yield c


class UvvisEmiQuantumYieldHeadingParser(BaseParser):
    """"""
    root = uvvis_emi_quantum_yield_heading

    def interpret(self, result, start, end):
        """"""
        # Yield an empty compound to signal that the Parser matched
        yield Compound()


class UvvisEmiQuantumYieldCellParser(BaseParser):
    """"""
    root = uvvis_emi_quantum_yield_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        qy = QuantumYield(
            value=first(result.xpath('./quantum_yield_value/text()'))
        )
        if qy.value:
            c.quantum_yields.append(qy)
            yield c


class FluorescenceLifetimeHeadingParser(BaseParser):
    """"""
    root = fluorescence_lifetime_heading

    def interpret(self, result, start, end):
        """"""
        fluorescence_lifetime_units = first(result.xpath('./fluorescence_lifetime_units/text()'))
        c = Compound()
        if fluorescence_lifetime_units:
            c.fluorescence_lifetimes.append(
                FluorescenceLifetime(units=fluorescence_lifetime_units)
            )
        yield c


class FluorescenceLifetimeCellParser(BaseParser):
    """"""
    root = fluorescence_lifetime_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        fl = FluorescenceLifetime(
            value=first(result.xpath('./fluorescence_lifetime_value/text()'))
        )
        if fl.value:
            c.fluorescence_lifetimes.append(fl)
            yield c


class MeltingPointHeadingParser(BaseParser):
    """"""
    root = melting_point_heading

    def interpret(self, result, start, end):
        """"""
        melting_point_units = first(result.xpath('./units/text()'))
        c = Compound()
        if melting_point_units:
            c.melting_points.append(
                MeltingPoint(units=melting_point_units)
            )
        yield c


class MeltingPointCellParser(BaseParser):
    """"""
    root = melting_point_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        for mp in result.xpath('./temp'):
            c.melting_points.append(
                MeltingPoint(
                    value=first(mp.xpath('./value/text()')),
                    units=first(mp.xpath('./units/text()'))
                )
            )
        if c.melting_points:
            yield c


class GlassTransitionHeadingParser(BaseParser):
    """"""
    root = glass_transition_heading

    def interpret(self, result, start, end):
        """"""
        glass_transition_units = first(result.xpath('./units/text()'))
        c = Compound()
        if glass_transition_units:
            c.glass_transitions.append(
                GlassTransition(units=glass_transition_units)
            )
        yield c


class GlassTransitionCellParser(BaseParser):
    """"""
    root = glass_transition_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        for tg in result.xpath('./temp'):
            c.glass_transitions.append(
                GlassTransition(
                    value=first(tg.xpath('./value/text()')),
                    units=first(tg.xpath('./units/text()'))
                )
            )
        if c.glass_transition:
            yield c


class ElectrochemicalPotentialHeadingParser(BaseParser):
    """"""
    root = electrochemical_potential_heading

    def interpret(self, result, start, end):
        """"""
        c = Compound(
            electrochemical_potentials=[
                ElectrochemicalPotential(
                    type=first(result.xpath('./electrochemical_potential_type/text()')),
                    units=first(result.xpath('./electrochemical_potential_units/text()'))
                )
            ]
        )
        yield c


class ElectrochemicalPotentialCellParser(BaseParser):
    """"""
    root = electrochemical_potential_cell

    def interpret(self, result, start, end):
        """"""
        c = Compound()
        for value in result.xpath('./electrochemical_potential_value/text()'):
            c.electrochemical_potentials.append(
                ElectrochemicalPotential(
                    value=value
                )
            )
        if c.electrochemical_potentials:
            yield c

"""These two parsers seem to run indefinitely, I'm not sure if it's just a consequence of their lax nature or if something is infinitely looping. They're currently unused"""
class McCompoundHeadingParser(BaseParser):
    """Checks for headings that indicate measured compounds. Currently unused"""
    root = mc_compound_heading

    def interpret(self, result, start, end):
        """"""
        yield Compound()


class McCompoundCellParser(BaseParser):
    """Checks for measured compound names and contextual units. Currently unused"""
    root = mc_compound_cell

    def interpret(self, result, start, end):
        c = Compound(
            names=result.xpath('./cem/name/text()'),
            labels=result.xpath('./cem/label/text()'),
        )
        mc_units = result.xpath('./units/text()')
        
        if mc_units:
            c.measured_concentrations.append(
                MeasuredConcentration(units=mc_units)
                )
            
        yield c

class McValueHeadingParser(BaseParser):
    """Checks for heading titles that indicate chemical concentration measurement values and contextual units"""
    root = mc_value_heading
    def interpret(self, result, start, end):
        mc_units = first(result.xpath('./units/text()'))
        mc_medium = first(result.xpath('./medium/text()'))
        c = Compound()
        if mc_units or mc_medium:
            c.measured_concentrations.append(
                MeasuredConcentration(units=mc_units, location=mc_medium)
            )
        yield c


class McValueCellParser(BaseParser):
    """Checks for chemical concentration values and units"""
    root = mc_value_cell
    
    def interpret(self, result, start, end):
        c = Compound()
        mc_units = first(result.xpath('./units/text()')) 
        for mc_value in result.xpath('./value/text()'):
            c.measured_concentrations.append(
                MeasuredConcentration(value=mc_value, units=mc_units)
            )
        yield c

class CaptionContextParser(BaseParser):
    """"""
    root = caption_context

    def __init__(self):
        pass

    def interpret(self, result, start, end):
        name = first(result.xpath('./subject_phrase/name/text()'))
        label = first(result.xpath('./subject_phrase/label/text()'))
        c = Compound()
        if name is not None:
            c.names = [name]
        if label is not None:
            c.labels = [label]
        context = {}
        # print(etree.tostring(result[0]))
        solvent = first(result.xpath('./solvent_phrase/name/text()'))
        if solvent is not None:
            context['solvent'] = solvent
        # Melting point shouldn't have contextual temperature
        if context:
            c.melting_points = [MeltingPoint(**context)]
        mc_units = first(result.xpath('./units/text()'))
        if mc_units is not None:
            context['units'] = mc_units
        if context:
            c.measured_concentrations = [MeasuredConcentration(**context)]
        temp = first(result.xpath('./temp_phrase'))
        if temp is not None:
            context['temperature'] = first(temp.xpath('./temp/value/text()'))
            context['temperature_units'] = first(temp.xpath('./temp/units/text()'))
        # Glass transition temperature shouldn't have contextual temperature
        if context:
            c.glass_transitions = [GlassTransition(**context)]
        temp = first(result.xpath('./temp_phrase'))
        if temp is not None:
            context['temperature'] = first(temp.xpath('./temp/value/text()'))
            context['temperature_units'] = first(temp.xpath('./temp/units/text()'))
        if context:
            c.quantum_yields = [QuantumYield(**context)]
            c.fluorescence_lifetimes = [FluorescenceLifetime(**context)]
            c.electrochemical_potentials = [ElectrochemicalPotential(**context)]
            c.uvvis_spectra = [UvvisSpectrum(**context)]
        if c.serialize():
            # print(c.to_primitive())
            yield c
