"""Table 3 of the article: h and s of the same states under different reference states.

CoolProp 6.4.1 gives water with the IAPWS convention (u = 0, s = 0 for saturated
liquid at the triple point). The other conventions are obtained by shifting h
and s so that the reference state of each convention takes its defined values:

  IIR     h = 200 kJ/kg, s = 1 kJ/(kg K) for saturated liquid at 0 degC
  ASHRAE  h = 0, s = 0 for saturated liquid at -40 degC
  NBP     h = 0, s = 0 for saturated liquid at the normal boiling point

Usage: python reference_states.py
"""
from CoolProp.CoolProp import PropsSI


def h_s(fluid, T, Q):
    """Saturated state at temperature T (K) and quality Q, in kJ/kg and kJ/(kg K)."""
    return PropsSI('H', 'T', T, 'Q', Q, fluid) / 1e3, PropsSI('S', 'T', T, 'Q', Q, fluid) / 1e3


def shifted(state, reference, h_ref, s_ref):
    """Values of state when the reference state is given h_ref and s_ref."""
    return state[0] - reference[0] + h_ref, state[1] - reference[1] + s_ref


STATES = [
    ('Water, saturated liquid, 100 degC', 'Water', 373.15, 0),
    ('R-134a, saturated vapour, -10 degC', 'R134a', 263.15, 1),
    ('Ammonia, saturated vapour, -10 degC', 'Ammonia', 263.15, 1),
]

for label, fluid, T, Q in STATES:
    state = h_s(fluid, T, Q)
    T_iir = 273.16 if fluid == 'Water' else 273.15  # water has no liquid below its triple point
    T_nbp = PropsSI('T', 'P', 101325, 'Q', 0, fluid)
    # Same column order as Table 2; water has no ASHRAE column, refrigerants no IAPWS one
    rows = {}
    if fluid == 'Water':
        rows['IAPWS'] = state
    rows['IIR'] = shifted(state, h_s(fluid, T_iir, 0), 200, 1)
    if fluid != 'Water':
        rows['ASHRAE'] = shifted(state, h_s(fluid, 233.15, 0), 0, 0)
    rows['NBP'] = shifted(state, h_s(fluid, T_nbp, 0), 0, 0)
    print(label)
    for convention, (h, s) in rows.items():
        print(f'  {convention:6s} h = {h:8.2f} kJ/kg   s = {s:5.2f} kJ/(kg K)')
