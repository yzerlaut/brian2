#!/usr/bin/env python
"""
Model of basal dendrites of L23 pyramidal cells
----------------------------------------

Rothman JS, Manis PB (2003) The roles potassium currents play in
regulating the electrical activity of ventral cochlear nucleus neurons.
J Neurophysiol 89:3097-113.

All model types differ only by the maximal conductances.

Adapted from their Neuron implementation by Romain Brette
"""
from brian2 import *

#defaultclock.dt=0.025*ms # for better precision


"""
------------------------------------
 Low-threshold Calcium current (iT)
------------------------------------

Genealogy (NEURON comment):

T-type Ca channel 
ca.mod to lead to thalamic ca current inspired by destexhe and huguenrd
Uses fixed eca instead of GHK eqn
changed from (AS Oct0899)
changed for use with Ri18  (B.Kampa 2005)
"""

gating_var = EquationTemplate('''d{name}/dt = q10*({name}__inf - {name})/tau_{name} : 1
                                 {name}_inf = {inf_rate_expression}                 : 1
                                 tau_{name} = {tau_expression}                      : second
                                 tau_{name} = {tau_scale}/({forward_rate} + {reverse_rate}) + {tau_base} : second''')


pos_sigmoid = ExpressionTemplate('1./(1+exp(-(vu - {midpoint}) / {scale}))')
sqrt_sigmoid = ExpressionTemplate('1./(1+exp(-(vu - {midpoint}) / {scale}))**0.5')
neg_sigmoid = ExpressionTemplate('1./(1+exp((vu - {midpoint}) / {scale}))')
exp_voltage_dep = ExpressionTemplate('{magnitude}*exp((vu-{midpoint})/{scale})')
neg_exp_voltage_dep = ExpressionTemplate('{magnitude}*exp(-(vu-{midpoint})/{scale})')

# Classical Na channel
m = gating_var(name='m',
               rate_expression=pos_sigmoid(midpoint=-38., scale=7.),
               tau_expression=pos_sigmoid(midpoint=-38., scale=7.),
               forward_rate=exp_voltage_dep(magnitude=5., midpoint=-60, scale=18.),
               reverse_rate=neg_exp_voltage_dep(magnitude=36., midpoint=-60, scale=25.),
               tau_base=0.04*ms, tau_scale=10*ms)


        I{name} = g{name} * (v - {E_Ca}*mV): amp/meter**2
        g{name} = gbar_{name} * m{name}**2 * h{name} : siemens/meter**2
        gbar_{name} : siemens/meter**2
        m{name}_inf = 1.0 / ( 1 + exp(-(v/mV+{v12m})/{vwm}) ) : 1 
	h{name}_inf = 1.0 / ( 1 + exp((v/mV+{v12h})/{vwh}) ) : 1 
	tau_m{name} = ( {am} + 1.0 / ( exp((v/mV+{vm1})/{wm1}) + exp(-(v/mV+{vm2})/{wm2}) ) ) * ms : second
	tau_h{name} = ( {ah} + 1.0 / ( exp((v/mV+{vh1})/{wh1}) + exp(-(v/mV+{vh2})/{wh2}) ) ) * ms : second
        dm{name}/dt = -(m{name} - m{name}_inf)/tau_m{name} : 1
        dh{name}/dt = -(h{name} - h{name}_inf)/tau_h{name} : 1
        """

        super().__init__(name, params)
        
    def default_params(self):
        return dict(E_Ca = 140, # mV
	            vshift = 0, #	(mV)		: voltage shift (affects all)
	            cao  = 2.5, #	(mM)	        : external ca concentration
	            v12m=50, #         	(mV)
	            v12h=78, #         	(mV)
	            vwm =7.4, #         	(mV)
	            vwh=5.0, #         	(mV)
	            am=3, #         	(mV)
	            ah=85, #         	(mV)
	            vm1=25, #         	(mV)
	            vm2=100, #         	(mV)
	            vh1=46, #         	(mV)
	            vh2=405, #         	(mV)
	            wm1=20, #         	(mV)
	            wm2=15, #         	(mV)
	            wh1=4, #         	(mV)
	            wh2=50, #         	(mV)
                    tadj = 2.5**((34-24)/10.0))
